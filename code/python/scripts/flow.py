#!/usr/bin/env python
"""
Find paths, by solving the min-cost flow problem, through a protein-protein interaction network that connect source proteins in the  <sources_file> with targets in the <targets_file>. The protein-protein interaction network 
given by <edges_file> has proteins identified by Ensembl protein identifiers (ENSP) and consists of edges or
interactions between proteins in the ABC format. That is, lines are three, space-delimited tokens:
<ensp_a> <ensp_b> <weight>
where the weight is the strength of the interaction. We typically use https://string-db.org/.
<min_sources> constrains the flow program to require it to include the number <min_sources> source proteins in the result. Higher values find pathways shared among more of the source proteins but at the expense of the strength of interaction evidence. A value of 1 imposes no constraint. Similar remarks can be made for <min_targets>.

The flow result is output to the file flow_result.graphml. If <flow_only> is not provided, the program then
performs Gene Set Enrichment Analysis on the connected components in the output network. Each connected
component can have multiple significant enrichments, so the first 10 significant enrichments are written
to files in the <output_dir>. For example the *_0_*.gv and *_0_*.png files report the most significant enrichment
for each connected component. The *_1_* files report the next most significant enrichment, and so on. The network
does not vary across these files, only the enrichments.

Because it is often convenient to use HGNC gene identifiers, you can provide a <mapping_file> to translate ENSP identifiers used by the program to HGNC identifiers to be used in the output files.

Authors: Anthony Gitter, Chris Magnano, Aaron Baker
"""
import argparse, sys
import os, os.path
import networkx as nx
import ppi.gsea
import ppi.plot
from ortools.graph import pywrapgraph
from gprofiler import GProfiler

def remove_node(G, node):
    """
    Wrapper around networkx.classes.graph.Graph.remove_node to return a dict mapping an incident node to its edge attributes
    """
    rv = {}
    for edge in G.edges(node, data=True):
        u = edge[0]
        v = edge[1]
        edge_attr = edge[2]
        t = None
        if u == node:
            t = v
        else:
            t = u
        rv[t] = edge_attr
    G.remove_node(node)
    return rv

def parse_nodes(node_file):
    ''' Parse a list of sources or targets and return a set '''
    with open(node_file) as node_f:
        lines = node_f.readlines()
        nodes = set(map(str.strip, lines))
    return nodes


def construct_digraph(edges_file, cap):
    ''' Parse a list of weighted undirected edges.  Construct a weighted
    directed graph in which an undirected edge is represented with a pair of
    directed edges.  Use the specified weight as the edge weight and a default
    capacity of 1.
    '''
    G = pywrapgraph.SimpleMinCostFlow()
    idDict = dict() #Hold names to number ids
    curID = 0
    default_capacity = int(cap)

    with open(edges_file) as edges_f:
        for line in edges_f:
            tokens = line.strip().split()
            node1 = tokens[0]
            if not node1 in idDict:
                idDict[node1] = curID
                curID += 1
            node2 = tokens[1]
            if not node2 in idDict:
                idDict[node2] = curID
                curID += 1
            # Google's solver can only handle int weights
            w = int((1-(float(tokens[2])))*100)
            G.AddArcWithCapacityAndUnitCost(idDict[node1],idDict[node2], default_capacity, int(w))
            G.AddArcWithCapacityAndUnitCost(idDict[node2],idDict[node1], default_capacity, int(w))
    idDict["maxID"] = curID
    return G,idDict


def print_graph(graph):
    ''' Print the edges in a graph '''
    print('\n'.join(sorted(map(str, graph.edges(data=True)))))


def add_sources_targets(G, sources, targets, idDict, source_capacity, default_target_capacity, target_capacity_dict=None):
    '''
    Similar to ResponseNet, add an artificial source node that is connected
    to the real source nodes with directed edges.  Unlike ResponseNet, these
    directed edges should have weight of 0 and Infinite capacity.  Also add an
    artificial target node that has directed edges from the real target nodes
    with the same weights and capacities as the source node edges.  The new
    nodes must be named "source" and "target".

    Parameters
    ==========
    G : ortools.graph.pywrapgraph.SimpleMinCostFlow
        graph object from ortools created by construct_digraph

    sources : list
        list of hashable node names which are used as sources in the min cost flow

    targets : list
        list of hashable node names which are used as targets in the min cost flow

    idDict : dict<string, int>
        mapping of node name to node id

    source_capacity : int
        positive integer to use as source node capacity

    default_target_capacity : int
        postitive integer to use as target node capacity

    target_capacity_dict: dict<string, int>
        mapping of node name to capacity for an edge from that node to the artifical target/sink node
    '''
    default_weight = 0
    curID = idDict["maxID"]
    idDict["source"] = curID
    curID += 1
    idDict["target"] = curID

    for source in sources:
        if source in idDict:
            G.AddArcWithCapacityAndUnitCost(idDict["source"],idDict[source], source_capacity, default_weight)

    for target in targets:
        if target in idDict:
            capacity = None
            if target_capacity_dict is not None and target in target_capacity_dict:
                capacity = target_capacity_dict[target]
            else:
                capacity = default_target_capacity
            G.AddArcWithCapacityAndUnitCost(idDict[target],idDict["target"], capacity, default_weight)

def write_output_to_sif(G,out_file_name,idDict):
    ''' Convert a flow dictionary from networkx.min_cost_flow into a list
    of directed edges with the flow.  Edges are represented as tuples.
    '''

    out_file = open(out_file_name,"w")
    names = {v: k for k, v in idDict.items()}
    numE = 0
    for i in range(G.NumArcs()):
        node1 = names[G.Tail(i)]
        node2 = names[G.Head(i)]
        flow = G.Flow(i)
        if flow <= 0:
            continue
        if node1 in ["source","target"]:
            continue
        if node2 in ["source","target"]:
            continue
        numE+=1
        out_file.write(node1+"\t"+node2+"\n")
    print("Final network had %d edges" %(numE))
    out_file.close()

    return

def min_cost_flow(G, flow, idDict, roots, targets):
    ''' Use the min cost flow algorithm to distribute the specified amount
    of flow from sources to targets.  The artificial source should have
    demand = -flow and the traget should have demand = flow.  output is the
    filename of the output file.  The graph should have artificial nodes
    named "source" and "target".
    '''
    G.SetNodeSupply(idDict['source'],int(flow))
    G.SetNodeSupply(idDict['target'],int(-1*flow))

    solved = None
    print("Computing min cost flow")
    if G.Solve() == G.OPTIMAL:
        print("Solved!")
        print(G.OptimalCost())
        solved = True
    else:
        print("There was an issue with the solver")
        solved = False

    return solved

def or2nx(G, idDict):
    H = nx.DiGraph()
    names = {v: k for k, v in idDict.items()} # node id to node name
    for i in range(G.NumArcs()):
        node1 = names[G.Tail(i)]
        node2 = names[G.Head(i)]
        flow = G.Flow(i)
        if flow <= 0:
            continue
        H.add_edge(node1, node2, **{'flow': flow})
    return H

def main(args):
    ''' Parse a weighted edge list, source list, and target list.  Run
    min cost flow or k-shortest paths on the graph to find source-target
    paths.  Write the solutions to a file.
    '''
    flow_meta_outfile = os.path.join(args.outdir, 'flow_meta.tsv')
    flow_outfile = os.path.join(args.outdir, 'flow_result.gv')
    comp_enrich_map_outfile = os.path.join(args.outdir, 'comp_enrich_map.txt')
    flow_graphml = os.path.join(args.outdir, 'flow_result.graphml')
    flow_target_outfile = os.path.join(args.outdir, 'target_flow.tsv')

    flow = args.min_sources * args.min_targets
    source_capacity = args.min_targets
    default_target_capacity = args.min_sources
    other_capacity = flow

    # parse optional arguments
    other_flow = None
    target_capacity_dict = {}
    if args.target_capacity_file is not None and args.flow_meta_file is not None:
        with open(args.flow_meta_file) as fh:
            for line in fh:
                line = line.rstrip()
                tokens = line.split()
                if tokens[0] == 'flow':
                    other_flow = int(tokens[1])

        with open(args.target_capacity_file) as fh:
            for line in fh:
                line = line.rstrip()
                tokens = line.split()
                target_capacity_dict[tokens[0]] = int(tokens[1])

    # rescale flow/capacity if optional arguments are present
    sources = parse_nodes(args.sources_file)
    targets = parse_nodes(args.targets_file)
    if args.verbose:
        sys.stdout.write('before construct_digraph\n')
        sys.stdout.flush()
    G,idDict = construct_digraph(args.edges_file, other_capacity)
    add_sources_targets(G, sources, targets, idDict, source_capacity, default_target_capacity)

    # update state of G with the solution
    if args.verbose:
        sys.stdout.write('before solve\n')
        sys.stdout.flush()
    solved = min_cost_flow(G, flow, idDict, sources, targets)
    if args.verbose:
        sys.stdout.write('after solve\n')
        sys.stdout.flush()

    if not solved:
        sys.stderr.write('Could not solve\n')
        if args.no_exit_on_fail:
            sys.exit(0)
        else:
            sys.exit(21)

    H = or2nx(G, idDict)
    source_edges_dict = remove_node(H, 'source')
    target_edges_dict = remove_node(H, 'target')

    # write flow result graphml
    nx.write_graphml(H, flow_graphml)

    if not args.flow_only:
      # perform GSEA on the connected components in the flow result graph
      enrich_dir = os.path.join(args.outdir, 'enrich')
      if not os.path.exists(enrich_dir):
          os.mkdir(enrich_dir)
      set_fp_pairs = ppi.gsea.gsea_connected_components(H, enrich_dir)

      # document which component is associated with which enrichment result
      # its a ragged csv where each line is a connected component
      # the first column is the file path to the enrichment
      # the 2..N-1 column is for gene 1, gene 2, ..., gene N in the connected component
      with open(comp_enrich_map_outfile, 'w') as fh:
       for gene_set, fp in set_fp_pairs:
         fh.write(",".join([fp] + list(gene_set)) + "\n")

      # write flow result graphviz
      weights = None
      if args.node_weights is not None:
        # parse --node-weights
        weights = {}
        with open(args.node_weights, 'r') as fh:
            for line in fh:
                line = line.rstrip()
                node, weight = line.split(',')
                weights[node] = float(weight)
      ppi.plot.vis_node_clusters_gv(H, open(flow_outfile, 'w'), sources, targets, weights=weights)

      # write enrichment graphviz
      if args.visualization == 'single':
          ppi.plot.vis_single_community(H, sources, targets, set_fp_pairs, args)

      elif args.visualization == 'multi':
          enrich_fps = list(map(lambda x: x[1], set_fp_pairs))
          ppi.plot.vis_multi_community(H, sources, targets, enrich_fps, args, weights=weights)

      with open(flow_meta_outfile, 'w') as fh:
          fh.write('{}\t{}\n'.format('flow', flow))

def add_flow_args(parser):
    parser.add_argument('--mapping-file')
    parser.add_argument('--edges-file',
                        help='edge file path with weights in [0,1]',
                        type=str,
                        required=True)
    parser.add_argument('--sources-file',
                        help='source node file path',
                        type=str,
                        required=True)
    parser.add_argument('--targets-file',
                        help='target node file path',
                        type=str,
                        required=True)
    parser.add_argument('--outdir',
                        help='Output file',
                        type=str,
                        required=True)
    parser.add_argument('--min-sources',
                        help='Minimum number of sources that flow must pass through',
                        type=int, 
                        required=True)
    parser.add_argument('--min-targets',
                        help='Minimum number of targets that flow must pass through',
                        type=int,
                        required=True)
    parser.add_argument('--target-capacity-file',
                        help='tsv file with lines as <node> <capacity> which specifies the capacity from that node to the sink',
                        type=str)
    parser.add_argument('--node-weights',
                        help='csv file with lines as <node>,<weight> where smaller weights result in larger nodes')
    parser.add_argument('--flow-meta-file',
                        help='tsv file with attribute value pairs of flow parameters including \'flow\'',
                        type=str)
    parser.add_argument('--verbose', '-v',
                        help='if set, produce verbose output',
                        action='store_true')
    parser.add_argument('--no-exit-on-fail',
                        help="If set, produce exit code zero even if min cost flow cannot be solved (useful for ignoring simulation failures)",
                        action='store_true')
    parser.add_argument('--visualization',
                        help='Visualization style. One of "single" or "multi"; default "multi".',
                        type=str,
                        default='multi')
    parser.add_argument("--image-format",
                        help="Output image file format for network images: either \"png\" or \"svg\". Default \"svg\".",
                        type=str,
                        default='svg')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    add_flow_args(parser)
    parser.add_argument('--flow-only', help='Do not perform GSEA and subsequent visualizations, just write the flow_result.graphml', action='store_true')
    args = parser.parse_args()
    main(args)
