#!/usr/bin/env python
import sys, argparse
import os, os.path
import networkx as nx
import ppi.parsers.abc

def main():
  parser = argparse.ArgumentParser(description="""
Compute node frequency of nodes in the resulting flow graph for many simulated runs.
""")
  parser.add_argument('--flow-results', nargs='+', help="Flow result graphml files", required=True)
  parser.add_argument('--edges-file', required=True)
  parser.add_argument('--outdir')
  args = parser.parse_args()

  node_to_count = {}
  G_ppi = ppi.parsers.abc.parse_abc(args.edges_file)
  for node in G_ppi.nodes():
    node_to_count[node] = 0

  for flow_result_fp in args.flow_results:
    if os.path.exists(flow_result_fp):
      G = nx.read_graphml(flow_result_fp)
      for node in G.nodes():
        node_to_count[node] += 1
    else:
      sys.stdout.write('[warning] missing flow result {}\n'.format(flow_result_fp))

  ofp = os.path.join(args.outdir, 'node_frequency.csv')
  with open(ofp, 'w') as ofh:
    for node, count in node_to_count.items():
      ofh.write('{},{}\n'.format(node, count))

if __name__ == "__main__":
  main()
