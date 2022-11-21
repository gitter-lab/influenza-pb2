from .gsea import enrichs_to_gv, read_enrich_line, score_enrichment
from . import SimPathException
import re
import subprocess
import sys
import os, os.path
import math
ROOT_COLOR = 'red'
TARGET_COLOR = 'gray'
ROOT_AND_TARGET_COLOR = 'red2'

# graph size and node width are in inches and are relative to each other
POINTS_PER_INCH = 72 # https://graphviz.gitlab.io/_pages/doc/info/attrs.html#points
GRAPH_SIZE = 50
MIN_WIDTH = 0.5
MAX_WIDTH = 2.5
MIN_WEIGHT = 0
MAX_WEIGHT = 1000 # TODO take this number as input --n-simulations
WIDTH_INTERP_MAX = -math.log(1/MAX_WEIGHT,10)
WIDTH_INTERP_SCALAR = (MAX_WIDTH - MIN_WIDTH) / WIDTH_INTERP_MAX
WEIGHT_LEN = MAX_WEIGHT - MIN_WEIGHT
WIDTH_LEN = MAX_WIDTH - MIN_WIDTH

def write_node_attrs(ofh, nodes, **attrs):
  attr_strs = []
  for key, value in attrs.items():
    attr_strs.append("{}={}".format(key, value))
  attr_str = ", ".join(attr_strs)
  for node in nodes:
    ofh.write("{} [{}]\n".format(node, attr_str))

def write_node_color(ofh, color, nodes):
  for node in nodes:
    if color in ['navyblue', 'indigo']:
      # use white fontcolor for dark fillcolors
      ofh.write("{} [fillcolor={}, style=filled, fontcolor={}]\n".format(node, color, 'white'))
    else:
      # use default black fontcolor
      ofh.write("{} [fillcolor={}, style=filled]\n".format(node, color))

def write_node_width(ofh, node_to_width):
  for node, width in node_to_width.items():
    ofh.write("{} [width={}, fontsize={}]\n".format(node, width, width*10))

def interpolate_weight_to_width(weight):
    weight_pos = weight / WEIGHT_LEN # fractional position on weight interval
    # smaller weight is more significant
    # adjust 0 -> 1 for log scale
    if weight_pos == 0:
      weight_pos = 1 / WEIGHT_LEN
    # the scalar below, for example, -log_10(weight) by design ranges from 0 to 3
    # 2/3 adjusts range from [0,3] to [0,2]
    width = -math.log(weight_pos,10) * WIDTH_INTERP_SCALAR + MIN_WIDTH
    return width

def compute_node_width(G, weights):
  """
  Scale node sizes (width of node shape) based on node significance scores in <weights>

  Parameters
  ----------
  G : nx.Graph
  
  weights : mapping of node in G to a numeric weight
  """
  node_to_width = {}
  for node in G.nodes():
    weight = weights.get(node)
    if weight is None:
      sys.stderr.write('[warning] no weight for node {}, using max\n'.format(node))
      weight = MAX_WEIGHT
    # linearly interpolate weight to width
    node_to_width[node] = interpolate_weight_to_width(weight)
  return node_to_width

def compute_node_width_legend():
  """
  Return a mapping of p-value to node width for a set of reference weights

  TODO details about hypothesis test here
  """
  reference_weights = [1, 10, 100, 500, 1000]
  rv = {}
  for weight in reference_weights:
    p_val = weight / MAX_WEIGHT
    width = interpolate_weight_to_width(weight)
    rv[p_val] = width
  return rv

def vis_legend_gv(ofh, node_to_width):
  """
  Prepare a Graphviz file to plot a legend in <node_to_width>, which is the output of compute_node_width_legend

  Parameters
  ----------
  ofh : output file handle, io-like
    location to write Graphviz data to
  node_to_width : dict
    mapping of node to its width
  """
  indent = "  "
  ofh.write(indent + "subgraph cluster_legend {\n")
  ofh.write(indent*2 + "label = \"Node sizes for \\nreference p-values\"\n")
  ofh.write(indent*2 + "fontsize = \"30\"\n")
  nodes = list(map(lambda x: str(x+1), range(len(node_to_width.keys()))))
  ofh.write(indent*2 + " ".join(nodes))
  ofh.write("\n" + indent + "}\n")
  # TODO order?
  i = 0
  for node, width in node_to_width.items():
    i = i+1
    ofh.write(indent + "{} [label=\"{:1.2f}\", width={}]\n".format(i, node, width))
  for i in range(len(node_to_width)-1):
    ofh.write(indent + "{} -> {} [style=invis]\n".format(i+1, i+2))

def vis_node_clusters_gv(G, ofh, roots=[], targets=[], cluster_label_pairs=[], weights=None):
  """
  Prepare a Graphviz file with special colors for <roots> and <targets> and with boxes identifying sets of nodes in <cluster_label_pairs>

  Parameters
  ----------
  cluster_label_pairs : (set, str)
    pairs of node members in cluster and the cluster label

  weights : None or dict
    if dict, a mapping of node name to a weight. a smaller weight results in a larger node.
  """
  def write_cluster(ofh, cluster_no, cluster_tpl, indent=2):
    ws = " "*indent
    cluster, label = cluster_tpl
    ofh.write("{}subgraph cluster_{} {{\n".format(ws, cluster_no))
    ofh.write("{}label = \"{}\"\n".format(ws*2, label))
    ofh.write("{}fontsize = \"{}\"\n".format(ws*2, MAX_WIDTH*10))
    ofh.write("{}{}\n".format(ws*2, " ".join(cluster)))
    ofh.write("{}}}\n".format(ws))

  roots = set(filter(lambda x: x in G, roots))
  targets = set(filter(lambda x: x in G, targets))
  roots_and_targets = roots.intersection(targets)
  roots = roots - roots_and_targets
  targets = targets - roots_and_targets

  # write graphviz
  if G.is_directed():
    ofh.write("digraph G {\n")
  else:
    ofh.write("graph G {\n")
  # standardize graph and node size
  ofh.write("graph [size={}]\n".format(GRAPH_SIZE))
  ofh.write("node [regular=true, fixedsize=true, width={}]\n".format(MIN_WIDTH))
  if len(cluster_label_pairs) != 0:
    # write clusters
    cluster_no = 0
    for cluster_label_pair in cluster_label_pairs:
      write_cluster(ofh, cluster_no, cluster_label_pair)
      cluster_no += 1
  # highlight sources and targets
  write_node_color(ofh, ROOT_COLOR, roots.union(roots_and_targets))
  write_node_color(ofh, TARGET_COLOR, targets)
  write_node_color(ofh, ROOT_AND_TARGET_COLOR, roots_and_targets)

  # write node widths
  if weights is not None:
    node_to_width = compute_node_width(G, weights)
    write_node_width(ofh, node_to_width)

  # write node ranks
  # e.g. { rank=same; 1; A;}
  if len(roots) > 0:
    ofh.write("{{ rank=min; {}; }}\n".format("; ".join(roots)))
  #ofh.write("{{ rank=max; {}; }}\n".format("; ".join(targets)))

  # write edges
  if G.is_directed():
    for edge in G.edges():
      ofh.write("{} -> {};\n".format(edge[0], edge[1]))
  else:
    for edge in G.edges():
      ofh.write("{} -- {};\n".format(edge[0], edge[1]))

  # write legend
  pval_to_width = compute_node_width_legend()
  vis_legend_gv(ofh, pval_to_width)

  # finalize gv graph
  ofh.write("}")

def vis_node_box_gv(G, ofh, roots=[], targets=[], cluster_label_pairs=[], weights=None):
  """
  Alternative to vis_node_clusters_gv. Instead of placing a box around nodes in the cluster,
  identify those nodes with a "box" node shape instead.
  """
  if len(cluster_label_pairs) != 1:
    raise SimPathException("cluster_label_pairs must contain only one element")

  roots = set(filter(lambda x: x in G, roots))
  targets = set(filter(lambda x: x in G, targets))
  roots_and_targets = roots.intersection(targets)
  roots = roots - roots_and_targets
  targets = targets - roots_and_targets

  # write graphviz
  if G.is_directed():
    ofh.write("digraph G {\n")
  else:
    ofh.write("graph G {\n")
  # standardize graph and node size
  ofh.write("graph [size={}]\n".format(GRAPH_SIZE))
  ofh.write("node [regular=true, fixedsize=true, width={}]\n".format(MIN_WIDTH))
  for nodes, label in cluster_label_pairs:
    for node in nodes:
      ofh.write("{} [shape=box]\n".format(node))

  # highlight sources and targets
  write_node_color(ofh, ROOT_COLOR, roots.union(roots_and_targets))
  write_node_color(ofh, TARGET_COLOR, targets)
  write_node_color(ofh, ROOT_AND_TARGET_COLOR, roots_and_targets)

  # write node widths
  if weights is None:
    weights = {}
    for node in G.nodes():
      weights[node] = (MAX_WEIGHT - MIN_WEIGHT) / 2
  node_to_width = compute_node_width(G, weights)
  write_node_width(ofh, node_to_width)

  # write node ranks
  # e.g. { rank=same; 1; A;}
  ofh.write("{{ rank=min; {}; }}\n".format("; ".join(roots)))
  #ofh.write("{{ rank=max; {}; }}\n".format("; ".join(targets)))

  # write edges
  for edge in G.edges():
    ofh.write("{} -> {};\n".format(edge[0], edge[1]))
  ofh.write("}")

def vis_single_community(G_prime, roots, targets, comp_enrich_pairs, args):
  """
  Create a network visualization for each pair in <comp_enrich_pairs>

  Parameters
  ----------
  G_prime: nx.Graph
    network to visualize

  roots : list
    list of nodes in G_prime which should be colored by ROOT_COLOR

  targets : list
    list of nodes in G_prime which should be colored by TARGET_COLOR

  comp_enrich_pairs : (set, str)
    pair of connected component in G_prime (as a node set) and a file with the G:Profiler enrichment results

  args : namespace
    namespace with the following elements:
    mapping_file - file path 
    outdir - directory to write results
    TODO change this to explicit arguments
  """
  n_plots = 3

  def get_valid_filename(s):
    s = str(s).strip().replace(' ', '_')
    return re.sub(r'(?u)[^-\w.]', '', s)

  gv_fps = []

  # prepare a set of enrichment figures for each node set/comp
  # each figure in the set corresponds to enrichments from GSEA
  i = 0
  for comp, enrich_fp in comp_enrich_pairs:
    i += 1
    G_sub = G_prime.subgraph(comp)
    comp_dir_fp = os.path.join(args.outdir, 'comp_{}'.format(i))
    if not os.path.exists(comp_dir_fp):
      os.mkdir(comp_dir_fp)

    with open(enrich_fp, 'r') as fh:
      # score enrichments and write enrichment scores to a file
      # skip header
      for line in fh:
        break
      csv_data = []
      tname_to_plot_elems = {}
      for line in fh:
        t_name, anno_hits, p_val, label = read_enrich_line(line)
        tname_to_plot_elems[t_name] = (anno_hits, label)
        csv_datum = (t_name, p_val, score_enrichment(line.rstrip().split('\t')))
        csv_data.append(csv_datum)
      comp_scores_fp = os.path.join(comp_dir_fp, "enrich_scores.csv")
      csv_data = sorted(csv_data, key=lambda x: x[2], reverse=True)
      with open(comp_scores_fp, 'w') as csv_ofh:
        for csv_datum in csv_data:
          csv_ofh.write(",".join(map(str, csv_datum)) + "\n")

      # plot the top enrichments
      for j in range(n_plots):
        csv_datum = csv_data[j]
        t_name = csv_datum[0]
        anno_hits, label = tname_to_plot_elems[t_name]

        gv_fp = os.path.join(comp_dir_fp, get_valid_filename(t_name) + ".gv")
        with open(gv_fp, 'w') as ofh:
          vis_node_box_gv(G_sub, ofh, roots=roots, targets=targets, cluster_label_pairs=[(anno_hits, label)])
        gv_fps.append(gv_fp)

  # generate figures from the graphviz source files
  process_gvs(gv_fps, args, gv_prog="dot")

def vis_multi_community(G_prime, roots, targets, enrich_fps, args, weights=None):
  max_enrich_gv = 10 # TODO other stopping condition

  # TODO how to best visualize different enrichment results? 
  # for now use n pngs to display n different annotation results for the same set of mcl communities
  enrich_fhs = []
  gv_out_fps = []
  for enrich_fp in enrich_fps:
    enrich_fh = open(enrich_fp, 'r')
    enrich_fh.__next__() # skip header
    enrich_fhs.append(enrich_fh)
  for i in range(max_enrich_gv):
    gv_out_fp = os.path.join(args.outdir, "graph_ensp_{}.gv".format(i))
    gv_out_fps.append(gv_out_fp)
    with open(gv_out_fp, 'w') as ofh:
      exhausted_inds, cluster_members = enrichs_to_gv(G_prime, enrich_fhs, ofh)
      vis_node_clusters_gv(G_prime, ofh, roots=roots, targets=targets, cluster_label_pairs=cluster_members, weights=weights)

      print("exhausted_inds: " + ",".join(map(str, exhausted_inds)))
      #sys.stderr.write("[warning] only wrote {} graphviz files, less than the requested {}\n".format(i, max_enrich_gv))
  for enrich_fh in enrich_fhs:
    enrich_fh.close()
  process_gvs(gv_out_fps, args)

def process_gvs(gv_fps, args, gv_prog="fdp"):
  """
  Translate ENSP to HGNC identifiers and process Graphviz input files with its layout program, fdp

  Parameters
  ----------
  args : namespace
    args.outdir
    args.mapping_file
    args.image_format

  gv_prog : str
    alternative program to use e.g. "dot"

  TODO
  ----
  probably dont want to have this subprocess stuff in here, not very portable for a library
  """
  if gv_prog not in ["fdp", "dot"]:
    raise SimPathException("Invalid gv_prog: {}".format(gv_prog))

  if args.image_format not in ["png", "svg"]:
    raise SimPathException("Invalid image_format: {}".format(args.image_format))

  # translate ensp identifiers to hgnc if mapping file is provided
  gv_infiles = []
  if(args.mapping_file is not None):
    for i in range(len(gv_fps)):
      gv_out_fp = gv_fps[i]
      gv_ensp_fp = gv_fps[i]
      gv_ensp_dn = os.path.dirname(gv_ensp_fp)
      gv_ensp_bn = os.path.basename(gv_ensp_fp)
      gv_ensp_bn_noext, ext = os.path.splitext(gv_ensp_bn)
      gv_hgnc_fp = os.path.join(gv_ensp_dn, gv_ensp_bn_noext + "_hgnc.gv")
      mapping_args = ["map_ensp_to_hgnc.py", "--for-graphviz", "--mapping-file", args.mapping_file, "--infile", 
        gv_out_fp, "--outfile", gv_hgnc_fp]
      stdout_fh = open(os.path.join(args.outdir, "map_gv_{}.out".format(i)), "w")
      stderr_fh = open(os.path.join(args.outdir, "map_gv_{}.err".format(i)), "w")
      sys.stdout.write("[STATUS] Launching {}\n".format(str(mapping_args)))
      process = subprocess.Popen(mapping_args, stdout=stdout_fh, stderr=stderr_fh)
      pid = process.pid
      sys.stdout.write("[STATUS] pid {}\n".format(pid))
      exit_code = process.wait()
      if(exit_code != 0):
        stdout_fh.flush()
        stderr_fh.flush()
        raise Exception("{} failed with exit code {}\n".format(mapping_args[0], exit_code))
      gv_infiles.append(gv_hgnc_fp)

  else:
    gv_infiles = gv_fps

  # process graphviz files
  for gv_infile in gv_infiles:
    # TODO dot or fdp?
    # TODO rank is not used in fdp but clusters are not used in dot
    gv_args = [gv_prog, "-T{}".format(args.image_format)]
    fp_no_ext, ext = os.path.splitext(gv_infile)
    png_fp = fp_no_ext + ".png"
    sys.stdout.write("[STATUS] Launching {} < {} > {}\n".format(str(gv_args), gv_infile, png_fp))
    process = subprocess.Popen(gv_args, stdin=open(gv_infile, 'r'), stdout=open(png_fp, 'wb'))
    pid = process.pid
    sys.stdout.write("[STATUS] pid {}\n".format(pid))
    exit_code = process.wait()
    if(exit_code != 0):
      raise Exception("{} failed with exit code {}\n".format(gv_args[0], exit_code))
