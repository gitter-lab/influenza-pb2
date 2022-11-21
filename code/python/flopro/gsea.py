from gprofiler import GProfiler
import math
import networkx as nx
import os, os.path

HEADER = [
  '#',
  'signf',
  'p-value',
  'T',
  'Q',
  'Q&T',
  'Q&T/Q',
  'Q&T/T',
  'term ID',
  't type',
  't group',
  't name',
  't depth',
  'Q&T list'
]
FIELD_TO_INDEX = {}
for i in range(len(HEADER)):
  FIELD_TO_INDEX[HEADER[i]] = i

# TODO handle empty enrichment result somewhere
def write_enrich(enrich, enrich_out_fp):
  """
  Write enrichment results from GProfiler to a file

  Parameters
  ----------
  enrich : list of list
    result of GProfiler.gprofile

  enrich_out_fp : str
    file path to write results to
  """
  enrich_sorted = sorted(enrich, key=score_enrichment, reverse=True)
  with open(enrich_out_fp, 'w') as efh:
    efh.write("#" + "\t".join(HEADER) + "\n")
    for rec in enrich_sorted:
      efh.write("\t".join(map(str, rec)) + "\n")

def read_enrich_line(line):
  line = line.rstrip()
  words = line.split('\t')
  p_val = float(words[FIELD_TO_INDEX['p-value']])
  t_name = words[FIELD_TO_INDEX['t name']]
  anno_hits = words[FIELD_TO_INDEX['Q&T list']].split(",")
  label = "{}\\np = {}".format(t_name, p_val)

  return t_name, anno_hits, p_val, label

def score_enrichment(datum):
  """
  Sorting function for enrichment results. Prioritize some combination of depth, p-value, 
  and percent of query in enrichment term (Q&T/T)

  Parameters
  ----------
  datum : list
    single record from GProfiler; fields are named in HEADER

  Returns
  -------
  vv : float
    sort key value (ascending)
  """
  # consider depth of 10, p-value of 10e-10, and fraction of 0.8 to be equivalently good
  # with default parameters
  p_val_comp = - math.log(float(datum[FIELD_TO_INDEX['p-value']]))
  depth_comp = float(datum[FIELD_TO_INDEX['t depth']])

  # goodness of enrichment fraction decays somewhat slowly but according to a quadratic function
  enrich_frac_comp = (float(datum[FIELD_TO_INDEX['Q&T/T']]) * 10.0/8.0) ** 2 * 10
  vv = p_val_comp + depth_comp + enrich_frac_comp
  return vv

def score_enrichments(enrich_recs, n_enrich=10):
  """
  Score top <n_enrich> records according to score_enrichment
  """
  total_score = 0.0
  enrich_recs = sorted(enrich_recs, key=score_enrichment, reverse=True)
  for enrich_rec in enrich_recs[:n_enrich]:
     score = score_enrichment(enrich_rec)
     total_score += score
  return total_score

def gsea_connected_components(G, outdir):
  """
  Perform Gene Set Enrichment Analysis on the connected components in G using GProfiler

  Returns
  -------
  rv : list of (set, str)
    tuples of gene set that was queried for enrichment and the enrichment output file
  """
  rv = []
  gp = GProfiler("FluPath/0.1")
  if nx.is_directed(G):
    G = G.to_undirected()
  comps = list(nx.connected_components(G))
  comp_no = 0
  for comp in comps:
    # TODO how are http errors handled?
    enrich_out_fp = os.path.join(outdir, "enrich_{}.tsv".format(comp_no))
    if not os.path.exists(enrich_out_fp):
      enrich = gp.gprofile(comp, src_filter=['GO:BP'])
      write_enrich(enrich, enrich_out_fp)
    rv.append((comp, enrich_out_fp))
    comp_no += 1 
  return rv

def enrichs_to_gv(G, enrich_fhs, ofh, **kwargs):
  """
  Read the next line in each file handle in <enrich_fhs> and parse the significant GO:BP
  annotation in that line. Then, write the GO:BP annotations, along with the original graph,
  to a graphviz open file handle <ofh>.

  Parameters
  ----------
  G : nx.Graph

  enrich_fhs : list of io-like
    GO:BP enrichment results from GProfiler, written by write_enrich

  ofh : io-like

  Returns
  -------
  exhausted

  TODO
  ----
  - handle case where any fh does not have next
  - this could have used the data from GProfiler python API directly to be faster but 
  this may be more reusable later
  """
  # parse annotation data to prepare graphviz "clusters"
  cluster_members = []
  enrich_ind = 0
  exhausted_inds = []
  for i in range(len(enrich_fhs)):
    enrich_fh = enrich_fhs[i]
    line = None
    try:
      line = enrich_fh.__next__().rstrip()
    except StopIteration as err:
      exhausted_inds.append(i)
    if line is not None:
      words = line.split('\t')
      p_val = float(words[FIELD_TO_INDEX['p-value']])
      t_name = words[FIELD_TO_INDEX['t name']]
      anno_hits = words[FIELD_TO_INDEX['Q&T list']].split(",")
      cluster = anno_hits
      label = "{}\\np = {}".format(t_name, p_val)
      tpl = (cluster, label)
      cluster_members.append(tpl)
    enrich_ind += 1

  return exhausted_inds, cluster_members
