#!/usr/bin/env python
import sys, argparse
import numpy as np
import ppi.parsers.abc
import flow
import os, os.path

def main():
  parser = argparse.ArgumentParser(description="""
Simulate genetic screening hits by uniform random sampling.
If --alt-sources-file is present, hits are sub-sampled from the alternative sources file.
Otherwise, hits are simulated from the provided network nodes in --edges-file.
""")
  parser.add_argument('--n-simulation', type=int, required=True)
  parser.add_argument('--sources-file', required=True)
  parser.add_argument('--edges-file', help='network file in abc format (node node weight)', required=True)
  parser.add_argument('--alt-sources-file', help='newline-delimited list of gene identifiers')
  parser.add_argument('--outdir', required=True)
  args = parser.parse_args()

  G = ppi.parsers.abc.parse_abc(args.edges_file)
  nodes = None # universe to sample from
  if args.alt_sources_file is not None:
    nodes = []
    with open(args.alt_sources_file) as fh:
      for line in fh:
        line = line.rstrip()
        nodes.append(line)
    n_nodes_pre = len(nodes)
    nodes = list(filter(lambda x: x in G, nodes))
    n_nodes_post = len(nodes)
    sys.stderr.write('[warning] {} nodes were removed from --alt-sources-file because they are not present in the --edges-file network\n'.format(n_nodes_pre - n_nodes_post))
  else:
    nodes = G.nodes()
  nodes = sorted(nodes)

  # note real sources file
  # TODO what if some of the input sources are not actually in the edges-file graph
  sources_real_fp = args.sources_file
  sources_real = flow.parse_nodes(sources_real_fp)

  # simulate screens
  for i in range(args.n_simulation):
    sample_inds = np.random.choice(len(nodes), size=len(sources_real), replace=False)
    sample_nodes = list(map(lambda x: nodes[x], sample_inds))
    ofp = os.path.join(args.outdir, 'sim{}.txt'.format(i))
    with open(ofp, 'w') as ofh:
      ofh.write('\n'.join(sample_nodes))

if __name__ == "__main__":
  main()
