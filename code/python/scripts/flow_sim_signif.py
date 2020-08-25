#!/usr/bin/env python
import sys, argparse
import os, os.path
import networkx as nx

def main():
  parser = argparse.ArgumentParser(description="""
Identify nodes in the resulting flow graph which are signifcant according to the simulated 
empirical distribution.
""")
  parser.add_argument('--flow-result', required=True)
  parser.add_argument('--node-frequency', required=True)
  parser.add_argument('--outdir', required=True)
  args = parser.parse_args()

  G = nx.read_graphml(args.flow_result)
  node_to_count = {}
  with open(args.node_frequency, 'r') as fh:
    for line in fh:
      line = line.rstrip()
      words = line.split(',')
      node = words[0]
      count = int(words[1])
      node_to_count[node] = count

  node_count_pairs = []
  for node in G.nodes():
    count = node_to_count[node] 
    node_count_pairs.append((node, count))

  # sort by count, ascending; most significant first
  node_count_pairs = sorted(node_count_pairs, key=lambda x: x[1])
  signif_ofp = os.path.join(args.outdir, 'signif.csv')
  with open(signif_ofp, 'w') as ofh:
    for node, count in node_count_pairs:
      ofh.write('{},{}\n'.format(node, count))

if __name__ == "__main__":
  main()
