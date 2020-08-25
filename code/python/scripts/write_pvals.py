#!/usr/bin/env python
import sys, argparse
import networkx as nx

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--node-frequency", "-q", help="2-column csv")
  parser.add_argument("--n-simulations", "-n", type=int, help="int")
  parser.add_argument("--flow-result", "-f", help="graphml")
  parser.add_argument("--outfile", "-o", help="output file")
  args = parser.parse_args()

  G = nx.read_graphml(args.flow_result)

  ofh = open(args.outfile, 'w')
  with open(args.node_frequency) as fh:
    for line in fh:
      line = line.rstrip()
      [id_v, freq_str] = line.split(",")
      if id_v in G:
        p_val = int(freq_str) / args.n_simulations
        ofh.write(",".join([id_v, "{:1.3f}".format(p_val)])+"\n")

if __name__ == "__main__":
  main()
