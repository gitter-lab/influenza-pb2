#!/usr/bin/env python
import argparse, sys
import ppi.parsers.abc
 
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--graph-type", type=str, default="abc")
  parser.add_argument("--graph", type=argparse.FileType("r"), required=True)
  parser.add_argument("--gene-list", type=argparse.FileType("r"), default=sys.stdin)
  parser.add_argument("--found", type=argparse.FileType("w"), default=sys.stdout)
  parser.add_argument("--missing", type=argparse.FileType("w"))
  args = parser.parse_args()

  if (args.graph_type is None):
    # infer graph type
    sys.stderr.write("Not implemented")
    sys.exit(1)
  else:
    graph_type = args.graph_type 

  # parse graph
  G = None
  if(graph_type == "abc"):
    G = ppi.parsers.abc.parse_abc(args.graph)
  else:
    sys.stderr.write("Not implemented")
    sys.exit(2)

  # parse gene list
  gene_list = []
  for line in args.gene_list:
    line = line.rstrip()
    gene_list.append(line)

  for gene in gene_list:
    if gene in G:
      args.found.write("{}\n".format(gene))
    else:
      if args.missing is not None:
        args.missing.write("{}\n".format(gene))

if __name__ == "__main__":
  main()
