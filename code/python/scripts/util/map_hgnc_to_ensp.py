#!/usr/bin/env python
import sys, argparse
from flopro import ensembl
from flopro.string_db import parse_string_fh

def main():
  parser = argparse.ArgumentParser(description="""
Convert HGNC gene identifiers to Ensembl protein identifiers. 
This is a one to many mapping. 
For each HGNC gene identifier, Include all mapped Ensembl protein identifiers.
""")
  # TODO use biomart python API instead?
  parser.add_argument("--network-file", "-n", help="ENSP protein interaction network e.g. STRING (many genes have evidence reported only for one of the proteins they encode)")
  parser.add_argument("--map-file", "-m", help="Output from dl_ensembl_map.R", type=argparse.FileType('r'))
  parser.add_argument("--gene-list", "-i", help="Newline delimited file of HGNC identifiers", type=argparse.FileType('r'))
  parser.add_argument("--outfile", "-o", help="File to write gene list of Ensembl protein IDs to", type=argparse.FileType('w'), default=sys.stdout)
  parser.add_argument("--out-type", help="If \"list\", print a list of identifiers to stdout. If \"map\", print a 2-column CSV file where the first column is HGNC and the second column is ENSP", default="list")
  args = parser.parse_args()

  G = None
  if(args.network_file is not None):
    G = parse_string_fh(open(args.network_file))

  hgnc_to_ensp_map = ensembl.map_hgnc_to_ensps(args.map_file)

  if(args.gene_list is not None):
    ensp_ids = ensembl.apply_map(hgnc_to_ensp_map, args.gene_list)
    for ensp_id in ensp_ids:
      args.outfile.write("{}\n".format(ensp_id))
  else:
    hgnc_to_ensp_map_filt = {}
    for hgnc, ensps in hgnc_to_ensp_map.iteritems():
      for ensp in ensps:
        if ensp in G:
          if hgnc in hgnc_to_ensp_map_filt:
            hgnc_to_ensp_map_filt[hgnc].append(ensp)
          else:
            hgnc_to_ensp_map_filt[hgnc] = [ensp]

    # TODO move this case up one level to include args.gene_list is not None case
    if(args.out_type == "list"):
      for hgnc, ensps in hgnc_to_ensp_map_filt.iteritems():
        args.outfile.write("\n".join(ensps))
    elif(args.out_type == "map"):
      for hgnc, ensps in hgnc_to_ensp_map_filt.iteritems():
        for ensp in ensps:
          args.outfile.write("{},{}\n".format(hgnc,ensp))

if __name__ == "__main__":
  main()
