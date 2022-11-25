#!/usr/bin/env python
import argparse
from flopro import irefindex as iref

def main():
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("-t", "--type", nargs=1, type=str, help="Identifier type. Default 'hgnc'.", default='hgnc')
  parser.add_argument("irefindex", nargs=1, type=argparse.FileType('r'), help="iRefIndex database file")
  parser.add_argument("id_file", nargs=1, type=argparse.FileType('r'), help="newline delimited file of identifiers")
  parser.add_argument("ids_match", nargs=1, type=argparse.FileType('w'), help="output file to write with identifiers found in irefindex")
  parser.add_argument("ids_missing", nargs=1, type=argparse.FileType('w'), help="output file to write with identifiers found in irefindex")
  args = parser.parse_args()

  all_hgncs = iref.get_hgnc_universe(args.irefindex[0])
  cand_hgncs = map(lambda x: x.rstrip(), args.id_file[0].readlines())
  founds, missings = iref.verify_hgnc(all_hgncs, cand_hgncs)

  for found in founds:
    args.ids_match[0].write(found + "\n")
  for missing in missings:
    args.ids_missing[0].write(missing + "\n")

if __name__ == "__main__":
  main()
