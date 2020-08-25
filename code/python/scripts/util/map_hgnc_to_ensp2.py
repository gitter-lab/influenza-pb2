#!/usr/bin/env python
import os, os.path
import argparse, sys
import re
from ppi import ensembl
from ppi import script_utils
from ppi.string_db import parse_string_fh

def main():
  parser = argparse.ArgumentParser(description="""
Translate HGNC words in <infile> to ENSP. If HGNC cannot be mapped, leave it.
""")
  parser.add_argument("--mapping-file", "-m", help="See dl_ensembl_map.R", required=True, type=argparse.FileType('r'))
  parser.add_argument("--network-file", "-n", help="ENSP protein interaction network e.g. STRING (many genes have evidence reported only for one of the proteins they encode)", required=True, type=argparse.FileType('r'))
  script_utils.add_file_and_dir_args(parser)
  args = parser.parse_args()
  io_pairs = script_utils.check_file_and_dir_args(args)

  hgnc_to_ensps_map = ensembl.map_hgnc_to_ensps(args.mapping_file)
  G = parse_string_fh(args.network_file)
  hgnc_to_ensp_map, copy_map = ensembl.filter_one_to_many(hgnc_to_ensps_map, G)
  for k, v in copy_map.items():
    sys.stderr.write("[warning] {} ENSPs associated to HGNC:{} found in <network-file>\n".format(v, k))

  sep_re = re.compile('\W')
  word = ""
  for infile, outfile in io_pairs:
    with open(infile, "r") as ifh:
      with open(outfile, "w") as ofh:
        for line in ifh:
          for char in line:
            match_data = sep_re.match(char)
            if(match_data is not None):
              # then process word
              ensp = hgnc_to_ensp_map.get(word)
              if ensp is None:
                ofh.write(word + match_data.group(0))
              else:
                ofh.write(ensp + match_data.group(0))
              word = ""
            else:
              word += char

if __name__ == "__main__":
  main()
