#!/usr/bin/env python
import argparse, sys
import re

# TODO see lib ensembl.py instead
def parse_mapping_file(fp):
  """
  """
  ensp_to_hgnc = {}
  with open(fp, "r") as fh:
    for line in fh:
      line = line.rstrip()
      words = line.split()
      if(len(words) >= 4):
        ensp_to_hgnc[words[3]] = words[0]

  return ensp_to_hgnc

# TODO specific to graphviz because some HGNC have "-" which this quotes
def main():
  parser = argparse.ArgumentParser(description="""
Translate ENSP words in <infile> to HGNC. If ENSP cannot be mapped, leave it.
""")
  parser.add_argument("--mapping-file", "-m", help="See dl_ensembl_map.R")
  parser.add_argument("--for-graphviz", action='store_true', help="If provided, will sanitize HGNC identifiers with the character \"-\" by quoting them")
  parser.add_argument("--infile", "-i")
  parser.add_argument("--outfile", "-o")
  args = parser.parse_args()

  ensp_to_hgnc = parse_mapping_file(args.mapping_file)

  sep_re = re.compile('\W')
  ensp_re = re.compile('(ENSP\d+)')
  word = ""
  with open(args.infile, "r") as ifh:
    with open(args.outfile, "w") as ofh:
      for line in ifh:
        for char in line:
          match_data = sep_re.match(char)
          if(match_data is not None):
            # then process word
            hgnc = ensp_to_hgnc.get(word)
            if hgnc is None:
              ofh.write(word + match_data.group(0))
            else:
              if(args.for_graphviz):
                if "-" in hgnc:
                  hgnc = "\"{}\"".format(hgnc)
              ofh.write(hgnc + match_data.group(0))
            word = ""
          else:
            word += char

if __name__ == "__main__":
  main()
