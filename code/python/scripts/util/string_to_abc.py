#!/usr/bin/env python
from ppi import string_db
import sys, argparse

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--infile", "-i", type=argparse.FileType('r'), default=sys.stdin)
  parser.add_argument("--outfile", "-o", type=argparse.FileType('w'), default=sys.stdout)
  args = parser.parse_args()
  string_db.string_to_abc(args.infile, args.outfile)

if __name__ == "__main__":
  main()
