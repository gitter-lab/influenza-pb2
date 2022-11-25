"""
The "abc" file format consists of lines each of which defines an edge and its confidence e.g.
A B 0
B C 1
"""
from ppi import SimPathException
import networkx as nx

class ParseAbcException(SimPathException):
  pass

def parse_abc(fh):
  if(type(fh) is str):
    fh = open(fh, 'r')
  G = nx.Graph()
  line_no = 0
  for line in fh:
    line_no += 1
    line = line.rstrip()
    words = line.split()
    if(len(words) != 3):
      raise ParseAbcException("Invalid line {}".format(line_no))
    G.add_edge(words[0], words[1], **{'weight': float(words[2])})
  return G
