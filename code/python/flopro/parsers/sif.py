"""
Parse .sif into networkx graph
http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats
"""
import networkx as nx

def parse_sif(fh):
  G = nx.Graph()
  for line in fh:
    words = line.split()
    src = words[0]
    edge_type = words[1]
    targets = words[2:]
    for target in targets:
      G.add_edge(src, target, **{'type': edge_type})
  return G
