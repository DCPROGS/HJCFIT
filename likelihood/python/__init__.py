""" DCProgs python library. """
import likelihood
import random

def network(qmatrix): 
  """ Creates networkx graph object from a QMatrix object.
  
      Vertices have an "open" attribute to indicate whether they are open or shut. Edges have a
      "k+" and "k-" attribute containing the transition rates for the node with smaller index to
      the node with larger index, and vice-versa. 
  """
  from networkx import Graph

  graph = Graph()
  for i in xrange(qmatrix.nopen): graph.add_node(i, open=True)
  for j in xrange(qmatrix.nshut): graph.add_node(i+j, open=False)

  for i in xrange(qmatrix.matrix.shape[0]):
    for j in xrange(i, qmatrix.matrix.shape[1]):
      if abs(qmatrix.matrix[i,j]) > 1e-8:
        graph.add_edge(i, j)
        graph[i][j]['k+'] = qmatrix.matrix[i, j]
        graph[i][j]['k-'] = qmatrix.matrix[j, i]
  return graph
