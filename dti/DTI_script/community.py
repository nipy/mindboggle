#!/usr/bin/python

# Implementation of Newman's community detection algorithm
# see Newman, MEJ (2009) PNAS 103(23):8577-8582

import networkx as nx
import numpy as np
from numpy import linalg as LA

def bg(g,i,j,nbunch):
  # returns the modularity matrix for a subset of nodes
  # as defined in Eq. 6 in the paper
  def a(i,j):
    # returns the adjacency matrix element for network g
    if j in g.neighbors(i): return 1
    else: return 0
  def b(i,j):
    # returns an element of the modularity matrix, defined in Eq. 3 of the paper
    m = len(g.edges())
    return (a(i,j) - 1.0*g.degree(i)*g.degree(j)/(2*m))
  if i == j and len(nbunch) < len(g):
    return b(i,j) - sum([b(i,k) for k in nbunch])
  else:
    return b(i,j)

def klRefine(s,B):
  # Kernighan-Lin refinement, described on pg. 8580
  def flip(v,pos):
    newv = v.copy()
    newv[pos] = -1*v[pos]
    return newv 
  def dQ(s2):
    return np.dot(np.dot(s2,B),s2) 
  sBest = s
  dQBest = dQ(sBest)
  while True:
    trials = [dQ(flip(sBest,i)) for i in range(len(sBest))]
    if max(trials) > dQBest:
      i = trials.index(max(trials))
      sBest = flip(sBest,i)
      dQBest = dQ(sBest)
    else:
      break
  return sBest

def split(g,nbunch,errorMargin=100,doKL=True):
  B = np.array([[bg(g,i,j,nbunch) for j in nbunch] for i in nbunch])
  w,v = LA.eigh(B) # returns eigenvalues, eigenvectors
  nb1 = []
  nb2 = []
  i = list(w).index(max(w))
  # s as defined on pg. 8579
  s = np.array([(1 if x > 0 else -1) for x in v[:,i]])
  # dQ as in eq. 2 and 5
  dQ = np.dot(np.dot(s,B),s)/(4*len(g.edges()))
  if dQ <= errorMargin*np.finfo(np.double).eps:
    return False 
  if doKL:
    s = klRefine(s,B)
    dQ = np.dot(np.dot(s,B),s)/(4*len(g.edges()))
  global Q
  Q += dQ
  for j in range(len(s)):
    if s[j] < 0:
      nb1.append(nbunch[j])
    else:
      nb2.append(nbunch[j])
  return (nb1,nb2)

def recursive(g,nbunch):
  nb = split(g,nbunch)
  if not nb:
    global resultList
    resultList.append(nbunch)
    return
  else:
    recursive(g,nb[0])
    recursive(g,nb[1])

def detect_communities(g):
  global resultList
  resultList = []
  global Q
  Q = 0
  recursive(g,g.nodes())  
  return Q,resultList
      
################################################################################
## Test code
################################################################################

if __name__ == '__main__':
  def draw_nodes(g,fileName):
    colorList = [0]*len(g.nodes())
    for nbunch in resultList:
      for n in nbunch:
        colorList[g.nodes().index(n)] = resultList.index(nbunch)
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8,8))
    pos = nx.graphviz_layout(g,prog='neato')
    nx.draw(g,pos,node_color=colorList,with_labels=False)
    plt.savefig(fileName)

  g = nx.heawood_graph()
  myQ = detect_communities(g)
  print "Q = ",myQ[0]
  draw_nodes(g,"heawood.png")


  



