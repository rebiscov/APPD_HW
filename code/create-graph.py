#!/usr/bin/env python2.7
import sys
import random

if (len(sys.argv) != 5):
    print >> sys.stderr, "Usage: ./create-graph [num-vertices] [num-edges] [max-weight] [file-name]"
    exit(1)

nVtx = int(sys.argv[1])
nEdge = int(sys.argv[2])
maxWeight = int(sys.argv[3])
fileName = sys.argv[4]

edgePerVtx = nEdge / nVtx
remEdge = nEdge % nVtx

component = dict()
size = dict()
for i in range(0,nVtx):
    component[i] = i
    size[i] = 1

def racine(i):
    j = i
    while j != component[j]:
        j = component[j]
    return j

connex = False

while not(connex):
    edges = dict()
    component = dict()
    for i in range(0, nVtx):
        component[i] = i

    for i in range(0, nVtx):
        j = 0
        while j < edgePerVtx:
            edge = random.randrange(0, nVtx)
            if ((i, edge) in edges or (edge, i) in edges):
                continue
            else:
                weight = random.randrange(1,maxWeight+1)
                edges[(i, edge)] = weight
                r1 = racine(i)
                r2 = racine(edge)
                if r1 != r2:
                    if (size[r1] < size[r2]):
                        component[r1] = r2
                        size[r2] += size[r1]
                    else:
                        component[r2] = r1
                        size[r1] += size[r2]
                j = j + 1

    for i in range(0, remEdge):
        l = random.randrange(0, nVtx);
        r = random.randrange(0, nVtx);
        while l == r or (l,r) in edges or (r,l) in edges:
            l = random.randrange(0, nVtx)
            r = random.randrange(0, nVtx)
        weight = random.randrange(1,maxWeight+1)
        edges[(l, r)] = weight
        r1 = racine(l)
        r2 = racine(r)
        if r1 != r2:
            if (size[r1] < size[r2]):
                component[r1] = r2
                size[r2] += size[r1]
            else:
                component[r2] = r1
                size[r1] += size[r2]

    connex = True
    r = racine(0)
    for i in range(1,nVtx):
        if racine(i) != r:
            connex = False
    if (not(connex)):
        print("non connexe")

f = open(fileName, 'w')
f.write('%d %d\n' % (nVtx, nEdge))
for i in edges:
    f.write('%d %d %d\n' % (i[0], i[1], edges[i]))

f.close()

print 'A graph with %d vertices and %d edges is written into the file %s.' % (nVtx, len(edges), fileName)
