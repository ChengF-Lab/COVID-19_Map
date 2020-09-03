# -*- coding: utf-8 -*-
# Author: Yadi Zhou

import sys
import random
import numpy as np
import numpy.ma as ma
import networkx as nx

"""
Examples:
python network_proximity.py example/Asthma.txt example/SARS2-DEG_lung.txt 1000 11096
python network_proximity.py DRUG example/SARS2-DEG.txt 1000 1024
"""


class DrugTarget(object):
    def __init__(self, dt="DrugTargetNetwork.txt"):
        self.d2t = {}
        with open(dt) as fi:
            for line in fi:
                d, t = line.strip("\n").split("\t")
                self.d2t.setdefault(d, []).append(t)
        self.drugs = sorted(self.d2t)

    def Screen(self, genes, repeat):
        interactome = Interactome(binSize=100)
        genes = [i for i in genes if i in interactome.G]
        for drug in self.drugs:
            targets = [t for t in self.d2t[drug] if t in interactome.G]
            if targets:
                d, z, p = interactome.ProximityZ(targets, genes, repeat=repeat)
                print(drug, " D: %.3f  Z: %s%.3f  P: %.3f" % (d, "" if z < 0 else " ", z, p))


class Interactome(object):
    def __init__(self, pathG="HumanInteractome.tsv", binSize=300):
        self.G = nx.read_edgelist(pathG, delimiter="\t", data=[("src", str), ("typ", str)])
        self.G.remove_edges_from(nx.selfloop_edges(self.G))
        self.nodes = sorted(self.G.nodes())
        self.n2d = {}
        self.d2n = {}
        for node, degree in self.G.degree():
            if degree in self.d2n:
                self.d2n[degree].append(node)
            else:
                self.d2n[degree] = [node]
            self.n2d[node] = degree
        self.d2b = {}
        self.b2i = {0: []}
        degrees = sorted(self.d2n.keys())
        self.dmin = degrees[0]
        self.dmax = degrees[-1]
        b = 0
        for curr, till in zip(degrees, degrees[1:] + [degrees[-1] + 1]):
            for d in range(curr, till):
                self.d2b[d] = b
            self.b2i[b].extend(self.d2n[curr])
            if curr != degrees[-1] and len(self.b2i[b]) >= binSize:
                b += 1
                self.b2i[b] = []
        if len(self.b2i[b]) < binSize and b > 0:
            for d in range(degrees[-1], -1, -1):
                if self.d2b[d] != b:
                    break
                self.d2b[d] = b - 1
            self.b2i[b - 1].extend(self.b2i[b])
            del self.b2i[b]
        self.binSize = binSize

    def DegreePreserveSampling(self, indexes):
        b2n = {}
        for index in indexes:
            b = self.d2b[max(self.dmin, min(self.dmax, self.n2d[index]))]
            if b not in b2n:
                b2n[b] = 0
            b2n[b] += 1
        for b in b2n:
            if b2n[b] > len(self.b2i[b]):
                raise ValueError("Number to sample > size of bin. Try to increase binSize")
        while True:
            yield sum([random.sample(self.b2i[b], b2n[b]) for b in sorted(b2n)], [])

    def Proximity(self, mod1, mod2):
        SD = np.zeros((len(mod1), len(mod2)), dtype=np.int)
        for row, source in enumerate(mod1):
            for col, target in enumerate(mod2):
                if source != target:
                    try:
                        SD[row, col] = nx.shortest_path_length(self.G, source, target)
                    except nx.NetworkXNoPath:
                        SD[row, col] = 9999
                else:
                    SD[row, col] = 0
        SD = ma.masked_values(SD, 9999)
        closest1 = SD.min(0)
        closest2 = SD.min(1)
        return (closest1.sum() + closest2.sum()) / (closest1.count() + closest2.count())

    def ProximityRandom(self, mod1, mod2, repeat):
        result = np.zeros(repeat)
        index = 0
        for mod1r, mod2r in zip(self.DegreePreserveSampling(mod1), self.DegreePreserveSampling(mod2)):
            v = self.Proximity(mod1r, mod2r)
            if not ma.is_masked(v):
                result[index] = v
                index += 1
                if index == repeat:
                    break
        return result

    def ProximityZ(self, mod1, mod2, repeat):
        d = self.Proximity(mod1, mod2)
        b = self.ProximityRandom(mod1, mod2, repeat=repeat)
        z, p = Z_Score(d, b)
        return d, z, p


def Z_Score(real, background):
    m = background.mean()
    s = background.std(ddof=1)
    z = (real - m) / s
    p = np.mean(background < real)
    return z, p


if __name__ == "__main__":
    input1 = sys.argv[1]
    input2 = sys.argv[2]
    repeat = int(sys.argv[3])
    random.seed(int(sys.argv[4]))

    if input1 == "DRUG":
        with open(input2) as fi:
            genes = fi.read().splitlines()
        dt = DrugTarget()
        dt.Screen(genes, repeat=repeat)
    else:
        with open(input1) as fi:
            genes1 = fi.read().splitlines()
        with open(input2) as fi:
            genes2 = fi.read().splitlines()
        interactome = Interactome()
        genes1 = [i for i in genes1 if i in interactome.G]
        genes2 = [i for i in genes2 if i in interactome.G]
        d, z, p = interactome.ProximityZ(genes1, genes2, repeat=repeat)
        print("D: %.3f  Z: %s%.3f  P: %.3f" % (d, "" if z < 0 else " ", z, p))
