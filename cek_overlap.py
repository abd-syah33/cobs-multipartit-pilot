# kenapa langkah overlap COBS (langkah 4) hampir tidak menambah apa2?
# hitung ulang similarity tiap simpul batas ke komunitas kandidatnya, pakai fungsi dari kode resmi
import os, sys, collections
import numpy as np
import pandas as pd

from graf import baca_graf
from cobs_asli import cobs_overlap as co, PARAM, HASIL

G = baca_graf()
tipe = {u: d["tipe"] for u, d in G.nodes(data=True)}
nama = {u: d["nama"] for u, d in G.nodes(data=True)}
deg = dict(G.degree())
N = {u: set(G.neighbors(u)) for u in G}
metrik = PARAM["sim_metric"]

for met in ["cobs_asli", "cobs_composite"]:
    d = pd.read_csv(f"{HASIL}/partisi_awal_{met}.csv")
    P = [set(g.simpul) for _, g in d.groupby("komunitas")]
    cid = {u: j for j, c in enumerate(P) for u in c}

    # langkah 2: simpul batas & komunitas kandidat
    K = collections.defaultdict(set)
    for u, v in G.edges():
        if cid[u] != cid[v]:
            K[u].add(cid[v]); K[v].add(cid[u])
    # langkah 3: profil komunitas = gabungan tetangga anggota
    NC = [set().union(*(N[x] for x in C)) for C in P]
    NCn = [co._norm_w2(s, deg, metrik) for s in NC]

    # langkah 4: similarity simpul -> komunitas kandidat
    sim = []
    for u in K:
        nu = co._norm_w2(N[u], deg, metrik)
        for j in K[u]:
            inter, dot = co._intersection_count_and_dot_w2(N[u], NC[j], deg, metrik)
            s = co._sim_from_stats(metrik, len(N[u]), len(NC[j]), inter, dot_w2=dot, normA_w2=nu, normB_w2=NCn[j])
            sim.append((s, u, j))
    a = np.array([x[0] for x in sim])

    print(f"\n{met}")
    tb = collections.Counter(tipe[u] for u in K)
    print(f"  simpul batas {len(K)} (obat {tb[0]}, target {tb[1]}, kanker {tb[2]}), pasangan dicek {len(a)}")
    print(f"  median derajat simpul batas {np.median([len(N[u]) for u in K])}, median ukuran profil N(C) {np.median([len(s) for s in NC])}")
    print(f"  similarity: median {np.median(a):.4f}, persentil 99 {np.percentile(a, 99):.4f}, maks {a.max():.4f}")
    for th in (0.1, 0.2, PARAM["theta_sim"]):
        print(f"  >= {th}: {(a >= th).sum()} pasangan")
    for s, u, j in sorted(sim, reverse=True)[:3]:
        print(f"   {nama[u]} (derajat {deg[u]}) -> komunitas {j} (|N(C)| = {len(NC[j])}), sim {s:.3f}")
