# COBS dg partisi awal composite modularity (modifikasi langkah 1 saja)
#
# Q = sum_r w_r * (1/m_r) * sum_c [ e_c^r - gamma * A_c^r * B_c^r / m_r ]
#   r      relasi (obat-target, target-kanker)
#   m_r    jumlah edge relasi r
#   e_c^r  edge relasi r di dalam komunitas c
#   A_c^r, B_c^r  jumlah derajat (di relasi r) simpul sisi A / sisi B yg ada di komunitas c
#   w_r = 1/jumlah relasi, gamma = 1
# tiap relasi pakai modularitas bipartit Barber, lalu dijumlah (composite, Liu dkk. 2014)
# optimasi pakai Louvain biasa, gain-nya diganti rumus di atas
import collections, random, time
import numpy as np
import networkx as nx

from graf import baca_graf
from cobs_asli import jalankan_cobs, simpan


def relasi(G, tipe):
    return sorted({tuple(sorted((tipe[u], tipe[v]))) for u, v in G.edges()})


def composite_modularity(G, komunitas, gamma=1.0):
    tipe = nx.get_node_attributes(G, "tipe")
    R = relasi(G, tipe)
    cid = {u: j for j, c in enumerate(komunitas) for u in c}
    m = collections.Counter(); e = collections.Counter()
    A = collections.Counter(); B = collections.Counter()
    for u, v in G.edges():
        r = tuple(sorted((tipe[u], tipe[v])))
        m[r] += 1
        if cid[u] == cid[v]:
            e[r, cid[u]] += 1
        a, b = (u, v) if tipe[u] == r[0] else (v, u)
        A[r, cid[a]] += 1
        B[r, cid[b]] += 1
    Q = 0
    for r in R:
        q = sum(e[r, c] - gamma * A[r, c] * B[r, c] / m[r] for c in range(len(komunitas))) / m[r]
        Q += q / len(R)
    return Q


def composite_louvain(G, gamma=1.0, seed=0):
    tipe = nx.get_node_attributes(G, "tipe")
    R = relasi(G, tipe); ri = {r: i for i, r in enumerate(R)}; nR = len(R)
    w = np.ones(nR) / nR
    m = np.zeros(nR)
    nodes = list(G); idx = {u: i for i, u in enumerate(nodes)}
    # derajat per relasi, sisi A dan sisi B
    As = np.zeros((len(nodes), nR)); Bs = np.zeros((len(nodes), nR))
    adj = [collections.defaultdict(lambda: np.zeros(nR)) for _ in nodes]
    for u, v in G.edges():
        r = tuple(sorted((tipe[u], tipe[v]))); k = ri[r]
        m[k] += 1
        a, b = (u, v) if tipe[u] == r[0] else (v, u)
        As[idx[a], k] += 1; Bs[idx[b], k] += 1
        adj[idx[u]][idx[v]][k] += 1; adj[idx[v]][idx[u]][k] += 1
    coef = w / m
    anggota = [[u] for u in nodes]
    rng = random.Random(seed)

    for level in range(10):
        n = len(anggota)
        kom = list(range(n))
        Ac = As.copy(); Bc = Bs.copy()
        ada_pindah = False
        pindah = True
        while pindah:
            pindah = False
            urutan = list(range(n)); rng.shuffle(urutan)
            for s in urutan:
                cs = kom[s]
                link = collections.defaultdict(lambda: np.zeros(nR))
                for t, wt in adj[s].items():
                    if t != s:
                        link[kom[t]] += wt
                Ac[cs] -= As[s]; Bc[cs] -= Bs[s]   # keluarkan s dulu

                def gain(c):
                    l = link.get(c, np.zeros(nR))
                    return float(np.dot(coef, l - gamma * (As[s] * Bc[c] + Bs[s] * Ac[c]) / m))

                best, bg = cs, gain(cs)
                for c in link:
                    g = gain(c)
                    if g > bg + 1e-12:
                        best, bg = c, g
                Ac[best] += As[s]; Bc[best] += Bs[s]
                if best != cs:
                    kom[s] = best; pindah = True; ada_pindah = True
        if not ada_pindah:
            break
        # agregasi jadi supernode
        lab = {c: i for i, c in enumerate(sorted(set(kom)))}
        nn = len(lab)
        baru = [[] for _ in range(nn)]
        nA = np.zeros((nn, nR)); nB = np.zeros((nn, nR))
        nadj = [collections.defaultdict(lambda: np.zeros(nR)) for _ in range(nn)]
        for s in range(n):
            c = lab[kom[s]]
            baru[c].extend(anggota[s]); nA[c] += As[s]; nB[c] += Bs[s]
            for t, wt in adj[s].items():
                nadj[c][lab[kom[t]]] += wt
        anggota, As, Bs, adj = baru, nA, nB, nadj
        if nn == n:
            break
    return [set(c) for c in anggota]


if __name__ == "__main__":
    G = baca_graf()
    print(G.number_of_nodes(), "simpul", G.number_of_edges(), "edge")

    t0 = time.time()
    # louvain hasilnya tergantung urutan, ambil yg Q-nya paling tinggi dari 5 seed
    P, qbest = None, -1
    for s in range(5):
        p = composite_louvain(G, seed=s)
        q = composite_modularity(G, p)
        print(" seed", s, len(p), "komunitas, Q =", round(q, 4))
        if q > qbest:
            P, qbest = p, q
    t1 = time.time()
    print("composite selesai", round(t1 - t0, 1), "detik,", len(P), "komunitas")

    hasil = jalankan_cobs(G, P)
    t2 = time.time()
    info = simpan("cobs_composite", G, P, hasil, t1 - t0, t2 - t0)
    print("komunitas akhir", info["final_k"], " simpul overlap", info["overlap_nodes"])
    print("shen modularity", round(info["shen_modularity"], 4))
