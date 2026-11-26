# COBS lengkap (langkah 1-5) asli vs composite di graf benchmark
# metrik: ONMI & omega thd ground truth, F1 simpul overlap
import os, sys, collections, time
import pandas as pd
import networkx as nx
from cdlib import NodeClustering, evaluation

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from cobs_asli import jalankan_cobs, PARAM
from cobs_composite import composite_louvain, composite_modularity
from data_benchmark import southern_women, tripartit_sintetis

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hasil_benchmark.csv")


def greedy(G):
    return [set(c) for c in nx.algorithms.community.greedy_modularity_communities(G)]


def composite(G):
    best, qb = None, -1
    for s in range(5):
        p = composite_louvain(G, seed=s)
        q = composite_modularity(G, p)
        if q > qb:
            best, qb = p, q
    return best


def nilai(G, kom, truth, pada):
    if pada is not None:   # southern women: potong ke perempuan saja
        kom = [c & pada for c in kom if c & pada]
        simpul = pada
    else:
        simpul = set(G)
    sisa = simpul - set().union(*kom)
    det = NodeClustering([list(c) for c in kom] + [[u] for u in sisa], G.subgraph(simpul), overlap=True)
    gt = NodeClustering([list(c) for c in truth], G.subgraph(simpul), overlap=True)
    onmi = evaluation.overlapping_normalized_mutual_information_MGH(det, gt).score
    omega = evaluation.omega(det, gt).score
    # F1 simpul overlap
    m_det = collections.Counter(u for c in kom for u in c)
    m_gt = collections.Counter(u for c in truth for u in c)
    od = {u for u, k in m_det.items() if k > 1}
    og = {u for u, k in m_gt.items() if k > 1}
    tp = len(od & og)
    f1 = 2 * tp / (len(od) + len(og)) if (od or og) else float("nan")
    return dict(K=len(kom), onmi=round(float(onmi), 3), omega=round(float(omega), 3),
                overlap_det=len(od), overlap_gt=len(og), f1_overlap=round(f1, 3))


baris = []
graf = [("southern women", 0) + southern_women()]
for s in range(5):
    graf.append(("tripartit sintetis", s) + tripartit_sintetis(seed=s))

for nama, seed, G, truth, pada in graf:
    for metode, fp in [("cobs_asli", greedy), ("cobs_composite", composite)]:
        t = time.time()
        P = fp(G)
        hasil = jalankan_cobs(G, P)
        akhir = [set(c) for c in hasil.communities]
        dt = time.time() - t
        r1 = nilai(G, P, truth, pada)       # sesudah langkah 1
        r2 = nilai(G, akhir, truth, pada)   # sesudah langkah 5
        baris.append(dict(graf=nama, seed=seed, metode=metode,
                          K_awal=r1["K"], onmi_awal=r1["onmi"],
                          K_akhir=r2["K"], onmi_akhir=r2["onmi"], omega_akhir=r2["omega"],
                          overlap_det=r2["overlap_det"], overlap_gt=r2["overlap_gt"], f1_overlap=r2["f1_overlap"],
                          detik=round(dt, 1)))
        print(baris[-1], flush=True)

df = pd.DataFrame(baris)
df.to_csv(OUT, index=False)
print()
print(df.groupby(["graf", "metode"]).mean(numeric_only=True).drop(columns="seed").round(3).to_string())
