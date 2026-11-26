# COBS lengkap (langkah 1-5) asli vs composite di graf benchmark
# metrik: ONMI & omega thd ground truth, F1 simpul overlap
import os, sys, collections, time
import pandas as pd
import networkx as nx
from cdlib import NodeClustering, evaluation
from sklearn.metrics import normalized_mutual_info_score as nmi

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from cobs_asli import jalankan_cobs, PARAM
from cobs_composite import composite_louvain, composite_modularity
from core.metrics import _objective_qov_penalized as shen
from data_benchmark import southern_women, tripartit_sintetis, baca_real, REAL

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


def tipe_kecil(G, kom):
    # tipe simpul paling sedikit (kanker / konferensi / sutradara dst): berapa komunitas memuatnya, maks per komunitas
    tipe = nx.get_node_attributes(G, "tipe")
    jml = collections.Counter(tipe.values())
    tk = min(jml, key=jml.get)
    isi = [sum(1 for u in c if tipe[u] == tk) for c in kom]
    isi = [x for x in isi if x]
    return len(isi), max(isi)


def nilai(G, kom, truth, pada):
    m_det = collections.Counter(u for c in kom for u in c)
    od = {u for u, k in m_det.items() if k > 1}          # simpul overlap, semua tipe
    r = dict(K=len(kom), overlap_det=len(od), onmi=float("nan"), omega=float("nan"), nmi=float("nan"),
             purity=float("nan"), overlap_gt=float("nan"), f1_overlap=float("nan"))
    r["kom_tipe_kecil"], r["maks_tipe_kecil"] = tipe_kecil(G, kom)
    if truth is None:   # tanpa ground truth
        return r
    simpul = pada if pada is not None else set(G)
    kom = [c & simpul for c in kom if c & simpul]      # southern women / dblp / imdb: potong ke simpul berlabel
    sisa = simpul - set().union(*kom)
    det = NodeClustering([list(c) for c in kom] + [[u] for u in sisa], G.subgraph(simpul), overlap=True)
    gt = NodeClustering([list(c) for c in truth], G.subgraph(simpul), overlap=True)
    r["onmi"] = round(float(evaluation.overlapping_normalized_mutual_information_MGH(det, gt).score), 3)
    # omega cdlib menghitung semua pasangan simpul, di DBLP/IMDB (4 ribu simpul) RAM 2 GB habis -> dilewati
    if len(simpul) <= 2000:
        r["omega"] = round(float(evaluation.omega(det, gt).score), 3)
    # NMI biasa & purity (keanggotaan pertama saja), buat pembanding ONMI
    lab = {u: k for k, c in enumerate(truth) for u in c}
    cid = {}
    for j, c in enumerate(kom):
        for u in c:
            cid.setdefault(u, j)
    urut = sorted(simpul, key=str)
    r["nmi"] = round(nmi([lab.get(u, -1) for u in urut], [cid.get(u, -1) for u in urut]), 3)
    r["purity"] = round(sum(collections.Counter(lab.get(u, -1) for u in c).most_common(1)[0][1] for c in kom) / len(simpul), 3)
    # F1 simpul overlap (di simpul berlabel)
    m_gt = collections.Counter(u for c in truth for u in c)
    og = {u for u, k in m_gt.items() if k > 1}
    odp = od & simpul
    r["overlap_gt"] = len(og)
    tp = len(odp & og)
    r["f1_overlap"] = round(2 * tp / (len(odp) + len(og)), 3) if (odp or og) else float("nan")
    return r


baris = []
graf = [("southern women", 0) + southern_women()]
for s in range(5):
    graf.append(("tripartit sintetis", s) + tripartit_sintetis(seed=s))
for nama in REAL:
    graf.append((nama, 0) + baca_real(nama))

for nama, seed, G, truth, pada in graf:
    for metode, fp in [("cobs_asli", greedy), ("cobs_composite", composite)]:
        t = time.time()
        P = fp(G)
        hasil = jalankan_cobs(G, P)
        akhir = [set(c) for c in hasil.communities]
        dt = time.time() - t
        r1 = nilai(G, P, truth, pada)       # sesudah langkah 1
        r2 = nilai(G, akhir, truth, pada)   # sesudah langkah 5
        baris.append(dict(graf=nama, seed=seed, metode=metode, n=G.number_of_nodes(), m=G.number_of_edges(),
                          K_awal=r1["K"], onmi_awal=r1["onmi"], nmi_awal=r1["nmi"], purity_awal=r1["purity"],
                          kom_tipe_kecil_awal=r1["kom_tipe_kecil"], maks_tipe_kecil_awal=r1["maks_tipe_kecil"],
                          q_newman_awal=round(nx.community.modularity(G, P), 4),
                          q_composite_awal=round(composite_modularity(G, P), 4),
                          K_akhir=r2["K"], onmi_akhir=r2["onmi"], omega_akhir=r2["omega"],
                          shen_akhir=round(shen(G, akhir), 4),
                          overlap_det=r2["overlap_det"], overlap_gt=r2["overlap_gt"], f1_overlap=r2["f1_overlap"],
                          detik=round(dt, 1)))
        print(baris[-1], flush=True)

df = pd.DataFrame(baris)
df.to_csv(OUT, index=False)
print()
ring = df.groupby(["graf", "metode"], sort=False).mean(numeric_only=True).drop(columns="seed").round(3)
ring.to_csv(OUT.replace(".csv", "_ringkas.csv"))
print(ring.to_string())
