# bandingkan COBS asli vs COBS + composite modularity
# baca hasil/ dari cobs_asli.py dan cobs_composite.py, tulis hasil/perbandingan.csv
import os, json
import pandas as pd
import networkx as nx

from graf import baca_graf, KANKER
from cobs_composite import composite_modularity

HASIL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hasil")
METODE = ["cobs_asli", "cobs_composite"]


def baca_komunitas(f):
    d = pd.read_csv(f)
    return [set(g.simpul) for _, g in d.groupby("komunitas")]


def kanker_per_komunitas(G, kom):
    tipe = nx.get_node_attributes(G, "tipe")
    jml = [sum(1 for u in c if tipe.get(u) == KANKER) for c in kom]
    jml = [x for x in jml if x > 0]
    return len(jml), max(jml), sorted(jml, reverse=True)


G = baca_graf()
nama = nx.get_node_attributes(G, "nama")
tipe = nx.get_node_attributes(G, "tipe")
baris = []
for met in METODE:
    info = json.load(open(f"{HASIL}/ringkasan_{met}.json"))
    P = baca_komunitas(f"{HASIL}/partisi_awal_{met}.csv")
    C = baca_komunitas(f"{HASIL}/komunitas_akhir_{met}.csv")
    k1, mx1, s1 = kanker_per_komunitas(G, P)
    k2, mx2, s2 = kanker_per_komunitas(G, C)
    baris.append({
        "metode": met,
        # langkah 1
        "K_awal": len(P),
        "Q_newman_awal": round(nx.community.modularity(G, P), 4),
        "Q_composite_awal": round(composite_modularity(G, P), 4),
        "kom_berkanker_awal": k1,
        "maks_kanker_per_kom_awal": mx1,
        "sebaran_kanker_awal": str(s1),
        # sesudah langkah 2-5
        "K_akhir": info["final_k"],
        "simpul_batas": info["boundary_nodes"],
        "simpul_overlap": info["overlap_nodes"],
        "overlap_node_ratio": round(info["overlap_node_ratio"], 4),
        "kom_berkanker_akhir": k2,
        "maks_kanker_per_kom_akhir": mx2,
        "ukuran_kom_terbesar": info["max_community_size"],
        "shen_modularity": round(info["shen_modularity"], 4),
        "internal_edge_density": round(info["internal_edge_density"], 4),
        "conductance": round(info["conductance"], 4),
        "normalized_cut": round(info["normalized_cut"], 4),
        "avg_transitivity": round(info["avg_transitivity"], 4),
        "detik_langkah1": info["detik_langkah1"],
        "detik_total": info["detik_total"],
    })

df = pd.DataFrame(baris).set_index("metode").T
df.to_csv(f"{HASIL}/perbandingan.csv")
print(df.to_string())

# komunitas yg berisi kanker, partisi awal
print()
for met in METODE:
    P = baca_komunitas(f"{HASIL}/partisi_awal_{met}.csv")
    print(met)
    for j, c in enumerate(P):
        kk = sorted(nama[u] for u in c if tipe[u] == KANKER)
        if kk:
            no = sum(1 for u in c if tipe[u] == 0); nt = sum(1 for u in c if tipe[u] == 1)
            print(f"  kom {j}: {no} obat, {nt} target, kanker: {', '.join(kk)}")
