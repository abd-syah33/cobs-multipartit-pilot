# COBS asli (kode resmi TropBRC-BioinfoLab/COBS) di graf obat-target-kanker
# langkah 1 = greedy modularity (CNM) dari networkx, sama seperti di kode resmi
import os, sys, json, time
import networkx as nx
import networkx.algorithms.community as nxc

COBS_DIR = os.environ.get("COBS_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "COBS"))
sys.path.insert(0, COBS_DIR)
from core.algos import cobs_overlap
from core.metrics import compute_evaluations

from graf import baca_graf, KANKER

# konfigurasi COBS yg dipakai di studi kasus kanker payudara
PARAM = {"sim_metric": "aa", "theta_sim": 0.30, "theta_overlap": 0.10}
HASIL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hasil")


def jalankan_cobs(G, partisi, param=PARAM):
    # kode resmi memanggil greedy_modularity_communities di langkah 1,
    # di sini diganti sementara supaya langkah 1 pakai partisi yg sudah dihitung.
    # langkah 2-5 tetap kode resmi tanpa diubah
    asli = nxc.greedy_modularity_communities
    nxc.greedy_modularity_communities = lambda *a, **k: [set(c) for c in partisi]
    try:
        hasil = cobs_overlap.SPEC.run(G, param)
    finally:
        nxc.greedy_modularity_communities = asli
    return hasil


def simpan(nama, G, partisi, hasil, detik_partisi, detik_total):
    os.makedirs(HASIL, exist_ok=True)
    ev = compute_evaluations(hasil)
    info = dict(hasil.method_parameters)
    info["membership_count_dist"] = {str(k): v for k, v in info["membership_count_dist"].items()}
    info.update(ev)
    info["detik_langkah1"] = round(detik_partisi, 1)
    info["detik_total"] = round(detik_total, 1)
    with open(f"{HASIL}/ringkasan_{nama}.json", "w") as f:
        json.dump(info, f, indent=1)

    # komunitas crisp dan akhir (format panjang: simpul, komunitas)
    with open(f"{HASIL}/partisi_awal_{nama}.csv", "w") as f:
        f.write("simpul,komunitas\n")
        for j, c in enumerate(partisi):
            for u in sorted(c):
                f.write(f"{u},{j}\n")
    with open(f"{HASIL}/komunitas_akhir_{nama}.csv", "w") as f:
        f.write("simpul,komunitas\n")
        for j, c in enumerate(hasil.communities):
            for u in c:
                f.write(f"{u},{j}\n")
    return info


if __name__ == "__main__":
    G = baca_graf()
    print(G.number_of_nodes(), "simpul", G.number_of_edges(), "edge")

    t0 = time.time()
    P = [set(c) for c in nx.algorithms.community.greedy_modularity_communities(G)]
    t1 = time.time()
    print("greedy selesai", round(t1 - t0, 1), "detik,", len(P), "komunitas")

    hasil = jalankan_cobs(G, P)
    t2 = time.time()
    info = simpan("cobs_asli", G, P, hasil, t1 - t0, t2 - t0)
    print("komunitas akhir", info["final_k"], " simpul overlap", info["overlap_nodes"])
    print("shen modularity", round(info["shen_modularity"], 4))
