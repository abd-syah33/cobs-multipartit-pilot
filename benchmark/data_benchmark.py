# graf benchmark yg ground truth-nya diketahui
#  - southern women (Davis 1941): bipartit perempuan x acara, 2 kelompok perempuan
#  - tripartit sintetis: obat-target-kanker buatan, komunitas + overlap ditanam
import os, glob, random
import pandas as pd
import networkx as nx

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# kelompok klasik Davis (dipakai di banyak paper, mis. Freeman 2003)
KEL1 = ["Evelyn Jefferson", "Laura Mandeville", "Theresa Anderson", "Brenda Rogers", "Charlotte McDowd",
        "Frances Anderson", "Eleanor Nye", "Pearl Oglethorpe", "Ruth DeSand"]
KEL2 = ["Verne Sanderson", "Myra Liddel", "Katherina Rogers", "Sylvia Avondale", "Nora Fayette",
        "Helen Lloyd", "Dorothy Murchison", "Olivia Carleton", "Flora Price"]


def southern_women():
    G = nx.davis_southern_women_graph()
    for u, d in G.nodes(data=True):
        d["tipe"] = d["bipartite"]   # 0 perempuan, 1 acara
    truth = [set(KEL1), set(KEL2)]
    return G, truth, set(KEL1) | set(KEL2)   # yg dievaluasi cuma perempuan


def tripartit_sintetis(n=(300, 300, 30), K=6, mu=0.2, ov=0.1, seed=0):
    # n = jumlah obat, target, kanker. K komunitas, tiap komunitas punya ketiga tipe
    # ov = proporsi simpul (tiap tipe) yg ikut 2 komunitas
    # mu = peluang edge keluar dari komunitas
    rng = random.Random(seed)
    anggota = {}
    for t, nt in enumerate(n):
        for i in range(nt):
            cs = {i % K}                       # dibagi rata ke K komunitas
            if rng.random() < ov:
                cs.add(rng.choice([c for c in range(K) if c != i % K]))
            anggota[f"{'otk'[t]}{i}"] = cs
    per_kom = {(t, c): [u for u in anggota if u[0] == 'otk'[t] and c in anggota[u]] for t in range(3) for c in range(K)}
    semua = {t: [u for u in anggota if u[0] == 'otk'[t]] for t in range(3)}

    G = nx.Graph()
    for u in anggota:
        G.add_node(u, tipe='otk'.index(u[0]))
    # obat -> target derajat 3-15, target -> kanker derajat 1-4
    for ta, tb, dmin, dmax in [(0, 1, 3, 15), (1, 2, 1, 4)]:
        for u in semua[ta]:
            for _ in range(rng.randint(dmin, dmax)):
                if rng.random() > mu:
                    c = rng.choice(sorted(anggota[u]))
                    v = rng.choice(per_kom[tb, c])
                else:
                    v = rng.choice(semua[tb])
                    while anggota[v] & anggota[u]:
                        v = rng.choice(semua[tb])
                G.add_edge(u, v)
    G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    truth = [set(u for u in G if c in anggota[u]) for c in range(K)]
    return G, truth, None


# dataset real world dari siapkan_data.py. urutan = urutan tipe simpul (0, 1, 2)
REAL = {
    "corporate_leadership": ["orang", "perusahaan"],
    "divorce": ["negara_bagian", "alasan"],
    "crime": ["orang", "kasus"],
    "protein_complex_drug": ["obat", "kompleks"],
    "dblp": ["penulis", "paper", "konferensi"],
    "imdb": ["sutradara", "film", "aktor"],
}


def baca_real(nama):
    d = os.path.join(DATA, nama)
    G = nx.Graph()
    label = {}
    for t, tipe in enumerate(REAL[nama]):
        df = pd.read_csv(f"{d}/node_{tipe}.csv")
        for i in df.id:
            G.add_node(i, tipe=t)
        if "label" in df:
            label.update(zip(df.id, df.label))
    for f in glob.glob(f"{d}/edge_*.csv"):
        df = pd.read_csv(f)
        G.add_edges_from(zip(df.iloc[:, 0], df.iloc[:, 1]))
    G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    if not label:
        return G, None, None       # tanpa ground truth
    pada = {u for u in label if u in G}
    truth = [set(u for u in pada if label[u] == k) for k in sorted(set(label.values()))]
    return G, truth, pada


if __name__ == "__main__":
    for nama in REAL:
        G, T, pada = baca_real(nama)
        print(nama, G.number_of_nodes(), "simpul", G.number_of_edges(), "edge", "" if T is None else f"label di {len(pada)} simpul, {len(T)} kelas")
    G, T, _ = southern_women()
    print("southern women", G.number_of_nodes(), "simpul", G.number_of_edges(), "edge")
    G, T, _ = tripartit_sintetis()
    ov = sum(1 for u in G if sum(u in c for c in T) > 1)
    print("tripartit sintetis", G.number_of_nodes(), "simpul", G.number_of_edges(), "edge,", ov, "simpul overlap")
