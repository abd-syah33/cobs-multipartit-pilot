# graf benchmark yg ground truth-nya diketahui
#  - southern women (Davis 1941): bipartit perempuan x acara, 2 kelompok perempuan
#  - tripartit sintetis: obat-target-kanker buatan, komunitas + overlap ditanam
import random
import networkx as nx

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


if __name__ == "__main__":
    G, T, _ = southern_women()
    print("southern women", G.number_of_nodes(), "simpul", G.number_of_edges(), "edge")
    G, T, _ = tripartit_sintetis()
    ov = sum(1 for u in G if sum(u in c for c in T) > 1)
    print("tripartit sintetis", G.number_of_nodes(), "simpul", G.number_of_edges(), "edge,", ov, "simpul overlap")
