import os
import pandas as pd
import networkx as nx

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# tipe simpul: 0 obat, 1 target, 2 kanker
OBAT, TARGET, KANKER = 0, 1, 2


def baca_graf(folder=DATA):
    obat = pd.read_csv(f"{folder}/node_obat.csv")
    target = pd.read_csv(f"{folder}/node_target.csv")
    kanker = pd.read_csv(f"{folder}/node_kanker.csv")
    e_ot = pd.read_csv(f"{folder}/edge_obat_target.csv")
    e_tk = pd.read_csv(f"{folder}/edge_target_kanker.csv")

    G = nx.Graph()
    for i, nama in zip(obat.id, obat.nama):
        G.add_node(i, tipe=OBAT, nama=nama)
    for i, s in zip(target.id, target.simbol):
        G.add_node(i, tipe=TARGET, nama=s)
    for i, nama in zip(kanker.id, kanker.nama):
        G.add_node(i, tipe=KANKER, nama=nama)

    # tanpa bobot (jangan kasih atribut weight/score, COBS otomatis pakai itu sbg bobot)
    G.add_edges_from(zip(e_ot.obat, e_ot.target))
    G.add_edges_from(zip(e_tk.target, e_tk.kanker))

    # komponen terbesar saja
    cc = max(nx.connected_components(G), key=len)
    return G.subgraph(cc).copy()


if __name__ == "__main__":
    G = baca_graf()
    tipe = nx.get_node_attributes(G, "tipe")
    n = [sum(1 for u in G if tipe[u] == t) for t in (OBAT, TARGET, KANKER)]
    print("obat", n[0], "target", n[1], "kanker", n[2])
    print("edge", G.number_of_edges())
    ot = sum(1 for u, v in G.edges() if {tipe[u], tipe[v]} == {OBAT, TARGET})
    print("  obat-target", ot, " target-kanker", G.number_of_edges() - ot)
