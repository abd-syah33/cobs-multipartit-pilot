# unduh dataset benchmark real world -> csv node/edge di benchmark/data/<nama>/
# bipartit : corporate leadership, crime (KONECT), divorce (UCI netdata), protein complex-drug (Nacher & Schwartz 2012)
# tripartit: DBLP (penulis-paper-konferensi) dan IMDB (sutradara-film-aktor), versi olahan MAGNN (Fu dkk. 2020)
# malaria gene-substring tidak ada: situs aslinya sudah mati, repo datanya cuma berisi jaringan gen-gen
# butuh: pip install rdata xlrd
import os, io, tarfile, zipfile, urllib.request, warnings
import numpy as np, pandas as pd, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
MENTAH = os.path.join(HERE, "mentah")      # file unduhan, tidak masuk git
DATA = os.path.join(HERE, "data")
os.makedirs(MENTAH, exist_ok=True)


def unduh(url, nama):
    f = os.path.join(MENTAH, nama)
    if not os.path.exists(f):
        print("unduh", url)
        urllib.request.urlretrieve(url, f)
    return f


def tulis(nama, nodes, edges):
    # nodes: {tipe: DataFrame(id, ...)}, edges: {(tipeA, tipeB): DataFrame(tipeA, tipeB)}
    d = os.path.join(DATA, nama); os.makedirs(d, exist_ok=True)
    for t, df in nodes.items():
        df.to_csv(f"{d}/node_{t}.csv", index=False)
    for (a, b), df in edges.items():
        df.drop_duplicates().to_csv(f"{d}/edge_{a}_{b}.csv", index=False)
    print(nama, {t: len(v) for t, v in nodes.items()}, {f"{a}-{b}": len(v.drop_duplicates()) for (a, b), v in edges.items()})


def konect(nama):
    f = unduh(f"http://konect.cc/files/download.tsv.{nama}.tar.bz2", nama + ".tar.bz2")
    tar = tarfile.open(f)
    baca = lambda akhiran: [m for m in tar.getmembers() if m.name.split("/")[-1].startswith(akhiran)]
    out = tar.extractfile(baca("out.")[0]).read().decode().splitlines()
    e = [tuple(map(int, l.split()[:2])) for l in out if not l.startswith("%")]
    return tar, baca, e


# --- corporate leadership: 20 orang x 24 perusahaan
tar, baca, e = konect("brunson_corporate-leadership")
e = pd.DataFrame([(f"o{a}", f"p{b}") for a, b in e], columns=["orang", "perusahaan"])
tulis("corporate_leadership",
      {"orang": pd.DataFrame({"id": sorted(set(e.orang))}), "perusahaan": pd.DataFrame({"id": sorted(set(e.perusahaan))})},
      {("orang", "perusahaan"): e})

# --- crime: orang x kasus, peran (suspect/victim/witness) disimpan di edge
tar, baca, e = konect("moreno_crime")
nm = tar.extractfile(baca("ent.moreno_crime_crime.person.name")[0]).read().decode().splitlines()
sx = tar.extractfile(baca("ent.moreno_crime_crime.person.sex")[0]).read().decode().splitlines()
role = tar.extractfile(baca("rel.moreno_crime_crime.person.role")[0]).read().decode().splitlines()
ed = pd.DataFrame([(f"o{a}", f"k{b}") for a, b in e], columns=["orang", "kasus"]); ed["peran"] = role
orang = pd.DataFrame({"id": [f"o{i+1}" for i in range(len(nm))], "nama": nm, "sex": sx})
tulis("crime", {"orang": orang[orang.id.isin(ed.orang)], "kasus": pd.DataFrame({"id": sorted(set(ed.kasus))})},
      {("orang", "kasus"): ed})

# --- divorce: 50 negara bagian x 9 alasan cerai (file R)
import rdata
warnings.filterwarnings("ignore")
d = rdata.read_rda(unduh("http://networkdata.ics.uci.edu/netdata/data/divorce.RData", "divorce.RData"))["divorce"]
nama = [str(v["vertex.names"][0]) for v in d["val"]]
nb = int(d["gal"]["bipartite"][0])     # 50 simpul pertama = negara bagian
e = []
for x in d["mel"]:
    a, b = int(x["outl"][0]) - 1, int(x["inl"][0]) - 1
    if a >= nb: a, b = b, a
    e.append((nama[a], nama[b]))
e = pd.DataFrame(e, columns=["negara_bagian", "alasan"])
tulis("divorce", {"negara_bagian": pd.DataFrame({"id": nama[:nb]}), "alasan": pd.DataFrame({"id": nama[nb:]})},
      {("negara_bagian", "alasan"): e})

# --- protein complex - drug (DrugBank x CORUM), suplemen PLoS One 7(1):e30028
url = "https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0030028.s00"
e = pd.read_csv(unduh(url + "1", "pcd_s001.txt"), sep="\t").rename(columns={"Drug_ID": "obat", "Complex_ID": "kompleks"})
cx = pd.read_excel(unduh(url + "2", "pcd_s002.xls"), engine="xlrd")
dr = pd.read_excel(unduh(url + "3", "pcd_s003.xls"), engine="xlrd")
dr = dr.rename(columns={"DrugBank ID": "id", "Drug name": "nama"})
cx["id"] = ["CORUM_%04d" % i for i in cx["CORUM ID"]]
cx = cx.rename(columns={"Complex name": "nama"})[["id", "nama"]]
tulis("protein_complex_drug", {"obat": dr[dr.id.isin(e.obat)], "kompleks": cx[cx.id.isin(e.kompleks)]},
      {("obat", "kompleks"): e})

# --- DBLP & IMDB versi MAGNN: adjM + node_types + labels
def magnn(nama, url):
    z = unduh(url, nama + "_processed.zip")
    zf = zipfile.ZipFile(z)
    A = sp.load_npz(io.BytesIO(zf.read("adjM.npz"))).tocoo()
    t = np.load(io.BytesIO(zf.read("node_types.npy")))
    y = np.load(io.BytesIO(zf.read("labels.npy")))
    return A, t, y

A, t, y = magnn("DBLP", "https://www.dropbox.com/s/yh4grpeks87ugr2/DBLP_processed.zip?dl=1")
# tipe: 0 penulis, 1 paper, 2 term (dibuang), 3 konferensi. label penulis: 4 bidang
tp = {0: "penulis", 1: "paper", 3: "konferensi"}
awalan = {0: "a", 1: "p", 3: "k"}
idn = lambda i: f"{awalan[t[i]]}{i}"
assert (t[:len(y)] == 0).all()   # label = urutan penulis di awal
m = A.row < A.col
e = {}
for a, b in zip(A.row[m], A.col[m]):
    if t[a] in tp and t[b] in tp:
        if t[a] > t[b]: a, b = b, a
        e.setdefault((tp[t[a]], tp[t[b]]), []).append((idn(a), idn(b)))
e = {k: pd.DataFrame(v, columns=list(k)) for k, v in e.items()}
pen = [i for i in range(len(t)) if t[i] == 0]
nodes = {"penulis": pd.DataFrame({"id": [idn(i) for i in pen], "label": y[:len(pen)]}),
         "paper": pd.DataFrame({"id": [idn(i) for i in range(len(t)) if t[i] == 1]}),
         "konferensi": pd.DataFrame({"id": [idn(i) for i in range(len(t)) if t[i] == 3]})}
tulis("dblp", nodes, e)

A, t, y = magnn("IMDB", "https://www.dropbox.com/s/g0btk9ctr1es39x/IMDB_processed.zip?dl=1")
# tipe: 0 film, 1 sutradara, 2 aktor. label film: 3 genre (action, comedy, drama)
tp = {0: "film", 1: "sutradara", 2: "aktor"}
awalan = {0: "f", 1: "s", 2: "a"}
idn = lambda i: f"{awalan[t[i]]}{i}"
assert (t[:len(y)] == 0).all()
m = A.row < A.col
e = {}
for a, b in zip(A.row[m], A.col[m]):
    if t[a] > t[b]: a, b = b, a
    e.setdefault((tp[t[a]], tp[t[b]]), []).append((idn(a), idn(b)))
e = {k: pd.DataFrame(v, columns=list(k)) for k, v in e.items()}
film = [i for i in range(len(t)) if t[i] == 0]
nodes = {"film": pd.DataFrame({"id": [idn(i) for i in film], "label": y[:len(film)]}),
         "sutradara": pd.DataFrame({"id": [idn(i) for i in range(len(t)) if t[i] == 1]}),
         "aktor": pd.DataFrame({"id": [idn(i) for i in range(len(t)) if t[i] == 2]})}
tulis("imdb", nodes, e)
