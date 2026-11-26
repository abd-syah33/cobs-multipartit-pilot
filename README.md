# cobs-multipartit-pilot

Percobaan COBS (Crisp-to-Overlap Boundary Similarity) di graf tripartit obat-target-kanker, lalu langkah 1 (partisi awal) diganti dengan composite modularity.

## Isi

```
data/
  node_obat.csv            id ChEMBL, nama, tipe, tahap_tertinggi
  node_target.csv          id Ensembl, simbol, nama
  node_kanker.csv          id MONDO, nama
  edge_obat_target.csv     mekanisme aksi (MoA)
  edge_target_kanker.csv   skor asosiasi >= 0,3
  validasi_obat_kanker.csv pasangan obat-kanker dari data klinis (tidak dipakai skrip, utk evaluasi nanti)
graf.py             baca csv -> graf
cobs_asli.py        COBS resmi (greedy modularity -> overlap)
cobs_composite.py   COBS dg langkah 1 composite modularity
bandingkan.py       tabel perbandingan
cek_overlap.py      cek similarity di langkah 4 (kenapa overlap hampir tidak ada)
hasil/
benchmark/
  data_benchmark.py       southern women + tripartit sintetis (ground truth diketahui)
  jalankan_benchmark.py   COBS lengkap asli vs composite, ONMI & F1 overlap
  hasil_benchmark.csv
```

## Data

Draf dataset (belum final), 11 kanker. Obat-target dari mekanisme aksi ChEMBL (lewat Open Targets 26.09). Target-kanker dari skor asosiasi Open Targets tanpa sumber yg sirkular (clinical precedence, cancer biomarkers, literatur), ambang 0,3.

Graf yg dipakai = komponen terbesar: 2.305 obat, 5.793 target, 11 kanker, 21.681 edge (7.678 obat-target, 14.003 target-kanker). Tanpa bobot.

Edge validasi obat-kanker sengaja tidak masuk graf supaya tidak bocor ke evaluasi.

## Cara jalan

Kode COBS tidak disalin ke sini, clone dulu repo resminya:

```
git clone https://github.com/TropBRC-BioinfoLab/COBS
pip install -r requirements.txt
export COBS_DIR=/path/ke/COBS      # default ../COBS

python graf.py
python cobs_asli.py        # +- 1 menit (greedy CNM)
python cobs_composite.py   # +- 10 detik
python bandingkan.py
python cek_overlap.py

cd benchmark && python jalankan_benchmark.py   # +- 15 detik
```

Parameter COBS sama utk keduanya: Adamic-Adar, theta_sim 0,30, theta_overlap 0,10. Yg beda cuma langkah 1. Langkah 2-5 kode resmi, tidak diubah; langkah 1 diganti dg menimpa sementara `greedy_modularity_communities` waktu COBS dipanggil (lihat `jalankan_cobs` di cobs_asli.py).

## Composite modularity

```
Q = sum_r w_r (1/m_r) sum_c [ e_c^r - gamma A_c^r B_c^r / m_r ]
```

r = relasi (obat-target, target-kanker), tiap relasi pakai modularitas bipartit Barber (2007), dijumlah seperti composite modularity Liu dkk. (2014). w_r = 1/2, gamma = 1. Dioptimasi dg Louvain, 5 seed, ambil Q tertinggi.

## Hasil (hasil/perbandingan.csv)

| | COBS asli | COBS + composite |
|---|---|---|
| komunitas awal | 131 | 49 |
| Q Newman awal | 0,579 | 0,540 |
| Q composite awal | 0,561 | 0,610 |
| komunitas yg berisi kanker | 4 | 9 |
| maks kanker per komunitas | 8 | 3 |
| komunitas terbesar | 2.859 | 1.315 |
| simpul batas | 1.770 | 2.066 |
| simpul overlap | 1 | 1 |
| Shen modularity | 0,580 | 0,525 |
| internal edge density | 0,442 | 0,119 |
| conductance | 0,189 | 0,194 |
| normalized cut | 0,191 | 0,202 |
| waktu langkah 1 | 62 s | 9 s |

Greedy menaruh 8 kanker di satu komunitas (145 obat, 2.706 target). Composite memisahkan: 8 kanker di komunitas sendiri-sendiri, 3 (glioblastoma, hati, ovarium) masih satu komunitas.

Catatan:
- Shen modularity & Q Newman pakai null model yg tidak membedakan tipe simpul, jadi condong ke greedy. Q composite condong ke composite.
- Internal edge density greedy tinggi krn banyak komunitas kecil (38 komunitas ukuran <= 3, median 8; composite median 56).
- Langkah overlap hampir tidak menambah apa2: cuma 1 simpul overlap di masing2 (CUDC-101 di COBS asli, kanker ovarium di composite). Ini masalah sendiri, belum dibahas di sini.
- Hasil composite tergantung gamma. gamma = 1 di sini; gamma lebih kecil -> komunitas lebih besar & kanker menumpuk lagi, gamma lebih besar -> lebih terpisah tapi komunitas makin banyak/kecil. Belum ada aturan pemilihan gamma.

## Kenapa overlap hampir tidak ada (cek_overlap.py)

Langkah 2-5 jalan semua, tapi hampir tidak ada yg lolos ambang:

| | COBS asli | COBS + composite |
|---|---|---|
| simpul batas | 1.770 | 2.066 |
| pasangan (simpul, komunitas kandidat) | 3.834 | 8.078 |
| median similarity | 0,008 | 0,009 |
| persentil 99 | 0,055 | 0,091 |
| lolos >= 0,30 | 1 | 2 |
| merge | 0 | 0 |

Simpul batas derajatnya kecil (median 4-8) sedangkan profil komunitas N(C) bisa ratusan-ribuan simpul, jadi similarity jenis cosine kira2 paling tinggi sqrt(|N(u)|/|N(C)|). Di graf tripartit juga tidak ada segitiga, irisan N(u) dan N(C) kecil. Yg lolos justru simpul derajat besar (kanker ovarium, derajat 877).

## Benchmark (benchmark/hasil_benchmark.csv)

COBS lengkap, parameter sama (AA, 0,30, 0,10), gamma = 1.

- Southern women: 18 perempuan x 14 acara, ground truth 2 kelompok klasik Davis, ONMI dihitung di perempuan saja.
- Tripartit sintetis: 300 obat, 300 target, 30 kanker, 6 komunitas (tiap komunitas punya 3 tipe), 10% simpul tiap tipe ikut 2 komunitas, mu = 0,2 (peluang edge keluar komunitas). 5 seed, dirata2.

| | COBS asli | COBS + composite |
|---|---|---|
| SW: komunitas akhir | 3 | 3 |
| SW: ONMI sesudah langkah 1 | 0,666 | 0,267 |
| SW: ONMI akhir | 0,666 | 0,366 |
| Sintetis: komunitas akhir | 5,6 | 6,0 |
| Sintetis: ONMI akhir | 0,564 | 0,634 |
| Sintetis: omega akhir | 0,660 | 0,716 |
| Sintetis: simpul overlap terdeteksi / tertanam | 0 / 65 | 0 / 65 |
| Sintetis: F1 simpul overlap | 0 | 0 |

- Di SW composite kalah: partisinya memecah jadi 4 (gamma = 1 terlalu halus utk graf sekecil ini), langkah merge COBS menggabungkan jadi 3 tapi tetap di bawah greedy.
- Di sintetis composite lebih baik di 4 dari 5 seed (seed 0 greedy sedikit lebih baik), di 2 seed greedy menggabungkan 2 komunitas (K = 5).
- Langkah overlap tidak menemukan satu pun dari 65 simpul overlap tertanam di kedua versi. Jadi masalah overlap bukan cuma di draf dataset.
