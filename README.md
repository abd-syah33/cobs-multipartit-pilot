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
  siapkan_data.py         unduh dataset real world -> data/<nama>/ (csv node/edge)
  data/                   corporate_leadership, divorce, crime, protein_complex_drug, dblp, imdb
  data_benchmark.py       southern women, tripartit sintetis, baca dataset real world
  jalankan_benchmark.py   COBS lengkap asli vs composite di semua graf benchmark
  hasil_benchmark.csv, hasil_benchmark_ringkas.csv
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

cd benchmark
python siapkan_data.py        # sekali saja, csv sudah ada di repo. DBLP mentah 258 MB
python jalankan_benchmark.py  # +- 3 menit (DBLP paling lama)
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

## Benchmark (benchmark/hasil_benchmark_ringkas.csv)

COBS lengkap (langkah 1-5), parameter sama (AA, 0,30, 0,10), gamma = 1. Semua graf diambil komponen terbesarnya.

Graf:

| graf | tipe simpul | simpul | edge | ground truth | sumber |
|---|---|---|---|---|---|
| southern women | perempuan, acara | 32 | 89 | 2 kelompok perempuan (Davis) | networkx |
| tripartit sintetis | obat, target, kanker | 630 | 3.145-3.341 | 6 komunitas tertanam, 10% overlap, mu 0,2, 5 seed | data_benchmark.py |
| corporate leadership | orang, perusahaan | 44 | 99 | - | KONECT |
| divorce | negara bagian, alasan cerai | 59 | 225 | - | UCI netdata |
| crime | orang, kasus | 1.263 | 1.377 | - | KONECT |
| protein complex-drug | obat (DrugBank), kompleks (CORUM) | 1.302 | 3.576 | - | Nacher & Schwartz 2012, PLoS One |
| DBLP | penulis, paper, konferensi | 18.405 | 33.973 | 4 bidang riset penulis | MAGNN (Fu dkk. 2020), term dibuang |
| IMDB | sutradara, film, aktor | 10.384 | 16.097 | 3 genre film | MAGNN (Fu dkk. 2020) |

Lisensi data ikut sumber masing2 (KONECT CC BY-SA, divorce CC BY-NC-ND, PLoS One CC BY, MAGNN tanpa keterangan), repo ini private.

Malaria gene-substring tidak bisa diambil: situs aslinya sudah mati dan repo data penulisnya cuma berisi jaringan gen-gen, bukan bipartit gen-substring.

Dengan ground truth (ONMI di simpul berlabel; NMI & purity dari partisi awal):

| graf | | COBS asli | COBS + composite |
|---|---|---|---|
| southern women | komunitas | 3 | 3 |
| | ONMI | **0,666** | 0,366 |
| | NMI / purity | **0,840 / 1,000** | 0,451 / 0,889 |
| tripartit sintetis (rata2 5 seed) | komunitas | 5,6 | 6,0 |
| | ONMI | 0,564 | **0,634** |
| | NMI / purity | 0,665 / 0,813 | **0,698 / 0,869** |
| DBLP | komunitas | 23 | 14 |
| | ONMI | **0,348** | 0,288 |
| | NMI / purity | 0,568 / 0,895 | 0,562 / 0,911 |
| | komunitas berisi konferensi / maks konferensi per komunitas | 12 / 4 | 14 / 3 |
| IMDB | komunitas | 71 | 85 |
| | ONMI | 0,000 | 0,000 |
| | NMI / purity | 0,027 / 0,469 | 0,042 / 0,519 |

Tanpa ground truth (cuma modularitas; Q composite = Barber utk graf bipartit):

| graf | komunitas (asli / composite) | Q Newman | Q composite | simpul overlap |
|---|---|---|---|---|
| corporate leadership | 6 / 4 | 0,388 / 0,424 | 0,399 / 0,425 | 8 / 3 |
| divorce | 4 / 3 | 0,170 / 0,178 | 0,177 / 0,185 | 10 / 8 |
| crime | 35 / 34 | 0,924 / 0,924 | 0,924 / 0,924 | 0 / 0 |
| protein complex-drug | 25 / 22 | 0,813 / 0,846 | 0,813 / 0,846 | 2 / 0 |

Catatan:
- Southern women: composite kalah, gamma = 1 memecah jadi 4, merge COBS menggabung jadi 3 tapi tetap di bawah greedy.
- Tripartit sintetis: composite lebih baik di 4 dari 5 seed.
- DBLP: hampir setara. NMI sama, ONMI greedy lebih tinggi, purity composite sedikit lebih tinggi. Mirip draf dataset: composite menyebar simpul tipe kecil (konferensi) ke lebih banyak komunitas (maks 3 vs 4 per komunitas).
- IMDB: dua2nya hampir tidak cocok dg genre (purity 0,47-0,52, padahal tebak genre terbanyak saja sudah 0,37). Genre bukan komunitas yg terbentuk dari sutradara/aktor, jadi IMDB kurang cocok jadi ground truth. ONMI 0 krn komunitasnya banyak (71-85) dan kecil.
- Tanpa ground truth: Q composite sama atau sedikit lebih tinggi (termasuk Q Newman), tapi ini cuma kualitas struktur, bukan benar/salah.
- Overlap: di semua graf langkah overlap hampir tidak menambah apa2 (0-10 simpul), di sintetis 0 dari 65 simpul overlap tertanam. Jadi masalah overlap bukan cuma di draf dataset.
- Omega tidak dihitung di DBLP & IMDB, fungsi cdlib-nya kehabisan RAM (2 GB).
