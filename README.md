# Simulasi Pendulum dengan Metode Integrasi Numerik

**Nama:** Jesaya David Gamalael N P  
**NPM:** 2306161965

## Deskripsi Program

Program ini mengimplementasikan simulasi pendulum menggunakan tiga metode integrasi numerik yang berbeda untuk membandingkan akurasi dan perilaku masing-masing metode. Simulasi dilakukan pada dua skenario: dengan gaya eksternal dan tanpa gaya eksternal.

### Metode yang Diimplementasikan

1. **Kuadratur Gauss (2-Titik)**
   - Menggunakan titik dan bobot optimal Gauss-Legendre
   - Memberikan akurasi tinggi untuk fungsi polynomial
   - Cocok untuk integral dengan bentuk yang smooth

2. **Metode Romberg**
   - Menggunakan ekstrapolasi Richardson pada aturan trapesium
   - Mencapai akurasi O(h⁴) dengan dua level iterasi
   - Efisien untuk fungsi yang dapat didekati dengan polynomial

3. **Kuadratur Adaptif**
   - Menyesuaikan subdivisi interval berdasarkan toleransi error
   - Menggunakan aturan Simpson dengan kontrol error otomatis
   - Optimal untuk fungsi dengan variasi kompleksitas yang tinggi

### Fitur Program

- **Simulasi Pendulum**: Model pendulum teredam dengan kemungkinan gaya eksternal
- **Perbandingan Metode**: Analisis komparatif ketiga metode integrasi
- **Visualisasi Data**: 
  - Plot sudut vs waktu untuk kedua skenario
  - Diagram fase (sudut vs kecepatan sudut)
  - Perbandingan visual antar metode
- **Export Data**: Hasil simulasi disimpan dalam format CSV
- **Parameter Fisik**: Dapat disesuaikan (panjang tali, massa, redaman, dll.)

### Struktur File

```
├── numInt.cpp          
├── plot.py             
├── csv/    
│   ├── pendulum_gauss_with_force.csv
│   ├── pendulum_gauss_no_force.csv
│   ├── pendulum_romberg_with_force.csv
│   ├── pendulum_romberg_no_force.csv
│   ├── pendulum_adaptive_with_force.csv
│   └── pendulum_adaptive_no_force.csv
└── README.md           
```

### Cara Menjalankan Program

1. **Kompilasi dan jalankan simulasi C++:**
   ```bash
   g++ -o numInt numInt.cpp
   ./numInt
   ```

2. **Visualisasi hasil dengan Python:**
   ```bash
   python plot.py
   ```

### Parameter Simulasi

- **Panjang tali (L)**: 1.0 m
- **Massa beban (m)**: 0.2 kg
- **Koefisien redaman (c)**: 0.1 N.s/m
- **Waktu simulasi**: 0 - 20 detik
- **Step size (h)**: 0.01 detik
- **Gaya eksternal**: A cos(ωt) dengan A = 0.5 rad/s², ω = 2.0 rad/s

### Output

Program menghasilkan:
- 6 file CSV dengan data simulasi untuk setiap kombinasi metode dan skenario
- Plot perbandingan visual yang disimpan sebagai PNG
- Analisis komparatif perilaku pendulum pada berbagai kondisi

Program ini berguna untuk memahami karakteristik berbagai metode integrasi numerik dalam konteks sistem dinamis dan dapat digunakan sebagai dasar untuk studi lebih lanjut tentang metode numerik dalam fisika pendulum.