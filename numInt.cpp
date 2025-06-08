// Jesaya David Gamalael N P
// 2306161965

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <functional>

using namespace std;

// =================================================================
// PARAMETER GLOBAL
// =================================================================

// Parameter fisik dari pendulum
const double L = 1.0;  // Panjang tali (m)
const double m = 0.2;  // Massa beban (kg)
const double c = 0.1;  // Koefisien redaman (N.s/m)
const double g = 9.81; // Percepatan gravitasi (m/s²)

// Parameter gaya eksternal (dapat diubah oleh fungsi simulasi)
double A = 0.5;       // Amplitudo gaya eksternal (rad/s²)
double omega_d = 2.0; // Frekuensi sudut gaya eksternal (rad/s)

// Kondisi awal
const double theta0 = 0.2; // Sudut awal (rad)
const double omega0 = 0.0; // Kecepatan sudut awal (rad/s)

// Parameter simulasi
const double t0 = 0.0;  // Waktu awal (s)
const double tf = 20.0; // Waktu akhir (s)
const double h = 0.01;  // Ukuran langkah waktu (s)

// =================================================================
// FUNGSI TURUNAN
// =================================================================

// Mendefinisikan turunan theta (dtheta/dt = omega)
double f_theta(double t, double theta, double omega)
{
    return omega;
}

// Mendefinisikan turunan omega (domega/dt = ...)
double f_omega(double t, double theta, double omega)
{
    return A * cos(omega_d * t) - (c / (m * L)) * omega - (g / L) * sin(theta);
}

// =================================================================
// METODE NUMERIK
// =================================================================

// --- 1. Metode Kuadratur Gauss (2-Titik) ---
// Metode ini menggunakan titik dan bobot optimal untuk akurasi tinggi.
void gauss_quadrature_method(double t, double &theta, double &omega)
{
    // Titik dan bobot Gauss-Legendre 2-titik
    const double w1 = 1.0, w2 = 1.0;
    const double x1 = -1.0 / sqrt(3.0);
    const double x2 = 1.0 / sqrt(3.0);

    // Transformasi titik dari interval [-1, 1] ke [t, t+h]
    double t1 = t + (h / 2.0) * (1.0 + x1);
    double t2 = t + (h / 2.0) * (1.0 + x2);

    // Evaluasi turunan pada titik-titik Gauss
    // Untuk ODE, kita asumsikan turunan dievaluasi berdasarkan keadaan di awal step (t, theta, omega)
    double dtheta_dt1 = f_theta(t1, theta, omega);
    double domega_dt1 = f_omega(t1, theta, omega);

    double dtheta_dt2 = f_theta(t2, theta, omega);
    double domega_dt2 = f_omega(t2, theta, omega);

    // Hitung integral menggunakan rumus kuadratur Gauss
    double integral_theta = (h / 2.0) * (w1 * dtheta_dt1 + w2 * dtheta_dt2);
    double integral_omega = (h / 2.0) * (w1 * domega_dt1 + w2 * domega_dt2);

    // Update nilai theta dan omega
    theta += integral_theta;
    omega += integral_omega;
}

// --- 2. Metode Integrasi Romberg ---
// Menggunakan ekstrapolasi Richardson pada aturan trapesium.
// Di sini kita implementasikan 2 level untuk mendapatkan hasil O(h^4).

// Helper: Aturan Trapesium
double trapezoidal_rule(const function<double(double)> &func, double a, double b, int n)
{
    double step = (b - a) / n;
    double sum = func(a) + func(b);
    for (int i = 1; i < n; ++i)
    {
        sum += 2 * func(a + i * step);
    }
    return step * sum / 2.0;
}

void romberg_method(double t, double &theta, double &omega)
{
    // Buat fungsi lambda untuk turunan, menangkap keadaan saat ini
    auto d_theta_func = [&](double tau)
    { return f_theta(tau, theta, omega); };
    auto d_omega_func = [&](double tau)
    { return f_omega(tau, theta, omega); };

    // Terapkan Romberg untuk theta
    double I1_theta = trapezoidal_rule(d_theta_func, t, t + h, 1); // h
    double I2_theta = trapezoidal_rule(d_theta_func, t, t + h, 2); // h/2
    double integral_theta = (4.0 * I2_theta - I1_theta) / 3.0;

    // Terapkan Romberg untuk omega
    double I1_omega = trapezoidal_rule(d_omega_func, t, t + h, 1); // h
    double I2_omega = trapezoidal_rule(d_omega_func, t, t + h, 2); // h/2
    double integral_omega = (4.0 * I2_omega - I1_omega) / 3.0;

    // Update nilai theta dan omega
    theta += integral_theta;
    omega += integral_omega;
}

// --- 3. Metode Kuadratur Adaptif ---
// Menyesuaikan step secara rekursif hingga toleransi error tercapai.
// CATATAN: Implementasi adaptif sejati akan mengubah 'h', yang memerlukan
// perubahan pada loop simulasi utama. Untuk menjaga struktur yang ada,
// fungsi ini menunjukkan prinsip adaptif dalam satu langkah 'h' yang tetap.

// Helper: Fungsi rekursif untuk kuadratur adaptif
double adaptive_quad_recursive(const function<double(double)> &func, double a, double b, double tolerance)
{
    double h_step = b - a;
    double c = (a + b) / 2.0;

    // Hitung integral dengan 1 langkah Simpson (kasar) dan 2 langkah Simpson (halus)
    double I_kasar = (h_step / 6.0) * (func(a) + 4 * func(c) + func(b));
    double I_halus = (h_step / 12.0) * (func(a) + 4 * func((a + c) / 2.0) + 2 * func(c) + 4 * func((c + b) / 2.0) + func(b));

    // Perkirakan error
    double error = abs(I_halus - I_kasar) / 15.0;

    if (error < tolerance)
    {
        return I_halus + (I_halus - I_kasar) / 15.0; // Hasil dengan koreksi error
    }
    else
    {
        // Jika error terlalu besar, bagi interval dan panggil secara rekursif
        return adaptive_quad_recursive(func, a, c, tolerance / 2.0) +
               adaptive_quad_recursive(func, c, b, tolerance / 2.0);
    }
}

void adaptive_method(double t, double &theta, double &omega)
{
    double tolerance = 1e-6; // Toleransi error untuk satu langkah

    // Buat fungsi lambda untuk turunan
    auto d_theta_func = [&](double tau)
    { return f_theta(tau, theta, omega); };
    auto d_omega_func = [&](double tau)
    { return f_omega(tau, theta, omega); };

    // Hitung integral adaptif
    double integral_theta = adaptive_quad_recursive(d_theta_func, t, t + h, tolerance);
    double integral_omega = adaptive_quad_recursive(d_omega_func, t, t + h, tolerance);

    // Update nilai theta dan omega
    theta += integral_theta;
    omega += integral_omega;
}

// =================================================================
// FUNGSI SIMULASI
// =================================================================

// Fungsi untuk menulis hasil ke file CSV
void write_to_csv(const string &filename, const vector<double> &time,
                  const vector<double> &theta, const vector<double> &omega)
{
    ofstream outfile(filename);
    if (!outfile.is_open())
    {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    outfile << "t,theta,omega" << endl;
    outfile << fixed << setprecision(8);

    for (size_t i = 0; i < time.size(); ++i)
    {
        outfile << time[i] << "," << theta[i] << "," << omega[i] << endl;
    }

    outfile.close();
    cout << "Data berhasil ditulis ke " << filename << endl;
}

// Fungsi utama simulasi
void run_simulation(bool with_external_force, const string &method_name,
                    void (*method_func)(double, double &, double &))
{

    // Set gaya eksternal berdasarkan kasus
    A = with_external_force ? 0.5 : 0.0;

    vector<double> time_vec, theta_vec, omega_vec;
    double t = t0;
    double theta = theta0;
    double omega = omega0;

    cout << "Metode: " << method_name << "..." << endl;

    while (t <= tf)
    {
        time_vec.push_back(t);
        theta_vec.push_back(theta);
        omega_vec.push_back(omega);

        method_func(t, theta, omega);
        t += h;
    }

    string scenario = with_external_force ? "with_force" : "no_force";
    string filename = "pendulum_" + method_name + "_" + scenario + ".csv";
    write_to_csv(filename, time_vec, theta_vec, omega_vec);
}

int main()
{

    // Kasus 1: Dengan gaya eksternal
    cout << "Memulai Skenario 1: Dengan Gaya Eksternal..." << endl;
    run_simulation(true, "gauss", gauss_quadrature_method);
    run_simulation(true, "romberg", romberg_method);
    run_simulation(true, "adaptive", adaptive_method);

    // Kasus 2: Tanpa gaya eksternal
    cout << "\nMemulai Skenario 2: Tanpa Gaya Eksternal..." << endl;
    run_simulation(false, "gauss", gauss_quadrature_method);
    run_simulation(false, "romberg", romberg_method);
    run_simulation(false, "adaptive", adaptive_method);

    cout << "\nSemua simulasi selesai." << endl;

    return 0;
}