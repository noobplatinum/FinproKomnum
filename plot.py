import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Konfigurasi plot
sns.set(style="whitegrid", palette="deep")
plt.rcParams['figure.figsize'] = (18, 22)
plt.rcParams['font.size'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 12

def load_data(method, scenario):
    """Fungsi untuk memuat data dari file CSV."""
    scenario_str = "with_force" if scenario else "no_force"
    filename = f"pendulum_{method}_{scenario_str}.csv"
    
    filepath = os.path.join("csv", filename)
    
    if not os.path.exists(filepath):
        print(f"PERINGATAN: File {filepath} tidak ditemukan.")
        return None
    
    df = pd.read_csv(filepath)
    method_map = {'gauss': 'Gauss', 'romberg': 'Romberg', 'adaptive': 'Adaptive'}
    df['Method'] = method_map.get(method, method.upper())
    df['Scenario'] = 'With Force' if scenario else 'No Force'
    return df

# Inisialisasi figure dan grid layout
fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(5, 3)

# Judul utama
fig.suptitle('Perbandingan Simulasi Pendulum: Metode Gauss vs Romberg vs Adaptif', 
            fontsize=20, y=1.03)

methods = ['gauss', 'romberg', 'adaptive']
scenarios = [True, False]
all_data = []

# Memuat semua data dari file CSV
for method in methods:
    for scenario in scenarios:
        data = load_data(method, scenario)
        if data is not None:
            all_data.append(data)

if not all_data:
    print("Tidak ada data yang dapat di-plot. Pastikan file CSV sudah ada.")
else:
    combined_data = pd.concat(all_data)

    # Plot 1: Sudut vs Waktu (Dengan Gaya Eksternal)
    ax1 = fig.add_subplot(gs[0:2, :])
    sns.lineplot(data=combined_data[combined_data['Scenario'] == 'With Force'], 
                x='t', y='theta', hue='Method', ax=ax1, linewidth=1.5)
    ax1.set_title('Sudut vs Waktu (Dengan Gaya Eksternal)')
    ax1.set_xlabel('Waktu (s)')
    ax1.set_ylabel('Sudut (rad)')
    ax1.legend(title='Metode')
    ax1.grid(True, linestyle='--', alpha=0.6)

    # Plot 2: Sudut vs Waktu (Tanpa Gaya Eksternal)
    ax2 = fig.add_subplot(gs[2:4, :])
    sns.lineplot(data=combined_data[combined_data['Scenario'] == 'No Force'], 
                x='t', y='theta', hue='Method', ax=ax2, linewidth=1.5)
    ax2.set_title('Sudut vs Waktu (Tanpa Gaya Eksternal)')
    ax2.set_xlabel('Waktu (s)')
    ax2.set_ylabel('Sudut (rad)')
    ax2.legend(title='Metode')
    ax2.grid(True, linestyle='--', alpha=0.6)

    # Diagram Fase untuk setiap metode
    plot_positions = {
        'Gauss (With Force)': gs[4, 0], 'Romberg (With Force)': gs[4, 1], 'Adaptive (With Force)': gs[4, 2]
    }

    for title, position in plot_positions.items():
        ax = fig.add_subplot(position)
        
        method_name = title.split(' ')[0]
        scenario_name = 'With Force'
        
        subset = combined_data[(combined_data['Method'] == method_name) & 
                            (combined_data['Scenario'] == scenario_name)]
        
        if not subset.empty:
            ax.plot(subset['theta'], subset['omega'], linewidth=1.0)
        
        ax.set_title(title)
        ax.set_xlabel('Sudut (rad)')
        ax.set_ylabel('Kecepatan Sudut (rad/s)')
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.tick_params(axis='x', rotation=45)

    # Simpan hasil plot
    output_filename = "pendulum_comparison_new_methods.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\nPlot perbandingan telah disimpan sebagai: {output_filename}")
    plt.close()