import os
import pandas as pd
import numpy as np
from scipy.optimize import nnls
from scipy.integrate import simpson
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

# =====================================================================
# 1. EXPERIMENTAL DATA LOADING & CLEANING
# =====================================================================
csv_filename = "dataeq.csv"
if not os.path.exists(csv_filename):
    raise FileNotFoundError(f"Error: The file '{csv_filename}' was not found in the directory.")

df_eq = pd.read_csv(csv_filename)
# Extracting numerical vectors from your dataset
Ce_raw = df_eq["Ce"].values.astype(float)
qe = df_eq["qe"].values.astype(float)

# Safety correction to avoid numerical indeterminacy (Ce <= 0)
Ce = np.where(Ce_raw <= 0, 1e-8, Ce_raw)

# Universal physical constants
R = 8.314e-3  # Gas constant in kJ/(mol·K)
T = 298.15    # Operational temperature in Kelvin (25 °C)

# =====================================================================
# 2. THERMODYNAMIC MODEL SELECTION (KERNEL CONFIGURATION)
# =====================================================================
KERNEL_TYPE = 'fowler-guggenheim'  
w_interaction = 1.0  # Lateral interaction parameter 'w' (kJ/mol)
K0 = 0.01            # Pre-exponential scale factor for chemical affinity

# =====================================================================
# 3. ENERGY SPECTRUM DISCRETIZATION (AED X-Axis)
# =====================================================================
E_min, E_max = 0.0, 40.0
num_E_points = 120  # Total sub-intervals for numerical integration
E_vector = np.linspace(E_min, E_max, num_E_points)

# =====================================================================
# 4. KERNEL MATRIX CONSTRUCTION (Fredholm Integral Equation)
# =====================================================================
print(f"Constructing local base matrix using model: {KERNEL_TYPE}...")
Kernel = np.zeros((len(Ce), num_E_points))

for j, E in enumerate(E_vector):
    # Calculate base affinity constant for this specific energy level
    Kl_base = K0 * np.exp(E / (R * T))
    
    if KERNEL_TYPE == 'langmuir':
        # Classic explicit Langmuir local model
        Kernel[:, j] = (Kl_base * Ce) / (1 + Kl_base * Ce)
        
    elif KERNEL_TYPE == 'fowler-guggenheim':
        # Non-linear implicit equation: solving local fractional coverage (theta)
        for i, c_val in enumerate(Ce):
            theta_space = np.linspace(0.001, 0.999, 1000)
            # Minimize error of the implicit Fowler-Guggenheim equality
            fg_error = np.abs(Kl_base * c_val - (theta_space / (1 - theta_space)) * np.exp(-2 * w_interaction * theta_space / (R * T)))
            optimal_theta = theta_space[np.argmin(fg_error)]
            Kernel[i, j] = optimal_theta

# =====================================================================
# 5. INVERSE PROBLEM RESOLUTION (Non-Negative Least Squares)
# =====================================================================
print("Solving the optimal adsorption energy distribution...")
# f_dist stores the calculated site density for each energy level
f_dist, residue = nnls(Kernel, qe)

# Matrix multiplication to obtain predicted qe values
qe_predicted = np.dot(Kernel, f_dist)

# =====================================================================
# 6. STATISTICAL GOODNESS-OF-FIT METRICS
# =====================================================================
# Coefficient of Determination (R²)
ss_res = np.sum((qe - qe_predicted) ** 2)
ss_tot = np.sum((qe - np.mean(qe)) ** 2)
r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

# Non-linear Chi-Squared Test (X²)
chi_squared = np.sum(((qe - qe_predicted) ** 2) / (qe_predicted + 1e-8))

# =====================================================================
# 7. QUANTITATIVE ANALYSIS VIA SIMPSON'S INTEGRATION
# =====================================================================
print("\n--- ACTIVE SITES REPOSITORY (RESULTS) ---")
total_area = simpson(y=f_dist, x=E_vector)
print(f"Total Integrated Capacity (Global Qm): {total_area:.2f} mg/g")
print(f"Goodness of Fit -> R²: {r_squared:.4f} | Chi-Squared (X²): {chi_squared:.4f}")

# Detect local maxima (peaks)
peaks, _ = find_peaks(f_dist, height=0.5)

# Partition at E = 8.0 kJ/mol (Physical vs Chemical Adsorption boundary)
idx_split = np.searchsorted(E_vector, 8.0)
area_physisorption = simpson(y=f_dist[:idx_split], x=E_vector[:idx_split])
area_chemisorption = simpson(y=f_dist[idx_split:], x=E_vector[idx_split:])

pct_phys = (area_physisorption / total_area) * 100 if total_area > 0 else 0
pct_chem = (area_chemisorption / total_area) * 100 if total_area > 0 else 0

print(f"-> Low-Energy Sites (E < 8 kJ/mol - Physisorption/Pores): {pct_phys:.1f}% of total capacity.")
print(f"-> High-Energy Sites (E >= 8 kJ/mol - Specific Chemisorption): {pct_chem:.1f}% of total capacity.")

if len(peaks) > 0:
    print("\nDetected Peak Maxima:")
    for p in peaks:
        print(f"   Site at E = {E_vector[p]:.1f} kJ/mol | Density f(E) = {f_dist[p]:.2f}")
print("-----------------------------------------\n")

# Export numerical raw profile data for Origin/Prism software
df_output = pd.DataFrame({"Energy_kJ_mol": E_vector, "Site_Density_fE": f_dist})
df_output.to_csv("aed_profile_data.csv", index=False)

# =====================================================================
# 8. RESULTS VISUALIZATION & HIGH-RESOLUTION SAVE
# =====================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))

# Graph 1: Adsorption Energy Distribution (AED) Profile Plot
ax1.plot(E_vector, f_dist, color='#16a085', linewidth=2.5, label=f'Kernel: {KERNEL_TYPE}')
ax1.fill_between(E_vector, f_dist, color='#16a085', alpha=0.15)
ax1.axvline(x=8.0, color='#7f8c8d', linestyle=':', label='Phys/Chem Boundary (8 kJ/mol)')
ax1.set_title('Adsorption Energy Distribution (AED)', fontsize=12, fontweight='bold')
ax1.set_xlabel('Interaction Energy $E$ (kJ/mol)', fontsize=11)
ax1.set_ylabel('Active Site Density $f(E)$ [(mg/g)/(kJ/mol)]', fontsize=11)
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.set_xlim(E_min, E_max)
ax1.legend()

# Graph 2: Statistical Fit Verification Plot
ax2.scatter(Ce, qe, color='#e74c3c', s=60, zorder=3, label='Experimental Data')
ax2.plot(Ce, qe_predicted, color='#2980b9', linestyle='--', linewidth=2.5, label=f'Model Fit (R² = {r_squared:.4f})')
ax2.set_title('Numerical Fit Verification', fontsize=12, fontweight='bold')
ax2.set_xlabel('$C_e$ (Equilibrium Concentration) [mg/L]', fontsize=11)
ax2.set_ylabel('$q_e$ (Adsorption Capacity) [mg/g]', fontsize=11)
ax2.legend()
ax2.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig("aed_results.png", dpi=300)
plt.show()
