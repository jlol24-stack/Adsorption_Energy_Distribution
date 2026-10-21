# Adsorption Energy Distribution of As(III) on Manganese Ferrite (MnFe₂O₄) Nanoparticles

This repository provides a python-based pipeline to resolve the continuous active site energy distribution f(E) of manganese ferrite nanoparticles during \(\text{As(III)}\) capture. It utilizes a non-linear **Fowler-Guggenheim local kernel** integrated into a **Fredholm equation of the first kind** and solved via **Non-Negative Least Squares (NNLS)** optimization.

---

## 📂 Repository Structure

*   `dataeq.csv`: Input experimental dataset with comma-separated `Ce` (mg/L) and `qe` (mg/g) columns.
*   `adsorption_energy_distribution.ipynb`: Main documented workflow with LaTeX equations and execution cells.
*   `aed_profile_data.csv`: Output raw distribution spectrum coordinates optimized for OriginLab or Prism.
*   `aed_results.png`: Exported high-resolution diagnostic and model fitting validation plot (300 DPI).

---

## 🚀 Quick Start

1. **Install requirements:**
   ```bash
   pip install pandas numpy scipy matplotlib notebook
   ```
2. **Run the workflow:** Open your terminal, launch `jupyter notebook`, and run all cells sequentially in `adsorption_energy_distribution.ipynb`.

---

## ✒️ Citation Guidelines

Please provide appropriate attribution when utilizing this pipeline in your research:

*   **This Specific Component:** López-Luna, J., Vázquez-Hipólito, V., Ramírez-Montes, L.E. et al. (2026). Adsorption Energy Distribution of As(III) on Manganese Ferrite (MnFe₂O₄) Nanoparticles [Jupyter Notebook Software]. Zenodo. https://doi.org
*   **General Software Repository:** López-Luna, J., Vázquez-Hipólito, V., Ramírez-Montes, L.E. et al. (2026). Adsorption Behavior of Trivalent Arsenic on Manganese Ferrite Nanoparticles: Kinetics, Isotherms, and Rate-Limiting Steps [Jupyter Notebook Software]. Zenodo. v1 10.5281/zenodo.22695744; v2.0.0 10.5281/zenodo.23172980.
*   **Foundational Study:** López-Luna, J., Ramírez-Montes, L.E. et al. Linear and nonlinear kinetic and isotherm adsorption models for arsenic removal by manganese ferrite nanoparticles. *SN Appl. Sci.* **1**, 950 (2019). https://doi.org
