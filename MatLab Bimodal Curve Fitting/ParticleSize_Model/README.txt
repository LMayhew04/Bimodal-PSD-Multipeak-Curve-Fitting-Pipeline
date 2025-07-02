# Particle Size Distribution Fitting Pipeline (Final, Improved Version)

## QUICK START GUIDE

1. **Requirements:**  
   - MATLAB R2021a or later  
   - Curve Fitting Toolbox  
   - Statistics and Machine Learning Toolbox  

2. **Input Data:**  
   - Place your Excel particle size data files in the `Data/` folder.  
   - Example: `Data/TSE001 75 min run ... .xlsx`  
   - Data columns should be labeled (e.g., `ParticleDiameter`, `SizeDistribution1`), or configure in the settings section.

3. **How to Run:**  
   - Open MATLAB and set your working directory to the root project folder.
   - Run the main driver script:  
     ```matlab
     runAllFits
     ```
   - By default, the script will:
     - Fit all sheets in all data files (except any marked as malformed or empty).
     - Save output plots to `Figures/` and the summary table as `Fit_Summary.csv`.

4. **Sheet/Column Selection:**  
   - To run on only certain sheets or columns, edit the “SETTINGS” section at the top of `runAllFits.m`.

5. **Results:**  
   - Each fit’s plot (PDF/CDF) is saved as a PNG in `Figures/`.
   - The table `Fit_Summary.csv` lists model type, residuals, peak locations, and visual notes.

## OVERVIEW

This pipeline fits particle size distributions from spreadsheet datasets using robust, visually-validated multimodal models (normal/lognormal mixtures).  
- The model is selected based on **visual fidelity to all peaks**, not just BIC.
- Improved minor/shoulder peak detection.
- Axis scaling is always data-driven.
- Fully beginner-friendly: all settings centralized, code heavily commented.

See CHANGELOG.txt for version-to-version details.  
See TROUBLESHOOTING.txt for help on common issues.

---

## FILE STRUCTURE

- `fitBestModel.m` — fits the optimal model for a single dataset
- `runAllFits.m` — runs the full pipeline on all files/sheets
- `plotFit.m` — generates PDF/CDF and saves PNGs
- `loadParticleData.m` — loads and preprocesses data from Excel
- `pdfHelpers/` — component PDFs (Normal, Lognormal, Rayleigh)
- `Data/` — input .xlsx files
- `Figures/` — output plots
- `Fit_Summary.csv` — fit quality/results summary