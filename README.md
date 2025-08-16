# NHS-Efficiency-Analysis-Data
This repository contains supplementary materials for the MSc dissertation on hospital efficiency analysis (University of Edinburgh, 2025).  
All data and scripts are provided for transparency, reproducibility, and further research.

## 📂 Contents
- `NHS_Scotland_Simulated_Dataset.csv`  
  Simulated hospital-level dataset (2019–2024), including key operational variables used in the ratio, regression, and DEA stages.  

- `Summary_Stats_all_years.csv`  
  Summary statistics of efficiency scores (TE, PTE, SE; EE, AE) with confidence intervals, grouped by year and hospital size.  

- `Efficiency_Scores_all_years_by_hospital.csv`  
  Hospital-level efficiency estimates across all years.  

- `Hospital_Slack_full_2024.csv`  
  Full slack diagnostics for all hospitals in 2024, including input and output slacks.  

- `dea_analysis.R`  
  R script used to run DEA analysis, bootstrapping, and summary table generation.  

## 🔎 How to Use
1. Download the CSV files or clone the repository.  
2. Open the R scripts in the “Code" folder. These include three parts:  
   - `4.1.1_RatioAnalysis.R`  
   - `4.1.2_Regression.R`  
   - `4.1.3_DEA.R`  
3. Update the file paths if needed, then run the scripts to reproduce the results.  


