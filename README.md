# R model using Partition survival structure

## Background
Non-small cell lung cancer (NSCLC) with MET exon 14 skipping mutations represents a significant clinical challenge. This project addresses the pharmacoeconomic evaluation of two MET inhibitors - gumarontinib versus savolitinib - from the perspective of the Chinese healthcare system, providing crucial cost-effectiveness insights for healthcare decision-makers. It performs survival extrapolation, calculates QALYs and costs, and conducts sensitivity analyses (DSA, PSA) to evaluate economic outcomes.

## Model Overview
This R implementation of a Partitioned Survival Model (PSM) evaluates the cost-effectiveness of gumarontinib compared to savolitinib for treating METex14 skipping NSCLC. The model follows CHEERS 2022 guidelines and China Pharmacoeconomic Evaluation Guidelines 2020. By leveraging clinical trial data and healthcare cost parameters, the model provides comprehensive insights into treatment costs, survival outcomes, quality-adjusted life years, and incremental cost-effectiveness ratios from the perspective of the Chinese healthcare system.

## Model Structure
- **Model Type**: 3-state Partitioned Survival Model
- **Health States**: Progression-Free Survival (PFS), Progressive Disease (PD), Death
- **Time Horizon**: Lifetime (until 99% patient mortality)
- **Cycle Length**: 1 month
- **Perspective**: Chinese Healthcare System
<img width="418" height="373" alt="image" src="https://github.com/user-attachments/assets/fc0376f9-8759-41fb-b076-cf752f25eddc" />

## Key Features
- **Survival Extrapolation**: Multiple distribution support (Exponential, Weibull, Lognormal, Gompertz, Loglogistic)
- **Matching-Adjusted Indirect Comparison**: Unanchored MAIC for cross-trial efficacy adjustment
- **Comprehensive Costing**: Drug costs, disease management, AE management, terminal care
- **Utility Integration**: Quality-adjusted life years calculation
- **Sensitivity Analysis**: Deterministic and probabilistic uncertainty assessment
<img width="409" height="317" alt="image" src="https://github.com/user-attachments/assets/7219a9b9-bf48-41d2-8601-279dcdf62a5e" />
<img width="408" height="301" alt="image" src="https://github.com/user-attachments/assets/83ee6571-4bc8-4b25-8149-4a8e47158686" />

## Pipeline
1. Data import and parameter initialization
2. Survival curve generation and extrapolation
3. Health state occupancy calculation
4. Cost and QALY accumulation
5. Incremental cost-effectiveness ratio calculation
6. Sensitivity and scenario analyses

## Model Parameters
### Clinical Inputs
- PFS and OS curves from GLORY and NCT02897479 trials
- Adverse event incidence and management
- Treatment duration and dosing schedules

### Economic Inputs
- Drug costs: Gumarontinib vs Savolitinib
- Disease management costs
- AE management costs
- Terminal care costs
- Discount rate: 5% annually

### Utility Parameters
- PFS utility: 0.804
- PD utility: 0.321
- AE disutilities

## Results
### Base Case Analysis
- **Incremental QALYs**: 0.977
- **Incremental Cost**: $1,131
- **ICER**: $11,582/QALY
- **WTP Threshold**: $35,007/QALY (3× China GDP per capita)

### Key Findings
Gumarontinib demonstrates cost-effectiveness compared to savolitinib, with ICER below the willingness-to-pay threshold.

## Uncertainty Analysis
### Deterministic Sensitivity Analysis
- Tornado diagrams for parameter impact assessment
- Key drivers: Drug costs and utility values

### Probabilistic Sensitivity Analysis
- 1,000 Monte Carlo simulations
- Cost-effectiveness acceptability curves
- Around 60.2% probability of cost-effectiveness at WTP threshold

### Scenario Analyses
- Dosage adjustment scenarios
- Alternative time horizons
- Parametric distribution variations

## Getting Started
### Prerequisites
```r
install.packages(c("readxl", "dplyr", "scales", "reshape2", "ggplot2"))
