# GOLD-R BioAge: A Residual Approach to Biological Age Estimation

## Overview

GOLD-R BioAge is a novel computational framework for estimating biological age using a residual-informed approach grounded in the Gompertz law of mortality. This method provides robust aging measurements across multiple data modalities including DNA methylation, proteomics, and clinical biomarkers.

## Background

Biological age (BA) quantifies an individual's aging status beyond chronological age. While conventional BA measures often lack robustness in heterogeneous populations, GOLD-R addresses these limitations by directly estimating BA residuals optimized for cross-sectional data - the most common scenario in clinical practice and research settings.

## Methodology

### Core Algorithm

The GOLD-R framework employs a multi-step, residual-informed approach:
https://jerryhaom.github.io/GOLDR-BioAge/Methods.html


## Applications

### 1. Epigenetic Aging (GOLD-R DNAmAge)
- **Data**: DNA methylation data from EWAS Data Hub (7,313 samples, 21 tissues)
- **Features**: 33 CpG sites (significantly fewer than GrimAge's 1,030 sites)
- **Performance**: Superior mortality prediction in pan-cancer cohorts (C-index: 0.623)

### 2. Proteomic Aging (GOLD-R ProtAge)
- **Data**: UK Biobank proteomics (53,014 participants, 2,923 proteins)
- **Scope**: Organismal and organ-specific aging measures
- **Performance**: Enhanced prediction of all-cause mortality (C-index: 0.740 vs 0.588 for conventional methods)

### 3. Clinical Biomarker Integration (GOLD-R DNAmCli Age)
- **Data**: NHANES and HRS cohorts
- **Approach**: Predicting epigenetic residuals using clinical biomarkers
- **Performance**: Competitive mortality prediction (C-index: 0.75-0.77)

## Key Features

- **Robust Performance**: Validated across multiple datasets and populations
- **Cross-sectional Optimization**: Designed for common clinical scenarios
- **Multi-modal Applicability**: Works with DNA methylation, proteomics, and clinical data
- **Organ-Specific Analysis**: Enables tissue-specific aging assessment
- **Clinical Translation**: Practical for research and healthcare applications

## Benchmark Results

GOLD-R demonstrates superior performance compared to established epigenetic clocks:
- Outperforms first-generation clocks (HorvathAge, HannumAge)
- Comparable to GrimAge with significantly fewer features
- Superior to phenotypic clocks (PhenoAge, KDM Age) in clinical applications

## Citation
If you use GOLD-R in your research, please cite:
1. Hao, M., Zhang, H., et al. (2025). Gompertz Law-Based Biological Age (GOLD BioAge): A Simple and Practical Measurement of Biological Ageing. *Advanced Science*.
2. Zhang, H., Zhang, S., et al. (2025). A Residual Approach to Estimate Biological Age from Gompertz Modeling. *Under Review*.
