# R Scripts for the manuscript entitled "Lifespan epigenetic clocks delineate organ-specific aging landscapes and centenarian resilience signatures"
This repository contains R scripts used to perform the core analyses described in the manuscript. The scripts cover differential methylation analysis, annotation of organ-enriched CpGs, construction of the lifespan and organ-specific DNA methylation clocks, calculation of established epigenetic clocks, and downstream association analyses.
## File Descriptions
| File | Description |
|------|-------------|
| `1-DMP.R` | Identification of differentially methylated positions (DMPs) associated with chronological age. |
| `2-Organ-enriched genes.R` | Identification of organ-enriched genes from GTEx data. |
| `3-CpG annotation.R` | Annotation of CpG sites mapped to organ-enriched genes. |
| `4-ML_clock.R` | Training and evaluation of the organismal lifespan DNA methylation clock. |
| `5-ML_clock_organ.R` | Training and evaluation of the eight organ-specific DNA methylation clocks. |
| `6-Other DNAm clock.R` | Calculation of established epigenetic clocks, including Horvath, Hannum, PhenoAge, and principal component versions. |
| `7-DNAmAgeDev.R` | Calculation of DNA methylation age deviation (DNAmAgeDev) using linear regression approach. |
| `8-Mortality association analyses.R` | Association analyses between DNAmAgeDev and mortality. |
| `9-Disease association analyses.R` | Association analyses between DNAmAgeDev and incident age-related and organ-specific diseases. |
| `10-Biomarker association analyses.R` | Age-stratified association analyses between the clocks and physiological, biochemical, and functional biomarkers. |
| `11-Analyses of centenarians.R` | Characterization of the Adipose–Kidney–Muscle slow-aging pattern and related analyses in centenarians. |
| `12-Lifestyle analyses.R` | Assessment of the relative importance of lifestyle factors and their associations with DNAmAgeDev across life stages. |
| `13-Aging peaks.R` | Identification of organismal and organ-specific aging peaks using the modified DE-SWAN framework. |
## Usage Notes
- Scripts are intended to be run in the order listed to reproduce the main analyses.
- Data paths and cohort-specific variables may need to be adjusted according to local file structures.
- Required R packages are loaded within each script.
