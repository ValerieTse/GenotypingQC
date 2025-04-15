# Genomic Data QC Pipeline

This repository outlines the quality control (QC) workflow for genomic data processing. The pipeline includes environmental setup, file updating, sample validation, sex checks, population structure analysis, and variant filtering. Follow these steps to ensure that the genotype data is clean and suitable for downstream analysis.

## 1. Preparing

### 1a. Creating a New Conda Environment
- **Objective:** Set up an isolated environment to manage dependencies.
- **Steps:**
  - Create a new conda environment.
  - Install the required packages:
    - `neural-admixture`
    - `illuminaio`

### 1b. Directory Setup
- **Objective:** Organize your work environment.
- **Steps:**
  - Create a dedicated working directory to manage the QC process.
  - Use subdirectories to separate raw data, intermediate files, and final results.

### 1c. Module Loading
- **Objective:** Ensure all necessary software modules are available.
- **Steps:**
  - Load required modules if using an HPC cluster or module system.
  - Verify that environment variables and paths for software like PLINK, Python, etc., are correctly configured.

---

## 2. Update the .fam File

- **Objective:** Enrich the `.fam` file with missing metadata.
- **Details:**
  - Most `.fam` entries are missing detailed annotations such as family ID, phenotype data, and annotated sex.
  - Use available metadata to update these fields.
- **Steps:**
  - Parse the metadata information.
  - Update the `.fam` file to include:
    - Family information
    - Phenotype
    - Annotated sex

---

## 3. Initial Check

- **Objective:** Evaluate the raw dataset before performing QC.
- **Steps:**
  - **Sample and Variant Counts:** Confirm the sample size and variant coverage.
  - **Missingness Analysis:** Check the missing data rates for both samples and variants.
  
---

## 4. Filtering Based on Missingness

- **Objective:** Remove low-quality samples and variants with high missingness.
- **Steps:**
  1. **Visualization:**
     - Plot the call rates for both samples and variants.
     - Identify outliers and potential quality issues.
  2. **Threshold Determination:**
     - Decide on a threshold for acceptable call rates.
  3. **Filtering:**
     - Remove samples or variants that fall below the set threshold.
  4. **Post-filtering Check:**
     - Re-assess sample size, variant count, and overall missingness statistics after filtering.

---

## 5. Sex Check

- **Objective:** Validate the sex of each sample using genetic data.
- **Steps:**
  1. **X Chromosome Heterozygosity:**
     - Determine the appropriate threshold for inferring genetic sex.
  2. **Comparative Analysis:**
     - Compare genetic sex assignments with the annotated (metadata) sex.
     - Identify and flag samples with discrepancies.
  3. **Y Chromosome Missingness:**
     - Evaluate the Y chromosome call rates in problematic samples.
  4. **Intensity Checks (Optional):**
     - If `.idat` files are available, extract and inspect X and Y chromosome intensities.
  5. **Addressing PLINK Limitations:**
     - The PLINK `--check-sex` command only provides an F-statistic.
     - To estimate the heterozygosity rate for the X chromosome:
       - Extract the X chromosome data.
       - Recode it as chromosome 1 (satisfying the requirement for the `--het` command).
       - Calculate the heterozygosity rate accordingly.

---

## 6. Population Structure Analysis

- **Objective:** Account for population stratification which can affect heterozygosity and HWE tests.
- **Steps:**
  1. **Reference Dataset:**
     - Download a reference dataset (e.g., 1000 Genomes Phase 3).
  2. **LD Pruning:**
     - Perform linkage disequilibrium (LD) pruning on both your data and the reference dataset.
  3. **Data Merging:**
     - Merge the pruned datasets.
     - Follow step-by-step checks to ensure data compatibility and correctness.
  4. **Principal Component Analysis (PCA):**
     - Conduct PCA to obtain 20 principal components.
  5. **Population Inference:**
     - Use the supervised version of neural-admixture to infer population structure.

---

## 7. Sample Heterozygosity Analysis

- **Objective:** Detect samples with aberrant heterozygosity which may indicate quality problems.
- **Steps:**
  - Use principal components to compute a population-corrected heterozygosity measure.
  - Flag samples with extreme heterozygosity values for further review.

---

## 8. Hardy-Weinberg Equilibrium (HWE) Testing

- **Objective:** Assess variant quality within subgroups by testing for deviations from HWE.
- **Steps:**
  - Conduct HWE tests separately for cases and controls.
  - Adjust the analysis based on the detected population structure.

---

## 9. Minor Allele Frequency (MAF) Filtering

- **Objective:** Filter variants based on their frequency to remove those that may be too rare for robust analysis.
- **Steps:**
  - Define a threshold for the minimum acceptable minor allele frequency.
  - Filter out variants falling below this threshold.
  - Re-assess variant counts post-filtering.

---

## Conclusion

This README details the steps required for robust genomic data quality control. Each stage—from initial preparation to detailed filtering—ensures that only high-quality data moves forward into further analysis. For complete command details and scripts, consult the supplementary documentation provided in each module of the pipeline.

---

## References & Acknowledgements

- Reference datasets (e.g., 1000 Genomes Phase 3)
- Software tools: PLINK, neural-admixture, illuminaio
- Additional documentation and tutorials available in the repository.

