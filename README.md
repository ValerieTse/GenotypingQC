# Genomic Data QC Pipeline

This repository outlines the quality control (QC) workflow for genomic data processing. The pipeline includes environmental setup, file updating, sample validation, sex checks, population structure analysis, and variant filtering. Follow these steps to ensure that the genotype data is clean and suitable for downstream analysis.

## 1. Preparing

### 1a. Creating a New Conda Environment
- **Objective:** Set up an isolated environment to manage dependencies.
- **Steps:**
  - Create a new conda environment.
  - Install the required packages:
    - `neural-admixture` -> For population inference.
    - `illuminaio` -> For chromosome X and Y intensities reading.
- **Commands:**
  ```bash
  conda create -n qc_env python=3.9
  conda activate qc_env
  pip install neural-admixture
  conda install bioconda::bioconductor-illuminaio

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
- **Commands:**
  ```bash
  plink --bfile your_data \
        --update_ids your_ids.txt \
        --update_sex your_sex.txt \
        --pheno your_pheno.txt \
        --make-bed --out your_data
  
---

## 3. Initial Check

- **Objective:** Evaluate the raw dataset before performing QC.
- **Steps:**
  - **Sample and Variant Counts:** Confirm the sample size and variant coverage.
  - **Missingness Analysis:** Check the missing data rates for both samples and variants.
- **Commands:**
  ```bash
  plink --bfile your_data \
        --missing \
        --out missingness
- **Outputs:**
  - If you used `plink1.9`, you will get `.lmiss` and `.imiss` files.
  - If you used `plink2`, you will get `.smiss` and `.vmiss` files.
  - Check the Plink website for the outputs.
  
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
     - **Commands:**
       ```bash
       plink --bfile your_data \
             --mind 0.1 \
             --geno 0.1 \
             --make-bed --out your_data
  4. **Post-filtering Check:**
     - Re-assess sample size, variant count, and overall missingness statistics after filtering.

---

## 5. Sex Check

- **Objective:** Validate the sex of each sample using genetic data.
- **Steps:**
  1. **X Chromosome Heterozygosity:**
     - Determine the appropriate threshold for inferring genetic sex.
     - **Commands:**
       ```bash
       plink --bfile your_data \
             --check-sex \
             --out your_data
  2. **Comparative Analysis:**
     - Compare genetic sex assignments with the annotated (metadata) sex.
     - Identify and flag samples with discrepancies.
  3. **Y Chromosome Missingness:**
     - Evaluate the Y chromosome call rates in problematic samples.
     - **Commands:**
       ```bash
       plink --bfile your_data \
             --keep problematic_list.txt \
             --chr Y --missing \
             --out your_data
  4. **Intensity Checks (Optional):**
     - If `.idat` files are available, extract and inspect X and Y chromosome intensities.
     - Using the `illuminaio` R package.
  5. **X Chromosome heterozygosity rate:**
     - The PLINK `--check-sex` command only provides an F-statistic.
     - To estimate the heterozygosity rate for the X chromosome:
       - Extract the X chromosome data. `--chr X`
       - Recode it as chromosome 1 (satisfying the requirement for the `--het` command).
       - Calculate the heterozygosity rate accordingly.

---

## 6. Population Structure Analysis

- **Objective:** Account for population stratification, which can affect heterozygosity and HWE tests.
- **Steps:**
  1. **Reference Dataset:**
     - Download a reference dataset (e.g., 1000 Genomes Phase 3).
  2. **LD Pruning:**
     - Perform linkage disequilibrium (LD) pruning on both your data and the reference dataset.
     - **Commands:**
     ```bash
     plink --bfile your_data \
           --indep-pairwise 50 5 0.2 \
           --out your_data
     
     plink --bfile your_data \
           --extract your_data.prune.in \
           --make-bed --out your_data
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
    - Raw heterozygosity rate = Corrected heterozygosity rate + PCs
  - Flag samples with extreme heterozygosity values for further review.
- **Command:**
  ```bash
  plink --bfile your_data \
        --het  \
        --out your_data
---

## 8. Hardy-Weinberg Equilibrium (HWE) Testing

- **Objective:** Assess variant quality within subgroups by testing for deviations from HWE.
- **Steps:**
  - Conduct HWE tests by race.
  - Adjust the threshold based on the detected population structure.
- **Command:**
  ```bash
  plink --bfile your_data \
        --hardy \
        --out your_data

---

## 9. Minor Allele Frequency (MAF) Filtering

- **Objective:** Filter variants based on their frequency to remove those that may be too rare for robust analysis.
- **Steps:**
  - Define a threshold for the minimum acceptable minor allele frequency.
  - Filter out variants falling below this threshold.
  - Re-assess variant counts post-filtering.
- **Command:**
  ```bash
  plink --bfile your_data \
        --freq  \
        --out your_data

---

## Conclusion

This README details the steps required for robust genomic data quality control. Each stage—from initial preparation to detailed filtering—ensures that only high-quality data moves forward into further analysis. For complete command details and scripts, consult the supplementary documentation provided in each module of the pipeline.

---

## References & Acknowledgements

- Reference datasets (e.g., 1000 Genomes Phase 3)
- Software tools: PLINK, neural-admixture, illuminaio
- Additional documentation and tutorials available in the repository.

