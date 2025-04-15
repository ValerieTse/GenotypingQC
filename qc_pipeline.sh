#!/bin/bash
# ==========================================================================================
# Script Name: qc_pipeline.sh
# Author: Jingyu Xie
# Date: 2025-04-14
# Version: v1.0
# Description: This script demonstrates the basic QC pipeline for genotypeing data in Bash.
# Data: The genotyping data are typical PLINK format.
# ==========================================================================================

# ----- Step 0: Create directory and module load -----
conda create -n geno_qc python=3.9
conda activate geno_qc 

module load plink2 # plink2 is still at the beta phase, if you download plink1.9, use it instead.
raw_data=/project/cshu0237_1136/Data/SPARK/iWES_v3/genotypes/genotyping_array/plink/SPARK.iWES_v3.2024_08.GSA-24_v2 # raw data, directory + prefix

mkdir DATA PROGRAM OUTPUTS Reference_data Release_Notes # creating the directory

# ----- Step 1: Update the .fam file -----
Rscript PROGRAM/update_info_extract.r /project/cshu0237_1136/Data/SPARK/iWES_v3/genotypes/genotyping_array/plink/SPARK.iWES_v3.2024_08.GSA-24_v2 /project/cshu0237_1136/Data/SPARK/iWES_v3/metadata/SPARK.iWES_v3.2024_08.sample_metadata.tsv

plink2 --bfile $raw_data --update-ids Release_Notes/update_id.txt --make-bed --out DATA/step1_update_ids
plink2 --bfile DATA/step1_update_ids --update-sex Release_Notes/update_sex.txt --make-bed --out DATA/step1_update_sex 
plink2 --bfile DATA/step1_update_sex --pheno Release_Notes/update_pheno.txt --make-bed --out DATA/step1_update_pheno

Rscript PROGRAM/counts.r DATA/step1_update_pheno # Check the number of samples and SNPs

# ----- Step 2: Initial check -----
plink2 --bfile DATA/step1_update_pheno --missing --out DATA/step2_overall_missing

Rscript PROGRAM/plot_call_rate.r DATA/step2_overall_missing.smiss DATA/step2_overall_missing.vmiss step2

plink2 --bfile DATA/step1_update_pheno --chr Y --missing --out DATA/step2_chrY_missing
plink2 --bfile DATA/step1_update_pheno --chr X --missing --out DATA/step2_chrX_missing
plink2 --bfile DATA/step1_update_pheno --chr MT --missing --out DATA/step2_chrMT_missing

# print the maximum missing rate for samples and SNPs
Rscript PROGRAM/call_rate.r DATA/step2_chrY_missing.smiss DATA/step2_chrY_missing.vmiss 
Rscript PROGRAM/call_rate.r DATA/step2_chrX_missing.smiss DATA/step2_chrX_missing.vmiss 
Rscript PROGRAM/call_rate.r DATA/step2_chrMT_missing.smiss DATA/step2_overall_missing.vmiss

# ----- Step 3: Missing call rate filtering -----
plink2 --bfile DATA/step1_update_pheno --mind 0.1 --make-bed --out DATA/step3_samples_filtered
plink2 --bfile DATA/step3_samples_filtered --geno 0.1 --make-bed --out DATA/step3_variants_filtered

plink2 --bfile DATA/step3_variants_filtered --missing --out DATA/step3_overall_missing

Rscript PROGRAM/plot_call_rate.r DATA/step3_overall_missing.smiss DATA/step3_overall_missing.vmiss step3

Rscript PROGRAM/counts.r DATA/step3_variants_filtered

# ----- Step 4: Sex check -----
## ---- Step 4.1: Chr X heterozygosity ----
plink1.9 --bfile DATA/step3_samples_filtered --check-sex --out DATA/step4_sexcheck # Hint: Using plink1.9 rather than plink2
grep 'PROBLEM' DATA/step4_sexcheck.sexcheck > DATA/step4_sex_discrepancy.txt 
wc -l DATA/step4_sex_discrepancy.txt  # count the problematic sample based on threshold 0.2 and 0.8

Rscript PROGRAM/plot_sexcheck.r DATA/step4_sexcheck.sexcheck step4

plink1.9 --bfile DATA/step3_samples_filtered --check-sex 0.4 0.6 --out DATA/step4_sexcheck_new # Based on previous results, re-select the threshold for genetic sex inference
grep 'PROBLEM' DATA/step4_sexcheck_new.sexcheck > DATA/step4_sex_discrepancy_new.txt
wc -l DATA/step4_sex_discrepancy_new.txt # count the problematic sample based on threshold 0.4 and 0.6
awk '{print $1, $2}' DATA/step4_sex_discrepancy_new.txt > DATA/step4_problematic_list.txt

plink1.9 --bfile DATA/step3_samples_filtered --recode AD --out DATA/step4_recode


## ---- Step 4.2: chr Y missingness ----
plink2 --bfile DATA/step3_samples_filtered --chr Y --keep DATA/step4_problematic_list.txt --missing --out DATA/step4_Y_missing

# If the annotated sex is female, and the missing rate of the Y chromosome is low, indicate the problem.
Rscript PROGRAM/sexcheck_chrY.r DATA/step4_Y_missing.smiss DATA/step4_sex_discrepancy_new.txt step4 # check the chromosome Y missing rate for this problematic sample

## ---- Step 4.3 Intensities of chr X and Y ----
Rscript PROGRAM/XY_Intensities.r /project/cshu0237_1136/Data/SPARK/iWES_v2/genotypes/genotyping_array/idat/ /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/release_notes/GSA-24v2-0_A2.csv /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA/step4_problematic_list.txt step4

# ----- Step 5: Population Structure -----
## Since the heterozygosity and HWE are sensitive to population structure, so we do this first.
## ---- Step 5.1: Downloading the reference data - 1000G phase 3 from plink2 website ----
refdir=Reference_data
cd $refdir 

# pgen=https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1
# pvar=https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1 
# sample=https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1

wget https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1 
mv 'all_hg38.pgen.zst?dl=1' all_phase3.pgen.zst 
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen 

wget https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1
mv 'all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz' all_phase3.pvar.zst 

wget https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1
mv 'hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x' all_phase3.psam

plink2 --pfile all_phase3 vzs\
       --allow-extra-chr \
       --max-alleles 2 \
       --make-bed \
       --out all_phase3 

cd .. # go back to the root dir

## ---- Step 5.2: Merge our data with reference data ----
### Extract 1-22 chr first  
plink2 --bfile DATA/step3_samples_filtered --chr 1-22 --make-bed --out DATA/step5_1-22
### Purning
plink2 --bfile DATA/step5_1-22 --indep-pairwise 50 5 0.2 --out DATA/step5_1-22
### Remove inds from pruned dataset
plink2 --bfile DATA/step5_1-22 --extract DATA/step5_1-22.prune.in --make-bed --out DATA/step5_purned
### Filter reference data for the same SNP set as in study
plink2 --bfile Reference_data/all_phase3 --allow-extra-chr --extract DATA/step5_1-22.prune.in --make-bed --out DATA/1000G_pruned

### Check duplicate
cut -f2 DATA/1000G_pruned.bim | sort | uniq -d > DATA/step5_dup_snps.txt
plink1.9 --bfile DATA/1000G_pruned --allow-extra-chr --exclude DATA/step5_dup_snps.txt --make-bed --out DATA/1000G_unique

### Check and correct chromosome mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} ($2 in a && a[$2] != $1) {print a[$2],$2}' DATA/step5_purned.bim DATA/1000G_unique.bim | sed -n '/^[XY]/!p' > DATA/1000G_toUpdateChr 

plink1.9 --bfile DATA/1000G_unique \
         --allow-extra-chr \
         --update-chr DATA/1000G_toUpdateChr 1 2 \
         --make-bed \
         --out DATA/1000G_updateChr 

### Position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4) {print a[$2],$2}' \
     DATA/step5_purned.bim DATA/1000G_unique.bim > \
     DATA/1000G_toUpdatePos

### Possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
     DATA/step5_purned.bim DATA/1000G_unique.bim > \
     DATA/1000G_toFlip

### Update positions and flip alleles
plink1.9 --bfile DATA/1000G_updateChr \
         --update-map DATA/1000G_toUpdatePos 1 2 \
         --flip DATA/1000G_toFlip \
         --make-bed \
         --out DATA/1000G_flipped

### Remove mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}'  DATA/step5_purned.bim DATA/1000G_flipped.bim > DATA/1000G_mismatch 
     
plink1.9 --bfile DATA/1000G_flipped \
         --exclude DATA/1000G_mismatch \
         --make-bed \
         --out DATA/1000G_finalcleaned

### Merged study genotypes and reference data
plink1.9 --bfile DATA/step5_purned.bim \
         --bmerge DATA/1000G_finalcleaned.bed DATA/1000G_finalcleaned.bim DATA/1000G_finalcleaned.fam \
         --make-bed \
         --out DATA/step5_merge_1000G

## ---- Step 5.3: PCA ----
plink2 --bfile DATA/step5_merge_1000G --threads 20 --pca approx --out DATA/step5_pca

Rscripts PROGRAM/plot_pca.r DATA/step5_merge_1000G.fam Reference_data/all_phase3.psam DATA/step5_pca.eigenvec step5

## ---- Step 5.4: Supervised population inference ----
### 1. Download the neural-admixture
pip install neural-admixture 

### 2. Prepare data
### When using 1 A40 GPU and 20 CPUs, our server can only work with less than 10000 samples
### Considering we will use 1000 Genomes phase 3 as reference (3202 samples), we need to separate our data by 5000 samples
awk '{print $1, $2}' DATA/step5_purned.fam > DATA/step5_sample_ids.txt
split -l 5000 DATA/step5_sample_ids.txt DATA/step5_sample_ids_

### Define the list of sample suffixes
samples=(aa ab ac ad ae af ag ah ai)

### Initialize counter for output file names
counter=1

### Loop through each sample suffix
for sample in "${samples[@]}"; do
    plink1.9 --bfile DATA/step5_purned --keep DATA/step5_sample_ids_"${sample}" --make-bed --out DATA/step5_subset"${counter}"
    counter=$((counter + 1))
done

### Merge our data with reference
for i in {1..9}; do
    plink1.9 --bfile DATA/step5_subset${i} \
             --bmerge DATA/1000G_finalcleaned.bed DATA/1000G_finalcleaned.bim DATA/1000G_finalcleaned.fam \
             --make-bed \
             --out DATA/step5_merge_subset${i}
done

### Generate population list file for supervised neural admixture
#### This file is a one column, no header, charater only file
#### Each row refers to the ancestry for corresponding sample 
#### For the reference samples, it should be the pure population they belong to
#### For our samples, it should be '-', or blank line, which indicates that the ancetsry need to be estimated.
#### Deatiled information in Admixture-manual-guide
for i in {1..9}; do
    PROGRAM/prepare_pop.r DATA/step5_merge_subset${i}.fam step5_merge_subset${i}
done

## 3. Request for GPU via discovery nodes
salloc --partition=gpu --ntasks=1 --gpus-per-task=a40:2 --mem=998G --cpus-per-task=64

module load cuda

## 4. Run Supervised Neural Admixture 
for i in {1..9}; do
CUDA_VISIBLE_DEVICES=0 neural-admixture train --k 5 \
                                              --supervised --populations_path  DATA/step5_merge_subset${i}.pop \
                                              --data_path  DATA/step5_merge_subset${i}.bed \
                                              --save_dir DATA  \
                                              --name  DATA/step5_merge_subset${i}\
                                              --seed 42 \
                                              --num_gpus 2 --num_cpus 64
done

## 5. Population Inference
Rscript PROGRAM/population_inference.r step5_merge_subset step5_

# ----- step 6: Sample heterozygosity by ethinity -----
## ---- Step 6.1: 
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --het --out DATA/step6_het
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --keep DATA/step5_AFR.txt --het --out DATA/step6_AFR_het
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --keep DATA/step5_EUR.txt --het --out DATA/step6_EUR_het
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --keep DATA/step5_AMR.txt --het --out DATA/step6_AMR_het
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --keep DATA/step5_EAS.txt --het --out DATA/step6_EAS_het
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --keep DATA/step5_SAS.txt --het --out DATA/step6_SAS_het

Rscript PROGRAM/check_heterozygosity.r DATA/step6_het.het step6_all
Rscript PROGRAM/check_heterozygosity.r DATA/step6_AFR_het.het step6_AFR
Rscript PROGRAM/check_heterozygosity.r DATA/step6_EUR_het.het step6_EUR
Rscript PROGRAM/check_heterozygosity.r DATA/step6_AMR_het.het step6_AMR
Rscript PROGRAM/check_heterozygosity.r DATA/step6_EAS_het.het step6_EAS
Rscript PROGRAM/check_heterozygosity.r DATA/step6_SAS_het.het step6_SAS

wc -l DATA/step6_all_heterozygosity_outliers.txt # 9
wc -l DATA/step6_AFR_heterozygosity_outliers.txt # 9
wc -l DATA/step6_EUR_heterozygosity_outliers.txt # 757
wc -l DATA/step6_AMR_heterozygosity_outliers.txt # 46
wc -l DATA/step6_EAS_heterozygosity_outliers.txt # 2
wc -l DATA/step6_SAS_heterozygosity_outliers.txt # 16

## ---- Step 6.2: Plotting heterozygosity rate with chromosome X and Y intensities ---
## ---- Step 6.3: Detecting long runs of homozygosity ----
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out DATA/step6_all_ROH
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --keep DATA/step5_AFR.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out DATA/step6_AFR_ROH
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --keep DATA/step5_EUR.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out DATA/step6_AMR_ROH
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --keep DATA/step5_AMR.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out DATA/step6_EAS_ROH
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --keep DATA/step5_EAS.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out DATA/step6_EAS_ROH
plink1.9 --bfile DATA/step3_samples_filtered --chr 1-22 --keep DATA/step5_SAS.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out DATA/step6_EAS_ROH

Rscript PROGRAM/plot_ROH.r DATA/step6_AFR_ROH.hom.indiv DATA/step6_AFR_het.het AFR
Rscript PROGRAM/plot_ROH.r DATA/step6_AMR_ROH.hom.indiv DATA/step6_AMR_het.het AMR
Rscript PROGRAM/plot_ROH.r DATA/step6_EAS_ROH.hom.indiv DATA/step6_EAS_het.het EAS
Rscript PROGRAM/plot_ROH.r DATA/step6_EUR_ROH.hom.indiv DATA/step6_EUR_het.het EUR
Rscript PROGRAM/plot_ROH.r DATA/step6_SAS_ROH.hom.indiv DATA/step6_SAS_het.het SAS

# Step7. Compute HWE in cases and controls by ancestry
plink1.9 --bfile step3_samples_filtered --filter-controls --keep step6_neural_ad/step6_AFR.txt --hardy --out step7_Controls_AFR_HWE
plink1.9 --bfile step3_samples_filtered --filter-controls --keep step6_neural_ad/step6_EUR.txt --hardy --out step7_Controls_EUR_HWE
plink1.9 --bfile step3_samples_filtered --filter-controls --keep step6_neural_ad/step6_AMR.txt --hardy --out step7_Controls_AMR_HWE
plink1.9 --bfile step3_samples_filtered --filter-controls --keep step6_neural_ad/step6_EAS.txt --hardy --out step7_Controls_EAS_HWE
plink1.9 --bfile step3_samples_filtered --filter-controls --keep step6_neural_ad/step6_SAS.txt --hardy --out step7_Controls_SAS_HWE

plink1.9 --bfile step3_samples_filtered --filter-cases --keep step6_neural_ad/step6_AFR.txt --hardy --out step7_Cases_AFR_HWE
plink1.9 --bfile step3_samples_filtered --filter-cases --keep step6_neural_ad/step6_EUR.txt --hardy --out step7_Cases_EUR_HWE
plink1.9 --bfile step3_samples_filtered --filter-cases --keep step6_neural_ad/step6_AMR.txt --hardy --out step7_Cases_AMR_HWE
plink1.9 --bfile step3_samples_filtered --filter-cases --keep step6_neural_ad/step6_EAS.txt --hardy --out step7_Cases_EAS_HWE
plink1.9 --bfile step3_samples_filtered --filter-cases --keep step6_neural_ad/step6_SAS.txt --hardy --out step7_Cases_SAS_HWE

## For controls, save the SNP list that P value < 1e-6
Rscript ../PROGRAM/plot_hwe.r step7_Controls_AFR_HWE.hwe AFR_Controls
Rscript ../PROGRAM/plot_hwe.r step7_Controls_EUR_HWE.hwe EUR_Controls
Rscript ../PROGRAM/plot_hwe.r step7_Controls_SAS_HWE.hwe SAS_Controls
Rscript ../PROGRAM/plot_hwe.r step7_Controls_AMR_HWE.hwe AMR_Controls
Rscript ../PROGRAM/plot_hwe.r step7_Controls_EAS_HWE.hwe EAS_Controls

## For cases, save the SNP list that P value < 1e-10
Rscript ../PROGRAM/plot_hwe.r step7_Cases_AFR_HWE.hwe AFR_Cases
Rscript ../PROGRAM/plot_hwe.r step7_Cases_EUR_HWE.hwe EUR_Cases
Rscript ../PROGRAM/plot_hwe.r step7_Cases_SAS_HWE.hwe SAS_Cases
Rscript ../PROGRAM/plot_hwe.r step7_Cases_AMR_HWE.hwe AMR_Cases
Rscript ../PROGRAM/plot_hwe.r step7_Cases_EAS_HWE.hwe EAS_Cases

wc -lc step7_AMR_Cases_HWE_list.txt # 978
wc -lc step7_AFR_Cases_HWE_list.txt # 1478
wc -lc step7_EUR_Cases_HWE_list.txt # 19682
wc -lc step7_SAS_Cases_HWE_list.txt # 158
wc -lc step7_EAS_Cases_HWE_list.txt # 894

wc -lc step7_AMR_Controls_HWE_list.txt # 2018
wc -lc step7_AFR_Controls_HWE_list.txt # 5040
wc -lc step7_EUR_Controls_HWE_list.txt # 30988
wc -lc step7_SAS_Controls_HWE_list.txt # 414
wc -lc step7_EAS_Controls_HWE_list.txt # 1878

# Step8. Compute MAF
plink1.9 --bfile step3_samples_filtered --freq --out step8_maf
awk 'NR==1 || $5 < 0.05' step8_maf.frq > step8_0.05.txt
wc -lc step8_0.05.txt # 335,658
awk 'NR==1 || $5 < 0.01' step8_maf.frq > step8_0.01.txt
wc -lc step8_0.01.txt # 127,314
