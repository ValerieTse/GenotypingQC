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

# plink2 is still at the beta phase, if you download plink1.9, use it instead.
# plink2 sometimes will have multi-thread issues, if it has an error, add `--thread 1` in your plink command.
# plink2 can do `-pca approx` but plink1.9 cannot, except this step and sex check (can only be done with plink1.9), others are the same.
# Carc has some issues with `module load plink1.9`, please download it by yourself.

module load plink2 

# For SPARK cohort, we have 3 datsets for the genotyping data (Array_v1, Array_v2, Sequencing)
mkdir Array_v1/{DATA,OUTPUTS} Array_v2/{DATA,OUTPUTS} Sequencing/{DATA,OUTPUTS} PROGRAM Reference_data Release_Notes # creating the directory

# Using Array_v2 as example
dir=Array_v2 
data_dir=/DATA
output_dir=/OUTPUTS
program_dir=PROGRAM
release_notes=Release_Notes
reference_data=Reference_data
raw_data=/project/cshu0237_1136/Data/SPARK/iWES_v3/genotypes/genotyping_array/plink/SPARK.iWES_v3.2024_08.GSA-24_v2 
# raw data = directory + prefix of the bfile.
meta_data=/project/cshu0237_1136/Data/SPARK/iWES_v3/metadata/SPARK.iWES_v3.2024_08.sample_metadata.tsv

# ----- Step 1: Update the .fam file -----
Rscript $program_dir/update_info_extract.r \
        $raw_data.fam \
        $meta_data

plink2 --bfile $raw_data --update-ids $release_notes/update_id.txt --make-bed --out $dir$data_dir/step1_update_ids
plink2 --bfile $dir$data_dir/step1_update_ids --update-sex $release_notes/update_sex.txt --make-bed --out $dir$data_dir/step1_update_sex 
plink2 --bfile $dir$data_dir/step1_update_sex --pheno $release_notes/update_pheno.txt --make-bed --out $dir$data_dir/step1_update_pheno

# Check the number of samples and SNPs
Rscript $program_dir/counts.r $dir$data_dir/step1_update_pheno 

# ----- Step 2: Initial check -----
plink2 --bfile $dir$data_dir/step1_update_pheno --missing --out $dir$data_dir/step2_overall_missing

Rscript $program_dir/plot_call_rate.r $dir$data_dir/step2_overall_missing.smiss $dir$data_dir/step2_overall_missing.vmiss step2 $dir

plink2 --bfile $dir$data_dir/step1_update_pheno --chr Y --missing --out $dir$data_dir/step2_chrY_missing
plink2 --bfile $dir$data_dir/step1_update_pheno --chr X --missing --out $dir$data_dir/step2_chrX_missing
plink2 --bfile $dir$data_dir/step1_update_pheno --chr MT --missing --out $dir$data_dir/step2_chrMT_missing

# print the maximum missing rate for samples and SNPs
Rscript $program_dir/call_rate.r $dir$data_dir/step2_chrY_missing.smiss $dir$data_dir/step2_chrY_missing.vmiss 
Rscript $program_dir/call_rate.r $dir$data_dir/step2_chrX_missing.smiss $dir$data_dir/step2_chrX_missing.vmiss 
Rscript $program_dir/call_rate.r $dir$data_dir/step2_chrMT_missing.smiss $dir$data_dir/step2_chrMT_missing.vmiss

# ----- Step 3: Missing call rate filtering -----
# Based on the results (hists) from step 2 to decide the threshold for missingness filtering
plink2 --bfile $dir$data_dir/step1_update_pheno --mind 0.1 --geno 0.1 --make-bed --out $dir$data_dir/step3_filtered

# Double check
plink2 --bfile $dir$data_dir/step3_filtered --missing --out $dir$data_dir/step3_overall_missing
Rscript $program_dir/plot_call_rate.r $dir$data_dir/step3_overall_missing.smiss $dir$data_dir/step3_overall_missing.vmiss step3 $dir
Rscript $program_dir/counts.r $dir$data_dir/step3_filtered

# ----- Step 4: Sex check -----
## ---- Step 4.1: Chr X heterozygosity ----
## Using `--check-sex` check for Chr X heterozygosity and detect the problematic samples.
plink --bfile $dir$data_dir/step3_filtered --check-sex --out $dir$data_dir/step4_sexcheck # Hint: Using plink rather than plink2
grep 'PROBLEM' $dir$data_dir/step4_sexcheck.sexcheck > $dir$data_dir/step4_sex_discrepancy.txt 
wc -l $dir$data_dir/step4_sex_discrepancy.txt  # count the problematic sample based on threshold 0.2 and 0.8

Rscript $program_dir/plot_sexcheck.r $dir$data_dir/step4_sexcheck.sexcheck step4 $dir

## Based on previous results, re-select the threshold for genetic sex inference
plink --bfile $dir$data_dir/step3_filtered --check-sex 0.4 0.6 --out $dir$data_dir/step4_sexcheck_new 
grep 'PROBLEM' $dir$data_dir/step4_sexcheck_new.sexcheck > $dir$data_dir/step4_sex_discrepancy_new.txt
wc -l $dir$data_dir/step4_sex_discrepancy_new.txt # count the problematic sample based on threshold 0.4 and 0.6
awk '{print $1, $2}' $dir$data_dir/step4_sex_discrepancy_new.txt > $dir$data_dir/step4_problematic_list.txt

## ---- Step 4.2: chr Y missingness ----
plink2 --bfile $dir$data_dir/step3_filtered --chr Y --keep $dir$data_dir/step4_problematic_list.txt --missing --out $dir$data_dir/step4_Y_missing

# If the annotated sex is female, and the missing rate of the Y chromosome is low, indicate the problem.
Rscript $program_dir/sexcheck_chrY.r \
        $dir$data_dir/step4_Y_missing.smiss \
        $dir$data_dir/step4_sex_discrepancy_new.txt \
        step4 \
        $dir # check the chromosome Y missing rate for this problematic sample

## ---- Step 4.3 Intensities of chr X and Y ----
Rscript $program_dir/XY_Intensities.r \
        /project/cshu0237_1136/Data/SPARK/iWES_v2/genotypes/genotyping_array/idat/ \
        $release_notes/GSA-24v2-0_A2.csv \
        $dir$data_dir/step4_problematic_list.txt \
        $dir$data_dir/step4_sexcheck_new.sexcheck
        step4 $dir

## ---- Step 4.4 Chr X heterozygosity rate ----
## `--check-sex` will only output the F-stat, recode chr X -> chr 1, and compute the heterozygosity rate for chr X only.
plink --bfile $dir$data_dir/step3_filtered --chr X --make-bed --out $dir$data_dir/step4_chrX

awk '{if ($1 == 23) $1=1; print}' $dir$data_dir/step4_chrX.bim > $dir$data_dir/step4_chrX.bim.tmp
mv $dir$data_dir/step4_chrX.bim.tmp $dir$data_dir/step4_chrX.bim

plink --bfile $dir$data_dir/step4_chrX --chr 1 --het --out $dir$data_dir/step4_chrx_het
Rscript $program_dir/het_rates.r $dir$data_dir/step4_chrx_het.het step4_chrx $dir

# ----- Step 5: Population Structure -----
## Since the heterozygosity and HWE are sensitive to population structure, we do this first.
## ---- Step 5.1: Downloading the reference data - 1000G phase 3 from the plink2 website ----
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
plink2 --bfile $dir$data_dir/step3_filtered --chr 1-22 --make-bed --out $dir$data_dir/step5_1-22
### Purning
plink2 --bfile $dir$data_dir/step5_1-22 --indep-pairwise 50 5 0.2 --out $dir$data_dir/step5_1-22
### Remove inds from pruned dataset
plink2 --bfile $dir$data_dir/step5_1-22 --extract $dir$data_dir/step5_1-22.prune.in --make-bed --out $dir$data_dir/step5_purned
### Filter reference data for the same SNP set as in study
plink2 --bfile $reference_data/all_phase3 --allow-extra-chr --extract $dir$data_dir/step5_1-22.prune.in --make-bed --out $dir$data_dir/1000G_pruned

### Check duplicate
cut -f2 $dir$data_dir/1000G_pruned.bim | sort | uniq -d > $dir$data_dir/step5_dup_snps.txt
plink --bfile $dir$data_dir/1000G_pruned --allow-extra-chr --exclude $dir$data_dir/step5_dup_snps.txt --make-bed --out $dir$data_dir/1000G_unique

### Check and correct chromosome mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} ($2 in a && a[$2] != $1) {print a[$2],$2}' $dir$data_dir/step5_purned.bim $dir$data_dir/1000G_unique.bim | sed -n '/^[XY]/!p' > $dir$data_dir/1000G_toUpdateChr 

plink --bfile $dir$data_dir/1000G_unique \
      --allow-extra-chr \
      --update-chr $dir$data_dir/1000G_toUpdateChr 1 2 \
      --make-bed \
      --out $dir$data_dir/1000G_updateChr 

### Position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4) {print a[$2],$2}' \
     $dir$data_dir/step5_purned.bim $dir$data_dir/1000G_unique.bim > \
     $dir$data_dir/1000G_toUpdatePos

### Possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
     $dir$data_dir/step5_purned.bim $dir$data_dir/1000G_unique.bim > \
     $dir$data_dir/1000G_toFlip

### Update positions and flip alleles
plink --bfile $dir$data_dir/1000G_updateChr \
      --update-map $dir$data_dir/1000G_toUpdatePos 1 2 \
      --flip $dir$data_dir/1000G_toFlip \
      --make-bed \
      --out $dir$data_dir/1000G_flipped

### Remove mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}'  $dir$data_dir/step5_purned.bim $dir$data_dir/1000G_flipped.bim > $dir$data_dir/1000G_mismatch 
     
plink --bfile $dir$data_dir/1000G_flipped \
      --exclude $dir$data_dir/1000G_mismatch \
      --make-bed \
      --out $dir$data_dir/1000G_finalcleaned

### Merged study genotypes and reference data
plink --bfile $dir$data_dir/step5_purned.bim \
      --bmerge $dir$data_dir/1000G_finalcleaned.bed $dir$data_dir/1000G_finalcleaned.bim $dir$data_dir/1000G_finalcleaned.fam \
      --make-bed \
      --out $dir$data_dir/step5_merge_1000G

## ---- Step 5.3: PCA ----
plink2 --bfile $dir$data_dir/step5_merge_1000G --threads 20 --pca approx --out $dir$data_dir/step5_pca

Rscripts $program_dir/plot_pca.r \
         $dir$data_dir/step5_merge_1000G.fam \
         $reference_data/all_phase3.psam \
         $dir$data_dir/step5_pca.eigenvec \
         step5 \
         $dir

## ---- Step 5.4: Supervised population inference ----
### 1. Download the neural-admixture
pip install neural-admixture 

### 2. Prepare data
### When using 1 A40 GPU and 20 CPUs, our server can only work with less than 10000 samples
### Considering we will use 1000 Genomes phase 3 as reference (3202 samples), we need to separate our data by 5000 samples
awk '{print $1, $2}' $dir$data_dir/step5_purned.fam > $dir$data_dir/step5_sample_ids.txt
split -l 5000 $dir$data_dir/step5_sample_ids.txt $dir$data_dir/step5_sample_ids_

### Define the list of sample suffixes
samples=(aa ab ac ad ae af ag ah ai)

### Initialize counter for output file names
counter=1

### Loop through each sample suffix
for sample in "${samples[@]}"; do
    plink --bfile $dir$data_dir/step5_purned --keep $dir$data_dir/step5_sample_ids_"${sample}" --make-bed --out $dir$data_dir/step5_subset"${counter}"
    counter=$((counter + 1))
done

### Merge our data with reference
for i in {1..9}; do
    plink --bfile $dir$data_dir/step5_subset${i} \
          --bmerge $dir$data_dir/1000G_finalcleaned.bed $dir$data_dir/1000G_finalcleaned.bim $dir$data_dir/1000G_finalcleaned.fam \
          --make-bed \
          --out $dir$data_dir/step5_merge_subset${i}
done

### Generate population list file for supervised neural admixture
#### This file is a one column, no header, charater only file
#### Each row refers to the ancestry for corresponding sample 
#### For the reference samples, it should be the pure population they belong to
#### For our samples, it should be '-', or blank line, which indicates that the ancetsry need to be estimated.
#### Deatiled information in Admixture-manual-guide
for i in {1..9}; do
    $program_dir/prepare_pop.r $dir$data_dir/step5_merge_subset${i}.fam step5_merge_subset${i} $dir
done

## 3. Request for GPU via discovery nodes
salloc --account=cshu0237_1136 --partition=gpu --ntasks=1 --gpus-per-task=a40:2 --cpus-per-task=32 --mem=248G # salloc in discovery can only require 1 hour

module load cuda

## 4. Run Supervised Neural Admixture 
for i in {1..9}; do
CUDA_VISIBLE_DEVICES=0 neural-admixture train --k 5 \
                                              --supervised --populations_path  $dir$data_dir/step5_merge_subset${i}.pop \
                                              --data_path  $dir$data_dir/step5_merge_subset${i}.bed \
                                              --save_dir DATA  \
                                              --name  $dir$data_dir/step5_merge_subset${i}\
                                              --seed 42 \
                                              --num_gpus 2 --num_cpus 64
done

## 5. Population Inference
Rscript $program_dir/population_inference.r step5_merge_subset step5_ $dir

# ----- step 6: Sample heterozygosity by ethinity -----
plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --het --out $dir$data_dir/step6_het
Rscript $program_dir/het_rates.r $dir$data_dir/step6_het.het step6 $dir
Rscript $program_dir/corrected_het.r $dir$data_dir/step6_het.het  /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA/PROCESSING/step6/step6_pca.eigenvec step6 $dir       

plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --keep $dir$data_dir/step5_AFR.txt --het --out $dir$data_dir/step6_AFR_het
plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --keep $dir$data_dir/step5_EUR.txt --het --out $dir$data_dir/step6_EUR_het
plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --keep $dir$data_dir/step5_AMR.txt --het --out $dir$data_dir/step6_AMR_het
plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --keep $dir$data_dir/step5_EAS.txt --het --out $dir$data_dir/step6_EAS_het
plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --keep $dir$data_dir/step5_SAS.txt --het --out $dir$data_dir/step6_SAS_het

Rscript $program_dir/check_heterozygosity.r $dir$data_dir/step6_het.het step6_all
Rscript $program_dir/check_heterozygosity.r $dir$data_dir/step6_AFR_het.het step6_AFR
Rscript $program_dir/check_heterozygosity.r $dir$data_dir/step6_EUR_het.het step6_EUR
Rscript $program_dir/check_heterozygosity.r $dir$data_dir/step6_AMR_het.het step6_AMR
Rscript $program_dir/check_heterozygosity.r $dir$data_dir/step6_EAS_het.het step6_EAS
Rscript $program_dir/check_heterozygosity.r $dir$data_dir/step6_SAS_het.het step6_SAS

wc -l $dir$data_dir/step6_all_heterozygosity_outliers.txt 
wc -l $dir$data_dir/step6_AFR_heterozygosity_outliers.txt 
wc -l $dir$data_dir/step6_EUR_heterozygosity_outliers.txt 
wc -l $dir$data_dir/step6_AMR_heterozygosity_outliers.txt 
wc -l $dir$data_dir/step6_EAS_heterozygosity_outliers.txt 
wc -l $dir$data_dir/step6_SAS_heterozygosity_outliers.txt 

## ---- Step 6.2: Plotting heterozygosity rate with chromosome X and Y intensities ---
## ---- Step 6.3: Detecting long runs of homozygosity ----
plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out $dir$data_dir/step6_all_ROH
plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --keep $dir$data_dir/step5_AFR.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out $dir$data_dir/step6_AFR_ROH
plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --keep $dir$data_dir/step5_EUR.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out $dir$data_dir/step6_AMR_ROH
plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --keep $dir$data_dir/step5_AMR.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out $dir$data_dir/step6_EAS_ROH
plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --keep $dir$data_dir/step5_EAS.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out $dir$data_dir/step6_EAS_ROH
plink --bfile $dir$data_dir/step3_filtered --chr 1-22 --keep $dir$data_dir/step5_SAS.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out $dir$data_dir/step6_EAS_ROH

Rscript $program_dir/plot_ROH.r $dir$data_dir/step6_AFR_ROH.hom.indiv $dir$data_dir/step6_AFR_het.het AFR
Rscript $program_dir/plot_ROH.r $dir$data_dir/step6_AMR_ROH.hom.indiv $dir$data_dir/step6_AMR_het.het AMR
Rscript $program_dir/plot_ROH.r $dir$data_dir/step6_EAS_ROH.hom.indiv $dir$data_dir/step6_EAS_het.het EAS
Rscript $program_dir/plot_ROH.r $dir$data_dir/step6_EUR_ROH.hom.indiv $dir$data_dir/step6_EUR_het.het EUR
Rscript $program_dir/plot_ROH.r $dir$data_dir/step6_SAS_ROH.hom.indiv $dir$data_dir/step6_SAS_het.het SAS

# ----- Step 7: Compute HWE in cases and controls by ancestry -----
plink --bfile $dir$data_dir/step3_filtered --filter-controls --keep $dir$data_dir/step5_AFR.txt --hardy --out $dir$data_dir/step7_Controls_AFR_HWE
plink --bfile $dir$data_dir/step3_filtered --filter-controls --keep $dir$data_dir/step5_EUR.txt --hardy --out $dir$data_dir/step7_Controls_EUR_HWE
plink --bfile $dir$data_dir/step3_filtered --filter-controls --keep $dir$data_dir/step5_AMR.txt --hardy --out $dir$data_dir/step7_Controls_AMR_HWE
plink --bfile $dir$data_dir/step3_filtered --filter-controls --keep $dir$data_dir/step5_EAS.txt --hardy --out $dir$data_dir/step7_Controls_EAS_HWE
plink --bfile $dir$data_dir/step3_filtered --filter-controls --keep $dir$data_dir/step5_SAS.txt --hardy --out $dir$data_dir/step7_Controls_SAS_HWE

plink --bfile $dir$data_dir/step3_filtered --filter-cases --keep $dir$data_dir/step5_AFR.txt --hardy --out $dir$data_dir/step7_Cases_AFR_HWE
plink --bfile $dir$data_dir/step3_filtered --filter-cases --keep $dir$data_dir/step5_EUR.txt --hardy --out $dir$data_dir/step7_Cases_EUR_HWE
plink --bfile $dir$data_dir/step3_filtered --filter-cases --keep $dir$data_dir/step5_AMR.txt --hardy --out $dir$data_dir/step7_Cases_AMR_HWE
plink --bfile $dir$data_dir/step3_filtered --filter-cases --keep $dir$data_dir/step5_EAS.txt --hardy --out $dir$data_dir/step7_Cases_EAS_HWE
plink --bfile $dir$data_dir/step3_filtered --filter-cases --keep $dir$data_dir/step5_SAS.txt --hardy --out $dir$data_dir/step7_Cases_SAS_HWE

## For controls, save the SNP list that P value < 1e-6
Rscript $program_dir/plot_hwe.r step7_Controls_AFR_HWE.hwe AFR_Controls
Rscript $program_dir/plot_hwe.r step7_Controls_EUR_HWE.hwe EUR_Controls
Rscript $program_dir/plot_hwe.r step7_Controls_SAS_HWE.hwe SAS_Controls
Rscript $program_dir/plot_hwe.r step7_Controls_AMR_HWE.hwe AMR_Controls
Rscript $program_dir/plot_hwe.r step7_Controls_EAS_HWE.hwe EAS_Controls

## For cases, save the SNP list that P value < 1e-10
Rscript $program_dir/plot_hwe.r $dir$data_dir/step7_Cases_AFR_HWE.hwe AFR_Cases
Rscript $program_dir/plot_hwe.r $dir$data_dir/step7_Cases_EUR_HWE.hwe EUR_Cases
Rscript $program_dir/plot_hwe.r $dir$data_dir/step7_Cases_SAS_HWE.hwe SAS_Cases
Rscript $program_dir/plot_hwe.r $dir$data_dir/step7_Cases_AMR_HWE.hwe AMR_Cases
Rscript $program_dir/plot_hwe.r $dir$data_dir/step7_Cases_EAS_HWE.hwe EAS_Cases

wc -l $dir$data_dir/step7_AMR_Cases_HWE_list.txt # 978
wc -l $dir$data_dir/step7_AFR_Cases_HWE_list.txt # 1478
wc -l $dir$data_dir/step7_EUR_Cases_HWE_list.txt # 19682
wc -l $dir$data_dir/step7_SAS_Cases_HWE_list.txt # 158
wc -l $dir$data_dir/step7_EAS_Cases_HWE_list.txt # 894

wc -l $dir$data_dir/step7_AMR_Controls_HWE_list.txt # 2018
wc -l $dir$data_dir/step7_AFR_Controls_HWE_list.txt # 5040
wc -l $dir$data_dir/step7_EUR_Controls_HWE_list.txt # 30988
wc -l $dir$data_dir/step7_SAS_Controls_HWE_list.txt # 414
wc -l $dir$data_dir/step7_EAS_Controls_HWE_list.txt # 1878

# ----- Step 8: Compute MAF -----
plink --bfile $dir$data_dir/step3_filtered --freq --out $dir$data_dir/step8_maf
awk 'NR==1 || $5 < 0.05' $dir$data_dir/step8_maf.frq > $dir$data_dir/step8_0.05.txt
wc -l $dir$data_dir/step8_0.05.txt # 335,658
awk 'NR==1 || $5 < 0.01' $dir$data_dir/step8_maf.frq > $dir$data_dir/step8_0.01.txt
wc -l $dir$data_dir/step8_0.01.txt # 127,314
