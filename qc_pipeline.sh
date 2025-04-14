# Step1: Update the .fam file
Rscript update_info_extract.r

plink2 --bfile /project/cshu0237_1136/Data/SPARK/iWES_v3/genotypes/genotyping_array/plink/SPARK.iWES_v3.2024_08.GSA-24_v2 --update-ids release_notes/update_ids.txt --make-bed --out DATA/step1_update_ids
plink2 --bfile DATA/step1_update_ids --update-sex release_notes/update_sex.txt --make-bed --out DATA/step1_update_sex 
plink2 --bfile DATA/step1_update_sex --pheno release_notes/pheno.txt --make-bed --out DATA/step1_update_pheno

cd DATA

Rscript ../PROGRAM/counts.r step1_update_pheno

# Step2: Initial check
plink2 --bfile step1_update_pheno --missing --out step2_overall_missing

Rscript ../PROGRAM/plot_call_rate.r step2_overall_missing.smiss step2_overall_missing.vmiss step2

plink2 --bfile step1_update_pheno --chr Y --missing --out step2_chrY_missing
plink2 --bfile step1_update_pheno --chr X --missing --out step2_chrX_missing
plink2 --bfile step1_update_pheno --chr MT --missing --out step2_chrMT_missing

Rscript ../PROGRAM/call_rate.r step2_chrY_missing.smiss step2_chrY_missing.vmiss step2
Rscript ../PROGRAM/call_rate.r step2_chrX_missing.smiss step2_chrX_missing.vmiss step2
Rscript ../PROGRAM/call_rate.r step2_chrMT_missing.smiss step2_overall_missing.vmiss step2

# Step3: Missing call rate filtering
plink2 --bfile step1_update_pheno --mind 0.1 --make-bed --out step3_samples_filtered
plink2 --bfile step3_samples_filtered --geno 0.1 --make-bed --out step3_variants_filtered

plink2 --bfile step3_variants_filtered --missing --out step3_overall_missing

Rscript ../PROGRAM/plot_call_rate.r step3_overall_missing.smiss step3_overall_missing.vmiss step3

Rscript ../PROGRAM/counts.r step3_variants_filtered

# Step4: Sex check
plink1.9  --bfile step3_samples_filtered --check-sex --out step4_sexcheck

Rscript ../PROGRAM/plot_sexcheck.r step4_sexcheck.sexcheck step4

plink1.9 --bfile step3_samples_filtered --check-sex 0.4 0.6 --out step4_sexcheck_new

grep 'PROBLEM' step4_sexcheck.sexcheck > step4_sex_discrepancy.txt
grep 'PROBLEM' step4_sexcheck_new.sexcheck > step4_sex_discrepancy_new.txt

awk '{print $1, $2}' step4_sex_discrepancy_new.txt > step4_problematic_list.txt

plink2 --bfile step3_samples_filtered --chr Y --keep step4_problematic_list.txt --missing --out step4_Y_missing
Rscript ../PROGRAM/sexcheck_chrY.r step4_Y_missing.smiss step4_sex_discrepancy_new.txt step4

plink2 --bfile step3_samples_filtered --keep step4_problematic_list.txt --make-bed --out step4_problematics

conda activate latest # contained the illuminaio package

# Step5: Sample heterozygosity
plink2 --bfile step3_samples_filtered --het --out step5_samples_het
plink2 --bfile step3_samples_filtered --chr 1-22 --het --out step5_samples_het_1-22

Rscript ../PROGRAM/check_heterozygosity.r step5_samples_het.het step5_all
Rscript ../PROGRAM/check_heterozygosity.r step5_samples_het_1-22.het step5_1-22

wc -lc step5_all_heterozygosity_outliers.txt # 989 
wc -lc step5_1-22_heterozygosity_outliers.txt # 989

wc -lc step5_1-22_heterozygosity_outliers.txt # 989



# Step6: Population Structure
## Download and decompress 1000 Genomes phase 3 data
refdir=1000G
mkdir $refdir
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

cd ..

## extract 1-22 chr first  
plink2 --bfile step3_samples_filtered --chr 1-22 --make-bed --out step6_1-22
## use pruned dataset 
plink2 --bfile step6_1-22 --indep-pairwise 50 5 0.2 --out step6_1-22_pruned
## remove inds from pruned dataset
plink2 --bfile step6_1-22 --extract step6_1-22_pruned.prune.in --make-bed --out step6_forpopcheck
## Filter reference data for the same SNP set as in study
plink2 --bfile 1000G/all_phase3 --allow-extra-chr --extract step6_1-22_pruned.prune.in --make-bed --out 1000G.pruned

## Check duplicate
cut -f2 1000G.pruned.bim | sort | uniq -d > dup_snps.txt
plink1.9 --bfile 1000G.pruned --allow-extra-chr --exclude dup_snps.txt --make-bed --out 1000G.cleaned

## Check and correct chromosome mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} ($2 in a && a[$2] != $1) {print a[$2],$2}' step6_forpopcheck.bim 1000G.cleaned.bim | sed -n '/^[XY]/!p' > 1000G.toUpdateChr 

plink1.9 --bfile 1000G.cleaned \
         --allow-extra-chr \
         --update-chr 1000G.toUpdateChr 1 2 \
         --make-bed \
         --out 1000G.updateChr 

## Position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4) {print a[$2],$2}' \
     step6_forpopcheck.bim 1000G.cleaned.bim > \
     1000G.toUpdatePos

## Possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
     step6_forpopcheck.bim 1000G.cleaned.bim > \
      1000G.toFlip

## Upate positions and flip alleles
plink1.9 --bfile 1000G.updateChr \
         --update-map 1000G.toUpdatePos 1 2 \
         --flip 1000G.toFlip \
         --make-bed \
         --out 1000G.flipped

## Remove mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}'  step6_forpopcheck.bim 1000G.flipped.bim > 1000G.mismatch 
     
plink1.9 --bfile 1000G.flipped \
         --exclude 1000G.mismatch \
         --make-bed \
         --out 1000G.finalclean

## Merged study genotypes and reference data
plink1.9 --bfile step6_forpopcheck \
         --bmerge 1000G.finalclean.bed 1000G.finalclean.bim 1000G.finalclean.fam \
         --make-bed \
         --out step6.merge.1000G

## PCA
plink2 --bfile step6.merge.1000G --threads 20 --pca approx --out step6_pca

## Ancestry inference by Rye
### https://github.com/healthdisparities/rye/tree/main
git clone https://github.com/healthdisparities/rye

Rscript -e 'install.packages(c('nnls','Hmisc','parallel', 'optparse', 'crayon'))'

### Download the repository
### Ensure permissions are correct
cd 
cd Software
git clone https://github.com/healthdisparities/rye

cd rye
### the uploader (me) wasn't smart enough to change the permissions the first time and refuses to reupload the file with right permissions
chmod +x rye.R

### Check if the install happened correctly
./rye.R -h

### Perform 
cd /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA
#### Required data are generate from ../PROGRAM/plot_pca_gendat_rye.r
#### Using 1000G as reference, there are 5 real ancestries (AFR, SAS, EAS, AMR, EUR)
/home1/jingyux/Software/rye-main/rye.R --eigenvec=step6_rye.eigenvec --eigenval=step6_pca.eigenval --pop2group=step6_rye_pop.txt --pcs=10 --out=step6_out 

### Using ../PROGRAM/ancestry_infer.r to infer ancestry and create ancetsry id list
#### step6_AFR.txt, etc.
wc -lc step6_AFR.txt # 28,858
wc -lc step6_EUR.txt # 3,950
wc -lc step6_AMR.txt # 9,439
wc -lc step6_SAS.txt # 695
wc -lc step6_EAS.txt # 0
wc -lc step6_forpopcheck.fam # 42,942

## Optional 
## Ancestry Inference by Admixture
### Skip this step due to resource limited
### Find the best K
for K in $(seq 3 10); do admixture --cv step6_forpopcheck.bed $K -j18 | tee step6_log${K}.out; done 
grep -h CV step6_log*.out

### run admixture
admixture step6_forpopcheck.bed k

## Ancestry Inference by neural-admixture
### This should run on discovery node
conda create -n nadmenv python=3.9
conda activate nadmenv

pip install neural-admixture

awk '{print $1, $2}' step6_forpopcheck.fam > step6_sample_ids.txt
split -l 5000 step6_sample_ids.txt step6_neural_ad/step6_sample_ids_
plink1.9 --bfile step6_forpopcheck --keep step6_neural_ad/step6_sample_ids_aa --make-bed --out step6_neural_ad/step6_subset1
plink1.9 --bfile step6_forpopcheck --keep step6_neural_ad/step6_sample_ids_ab --make-bed --out step6_neural_ad/step6_subset2
plink1.9 --bfile step6_forpopcheck --keep step6_neural_ad/step6_sample_ids_ac --make-bed --out step6_neural_ad/step6_subset3
plink1.9 --bfile step6_forpopcheck --keep step6_neural_ad/step6_sample_ids_ad --make-bed --out step6_neural_ad/step6_subset4
plink1.9 --bfile step6_forpopcheck --keep step6_neural_ad/step6_sample_ids_ae --make-bed --out step6_neural_ad/step6_subset5
plink1.9 --bfile step6_forpopcheck --keep step6_neural_ad/step6_sample_ids_af --make-bed --out step6_neural_ad/step6_subset6
plink1.9 --bfile step6_forpopcheck --keep step6_neural_ad/step6_sample_ids_ag --make-bed --out step6_neural_ad/step6_subset7
plink1.9 --bfile step6_forpopcheck --keep step6_neural_ad/step6_sample_ids_ah --make-bed --out step6_neural_ad/step6_subset8
plink1.9 --bfile step6_forpopcheck --keep step6_neural_ad/step6_sample_ids_ai --make-bed --out step6_neural_ad/step6_subset9

plink1.9 --bfile step6_neural_ad/step6_subset1 \
         --bmerge 1000G.finalclean.bed 1000G.finalclean.bim 1000G.finalclean.fam \
         --make-bed \
         --out step6_neural_ad/step6_merge_subset1

plink1.9 --bfile step6_neural_ad/step6_subset2 \
         --bmerge 1000G.finalclean.bed 1000G.finalclean.bim 1000G.finalclean.fam \
         --make-bed \
         --out step6_neural_ad/step6_merge_subset2

plink1.9 --bfile step6_neural_ad/step6_subset3 \
         --bmerge 1000G.finalclean.bed 1000G.finalclean.bim 1000G.finalclean.fam \
         --make-bed \
         --out step6_neural_ad/step6_merge_subset3

plink1.9 --bfile step6_neural_ad/step6_subset4 \
         --bmerge 1000G.finalclean.bed 1000G.finalclean.bim 1000G.finalclean.fam \
         --make-bed \
         --out step6_neural_ad/step6_merge_subset4

plink1.9 --bfile step6_neural_ad/step6_subset5 \
         --bmerge 1000G.finalclean.bed 1000G.finalclean.bim 1000G.finalclean.fam \
         --make-bed \
         --out step6_neural_ad/step6_merge_subset5

plink1.9 --bfile step6_neural_ad/step6_subset6 \
         --bmerge 1000G.finalclean.bed 1000G.finalclean.bim 1000G.finalclean.fam \
         --make-bed \
         --out step6_neural_ad/step6_merge_subset6

plink1.9 --bfile step6_neural_ad/step6_subset7 \
         --bmerge 1000G.finalclean.bed 1000G.finalclean.bim 1000G.finalclean.fam \
         --make-bed \
         --out step6_neural_ad/step6_merge_subset7

plink1.9 --bfile step6_neural_ad/step6_subset8 \
         --bmerge 1000G.finalclean.bed 1000G.finalclean.bim 1000G.finalclean.fam \
         --make-bed \
         --out step6_neural_ad/step6_merge_subset8

plink1.9 --bfile step6_neural_ad/step6_subset9 \
         --bmerge 1000G.finalclean.bed 1000G.finalclean.bim 1000G.finalclean.fam \
         --make-bed \
         --out step6_neural_ad/step6_merge_subset9

cd step6_neural_ad
Rscript ../../PROGRAM/prepare_pop.r step6_merge_subset1.fam subset1
Rscript ../../PROGRAM/prepare_pop.r step6_merge_subset2.fam subset2
Rscript ../../PROGRAM/prepare_pop.r step6_merge_subset3.fam subset3
Rscript ../../PROGRAM/prepare_pop.r step6_merge_subset4.fam subset4
Rscript ../../PROGRAM/prepare_pop.r step6_merge_subset5.fam subset5
Rscript ../../PROGRAM/prepare_pop.r step6_merge_subset6.fam subset6
Rscript ../../PROGRAM/prepare_pop.r step6_merge_subset7.fam subset7
Rscript ../../PROGRAM/prepare_pop.r step6_merge_subset8.fam subset8
Rscript ../../PROGRAM/prepare_pop.r step6_merge_subset9.fam subset9

### Request for the GPU
salloc --partition=gpu --ntasks=1 --gpus-per-task=a40:2 --mem=64G --cpus-per-task=20

module load cuda

CUDA_VISIBLE_DEVICES=0 neural-admixture train --k 5 --supervised --populations_path step6_refpop_subset1.pop --data_path step6_merge_subset1.bed --save_dir /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA/step6_neural_ad  --name step6_merge_subset1 --seed 42 --num_gpus 2 --num_cpus 20
CUDA_VISIBLE_DEVICES=0 neural-admixture train --k 5 --supervised --populations_path step6_refpop_subset2.pop --data_path step6_merge_subset2.bed --save_dir /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA/step6_neural_ad  --name step6_merge_subset2 --seed 42 --num_gpus 2 --num_cpus 20
CUDA_VISIBLE_DEVICES=0 neural-admixture train --k 5 --supervised --populations_path step6_refpop_subset3.pop --data_path step6_merge_subset3.bed --save_dir /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA/step6_neural_ad  --name step6_merge_subset3 --seed 42 --num_gpus 2 --num_cpus 20
CUDA_VISIBLE_DEVICES=0 neural-admixture train --k 5 --supervised --populations_path step6_refpop_subset4.pop --data_path step6_merge_subset4.bed --save_dir /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA/step6_neural_ad  --name step6_merge_subset4 --seed 42 --num_gpus 2 --num_cpus 20
CUDA_VISIBLE_DEVICES=0 neural-admixture train --k 5 --supervised --populations_path step6_refpop_subset5.pop --data_path step6_merge_subset5.bed --save_dir /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA/step6_neural_ad  --name step6_merge_subset5 --seed 42 --num_gpus 2 --num_cpus 20
CUDA_VISIBLE_DEVICES=0 neural-admixture train --k 5 --supervised --populations_path step6_refpop_subset6.pop --data_path step6_merge_subset6.bed --save_dir /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA/step6_neural_ad  --name step6_merge_subset6 --seed 42 --num_gpus 2 --num_cpus 20
CUDA_VISIBLE_DEVICES=0 neural-admixture train --k 5 --supervised --populations_path step6_refpop_subset7.pop --data_path step6_merge_subset7.bed --save_dir /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA/step6_neural_ad  --name step6_merge_subset7 --seed 42 --num_gpus 2 --num_cpus 20
CUDA_VISIBLE_DEVICES=0 neural-admixture train --k 5 --supervised --populations_path step6_refpop_subset8.pop --data_path step6_merge_subset8.bed --save_dir /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA/step6_neural_ad  --name step6_merge_subset8 --seed 42 --num_gpus 2 --num_cpus 20
CUDA_VISIBLE_DEVICES=0 neural-admixture train --k 5 --supervised --populations_path step6_refpop_subset9.pop --data_path step6_merge_subset9.bed --save_dir /project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA/step6_neural_ad  --name step6_merge_subset9 --seed 42 --num_gpus 2 --num_cpus 20

# step5. additional sample heterozygosity by ethinity
plink1.9 --bfile step3_samples_filtered --keep step6_neural_ad/step6_AFR.txt --het --out step5_AFR_het
plink1.9 --bfile step3_samples_filtered --keep step6_neural_ad/step6_EUR.txt --het --out step5_EUR_het
plink1.9 --bfile step3_samples_filtered --keep step6_neural_ad/step6_AMR.txt --het --out step5_AMR_het
plink1.9 --bfile step3_samples_filtered --keep step6_neural_ad/step6_EAS.txt --het --out step5_EAS_het
plink1.9 --bfile step3_samples_filtered --keep step6_neural_ad/step6_SAS.txt --het --out step5_SAS_het

plink1.9 --bfile step3_samples_filtered --chr 1-22 --keep step6_neural_ad/step6_AFR.txt --het --out step5_AFR_1-22_het
plink1.9 --bfile step3_samples_filtered --chr 1-22 --keep step6_neural_ad/step6_EUR.txt --het --out step5_EUR_1-22_het
plink1.9 --bfile step3_samples_filtered --chr 1-22 --keep step6_neural_ad/step6_AMR.txt --het --out step5_AMR_1-22_het
plink1.9 --bfile step3_samples_filtered --chr 1-22 --keep step6_neural_ad/step6_EAS.txt --het --out step5_EAS_1-22_het
plink1.9 --bfile step3_samples_filtered --chr 1-22 --keep step6_neural_ad/step6_SAS.txt --het --out step5_SAS_1-22_het

Rscript ../PROGRAM/check_heterozygosity.r step5_AFR_1-22_het.het step5_AFR_1-22
Rscript ../PROGRAM/check_heterozygosity.r step5_EUR_1-22_het.het step5_EUR_1-22
Rscript ../PROGRAM/check_heterozygosity.r step5_AMR_1-22_het.het step5_AMR_1-22
Rscript ../PROGRAM/check_heterozygosity.r step5_EAS_1-22_het.het step5_EAS_1-22
Rscript ../PROGRAM/check_heterozygosity.r step5_SAS_1-22_het.het step5_SAS_1-22


wc -lc step5_AFR_1-22_heterozygosity_outliers.txt # 9
wc -lc step5_EUR_1-22_heterozygosity_outliers.txt # 757
wc -lc step5_AMR_1-22_heterozygosity_outliers.txt # 46
wc -lc step5_EAS_1-22_heterozygosity_outliers.txt # 2
wc -lc step5_SAS_1-22_heterozygosity_outliers.txt # 16

# step5. ROH
plink1.9 --bfile step3_samples_filtered --chr 1-22 --keep step6_neural_ad/step6_AFR.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out step5_AFR_ROH
plink1.9 --bfile step3_samples_filtered --chr 1-22 --keep step6_neural_ad/step6_AMR.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out step5_AMR_ROH
plink1.9 --bfile step3_samples_filtered --chr 1-22 --keep step6_neural_ad/step6_EAS.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out step5_EAS_ROH
plink1.9 --bfile step3_samples_filtered --chr 1-22 --keep step6_neural_ad/step6_EUR.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out step5_EUR_ROH
plink1.9 --bfile step3_samples_filtered --chr 1-22 --keep step6_neural_ad/step6_SAS.txt --homozyg --homozyg-kb 1000 --homozyg-snp 50 --homozyg-gap 100 --out step5_SAS_ROH

Rscript ../PROGRAM/plot_ROH.r step5_AFR_ROH.hom.indiv step5_AFR_het.het AFR
Rscript ../PROGRAM/plot_ROH.r step5_AMR_ROH.hom.indiv step5_AMR_het.het AMR
Rscript ../PROGRAM/plot_ROH.r step5_EAS_ROH.hom.indiv step5_EAS_het.het EAS
Rscript ../PROGRAM/plot_ROH.r step5_EUR_ROH.hom.indiv step5_EUR_het.het EUR
Rscript ../PROGRAM/plot_ROH.r step5_SAS_ROH.hom.indiv step5_SAS_het.het SAS

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
