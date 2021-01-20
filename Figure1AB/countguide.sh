#!/bin/bash
#SBATCH -n 8                              # 1 core
#SBATCH -t 1-00:00                         # in D-HH:MM format
#SBATCH -p medium                           # Run in short partition
#SBATCH --mem=64G
#SBATCH -o countguide.out
#SBATCH -e countguide.err

### Above is important when submitting as a job

### Please download the fastq files from GSE160456, and put the fastq files into the "data" folder
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160456

folder=./data/


python guide_counter_general.py -l \
C9B_1_high \
C9B_1_presort \
C9B_2_high \
C9B_2_presort \
C9B_3_high \
C9B_3_presort \
dCas9_1_high \
dCas9_1_presort \
dCas9_2_high \
dCas9_2_presort \
dCas9_3_high \
dCas9_3_presort \
KRAB_1_high \
KRAB_1_presort \
KRAB_2_high \
KRAB_2_presort \
KRAB_3_high \
KRAB_3_presort \
VP64_1_high \
VP64_1_presort \
VP64_2_high \
VP64_2_presort \
VP64_3_high \
VP64_3_presort \
-n hbb_screen -o ./ -p 8 sgRNA_lib.txt \
${folder}/C9B_1_high_S2_R1_001.fastq.gz,${folder}/C9B_1_high_S2_R2_001.fastq.gz \
${folder}/C9B_1_presort_S1_R1_001.fastq.gz,${folder}/C9B_1_presort_S1_R2_001.fastq.gz \
${folder}/C9B_2_high_S5_R1_001.fastq.gz,${folder}/C9B_2_high_S5_R2_001.fastq.gz \
${folder}/C9B_2_presort_S4_R1_001.fastq.gz,${folder}/C9B_2_presort_S4_R2_001.fastq.gz \
${folder}/C9B_3_high_S7_R1_001.fastq.gz,${folder}/C9B_3_high_S7_R2_001.fastq.gz \
${folder}/C9B_3_presort_S6_R1_001.fastq.gz,${folder}/C9B_3_presort_S6_R2_001.fastq.gz \
${folder}/dCas9_1_high_S16_R1_001.fastq.gz,${folder}/dCas9_1_high_S16_R2_001.fastq.gz \
${folder}/dCas9_1_presort_S15_R1_001.fastq.gz,${folder}/dCas9_1_presort_S15_R2_001.fastq.gz \
${folder}/dCas9_2_high_S19_R1_001.fastq.gz,${folder}/dCas9_2_high_S19_R2_001.fastq.gz \
${folder}/dCas9_2_presort_S18_R1_001.fastq.gz,${folder}/dCas9_2_presort_S18_R2_001.fastq.gz \
${folder}/dCas9_3_high_S21_R1_001.fastq.gz,${folder}/dCas9_3_high_S21_R2_001.fastq.gz \
${folder}/dCas9_3_presort_S20_R1_001.fastq.gz,${folder}/dCas9_3_presort_S20_R2_001.fastq.gz \
${folder}/KRAB_1_high_S9_R1_001.fastq.gz,${folder}/KRAB_1_high_S9_R2_001.fastq.gz \
${folder}/KRAB_1_presort_S8_R1_001.fastq.gz,${folder}/KRAB_1_presort_S8_R2_001.fastq.gz \
${folder}/KRAB_2_high_S12_R1_001.fastq.gz,${folder}/KRAB_2_high_S12_R2_001.fastq.gz \
${folder}/KRAB_2_presort_S11_R1_001.fastq.gz,${folder}/KRAB_2_presort_S11_R2_001.fastq.gz \
${folder}/KRAB_3_high_S14_R1_001.fastq.gz,${folder}/KRAB_3_high_S14_R2_001.fastq.gz \
${folder}/KRAB_3_presort_S13_R1_001.fastq.gz,${folder}/KRAB_3_presort_S13_R2_001.fastq.gz \
${folder}/VP64_1_high_S23_R1_001.fastq.gz,${folder}/VP64_1_high_S23_R2_001.fastq.gz \
${folder}/VP64_1_presort_S22_R1_001.fastq.gz,${folder}/VP64_1_presort_S22_R2_001.fastq.gz \
${folder}/VP64_2_high_S26_R1_001.fastq.gz,${folder}/VP64_2_high_S26_R2_001.fastq.gz \
${folder}/VP64_2_presort_S25_R1_001.fastq.gz,${folder}/VP64_2_presort_S25_R1_001.fastq.gz \
${folder}/VP64_3_high_S28_R1_001.fastq.gz,${folder}/VP64_3_high_S28_R1_001.fastq.gz \
${folder}/VP64_3_presort_S27_R1_001.fastq.gz,${folder}/VP64_3_presort_S27_R1_001.fastq.gz
