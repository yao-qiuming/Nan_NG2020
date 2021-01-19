#From Nan Liu 

#run CUT&RUN tools, locate to workdir

sbatch ./integrated.sh CR_NFYA_cloneB6p_Day3_r1_S3_R1_001.fastq.gz
cd aligned.aug10
sbatch ./integrated.step2.sh CR_NFYA_cloneB6p_Day3_r1_S3_aligned_reads.bam
cd ../macs2.narrow.aug18
sbatch ./integrate.motif.find.sh CR_NFYA_cloneB6p_Day3_r1_S3_aligned_reads_peaks.narrowPeak
sbatch ./integrate.footprinting.sh CR_NFYA_cloneB6p_Day3_r1_S3_aligned_reads_peaks.narrowPeak



#obtain cut frequency of the HBG promoter
#locate to workdir
cd macs2.narrow.aug18/
./get_cuts_single_locus.sh chr11:5,267,523-5,277,680 ../aligned.aug10/dup.marked.120bp/CR_NFYA_cloneB6p_Day3_r1_S3_aligned_reads.bam single.locus

