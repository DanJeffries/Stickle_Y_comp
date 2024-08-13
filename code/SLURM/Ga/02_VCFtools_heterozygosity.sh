#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --job-name=VCFtools_Hardy
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err

module load vital-it
module add UHTS/Analysis/vcftools/0.1.15;

INVCF=/storage/research/iee_evol/StickleYproject_DJandZY_2024/99_variant_calls/Ga/vcfs/Ga.SNP.Q30.GQ20.meanDP10-50.minDP8.vcf.gz
OUT_PREFIX=/storage/scratch/iee/dj20y461/Stickleback/Y_comp/Sex_chrom_stats/Heterozygosity/Ga/Ga.SNP.Q30.GQ20.meanDP10-50.minDP8
MALES=/storage/homefs/dj20y461/Stickleback/Stickle_Y_comp/code/SLURM/Ga/males.txt
FEMALES=/storage/homefs/dj20y461/Stickleback/Stickle_Y_comp/code/SLURM/Ga/females.txt

vcftools \
    --gzvcf $INVCF \
    --keep $MALES \
    --max-missing-count 10 \
    --hardy \
    --out ${OUT_PREFIX}_MALES

vcftools \
    --gzvcf $INVCF \
    --keep $FEMALES \
    --max-missing-count 10 \
    --hardy \
    --out ${OUT_PREFIX}_FEMALES
