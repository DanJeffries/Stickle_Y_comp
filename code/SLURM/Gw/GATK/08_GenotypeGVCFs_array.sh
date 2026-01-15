#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --job-name=GT_GVCFs
#SBATCH --array=2-22
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
export GATK_LOCAL_JAR=/software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/

INTERVALS=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/code/array_files/Intervals_GVCF.txt
INTERVAL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $INTERVALS)

ASSEMBLY=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
GENDB=$WD/GenDBs/GenDB_${INTERVAL}
OUTDIR=$WD/Joint_VCFs/Joint_${INTERVAL}/

if [ ! -d $OUTDIR ]; then
   mkdir $OUTDIR
fi

gatk --java-options "-Xmx24g -Xms4g" GenotypeGVCFs \
   -R ${ASSEMBLY} \
   -V gendb://${GENDB} \
   -O ${OUTDIR}/${INTERVAL}.vcf 



