#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --array=1-18
#SBATCH --job-name=SAMSTATS
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/samtools/1.10

ITERFILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/code/samples.txt

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITERFILE)

samtools stats $WD/bams/${SAMPLE_NAME}.fixmate.coordsorted.bam > $WD/bams/${SAMPLE_NAME}.stats
samtools coverage $WD/bams/${SAMPLE_NAME}.fixmate.coordsorted.bam > $WD/bams/${SAMPLE_NAME}.depths

## Calculate some averages and add to the end of each file
egrep -v '#|chrUn|chrM' $WD/bams/${SAMPLE_NAME}.depths | awk '{ sum += $6; n++ } END { if (n > 0) print "\n## Avg. perc. bases. covered = " sum / n; }' >> $WD/bams/${SAMPLE_NAME}.depths
egrep -v '#|chrUn|chrM' $WD/bams/${SAMPLE_NAME}.depths | awk '{ sum += $7; n++ } END { if (n > 0) print "## Avg. coverage = " sum / n; }' >> $WD/bams/${SAMPLE_NAME}.depths

~
