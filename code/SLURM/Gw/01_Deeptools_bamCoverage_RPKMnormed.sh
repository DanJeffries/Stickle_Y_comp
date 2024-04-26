#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --array=1-36
#SBATCH --job-name=Ga_DEEP_1000
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module add UHTS/Analysis/deepTools/2.5.4;

ID=$SLURM_ARRAY_TASK_ID

WD=/storage/scratch/iee/dj20y461/Stickleback/Y_comp/Sex_chrom_stats/Read_depths/Gw
ITER_FILE=/storage/homefs/dj20y461/Stickleback/Y_comp/Find_2n_windows/Gw/code/samples.txt
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)

BAM_DIR=$WD/bam
BAM=${BAM_DIR}/${SAMPLE_NAME}.markdup.bam

OUTDIR=$WD/depths

if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
fi


bamCoverage --bam $BAM \
            --binSize 1000 \
            --numberOfProcessors 16 \
            --verbose \
	    --normalizeUsingRPKM \
            --outFileName ${OUTDIR}/${SAMPLE_NAME}.1kb.RPKM.depth \
            --outFileFormat bedgraph




