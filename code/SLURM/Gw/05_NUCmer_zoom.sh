#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --export=NONE
#SBATCH --job-name=Gw_NUCmer_XY
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load SAMtools/1.13-GCC-10.3.0
module load gnuplot/5.4.2-GCCcore-10.3.0


MUMmer=~/Software/MUMmer/bin/mummer
NUCmer=~/Software/MUMmer/bin/nucmer
MUMmerplot=~/Software/MUMmer/bin/mummerplot
SHOW_COORDS=~/Software/MUMmer/bin/show-coords

GwX=/storage/scratch/iee_evol/dj20y461/Stickle_Y_comp/Gw_MUMmer/Gw_X.fa
GwY=/storage/scratch/iee_evol/dj20y461/Stickle_Y_comp/Gw_MUMmer/Gw_Y.fa
OUTDIR=/storage/scratch/iee_evol/dj20y461/Stickle_Y_comp/Gw_MUMmer/

### Make fa files for the comparison (this is better than doing the whole aligment and zooming in with mummerplot as it avoids many alignments across the genome) 

## REF and QUERY start and stops are 1kb up and downstream of the breakpoint inferred from the MUMmer alignment for the whole inversion. 

WINDOW_SIZE=2000
SPACER=$(($WINDOW_SIZE/2)) ## dist to add or subtract around the breakpoint to acheive window size

REF_BREAKPOINT=19389071

REF_CHROM='chr12'
REF_START=$(($REF_BREAKPOINT-$SPACER))
REF_STOP=$(($REF_BREAKPOINT+$SPACER))
REF=$OUTDIR/Gw_X_${REF_CHROM}_${REF_START}_${REF_STOP}.fa

QUERY_BREAKPOINT=21306767

QUERY_CHROM='chr12'
QUERY_START=$(($QUERY_BREAKPOINT-$SPACER))
QUERY_STOP=$(($QUERY_BREAKPOINT+$SPACER))
QUERY=$OUTDIR/Gw_X_${QUERY_CHROM}_${QUERY_START}_${QUERY_STOP}.fa

samtools faidx $GwX $REF_CHROM:$REF_START-$REF_STOP > $REF
samtools faidx $GwX $QUERY_CHROM:$QUERY_START-$QUERY_STOP > $QUERY

LEN=20

NUCmer_OUT=$OUTDIR/${REF_CHROM}_${REF_START}_${REF_STOP}_vs_${QUERY_CHROM}_${QUERY_START}_${QUERY_STOP}_LEN${LEN}_maxmatch_NUC_c100

$NUCmer --maxmatch \
       	-p $NUCmer_OUT \
	-l $LEN \
	-c 100 \
	$REF \
        $QUERY

$SHOW_COORDS -r -c -l $NUCmer_OUT.delta > $NUCmer_OUT.coords

$MUMmerplot  -t png \
	     -p $OUTDIR/${REF_CHROM}_${REF_START}_${REF_STOP}_vs_${QUERY_CHROM}_${QUERY_START}_${QUERY_STOP}_LEN${LEN}_maxmatch_NUC_c100 \
	     $NUCmer_OUT.delta



#MUMmer_OUT=$OUTDIR/${REF_CHROM}_${REF_START}_${REF_STOP}_vs_${QUERY_CHROM}_${QUERY_START}_${QUERY_STOP}_LEN${LEN}_maxmatch.mums

#$MUMmer -mum \
#       -b \
#       -c \
#       -l $LEN \
#       --maxmatch \
#        $REF \
#       $QUERY \
#       > $MUMmer_OUT
