#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --job-name=Gw_MUMmer_XY
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load gnuplot/5.4.2-GCCcore-10.3.0

GwX=/storage/scratch/iee_evol/dj20y461/Stickle_Y_comp/Gw_MUMmer/Gw_X.fa
GwY=/storage/scratch/iee_evol/dj20y461/Stickle_Y_comp/Gw_MUMmer/Gw_Y.fa
OUTDIR=/storage/scratch/iee_evol/dj20y461/Stickle_Y_comp/Gw_MUMmer/

MUMmer=~/Software/MUMmer/bin/mummer
MUMmerplot=~/Software/MUMmer/bin/mummerplot

LEN=500

$MUMmer -mum \
	-b \
	-c \
	-l $LEN \
        $GwX \
	$GwY \
	> $OUTDIR/Gw_XY_l${LEN}.mums	


$MUMmerplot  -t png \
	     -p $OUTDIR/Gw_XY_l${LEN} \
	     -R $GwX \
	     --filter \
	     --layout \
	     $OUTDIR/Gw_XY_l${LEN}.mums


