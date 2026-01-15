#!/bin/bash
#----------------------
#SBATCH --account=gratis
#----------------------
#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --export=NONE
#SBATCH --array=2-9326
#SBATCH --job-name=Trim_and_align_all_virus
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module add Trimmomatic/0.39-Java-11
module add BWA/0.7.17-GCC-10.3.0
module add samblaster/0.1.26-GCC-10.3.0
module add SAMtools/1.13-GCC-10.3.0
module add Sambamba/0.8.2-GCC-10.3.0

echo "start time"
date

#WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_calling

WD=/storage/research/iee_temp_dj20y461/DV_calling/virus_screening_full

ID=$SLURM_ARRAY_TASK_ID
#ITER_FILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/sample_metadata/random_125perpop_paths.txt
## positive infection samples
#ITER_FILE=/storage/research/iee_temp_dj20y461/DV_calling/virus_screening/raw_positives/positives_path.txt
## All samples
ITER_FILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/sample_metadata/sample_paths_rmdup.txt

SAMPLE_FASTQ_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)
SAMPLE_NAME=$(echo $SAMPLE_FASTQ_PATH | rev | cut -f1 -d'/' | rev)

######### Trimming ##############################################################################################################################################################################################

ADAPTERS=/software.9/software/Trimmomatic/0.39-Java-11/adapters/NexteraPE-PE.fa

TRIMMED_OUTDIR=$WD/trimmed

if [ ! -d "$TRIMMED_OUTDIR" ]
then
    mkdir -p $TRIMMED_OUTDIR
fi

## OUT FILES

SAMPLE_R1_trimmed=$TRIMMED_OUTDIR/${SAMPLE_NAME}.R1.trimmed.fastq.gz
SAMPLE_R1_unpaired=$TRIMMED_OUTDIR/${SAMPLE_NAME}.R1.unpaired.fastq.gz
SAMPLE_R2_trimmed=$TRIMMED_OUTDIR/${SAMPLE_NAME}.R2.trimmed.fastq.gz
SAMPLE_R2_unpaired=$TRIMMED_OUTDIR/${SAMPLE_NAME}.R2.unpaired.fastq.gz

if [ ! -f "$SAMPLE_R1_trimmed" ]
then
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
         -threads 8 \
      	 $SAMPLE_FASTQ_PATH*_R1*fastq.gz $SAMPLE_FASTQ_PATH*_R2*fastq.gz \
       	 $SAMPLE_R1_trimmed $SAMPLE_R1_unpaired \
       	 $SAMPLE_R2_trimmed $SAMPLE_R2_unpaired \
       	 ILLUMINACLIP:${ADAPTERS}:3:30:10 SLIDINGWINDOW:4:20 MINLEN:25 HEADCROP:7
fi


######### Aligning ##############################################################################################################################################################################################

DATA_DIR=$WD/trimmed

GENOME_IDX=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/ref/Gacu_Herpes_Irido_ref/bwa/Gacu_AXYMtUn_Herpes_Irido.fa

## Get the read group info to add to bams from the first read of every file

HEADER=$(zcat $SAMPLE_R1_trimmed | head -n1) ## get first read ID line 

INSTRUMENT=$(echo $HEADER | cut -f1 -d':'| sed 's/@//g')
RUN=$(echo $HEADER | cut -f2 -d':')
FLOWCELL=$(echo $HEADER | cut -f3 -d':')
LANE=$(echo $HEADER | cut -f4 -d':')
BARCODES=$(echo $HEADER | cut -f2 -d' '| cut -f4 -d':' | sed 's/+/-/g')

echo "@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tPU:${FLOWCELL}.${LANE}\tSM:${SAMPLE_NAME}\tCN:MCGILL\tBC:${BARCODES}"

echo $READGROUP_STRING

BAMDIR=$WD/bams/

if [ ! -d "$BAMDIR" ]; then
   mkdir -p $BAMDIR
fi

BAM_FIX_COORDSORT=${BAMDIR}/${SAMPLE_NAME}.fixmate.coordsorted.bam

echo "bwa mem -t 20 \
        $GENOME_IDX \
        $SAMPLE_R1_trimmed \
        $SAMPLE_R2_trimmed"

#-R $(echo "@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tPU:${FLOWCELL}.${LANE}\tSM:${SAMPLE_NAME}\tCN:MCGILL\tBC:${BARCODES}") \

#bwa mem -t 20 \
#        $GENOME_IDX \
#        $SAMPLE_R1_trimmed \
#        $SAMPLE_R2_trimmed | \
#samblaster | \
#samtools fixmate \
#	-m \
#        -@ 20 \
#        - \
#	- | \
#samtools sort \
#	-@ 20 \
#	-O BAM \
#	-o $BAM_FIX_COORDSORT \
#	-
#
#samtools index -c -@ 20 $BAM_FIX_COORDSORT

## End time
#echo "end time"
#date

##################################################################
####### calculate some stats ####################################
#################################################################

#samtools stats $BAM_FIX_COORDSORT > ${BAMDIR}/${SAMPLE_NAME}.stats
#samtools coverage $BAM_FIX_COORDSORT  > ${BAMDIR}/${SAMPLE_NAME}.depths

## Calculate some averages and add to the end of each file
#grep 'chromosome' ${BAMDIR}/${SAMPLE_NAME}.depths | awk '{ sum += $6; n++ } END { if (n > 0) print "\n## Avg. perc. bases. covered (chroms) = " sum / n; }' >> ${BAMDIR}/${SAMPLE_NAME}.depths
#grep 'chromosome' ${BAMDIR}/${SAMPLE_NAME}.depths | awk '{ sum += $7; n++ } END { if (n > 0) print "## Avg. coverage (chroms) = " sum / n; }' >> ${BAMDIR}/${SAMPLE_NAME}.depths


######### Cleanup ##############################################################################################################################################################################################

## the raw bam files probably don't need to be kept. . . . lets see how it goes for now. 

## Can I store bams as CRAM? 



