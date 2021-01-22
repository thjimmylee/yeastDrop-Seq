#!/bin/bash

### Check fastq file quality with fastqc
fastqc DR1.fastq DR2.fastq.gz GMR1.fastq.gz GMR2.fastq.gz GR1.fastq.gz GR2.fastq.gz MR1.fastq.gz MR2.fastq.gz  -o ../../fastqc/srr_fastqc/


file_head=(DR GMR GR MR)
SAMPLE=(DMSO Guanine-MPA Guanine MPA)

### Sorting fastq.gz

WORKINGDIR=/lustre/scratch117/cellgen/team218/JL/yeastdropseq/GSE144636/fastq/00_srr_fastq
OUTPUTDIR=/lustre/scratch117/cellgen/team218/JL/yeastdropseq/GSE144636/fastq/00_srr_fastq_sorted

INPUTFILES=$WORKINGDIR/*.fastq.gz

REPORTPREFFIX=../../reports/

mkdir -p $OUTPUTDIR

cd $WORKINGDIR

for FILE in $INPUTFILES ; do
	j=`echo $FILE | awk '{ gsub(".fastq.gz", "_sorted.fastq"); print $0 }'`

	jobname=`basename $FILE`

	bsub -R"select[mem>10000] rusage[mem=10000]" -M10000 -q normal -J $jobname  \
    in_house_scripts/DO_FASTQGZ_SORTING.sh $FILE $j
done

### Trimming adaptors

WORKINGDIR=/lustre/scratch117/cellgen/team218/JL/yeastdropseq/GSE144636
INPUTDIR=00_srr_fastq
OUTPUTDIR=01_fastq_trimmed

REPORTPREFFIX=$WORKINGDIR"/reports/"$INPUTDIR"-to-"$OUTPUTDIR

mkdir -p $WORKINGDIR"/fastq/"$OUTPUTDIR

cd $WORKINGDIR

for i in {0..3}; do
	bsub -R"select[mem>10000] rusage[mem=10000]" -M10000 -q normal \
    -J ${SAMPLE[i]} \
    in_house_scripts/DO_TRIM_PE_ADAPTERS.sh ${SAMPLE[i]} $WORKINGDIR"/fastq/"$INPUTDIR/${file_head[i]}"1_sorted.fastq.gz" $WORKINGDIR"/fastq/"$INPUTDIR/${file_head[i]}"2_sorted.fastq.gz" $WORKINGDIR/$OUTPUTDIR/
done


### Create barcode whitelist

INPUTDIR=01_fastq_trimmed
OUTPUTDIR=whitelist

REPORTPREFFIX=$WORKINGDIR/reports/$INPUTDIR"-to-"$OUTPUTDIR

mkdir -p $WORKINGDIR/$OUTPUTDIR

cell_no=(527 794 384 857)

cd $WORKINGDIR

for i in {0..3}; do
	bsub -R"select[mem>10000] rusage[mem=10000]" -M10000 -q normal -J ${SAMPLE[i]}  \
    in_house_scripts/DO_UMITOOLS_WHITELIST.sh $WORKINGDIR"/fastq/"$INPUTDIR"/TRIMMED-"${file_head[i]}"1_sorted.fastq.gz" CCCCCCCCCCCCNNNNNNNN ${cell_no[i]} $WORKINGDIR/$OUTPUTDIR"/TRIMMED_"${SAMPLE[i]}"_whitelist.txt" 
done

### umi-tools extract

INPUTDIR=01_fastq_trimmed
OUTPUTDIR=02_fastq_extract

REPORTPREFFIX=$WORKINGDIR"/reports/"$INPUTDIR"-to-"$OUTPUTDIR

mkdir -p $WORKINGDIR"/fastq/"$OUTPUTDIR

cd $WORKINGDIR

for i in {0..3}; do
	bsub -R"select[mem>10000] rusage[mem=10000]" -M10000 -q normal -J ${SAMPLE[i]} \
    in_house_scripts/DO_UMITOOLS_EXTRACT.sh $WORKINGDIR"/fastq/"$INPUTDIR"/TRIMMED-"${file_head[i]}"1_sorted.fastq.gz" $WORKINGDIR"/fastq/"$INPUTDIR"/TRIMMED-"${file_head[i]}"2_sorted.fastq.gz" CCCCCCCCCCCCNNNNNNNN ./whitelist/TRIMMED_${SAMPLE[i]}"_whitelist.txt" $WORKINGDIR"/fastq/"$OUTPUTDIR
done
