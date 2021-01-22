#!/bin/bash

WORKINGDIR=/lustre/scratch117/cellgen/team218/JL/yeastdropseq/GSE144636

SAMPLE=(DMSO Guanine-MPA Guanine MPA)

### featureCounts

INPUTDIR=03_bam_map-to-SC
OUTPUTDIR=04_bam_featureCounts

REPORTPREFFIX=$WORKINGDIR"/reports/"$INPUTDIR"-to-"$OUTPUTDIR

cd $WORKINGDIR

mkdir -p $WORKINGDIR"/bam/"$OUTPUTDIR

for i in {0..3}; do
    bsub -n16 -R"span[hosts=1] select[mem>60000] rusage[mem=60000]" -M60000 -q normal -J ${SAMPLE[i]}"_fC" \
    featureCounts -T 16 -a /lustre/scratch117/cellgen/team218/JL/GENOMES/yeast/Genes.gtf -o $WORKINGDIR"/bam/"$OUTPUTDIR/ -R BAM $WORKINGDIR"/bam/"$INPUTDIR/${SAMPLE[i]}"_Aligned.sortedByCoord.out.bam"

done

### sorting

INPUTDIR=04_bam_featureCounts
OUTPUTDIR=05_sortedfeatureCounts

REPORTPREFFIX=$WORKINGDIR"/reports/"$INPUTDIR"-to-"$OUTPUTDIR

cd $WORKINGDIR

mkdir -p $WORKINGDIR"/bam/"$OUTPUTDIR

for i in {0..3}; do

	bsub -R"select[mem>10000] rusage[mem=10000]" -M10000 -q normal -J ${SAMPLE[i]}"_sortfC"  \
    samtools sort $WORKINGDIR"/bam/"$INPUTDIR/${SAMPLE[i]}"_Aligned.sortedByCoord.out.bam.featureCounts.bam" -o $WORKINGDIR"/bam/"$OUTPUTDIR/${SAMPLE[i]}"_assigned_sorted.bam"
    
done

## index sorted bam

OUTPUTDIR=05_sortedfeatureCounts

cd $WORKINGDIR

for i in {0..3}; do
	samtools index $WORKINGDIR"/bam/"$OUTPUTDIR/${SAMPLE[i]}"_assigned_sorted.bam"
done


### creating raw count matrix 

INPUTDIR=05_sortedfeatureCounts
OUTPUTDIR=count_matrix

REPORTPREFFIX=$WORKINGDIR"/reports/"$INPUTDIR"-to-"$OUTPUTDIR

cd $WORKINGDIR

mkdir -p $WORKINGDIR/$OUTPUTDIR

for i in {0..3}; do
	bsub -R"select[mem>10000] rusage[mem=10000]" -M10000 -q normal -J ${SAMPLE[i]}"_countMX" \
    umi_tools count --per-gene --wide-format-cell-counts --gene-tag=XT --assigned-status-tag=XS --per-cell -I $WORKINGDIR"/bam/"$INPUTDIR/${SAMPLE[i]}"_assigned_sorted.bam" -S $WORKINGDIR/$OUTPUTDIR/${SAMPLE[i]}"_counts.tsv.gz"
done
