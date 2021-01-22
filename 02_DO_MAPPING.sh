#!/bin/bash

WORKINGDIR=/lustre/scratch117/cellgen/team218/JL/yeastdropseq/GSE144636

file_head=(DR GMR GR MR)
SAMPLE=(DMSO Guanine-MPA Guanine MPA)

## Genome Indexing

# For SCere
bsub -n24 -R"span[hosts=1] select[mem>35000] rusage[mem=35000]" -M35000 -q long -J IndexingSC \
STAR --runThreadN 24 \
        --runMode genomeGenerate \
        --genomeDir /lustre/scratch117/cellgen/team218/JL/GENOMES/yeast/SCindex \
        --genomeFastaFiles /lustre/scratch117/cellgen/team218/JL/GENOMES/yeast/Scer.fa \
        --sjdbGTFfile /lustre/scratch117/cellgen/team218/JL/GENOMES/yeast/Genes.gtf \
        --sjdbOverhang 64 \
        --limitGenomeGenerateRAM 31000000000

## Mapping to SCere index

INPUTDIR=02_fastq_extract
OUTPUTDIR=03_bam_map-to-SC

REPORTPREFFIX=$WORKINGDIR"/reports/"$INPUTDIR"-to-"$OUTPUTDIR

cd $WORKINGDIR

mkdir -p $WORKINGDIR"/bam/"$OUTPUTDIR

for i in {0..3}; do
	bsub -n24 -R"span[hosts=1] select[mem>90000] rusage[mem=90000]" -M90000 -q normal -J ${SAMPLE[i]}"_Map" \
STAR --runThreadN 24 \
       --genomeDir /lustre/scratch117/cellgen/team218/JL/GENOMES/yeast/SCindex \
       --readFilesIn $WORKINGDIR"/fastq/"$INPUTDIR"/EXTRACT-"${file_head[i]}"2_sorted.fastq.gz" \
       --readFilesCommand zcat \
       --outFilterMultimapNmax 1 \
       --runMode alignReads \
       --parametersFiles in_house_scripts/STAR_Parameters.txt \
       --outFileNamePrefix $WORKINGDIR"/bam/"$OUTPUTDIR/${SAMPLE[i]}"_"
done
