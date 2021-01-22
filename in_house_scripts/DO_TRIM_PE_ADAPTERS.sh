#!/bin/bash

echo Running $0 $1 $2 $3 $4

SAMPLE=$1 # 1
fastqR1=$2
fastqR2=$3
OUTPUTDIR=$4
TRIMMER=/software/team218/JL/bin/Trimmomatic-0.39/trimmomatic-0.39.jar

if [ -z "$SAMPLE" ] ; then
  echo "Please provide sample ID (Argument 1/4)"
  exit 1
fi

if [ -z "$fastqR1" ] ; then
	echo "Please provide R1.fastq.gz (Argument 2/4)"
	exit 1
fi

if [ -z "$fastqR2" ] ; then
	echo "Please provide R2.fastq.gz (Argument 3/4)"

fi

if [ -z "$OUTPUTDIR" ] ; then
  echo "Please provide a directory for outputfiles (Argument 4/4)"
  exit 1
fi

if [ ! -f "$TRIMMER" ] ; then
  echo "Sorry $TRIMMER not available "
  exit 1
fi

mkdir -p $OUTPUTDIR
mkdir -p $OUTPUTDIR""output_unpaired/


export _JAVA_OPTIONS="-Xmx1000M -XX:MaxHeapSize=1000m"

java -jar $TRIMMER PE -phred33 $fastqR1 $fastqR2 $OUTPUTDIR""TRIMMED-$(basename $fastqR1) $OUTPUTDIR""output_unpaired/$(basename $fastqR1) $OUTPUTDIR""TRIMMED-$(basename $fastqR2) $OUTPUTDIR""output_unpaired/$(basename $fastqR2) ILLUMINACLIP:/software/team218/JL/bin/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30

echo $SAMPLE"",$((`zcat $fastqR1 | wc -l`/4))"",$((`zcat $fastqR2 | wc -l`/4))"",$((`zcat $OUTPUTDIR""TRIMMED-$(basename $fastqR1) | wc -l`/4))"",$((`zcat $OUTPUTDIR""TRIMMED-$(basename $fastqR2) | wc -l`/4))"",$((`zcat $OUTPUTDIR""output_unpaired/$(basename $fastqR1) | wc -l`/4))"",$((`zcat $OUTPUTDIR""output_unpaired/$(basename $fastqR2) | wc -l`/4)) >> $OUTPUTDIR""$SAMPLE.csv


echo DONE $0 $1 $2 $3
echo debug $1
