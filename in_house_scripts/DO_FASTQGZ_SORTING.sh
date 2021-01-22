#!/bin/bash

echo Running $0 $1 $2

INPUTFILE=$1
OUTPUTFILE=$2

if [ -z "$INPUTFILE" ] ; then
  echo "Please provide file path to unsorted fastq.gz (Argument 1/2)"
  exit 1
fi

if [ -z "$OUTPUTFILE" ] ; then
  echo "Please provide file path to output fastq.gz (Argument 2/2)"
  exit 1
fi

zcat $INPUTFILE | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > $OUTPUTFILE

gzip $OUTPUTFILE

echo DONE $0 $1 $2
echo debug $1

