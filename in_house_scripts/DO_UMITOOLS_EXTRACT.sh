#!/bin/bash

echo Running $1 $2 $3 $4 $5

INPUTR1FASTQGZ=$1
INPUTR2FASTQGZ=$2
BARCODEPATTERN=$3
WHITELISTFILE=$4
OUTPUTDIR=$5

if [ -z "$INPUTR1FASTQGZ" ] ; then
	echo "Please provide file path to R1.fastq.gz (Argument 1/5)"
  exit 1
fi

if [ -z "$INPUTR2FASTQGZ" ] ; then
        echo "Please provide file path to R2.fastq.gz (Argument 2/5)"
  exit 1
fi

if [ -z "$BARCODEPATTERN" ] ; then
  echo "Please provide a valid barcode sequence e.g. CCCCCCCCCCCCCCCCNNNNNNNNNNi which C for barcode and N for umi (Argument 3/5)"
  exit 1
fi

if [ -z "$WHITELISTFILE" ] ; then
  echo "Please provide path to the whitelist. To create one, use DO_UMITOOLS_WHITELIST.sh (Argument 4/5)"
  exit 1
fi

if [ -z "$OUTPUTDIR" ] ; then
  echo "Please provide directory path for the output file (Argument 5/5)"
  exit 1
fi

mkdir -p $OUTPUTDIR

	OUTPUTR1FQGZ=`basename $INPUTR1FASTQGZ`
	OUTPUTR1FQGZ=$OUTPUTDIR"/"$OUTPUTR1FQGZ
	OUTPUTR1FQGZ=${OUTPUTR1FQGZ/TRIMMED/EXTRACT}
	OUTPUTR2FQGZ=`basename $INPUTR2FASTQGZ`
        OUTPUTR2FQGZ=$OUTPUTDIR"/"$OUTPUTR2FQGZ
        OUTPUTR2FQGZ=${OUTPUTR2FQGZ/TRIMMED/EXTRACT}

umi_tools extract --bc-pattern=$BARCODEPATTERN --stdin $INPUTR1FASTQGZ --stdout $OUTPUTR1FQGZ --read2-in $INPUTR2FASTQGZ --read2-out=$OUTPUTR2FQGZ --filter-cell-barcode --whitelist=$WHITELISTFILE 

echo DONE $0 $1 $2 $3 $4 $5
echo debug $1

