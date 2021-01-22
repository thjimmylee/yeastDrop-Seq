#!/bin/bash

echo Running $0 $1 $2 $3 $4

INPUTR1FASTQGZ=$1
BARCODEPATTERN=$2
CELLNUMBER=$3
OUTPUTFILE=$4

if [ -z "$INPUTR1FASTQGZ" ] ; then
  echo "Please provide file path to R1 (Argument 1/4)"
  exit 1
fi

if [ -z "$BARCODEPATTERN" ] ; then
  echo "Please provide a valid barcode sequence e.g. CCCCCCCCCCCCCCCCNNNNNNNNNNi which C for barcode and N for umi (Argument 2/4)"
	exit 1
fi

if [ -z "$CELLNUMBER" ] ; then
  echo "Please provide an approximate cell number of your sample (Argument 3/4)"
  exit 1
fi

if [ -z "$OUTPUTFILE" ] ; then
  echo "Please provide file path to the output file (Argument 4/4)"
  exit 1
fi

umi_tools whitelist --stdin $INPUTR1FASTQGZ --bc-pattern=$BARCODEPATTERN --set-cell-number=$CELLNUMBER --log2stderr > $OUTPUTFILE

echo DONE $0 $1 $2 $3
echo debug $1

