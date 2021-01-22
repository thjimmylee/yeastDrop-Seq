#!/bin/bash

# samtools
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar -xjf samtools-1.10.tar.bz2
./configure --prefix=/software/team218/JL/bin
make
make install

# Fastqc
unzip fastqc_v0.11.9.zip 

# STAR
cd /software/team218/JL/app
wget https://github.com/alexdobin/STAR/archive/2.7.5a.tar.gz
tar -xzf 2.7.5a.tar.gz


# Trimmomatic
cd /software/team218/JL/bin
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip

# featureCounts from the subread package
cd /software/team218/JL/app
tar zxvf subread-1.x.x.tar.gz

# UMI-tools
pip install umi_tools
