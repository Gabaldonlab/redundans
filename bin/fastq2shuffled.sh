#!/bin/bash
# Shuffle two FastQ files
# USAGE: fastq2shuffle.sh fastq1.gz fastq2.gz > shuffled.fq

paste <(zcat $1 | paste - - - -) <(zcat $2 | paste - - - -) | tr '\t' '\n'
