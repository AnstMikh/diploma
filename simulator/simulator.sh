#!/bin/sh

while getopts ":h" option; do
   case $option in
      r) python simulator_readerror.py ancient.bam reference.fasta modern.bam modern.vcf
         exit;;
     \?) python simulator_noerror.py ancient.bam reference.fasta modern.bam modern.vcf
         exit;;
   esac
done
