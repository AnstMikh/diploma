#!/bin/sh

while getopts ":h" option; do
   case $option in
      r) python simulator.py ancient.bam reference.fasta modern.bam modern.vcf
         exit;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done
