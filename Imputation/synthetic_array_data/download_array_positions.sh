#!/bin/bash

# Author: Peyton McClelland <peyton.mcclelland@mail.mcgill.ca>
# Year: 2024

wget https://www.well.ox.ac.uk/~wrayner/strand/GSAMD-24v2-0_20024620_A1-b38-strand.zip
unzip GSAMD-24v2-0_20024620_A1-b38-strand.zip
awk '{print "chr"$2,$3}' OFS='\t' GSAMD-24v2-0_20024620_A1-b38.strand | sort -k1,1 -k2n,2 | uniq > Illumina_GSAv2.positions.txt

rm *.miss *.multiple *.strand
rm GSAMD-24v2-0_20024620_A1-b38-strand.zip
