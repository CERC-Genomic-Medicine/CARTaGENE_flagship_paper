#!/bin/bash

path=“/path/to/plink_files/*bim”  # path to CAG PLINK binary format's set, each set should contain a .bed .bim .fam file.
remove_samples_file='/path/to/list_of_samples_to_remove/file.samples' # path to the list of samples to remove

for bim in ${path}; do 
        plink -bfile ${bim%.*} --remove ${remove_samples_file} --make-bed --output-chr MT --out ${bim%.*}_consented
