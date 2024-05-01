#!/bin/bash

declare -a paths=("/path/to/plink_files/[prefix]" ...)
# list of path to CAG PLINK binary format's set, each set should contain a .bed .bim .fam file.
remove_samples_file='/path/to/list_of_samples_to_remove/file.samples'
# path to the list of samples to remove

for set in ${paths[@]}; do 
        plink -bfile ${set} --remove ${remove_samples_file} --make-bed --output-chr MT --out ${set}_consented
