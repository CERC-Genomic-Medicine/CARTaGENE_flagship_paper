#!/bin/bash

# Script
## Declare Variables

path=“path/to/file/*.bim”  # path to CAG PLINK binary format's bim file, each set should contain a .bed .bim .fam file.
output_name="CaG"
threads=5

mkdir Merger; cd Merger

# Find List of overlapping individuals
awk '{print $1}' $(ls ${filenames}) | sort | uniq -c | awk '$1 > 1 n {print $2,$2}' > Samples_to_exclude.txt
sed -i '1i#FID IID' Samples_to_exclude.txt
# Order files per nb of variants

list_files_order=$(wc -l ${path} | sort -n | awk '{print $2}' | head -n -1)

# Iter over list to remove duplicate individual (removing individuals form list after each)
for file in $list_files_order
do
  filename=${file##*/}
  plink --bfile ${file%.*} --remove Samples_to_exclude.txt --keep-allele-order --make-bed --output-chr chrMT --threads ${threads} --out ${filename%.*}_temp
  grep -v ${file%.*}.fam Samples_to_exclude.txt > tmp ## Remove any duplicate just removed from list
  mv tmp Samples_to_exclude.txt
  if cmp -s ${file%.*}.fam ${filename%.*}_temp.fam; then
    rm *temp*
  else
    mv ${file%.*}_temp.fam ${file%.*}.fam
    mv ${file%.*}_temp.bed ${file%.*}.bed
    mv ${file%.*}_temp.bim ${file%.*}.bim
  fi
  if [$(wc -l exclude_dup_sample.txt ) -eq 1]; then # If there is no longer any duplicates (1st line FID IID)
    break
  fi
done

## Set the working files to current directory (updated files)
for i in "${!list_files_order[@]}"; do
    # This loop removes everything after the first dot encountered
    filenames[i]="${list_files_order[i]##*/}"
done


# Find variant present in all arrays
nfile=$(ls -lh ${filenames} | wc -l)
awk '{print $2}' $(ls ${filenames}) | sort | uniq -c | awk -v n="$nfile" '$1 == n {print $2}' | grep -v 'chrY' - > shared.txt
# Extract variant
for bim in ${filenames} ; do
  plink --bfile ${bim%.*} --extract shared.txt  --keep-allele-order --make-bed --output-chr chrMT --out ${bim%.*}_shared --threads ${threads}
done

files=($(ls ${filenames}))
base_file="${files[0]%.*}_shared"
unset files[0] # pop out the first file

output="merge_temp"

# Loop through the remaining files and merge them one by one
for file in "${files[@]}"
do
    # Determine the next output name
    next_output="${output}_next"

    # Run plink to merge the current base with the next file
    plink --bfile $base_file --keep-allele-order --output-chr chrMT --bmerge ${file%.*}_shared --threads ${threads} --out $next_output

    # Update base_file to the new merged file for the next iteration
    base_file=$next_output
    output=$next_output
done

mv ${output}.bed ${output_name}.bed
mv ${output}.bim ${output_name}.bim
mv ${output}.fam ${output_name}.fam

rm *temp*
