#!/usr/bin/python3

import pandas as pd
import glob
import os

directory = 'path/to/genotype_files/'
pattern = '*.fam,'  # Adjust the pattern according to your file naming convention

# Construct the full path with the pattern
search_pattern = os.path.join(directory, pattern)

# Find all files matching the pattern
files = glob.glob(search_pattern)

df = pd.concat([
pd.read_csv(fam, sep = "\s+",header=None, low_memory = False,usecols = [0,1], names = ['FID','IID']) for fam in files])
df.groupby(['FID', 'IID']).filter(lambda x: len(x) > 1).drop_duplicates().to_csv('exclude_dup_sample_test.txt',index=False,header=True, sep='\t')
