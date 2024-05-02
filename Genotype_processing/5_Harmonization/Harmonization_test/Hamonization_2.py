#!/usr/bin/python3
import pandas as pd

Pval_threshold = 0.05
AF_ftest_4PC = "path/to/[file]" # concatenated version of the summary of the tests performed

df= pd.read_csv(AF_ftest_4PC, sep = "\s+", low_memory = False, na_values = ['None'])
df = df[~df.LRT_PVALUE.isna()]
n = len(df)
pval_threshold = threshold / n
print(n, pval_threshold) #483001 1.0351945441106748e-07
df[df.LRT_PVALUE<=(0.05/len(df.LRT_PVALUE))].drop_duplicates(subset='VARIANT').to_csv('Variants_failing_ancestry_LRTtest.txt',index=False,header=False, columns=["VARIANT"])
