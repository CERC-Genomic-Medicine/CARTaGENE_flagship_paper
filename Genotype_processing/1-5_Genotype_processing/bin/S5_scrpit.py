#!/usr/bin/env python3
import subprocess
import pandas as pd
import argparse

argparser = argparse.ArgumentParser(description = 'Filter script for Genotype processing')
argparser.add_argument('-i', metavar = 'name', dest = 'input', type = str, required = True, help = 'Result from test script')
argparser.add_argument('-p', metavar = 'name', dest = 'threshold', type = str, required = True, help = 'unajusted pval threshold')




if __name__ == '__main__':
    args = argparser.parse_args()
    df= pd.read_csv(args.input, sep = "\t", low_memory = False, na_values = ['None'])
    df = df[~df.LRT_PVALUE.isna()]
    n = len(df)
    pval_threshold = pd.to_numeric(args.threshold) / n
    df.loc[df.LRT_PVALUE<=pval_threshold,:].drop_duplicates(subset='VARIANT').to_csv('Variants_failing_ancestry_LRTtest.txt',index=False,header=False, columns=["VARIANT"])


