import numpy as np
from scipy.stats import fisher_exact
import pandas as pd
import argparse

argparser = argparse.ArgumentParser(description = 'Develop a Python script for performing a two-sided Fisher\'s Exact Test. This script will be specifically designed to compare allele frequencies obtained from genotyping array data with those from a reference panel.')
argparser.add_argument('-ig', '--input_genotyping_array', dest = 'genotyping_array_file_path', required = True, help = 'Path to the input txt file containing following columns for the genotyping array data:[chromosome, position, reference allele, alternate allele, allele frequency].')
argparser.add_argument('-ir', '--input_reference_panel', dest = 'reference_panel_file_path', required = True, help = 'Path to the input txt file containing following columns for the reference panel data:[chromosome, position, reference allele, alternate allele, allele frequency].')
argparser.add_argument('-ng', '--sample_size_genotyping_data', dest = 'genotyping_data_sample_size', required = True, help = 'sample size of genotyping data.')
argparser.add_argument('-nr', '--sample_size_reference_data', dest = 'reference_data_sample_size', required = True, help = 'sample size of genotyping data.')
argparser.add_argument('-o', '--output', dest = 'output_file_path', required = True, help = 'Path to the output txt file containing following columns for the reference panel data:[chromosome, position, reference allele, alternate allele, reference panel allele frequency, genotyping array allele frequency, fisher exact test p_values].')

if __name__ == "__main__":
    args = argparser.parse_args()

    genotyping_array_df = pd.read_csv(args.genotyping_array_file_path, sep='\t', header = None, names = ['CHR', 'POS', 'REF', 'ALT', 'AF_GT'])
    reference_panel_df = pd.read_csv(args.reference_panel_file_path, sep='\t', header = None, names = ['CHR', 'POS', 'REF', 'ALT', 'AF_REF'])

    df_merge = genotyping_array_df.merge(reference_panel_df , on = ['CHR', 'POS', 'REF', 'ALT'])
    p_values = []
    AF_array = list(df_merge['AF_GT'])
    AF_reference = list(df_merge['AF_REF'])
    for i in range(len(AF_array)):
        table = [[args.genotyping_data_sample_size*2*(AF_array[i]), args.genotyping_data_sample_size*2*(1 - AF_array[i])], [args.reference_data_sample_size*2*AF_reference[i], args.reference_data_sample_size*2*(1 - AF_reference[i])]]
        res = fisher_exact(table, alternative='two-sided')
        p_values.append(res[1])
    df_merge['p_values'] = p_values
    df_merge.to_csv(args.output_file_path, index = None, sep='\t')
