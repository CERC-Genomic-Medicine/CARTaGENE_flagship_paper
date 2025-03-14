#!/usr/bin/env python

import argparse
import gzip
import numpy as np
from scipy.stats import norm
import pandas as pd


CHROMOSOME_CODES = set(
        [str(i) for i in range(1, 24)] + 
        [f'chr{i}'for i in range(1, 23)] + 
        ['X', 'chrX']
    )


argparser = argparse.ArgumentParser(description = 'This script compute credible sets based in P-values as described in doi:10.1038/ng.2435.')
argparser.add_argument('-g', '--gwas', metavar = 'file', dest = 'in_gwas_file', type = str, required = True, help = 'Compressed (gzip) GWAS result file in Regenie format.')
argparser.add_argument('-f', '--min-af', metavar = 'float', dest = 'min_af', type = float, required = False, default = 0.0, help = 'Threshold for the minimal alternate allele frequency (AF). Default: 0.0')
argparser.add_argument('-q', '--min-imp-quality', metavar = 'float', dest = 'min_info', type = float, required = False, default = 0.0, help = 'Threshold for the minimal imputation quality (INFO or Rsq). Defailt: 0.0') 
argparser.add_argument('-s', '--set-probability', metavar = 'float', dest = 'set_probability', type = float, required = False, default = 0.99, help = 'Probability of causal variant being included in the credible set. Default: 0.99')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_filename', type = str, required = True, help = 'Output file name.')


def read_regenie(gwas_filename, min_af, min_info):
    required_columns = ['CHROM', 'ID', 'A1FREQ', 'LOG10P']
    variant_log10p = dict()
    n_loaded = 0
    with gzip.open(gwas_filename, 'rt') as ifile:
        print(f'Scanning {gwas_filename}...')
        header = ifile.readline().rstrip().split()
        if any(c not in header for  c in required_columns):
            raise Exception('Required CHROM, ID, A1FREQ, and LOG10P columns are missing from the GWAS header.')
        chrom_idx = header.index('CHROM')
        id_idx = header.index('ID')
        a1freq_idx = header.index('A1FREQ')
        log10p_idx = header.index('LOG10P')
        if 'INFO' in header:
            info_idx = header.index('INFO')
        else:
            info_idx = None
        for n, line in enumerate(ifile, 1):
            fields = line.rstrip().split()
            if n % 1000000 == 0:
                print(f'\r{n} records scanned...', end = '', flush = True)
            chrom = fields[chrom_idx]
            if chrom not in CHROMOSOME_CODES:
                continue
            if float(fields[a1freq_idx]) < min_af:
                continue
            if info_idx is not None and float(fields[info_idx]) < min_info:
                continue
            log10p = float(fields[log10p_idx])
            variant_log10p[fields[id_idx]] = log10p
            n_loaded += 1
        print(f'\rDone. {n} records scanned. {n_loaded} records loaded.\t\t')
    return variant_log10p


# Here log10p is a negative log_10 of a pvalue.
def log10p_to_bf(log10p):
    pvalue = np.power(10.0, -1.0 * log10p)
    z = norm.ppf(pvalue / 2)
    bf = np.exp(z**2 / 2)
    return(bf)


if __name__ == '__main__':
    args = argparser.parse_args()

    variant_log10p = read_regenie(args.in_gwas_file, args.min_af, args.min_info)
    print('Computing credible sets... ', end = '', flush = True)
    bf_sum = 0
    variant_bf = dict()
    for variant_id, log10p in variant_log10p.items():
        bf = log10p_to_bf(log10p)
        bf_sum += bf
        variant_bf[variant_id] = bf

    variant_pip = []
    for variant_id, bf in variant_bf.items():
        variant_pip.append({'ID': variant_id, 'PIP': (bf / bf_sum) * 100})
    
    df = pd.DataFrame.from_records(variant_pip).sort_values(by = ['PIP'], ascending = False)
    df['PIP_CUMSUM'] = df['PIP'].cumsum()
    df = df[(df['PIP_CUMSUM'] <= args.set_probability * 100) | (df['PIP_CUMSUM'].shift(1).isna()) | (df['PIP_CUMSUM'].shift(1) < args.set_probability * 100)] # using a special trick to include the first row exceeding the threshold if the previous was still below the threshold.
    df.to_csv(args.out_filename, sep = '\t', header = True, index = False)
    print('Done.')

