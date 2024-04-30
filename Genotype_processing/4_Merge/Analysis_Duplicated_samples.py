#!/usr/bin/python3

import pandas as pd
importnumpy as np
path=''
df = pd.concat([
pd.read_csv(path + "5300_hg38.fam", sep = "\s+",header=None, low_memory = False,usecols = [0,1], names = ['FID','IID']),
pd.read_csv(path +"4224_hg38.fam", sep = "\s+", header=None,low_memory = False,usecols = [0,1], names = ['FID','IID']),
pd.read_csv(path +"archi_hg38.fam", sep = "\s+", header=None,low_memory = False,usecols = [0,1], names = ['FID','IID']),
pd.read_csv(path +"17k_hg38.fam", sep = "\s+", header=None,low_memory = False,usecols = [0,1], names = ['FID','IID']),
pd.read_csv(path +"760_hg38.fam", sep = "\s+", header=None,low_memory = False,usecols = [0,1], names = ['FID','IID'])])
df.groupby(['FID', 'IID']).filter(lambda x: len(x) > 1).drop_duplicates().to_csv('exclude_dup_sample_test.txt',index=False,header=True, sep='\t')
