#!/usr/bin/python3

import pandas as pd
SEX= pd.read_csv("SEX_ID", sep = "\s+", low_memory = False)
PCA= pd.read_csv("CARTaGENE.eigenvec", sep = "\s+", low_memory = False)
PCA.rename(columns = {'#IID': 'IID'},inplace=True)
pd.merge(SEX,PCA,on="IID").to_csv('COVARIANT.file',sep='\t',index=False,header=True)
