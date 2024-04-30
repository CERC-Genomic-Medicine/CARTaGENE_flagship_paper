#!/usr/bin/python3
import pandas as pd
df= pd.read_table("CARTaGENE.ProPC.coord",header=0,sep="\t",low_memory=False, usecols = ['indivID'] + [f'PC{i}' for i in range(1, 21)])
df.rename(columns = {'indivID': "IID"}, inplace=True)
df.to_csv('CARTaGENE_projection.PCA', sep= '\t', index=False)
