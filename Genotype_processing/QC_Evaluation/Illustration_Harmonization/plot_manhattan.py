#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list
      

def manhatten(Fig, file,file_test,case,ncases,control,ncontrols):
	df= pd.read_table(file, header=0,sep="\t",low_memory=False)
	df['minuslog10pvalue'] = -np.log10(df["P"])
	df.loc[df.minuslog10pvalue>11,'minuslog10pvalue'] = 11
	df["chromosome"]=df["#CHROM"]
	df=df[~df.chromosome.isin(["Y"])]
	df.chromosome = df.chromosome.astype('category')
	df.chromosome = df.chromosome.cat.set_categories([ i for i in unique(df.chromosome)], ordered=True)
	df = df.sort_values('chromosome')
	df['ind'] = range(len(df))
	df_grouped = df.groupby('chromosome')
	ax = Fig
	colors = ['black','gray']
	x_labels = []
	x_labels_pos = []
	for num, (name, group) in enumerate(df_grouped):
		group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=colors[num % len(colors)], ax=ax)
		x_labels.append(name)
		x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
	ax.set_xticks(x_labels_pos)
	ax.set_xticklabels(x_labels)
	# set axis limits
	ax.set_xlim([0, len(df)])
	ax.set_ylim([0, 12])
	ax.axhline(y=7.30102999566 , color='r', linestyle='-')
	ax.set_xlabel('Chromosome')
	ax.set_ylabel('-'+r'$\log_{10}$' +'(P-value)')
	ax.title.set_text("Unfiltered")



def manhatten_f(Fig,file,file_test,test,case,ncases,control,ncontrols):
	df= pd.read_table(file,header=0,sep="\t",low_memory=False)
	filtered= pd.read_table(file_test,header=0,sep=" ",low_memory=False)
	df['minuslog10pvalue'] = -np.log10(df["P"])
	df.loc[df.minuslog10pvalue>11,'minuslog10pvalue'] = 11
	df=df[~df.ID.isin(filtered.VARIANT)]
	df["chromosome"]=df["#CHROM"]
	df=df[~df.chromosome.isin(["Y"])]
	df.chromosome = df.chromosome.astype('category')
	df.chromosome = df.chromosome.cat.set_categories([ i for i in unique(df.chromosome)], ordered=True)
	df = df.sort_values('chromosome')
	df['ind'] = range(len(df))
	df_grouped = df.groupby('chromosome')
	ax = Fig
	colors = ['black','gray']
	x_labels = []
	x_labels_pos = []
	for num, (name, group) in enumerate(df_grouped):
		group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=colors[num % len(colors)], ax=ax)
		x_labels.append(name)
		x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
	ax.set_xticks(x_labels_pos)
	ax.set_xticklabels(x_labels)
	# set axis limits
	ax.set_xlim([0, len(df)])
	ax.set_ylim([0, 12])
	ax.axhline(y=7.30102999566 , color='r', linestyle='-')
	ax.axes.get_yaxis().set_visible(False)
	ax.set_xlabel('Chromosome')
	ax.title.set_text("Filtered")

fig, axs = plt.subplots(1, 2, constrained_layout=True,figsize=(16,5))
manhatten(axs[0],"17Kphase1vs5Kphase1.17Kphase1vs5Kphase1.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","17K Phase 1","8,860","5K Phase 1","2,189")
manhatten_f(axs[1],"17Kphase1vs5Kphase1.17Kphase1vs5Kphase1.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","LRT_PVALUE","17K Phase 1","8,860","5K Phase 1","2,189")
plt.suptitle("17K Phase 1 (N=8,904) vs 5K Phase 1 (N=2,199)", fontsize=12)
fig.savefig('17kvs5k.png', dpi=fig.dpi)



fig, axs = plt.subplots(1, 2, constrained_layout=True,figsize=(16,5))
manhatten(axs[0],"4Kphase1vs5Kphase1.4Kphase1vs5Kphase1.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","4K Phase 1","3,549","5K Phase 1","2,189")
manhatten_f(axs[1],"4Kphase1vs5Kphase1.4Kphase1vs5Kphase1.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","LRT_PVALUE","4K Phase 1","3,549","5K Phase 1","2,189")
plt.suptitle("4K Phase 1 (N=3,540) vs 5K Phase 1 (N=2,199)", fontsize=12)
fig.savefig('4kvs5k.png', dpi=fig.dpi)


fig, axs = plt.subplots(1, 2, constrained_layout=True,figsize=(16,5))
manhatten(axs[0],"4Kphase1vsarchiphase1.4Kphase1vsarchiphase1.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","4K","3,549","Archi","1,662")
manhatten_f(axs[1],"4Kphase1vsarchiphase1.4Kphase1vsarchiphase1.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","LRT_PVALUE","4K","3,549","Archi","1,662")
plt.suptitle("4K (N=3,540) vs Archi (N=1,655)", fontsize=12)
fig.savefig('4kvsArchi.png', dpi=fig.dpi)

fig, axs = plt.subplots(1, 2, constrained_layout=True,figsize=(16,5))
manhatten(axs[0],"archiphase1vs760phase1.archiphase1vs760phase1.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","Archi","1,662","760","589")
manhatten_f(axs[1],"archiphase1vs760phase1.archiphase1vs760phase1.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","LRT_PVALUE","Archi","1,662","760","589")
plt.suptitle("Archi (N=1,655) vs 760 (N=596)", fontsize=12)
fig.savefig('Archivs760.png', dpi=fig.dpi)


fig, axs = plt.subplots(1, 2, constrained_layout=True,figsize=(16,5))
manhatten(axs[0],"17Kphase1vs17Kphase2.17Kphase1vs17Kphase2.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","17K Phase 1","8,904","17K Phase 2","6,486")
manhatten_f(axs[1],"17Kphase1vs17Kphase2.17Kphase1vs17Kphase2.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","LRT_PVALUE","17K Phase 1","8,904","17K Phase 2","6,486")
plt.suptitle("17K Phase 1 (N=8,904) vs 17K Phase 2 (N=6,486)", fontsize=12)
fig.savefig('17kvs17k.png', dpi=fig.dpi)

fig, axs = plt.subplots(1, 2, constrained_layout=True,figsize=(16,5))
manhatten(axs[0],"5Kphase1vs5Kphase2.5Kphase1vs5Kphase2.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","5K Phase 1","2,199","5K Phase 2","2,353")
manhatten_f(axs[1],"5Kphase1vs5Kphase2.5Kphase1vs5Kphase2.glm.logistic.hybrid","../Variants_failling_ancestry_LRTtest.txt","LRT_PVALUE","5K Phase 1","2,199","5K Phase 2","2,353")
plt.suptitle("5K Phase 1 (N=2,199) vs 5K Phase 2 (N=2,353)", fontsize=12)
fig.savefig('5kvs5k.png', dpi=fig.dpi)
