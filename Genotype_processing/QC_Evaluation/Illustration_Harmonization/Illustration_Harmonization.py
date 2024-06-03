#!/usr/bin/python3

### Input file
phase_file = ""                             # File detailing the recruitement phase of each sample
cag_unmerged = ""                                     # path to **FOLDER** containing the genotyping arrays PLINK binary format
variants = ""
cag_merged_unfiltered = ''   #path to CaG unfiltered merged PLINK binary format (with bed bim fam files)
unrelated_individuals = ""                  # path to file containing the list of unrelated individuals (obtained in step 5 Harmonization
array_list=["archi","760","5300","4224","17k"]
PC = '5' 																# number of PCS as covariates
threads = '5'


import pandas as pd
import numpy as np
import glob
import os
import subprocess
import matplotlib.pyplot as plt

def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list
      
# Plot Manhattan (basic)
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


# Plot Manhatten without variant in {file_test}
def manhatten_f(Fig,file,file_test,test,case,ncases,control,ncontrols):
	df= pd.read_table(file,header=0,sep="\t",low_memory=False)
	filtered= pd.read_table(file_test,header=None,sep=" ",low_memory=False)
	df['minuslog10pvalue'] = -np.log10(df["P"])
	df.loc[df.minuslog10pvalue>11,'minuslog10pvalue'] = 11
	df=df[~df.ID.isin(filtered[0])]
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

# Create Binairy + Categorical full matrix
def create_comparison_df(df_categorical, df_binary):
    # Merge the categorical and binary dataframes on 'IID'
    df_merged = pd.merge(df_categorical, df_binary, on='IID')
    
    # Create a new column that combines the category and binary value
    df_merged['group'] = df_merged['Source']+ ' phase ' + df_merged['PHASE'].astype(str)
    
    # Find all unique group combinations
    unique_groups = df_merged['group'].unique()
    print(df_merged)
    # Prepare a new DataFrame to store the results
    results_df = pd.DataFrame(index=df_merged['IID'].unique())
    
    # Iterate over each pair of groups to create columns
    for i in range(len(unique_groups)):
        for j in range(i + 1, len(unique_groups)):
            group1 = unique_groups[i]
            group2 = unique_groups[j]
            
            # Create a column for each pair
            col_name = f"{group1} vs {group2}"
            results_df[col_name] = np.nan  # Initialize with NaN
            
            # Fill the column based on group association
            mask1 = df_merged['group'] == group1
            mask2 = df_merged['group'] == group2
            
            results_df.loc[df_merged[mask1]['IID'], col_name] = 0  # Group1 association
            results_df.loc[df_merged[mask2]['IID'], col_name] = 1  # Group2 association
            
    return results_df

## Create Array Dataframe # Assume *{pattern}*shared.fam input
def aggregate_ids_with_glob(path, list_pattern):
	# Initialize an empty DataFrame to hold all the data
	aggregated_df = pd.DataFrame()
	# Construct the glob pattern
	for pattern in list_pattern :
		file_pattern = '*' + pattern + '*shared.fam'
		glob_pattern = os.path.join(path, file_pattern)
		# Use glob to find files matching the pattern
		for file_path in glob.glob(glob_pattern):
			file_name = os.path.basename(file_path)
			df = pd.read_csv(file_path,header=None,sep='\s+')
			df.columns=['FID','IID','x','y','Sex','NA']
			temp_df = pd.DataFrame()
			temp_df['IID'] = df['IID']
			temp_df['Source'] = pattern
			print(temp_df)
			aggregated_df = pd.concat([aggregated_df, temp_df], ignore_index=True)
	return aggregated_df


if __name__ == '__main__':
	# READ files
	## Step 1
	phases = pd.read_csv(phase_file,usecols=['IID','PHASE'])
	phase_list = phases['PHASE'].unique()
	unrelated=pd.read_csv('/scratch/vct/CARTaGENE_v1.1/Mooser_433651_genotypes_hg38_Harmonized_dataset/CARTaGENE_hg38_shared_unrelated.king.cutoff.in.id')

	# Read all Arrays
	Array =  aggregate_ids_with_glob(cag_unmerged, array_list)
	matrix = create_comparison_df(Array,phases)

	matrix=matrix.loc[unrelated.iloc[:,0],]
	filter=[matrix[i].sum()>20 and sum(matrix[i]==0)>20 for i in matrix.columns]
	final_matrix=matrix.loc[:,filter]
	final_matrix.columns= [i.replace(' ','_') for i in final_matrix.columns]
	columns_filtered = [col.split("_")[0]==col.split("_")[4] or col.split("_")[2]==col.split("_")[6] for col in final_matrix.columns]
	final_matrix = final_matrix.loc[:,columns_filtered]
	final_matrix.index.name = "#FID"
	final_matrix.insert(loc=0, column='IID', value=final_matrix.index)
	final_matrix.insert(loc=0, column='FID', value=final_matrix.index)
	final_matrix.to_csv('Array_Phase.pheno', sep = '\t',na_rep='NA', index=False)
	### Create PCA covariates
	subprocess.run(f"plink2 --bfile {cag_merged_unfiltered[:-4]} --keep {unrelated_individuals} --pca {PC} --out PCA_covar --threads {threads}", shell=True, check=True, stderr=subprocess.PIPE, text=True)
	### Association Testing
	subprocess.run(f"plink2 --bfile ' + cag_merged_unfiltered[:-4]} --keep {Unrelated_individuals} --pheno Array_Phase.pheno --1 --covar PCA_covar.eigenvec --glm hide-covar --threads {threads}", shell=True, check=True, stderr=subprocess.PIPE, text=True)
	
	### Process each association file into a comparison figure
	file_pattern = "plink2.*.glm.logistic.hybrid"
	glob_pattern = os.path.join(os.getcwd(), file_pattern)
	for file_path in glob.glob(glob_pattern):
		comparison = file_path.replace(os.getcwd() + "/plink2.","").replace(".glm.logistic.hybrid","")
		Array1 = comparison.split("_")[0]
		Array2 = comparison.split("_")[4]
		Phase1 = comparison.split("_")[2]
		Phase2 = comparison.split("_")[6]
		n1 = int(sum(final_matrix[comparison]==0)) ## N in comparison
		n2 = int(sum(final_matrix[comparison]==1)) ## N in comparison
		fig, axs = plt.subplots(1, 2, constrained_layout=True,figsize=(16,5))
		manhatten(axs[0],file_path,variants,Array1,f"{n1:,}", Array2 ,f"{n2:,}")
		manhatten_f(axs[1],file_path,variants,"LRT_PVALUE", Array1 ,f"{n1:,}", Array2 ,f"{n2:,}")
		plt.suptitle(f"{Array1} (N={n1:,}) vs {Array2} (N={n2:,})", fontsize=12)
		fig.savefig(f"{comparison}.png", dpi=fig.dpi)
