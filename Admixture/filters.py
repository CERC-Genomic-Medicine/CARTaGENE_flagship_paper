#!/usr/bin/python3.8

# modules used
import pandas as pd
import numpy as np
import re
import os
import csv

# Define functions
def parse_metadata(metfile, catalog, fields="file111,PHASE,ethnic,country", code='936',drop_other=True):
    '''
    Takes the metadata file and the catalog, returns the parsed data 
    '''
    #data = read_metadata(metfile, fields, code)
    data = read_metadata_large(metfile, fields)
    meta = read_catalog(catalog)

    data_filtered = code_to_value(data,meta,drop_other=drop_other)

    return data_filtered

def read_metadata(metfile, fields="file111,country", code='936'):
    '''
    Reads cartagene metadata file and parses it
    Extracts relevant fields when stated
    '''
    print("Reading metadata")
    rows=[] # List of lists, to save the dataframe rows
    with open(metfile,'r') as f:
        file = csv.reader(f,delimiter=",") # Read with csv as it already considers commas within open quotes
        header = next(file) # Get header 
        for row in file:
            line=",".join(row) # Join list to filter as string
            # Remove lines that do not contain metadata for individuals (starting with 937) and that have attention notes 
            #if re.match('936',line) and not re.search('\*\*\* Attention',line): # removed attention notes filter as it was filtering out individuals
            if re.match(code,line):
                assert len(row) == len(header) # Check that we are retrieving the correct number of columns
                rows.append(row) # Add list to rows
        
    data = pd.DataFrame(rows,columns=header)
    print("Finished reading data")
    # Extract relevant fields from metadata 
    expression=str.replace(fields,",","|")
    data_filtered = data.filter(regex=re.compile(expression, re.IGNORECASE))
    return data_filtered

def read_metadata_large(metfile, fields="file111,country"):
    expression=str.replace(fields,",","|")
    data = pd.read_csv(metfile, sep=",",dtype='object')
    data = data.filter(regex=re.compile(expression, re.IGNORECASE))
    return(data)

def read_catalog(catalog_file):
    # Read metadata info and extract the categories field in the categories sheet
    meta = pd.read_excel(catalog_file,sheet_name=['Categories'])['Categories']
    # Retain only fields that matter: Variable (name of field in metadata), Code (Code in metadata), Category (Code meaning)
    meta = meta.filter(items=["VARIABLE","CODE","CATEGORY"])

    return meta

def code_to_value(data,meta,drop_other=True):
    # Iterate over relevant fields and change code value to category value
    lastkey = ""
    for key in data:
        code_dic = dict(zip(meta.loc[meta["VARIABLE"]==key, "CODE"].copy(), meta.loc[meta["VARIABLE"]==key,"CATEGORY"].copy()))
        if len(code_dic.values()) > 1:
            #if ("COUNTRY" in key) or ("ETHNIC" in key): # If dictionary is empty, that category is not part of the metadata info, gets skipped
            try:
                assert len(code_dic) > 0
                data.loc[:, key] = [code_dic[int(item)] if pd.notna(item) else np.nan for item in data[key].values]
            except:
                print(key)
                print(code_dic)
                next
    return data

def correct_iso(dataog,phase2_file,iso_file,drop_other=True):
    '''
    Takes the parsed, filtered data and corrects the country issues in phase B by using the iso country codes.
    '''
    data_cp = dataog.copy()
    iso = {}
    with open(iso_file,"r") as file:
        for line in file:
            code=line.strip().split(" ")[0]
            country=line.strip().split(" ")[1]
            iso[code]=country

    phase2 = {}
    with open(phase2_file,"r") as file:
        for line in file:
            code=line.strip().split(" ")[0]
            country=line.strip().split(" ")[1]
            phase2[country]=code
    phase2['LIBYA'] = 'LY'

    last_key = ""
    for key in data_cp:
        if "COUNTRY" and "_OTHER" in key:        
            data_cp.loc[:,key][data_cp.phase == "PHASE B - SECOND WAVE"] = [ iso[phase2[wrong_country]] if (pd.notna(wrong_country) and wrong_country != "NO ANSWER") else wrong_country for wrong_country in data_cp[key][data_cp.phase == "PHASE B - SECOND WAVE"] ]
            data_cp[last_key][pd.notna(data_cp[key])] = data_cp[key][pd.notna(data_cp[key])]
            if drop_other:
                data_cp = data_cp.drop(str(key), axis=1) # Remove the other column
        last_key = key
    return data_cp

def filter_by_value(data,filters):
    ''' 
    Filter rows based on criteria
    '''
    filter_data = data[data.isin(filters).any(axis=1)].copy()
    return filter_data

def sample_pop(data,samplepops,samplen=30,avoid=None,strict=False):
    '''
    sample population based on criteria and sample length
    '''
    sample = pd.DataFrame()
    if avoid is None:
        avoid = []
    else:
        avoid = avoid.copy() # Copy because it passes as pointer and we do not want to change the original filters
    [avoid.append(f) for f in samplepops]
    for pop in samplepops:
        avoid2 = avoid.copy()
        avoid2.remove(pop)
        sampled = data[(data.isin([pop]).any(axis=1)) & (-data.isin(avoid2).any(axis=1))]
        if(strict & (pop in list(sampled.COUNTRY_BIRTH))):
            expression='MOTHERS|FATHERS'
            sampled.loc[:,"NO_GRANDPARENTS"] = (sampled.filter(regex=re.compile(expression, re.IGNORECASE)).isin([pop])).sum(1) # Number of grandparents in that country
            sampled = sampled[sampled.NO_GRANDPARENTS == 4]
        sampled = sampled.sample(samplen).copy()
        sampled.loc[:,"Population"] = np.repeat(pop, len(sampled.file111), axis=0)
        sample = pd.concat([sample,sampled])

    return sample

def write_keep(outfile,*data):
    '''
    Writes keep file in the format Population IID IID
    '''
    if os.path.exists(outfile):
        os.remove(outfile)
        print(outfile+" already exists, overwritting")
    for pop in data:
        pop[["Population","file111","file111"]].to_csv(outfile,mode="a", sep="\t", quoting=csv.QUOTE_NONE,header=False,index=False)
