# Python script to carry out UMAP on PC data

import argparse
import numpy as np
# import logging
import os
import sys
import time
import timeit

import umap

# Desired inputs with argparse
# -in (filename)
# -pc (# PCs)
# -nn 
# -md
# -nc
# -dist
# -gt (using genotype rather than PC data)
# -outdir 


def str2bool(v):
    # Define a str2bool function to intake the -head argument
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# Define parser
parser = argparse.ArgumentParser(description='Run UMAP on specified datasets.')

# Input variables
parser.add_argument('-dset', metavar='dataset', type=str,
                    help='Input dataset (assume PCs unless -gt specified)')
parser.add_argument('-pc', metavar='PCs', type=int,
                    default='10',
                    help='Number of top PCs to use (default 10)')
parser.add_argument('-nn', metavar='neighbours', type=int,
                    default=15,
                    help='Number of neighbours for UMAP')
parser.add_argument('-md', metavar='min_dist', type=float,
                    default=0.1,
                    help='Minimum distance for UMAP (default 0.1)')
parser.add_argument('-nc', metavar='components', type=int,
                    default=2,
                    help='Low dimensional components to project to (default 2D)')
parser.add_argument('-met', metavar='distance', type=str,
                    default='euclidean',
                    help='Type of distance metric to use (default euclidean)')
parser.add_argument('-outdir', type=str,
                    help='Output directory')
parser.add_argument('-head', metavar='headers', type=str2bool,
                    help='Indicate whether the file has headers')
parser.add_argument('-n_id', type=int,
                    default=2,
                    help='Number of ID columns')
parser.add_argument('-log', type=str,
                    help='Log directory')

args = parser.parse_args()
tstamp = time.strftime('%Y%m%d_%H%M%S',time.localtime(time.time()))

# Import arguments
dset = args.dset
pcs = args.pc
nn = args.nn
md = args.md
nc = args.nc
met = args.met.lower()
out_dir = args.outdir
has_headers = args.head
n_id = args.n_id
log_dir = args.log

# Check if important parameters have been left empty
if dset is None:
    print('ERROR: No input dataset specified.')
    sys.exit(1)
elif out_dir is None:
    print('ERROR: No output directory specified.')
    sys.exit(1)
elif has_headers is None:
    print('ERROR: Headers in file not specified.')
    sys.exit(1)

# Make sure the number of components is >= the number of PCs
if pcs < nc:
    print('ERROR: Number of PCs is less than request dimensions.')
    sys.exit(1)

# Print the parameters
param_str = dset.split('/')[-1].split('.txt')[0] + '_UMAP_PC' + str(pcs) + '_NC' + str(nc) + '_NN' \
            + str(nn) + '_MD' + str(md) + '_' + met

log_file = os.path.join(log_dir, 'log_umap_' + param_str + '_' + tstamp + '.txt')

print('Beginning import of data:', dset)
print('Parameters: ', '\n PCs:', str(pcs), '\n NC:', str(nc), '\n NN:', str(nn), '\n MD:', str(md),
      '\n Metric:', met, '\n Has headers:', str(has_headers))

# set up logging
orig_stdout = sys.stdout  # print() statements
orig_stderr = sys.stderr  # terminal statements
f = open(log_file, 'w')
sys.stdout = f
sys.stderr = f

print('Parameters: ', '\n PCs:', str(pcs), '\n NC:', str(nc), '\n NN:', str(nn), '\n MD:', str(md),
      '\n Metric:', met, '\n Has headers:', str(has_headers))

try:
    with open(dset) as data:
        data_contents = data.readlines()

        pca_data = []

        # import top PCs
        if has_headers:
            # If it has headers, skip the first row
            for pc in data_contents[1:]:
                # pca_data.append(pc.split()[2:len(pc)])
                pca_data.append(pc.split()[n_id:len(pc)])  # Skip the ID columns
        else:
            # If there are no headers, just strip the ID columns
            for pc in data_contents:
                pca_data.append(pc.split()[n_id:len(pc)])  # Skip the ID columns

        pca_data_array = np.array(pca_data).astype(np.float)

        print(pca_data_array.shape)

        # Cleanup
        del pca_data
        del pc
        del data_contents
except Exception as e:
    print(e)
    print('Error during data import')

    f.close()
    print('Error during data import')

    sys.exit(1)

# fname = dset.split('.txt')[0] + '_UMAP_PC' + str(pcs) + '_NC' + str(nc) + '_NN' + str(nn) + '_MD' + str(md) + '_' \
# + met + "_" + tstamp + ".txt"

fname = param_str + '_' + tstamp + '.txt'

# preamble for log
print()
print("Using UMAP version: " + umap.__version__)
print("Reducing to " + str(nc) + " components")
print("Using " + str(nn) + " neighbours")
print("Using minimum distance of " + str(md))
print("Using metric: " + met)
print("Using " + str(pcs) + " PCs")
print()
print("Input data shape: ", pca_data_array.shape)

try:
    # Carry out UMAP
    start = timeit.default_timer()
    umap_proj = umap.UMAP(n_components=nc, n_neighbors=nn,min_dist=md,metric=met,
                          verbose=True).fit_transform(pca_data_array[:, :pcs])
    stop = timeit.default_timer()
except Exception as e:
    print(e)
    print('Error during UMAP')

    f.close()
    print('Error during UMAP')
    sys.exit(1)

print()
print("UMAP runtime: ", stop - start)

out_file = os.path.join(out_dir, fname)

print()
print("Output file:", out_file)
print("Output data shape:", umap_proj.shape)

np.savetxt(out_file, umap_proj)

del umap_proj
del pca_data_array

# restore print statements to terminal
sys.stdout = orig_stdout
sys.stderr = orig_stderr
f.close()

# print runtime to terminal.
print("Finished successfully! UMAP runtime: ", stop - start)
