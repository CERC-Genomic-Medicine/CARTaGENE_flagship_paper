# This program is inteded to select and parameterize the HDBSCAN script

# The HDBSCAN script is already set up to create the clusters but needs
# input parameters (as well as directory information)

# Rather than modify that already-functional program, this one will call it
import os

num_dim = 3
# We will loop this over every N+ dimensional projection that we have.
# This is specified in the loop

# where your UMAP projections are stored
proj_dir = 'cartagene/projections/hla_removed_ld_thinned_0.1'

# where to output the cluster labels
out_dir = 'cartagene/hdbscan_clusters/hla_removed_ld_thinned_0.1'

# where to output the cluster logs
log_dir = 'cartagene/logs'

# location of the python script that does the clustering
pgm_path = 'cartagene/code/clustering/general_hdbscan_clustering.py'

head = 'F'

eps = 0.5

# minimum numbers of points to use in clusters
min_points = [25, 50, 100]
n = 1

# Find the number of total iterations:
number_iterations = 0
for min_point in min_points:
    for fname in os.listdir(proj_dir):
        if 'NC' in fname:
            if int(fname.split('NC')[1][0]) >= num_dim:
                number_iterations += 1

cur_iteration = 0

for min_point in min_points:
    for fname in os.listdir(proj_dir):
        if 'NC' in fname:
            if int(fname.split('NC')[1][0]) >= num_dim:
                cur_iteration += 1
                dset = os.path.join(proj_dir, fname)

                command_str = '''\
python {_pgm_path} \
-dset {_dset} \
-min_points {_min_points} \
-eps {_eps} \
-n {_nid} \
-outdir {_outdir} \
-head {_head} \
-log {_logdir} \
'''.format(_pgm_path = pgm_path, _dset=dset, _min_points=min_point,
                    _eps=eps, _nid=n, _outdir=out_dir, _head=head, _logdir=log_dir)

                print('Beginning iteration', cur_iteration, 'of', number_iterations)
                print()
                os.system(command_str)
