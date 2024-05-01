import pandas as pd
import glob
import os

### Step 1 - find all shared variants

# Define the directory and pattern for file matching
directory = 'path/to/genotype_files'
pattern = '*.bim'  

# Construct the full path with the pattern
search_pattern = os.path.join(directory, pattern)

# Find all files matching the pattern
files = glob.glob(search_pattern)

# Load the first file into a DataFrame to start the comparison
shared_elements = pd.read_csv(files[0], delimiter='\t', header=None, usecols=[1], squeeze=True).unique()

# Process each subsequent file and intersect the sets
for file_path in files[1:]:
    current_elements = pd.read_csv(file_path, delimiter='\t', header=None, usecols=[1], squeeze=True).unique()
    shared_elements = set(shared_elements).intersection(current_elements)

# Output the results
output_file = 'shared.txt'
with open(output_file, 'w') as f:
    for element in shared_elements:
        f.write(element + "\n")
