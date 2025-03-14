import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from fuzzywuzzy import process
# Read the region dictionary file with space delimiter
region_dict = pd.read_csv("/Users/amgaricia/Downloads/dictionary_regionsISGen.txt", delim_whitespace=True, header=None, names=["RegionName", "ID"])

# Display the first few rows to verify the data
print(region_dict.head())
genetic_data = pd.read_csv("/Users/amgaricia/Downloads/finalfreqSPG7", header=None, names=["ID", "2.5%", "mean", "97.5%"])
# Remove the first row
genetic_data = genetic_data.iloc[1:]

# Display the modified DataFrame
genetic_data

genetic_data["ID"] = genetic_data["ID"].astype(str)
region_dict["ID"] = region_dict["ID"].astype(str)
merged_data2 = pd.merge(genetic_data, region_dict, left_on="ID", right_on="ID")
# Convert the ID column to integer for plotting
merged_data2["ID"] = merged_data2["ID"].astype(int)


# Sort the merged DataFrame by the mean frequency in descending order
sorted_data = merged_data2.sort_values(by="mean", ascending=False)

# Load the GeoJSON file with geographic data
geo_data = gpd.read_file('BALSAC.geojson')

# Load the CSV file with frequency values
freq_data = merged_data2
    
# List of regions in geographic data
geo_regions = geo_data['Name'].unique()

# Function to match regions using fuzzy matching
def match_region(region_name, choices, threshold=80):
    match, score = process.extractOne(region_name, choices)
    if score >= threshold:
        return match
    else:
        return None

# Apply fuzzy matching to freq_data regions
freq_data['matched_region'] = freq_data['RegionName'].apply(match_region, args=(geo_regions,))

# Merge frequency data with geographic data using the matched regions
merged_data = geo_data.merge(freq_data, left_on='Name', right_on='matched_region', how='left')
# Locate the row with "Québec (agglomération)" in the Name column and extract its mean value
mean_value_agglomeration = merged_data.loc[merged_data['Name'] == 'Québec (agglomération)', 'mean'].iloc[0]

# Locate the row with "Québec (région de)" in the Name column and update its mean value
merged_data['mean'] = merged_data['mean'].fillna(mean_value_agglomeration)

# Now, the mean value for "Québec (région de)" should include the mean value of "Québec (agglomération)"
merged_data

merged_data['mean'] = merged_data['mean'].fillna(mean_value_agglomeration)
merged_data["mean"]

# Convert 'mean' column to numeric
merged_data['mean'] = pd.to_numeric(merged_data['mean'], errors='coerce')

# Check the data type of the 'mean' column after conversion
print(merged_data['mean'].dtype)

# Example DataFrame (replace with your own)
# merged_data = gpd.read_file('your_shapefile.shp')

# Creating a plot with higher DPI for better resolution
fig, ax = plt.subplots(1, 1, figsize=(30, 10), dpi=600)  # Increase the height of the figure

# Plotting with vivid colors and better borders
merged_data.plot(
    column="mean",
    cmap='viridis',  # Using 'viridis' colormap for vivid colors
    legend=False,
    linewidth=0.3,  # Thicker borders for better visibility
    edgecolor='black',  # Black borders for strong contrast
    ax=ax,
    alpha=1,
    legend_kwds={
        'label': "Mean Values",
        'orientation': "horizontal",
        'shrink': 0.5,  # Adjust the size of the legend bar
        'pad': 0.05,    # Adjust the padding around the legend bar
        'aspect': 40,   # Adjust the aspect ratio of the legend bar
        'anchor': (0.4, 0.5)  # Center the legend bar
    }
)

# Remove X and Y axes
ax.set_axis_off()

# Retrieve the colorbar and set the font size for the tick labels
cbar = ax.get_figure().colorbar(ax.collections[0], ax=ax, orientation="horizontal", shrink=0.5, pad=0.05, aspect=40, anchor=(0.4, 0.5))
cbar.ax.tick_params(labelsize=45, pad=15)  # Set the font size for colorbar tick labels
#cbar.set_label("Carrier frequency", fontsize=90)  # Set the font size for the colorbar label

# Save the figure with high resolution
plt.savefig('SPG7.png', dpi=600, bbox_inches='tight')  # Save with high DPI
    
plt.show()
