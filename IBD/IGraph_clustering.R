library(igraph)       # For graph-based clustering
library(visNetwork)   # For network visualization
library(dplyr)        # For data manipulation
library(stringr)      # For string operations

# Load clustering results if needed
# load("Clustering.RData")

# Read input from IBD detection pipeline
# The file contains pairwise IBD segments detected between individuals
df <- read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/CaG_Result/all_chr.per_pair_unique.hapIBD_phaseIBD.ibd", header = TRUE)

# Filter data: Keep all segments (here 0.00 quantile means no filtering is applied)
data <- df[which(df$size_cM_hapIBD >= quantile(df$size_cM_hapIBD, 0.00, na.rm = TRUE)),]

# Create a graph where:
#   - Nodes are individuals
#   - Edges represent the total IBD segment length shared between individuals
nodes <- data.frame(id1 = unique(c(df$id1, df$id2)))
nodes$group <- 0  # Initialize all nodes with group 0
edges <- data.frame(from = data$id1, to = data$id2, width = data$size_cM_hapIBD)

# Perform Louvain clustering iteratively
j <- 1  # Iteration counter
while (j <= 4) {
  print(paste0("ITERATION: ", j))
  
  # Within each iteration, refine clustering for each existing cluster
  i <- 0  # Sub-cluster counter
  new_nodes <- data.frame()  # Dataframe to store updated node assignments
  
  for (cluster in unique(nodes$group)) {
    print(paste0("Cluster: ", cluster, " ---> ", nrow(nodes[nodes$group == cluster, ]), " samples"))
    
    # Extract individuals in the current cluster
    g <- nodes[nodes$group == cluster, ]$id1
    tmp_data <- data[(data$id1 %in% g) & (data$id2 %in% g), ]
    tmp_edges <- data.frame(from = tmp_data$id1, to = tmp_data$id2, width = tmp_data$size_cM_hapIBD)
    tmp_nodes <- data.frame(id1 = unique(c(tmp_data$id1, tmp_data$id2)))
    
    # Create graph for Louvain clustering
    tmp_graph <- graph_from_data_frame(tmp_edges, directed = FALSE)
    print("Graph constructed")
    
    # Perform Louvain community detection
    tmp_louvain <- cluster_louvain(tmp_graph)
    print("Louvain clustering completed")
    
    # Extract cluster memberships
    cluster_df <- data.frame(as.list(membership(tmp_louvain)))
    cluster_df <- as.data.frame(t(cluster_df))
    cluster_df$id2 <- rownames(cluster_df)
    id_tmp <- as.data.frame(str_split_fixed(cluster_df$id2, "[X]", 2))
    cluster_df$id1 <- as.numeric(id_tmp$V2)  # Extract numeric IDs
    print(head(cluster_df))
    
    # Assign new cluster groups
    tmp_nodes <- left_join(tmp_nodes, cluster_df, by = "id1")
    colnames(tmp_nodes)[2] <- "group"  # Rename column for clarity
    print(table(tmp_nodes$group))
    
    # Offset group numbering to ensure uniqueness across iterations
    tmp_nodes$group <- tmp_nodes$group + i
    new_nodes <- rbind(new_nodes, tmp_nodes)  # Append new assignments
    
    print(paste0("New cluster IDs: ", unique(tmp_nodes$group)))
    
    # Update the offset for the next cluster
    i <- i + length(unique(tmp_nodes$group))
    print(i)
  }
  
  # Output the updated clustering assignments
  table(new_nodes$group)
  write.table(new_nodes, paste0("v2_clusters_hap_IBD_louvain_iteration", j, ".txt"), row.names = FALSE, quote = FALSE)
  
  # Update node assignments for the next iteration
  nodes <- new_nodes
  j <- j + 1  # Move to next iteration
}
