# Load cluster information
cluster <- read.table("group3_cluster_info.txt", header = FALSE)
names(cluster) <- c("FID", "IID", "group")

# Load cluster size information
nb_cluster <- read.table("group3_cluster_nb_per_cluster.txt", header = FALSE)
names(nb_cluster) <- c("cluster", "nb")
small_clusters <- nb_cluster[which(nb_cluster$nb <= 10), ]  # Identify small clusters

# Load pairwise Fst data
fst <- read.table("pairwise_FST_CAG_cluster3.fst.summary", header = FALSE)
names(fst) <- c("group1", "group2", "Fst")

# Sort Fst values for better analysis
df2 <- fst[order(fst$Fst), ]

#---------------- Merge clusters based on Fst values ----------------#

# Identify clusters to merge (Fst < 0.0005)
to_merge <- fst[which(fst$Fst < 0.0005), ]

for (lines in 1:nrow(to_merge)) {
    from <- to_merge[lines, ]$group2    # Source cluster
    to <- to_merge[lines, ]$group1      # Target cluster
    
    print(paste0("Merging from: ", from, " to: ", to))
    
    if (nrow(cluster[which(cluster$group == from), ]) > 0) {
        print("Merging confirmed")
        cluster[which(cluster$group == from), ]$group <- to  # Assign new cluster
    }
}

# Output number of remaining unique clusters
length(unique(cluster$group))

# Write pruned cluster assignments to file
write.table(data.frame(cluster$IID, cluster$group), "FST_pruned_clusters3.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#---------------- Plot Fst Distribution ----------------#
library(ggplot2)

# Plot density distribution with thresholds
p <- ggplot(fst, aes(x = Fst)) + 
    geom_density() + 
    geom_vline(xintercept = 0.0005, color = "red") + 
    geom_vline(xintercept = 0.001, color = "blue")

ggsave("pairwise_Fst_CAG_cluster3_distribution.png", p)

#---------------- Generate Fst Heatmap ----------------#
p <- ggplot(fst, aes(x = group2, y = group1, fill = Fst)) +  
    geom_tile(color = "black") + 
    scale_fill_gradient2(low = "#075AFF", mid = "#FFFFCC", high = "#FF0000") +  
    coord_fixed() + 
    scale_y_continuous(breaks = seq(1, 56, 1)) + 
    scale_x_continuous(breaks = seq(1, 56, 1)) +
    xlab("Cluster") + ylab("Cluster") + 
    theme_classic()

ggsave("pairwise_Fst_CAG_cluster3_heatmap.png", width = 11, height = 11, p)
