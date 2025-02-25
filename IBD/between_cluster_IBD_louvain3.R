library(ggplot2)
library(dplyr)


#LOUVAIN iteration
louvain=3

#PCA
pca=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/PCA/CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all.vcf.gz.LDpruned.no_outliers.PC", header=T)
#louvain clustering
louvain1=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/Clustering/v2_clusters_hap_IBD_louvain_iteration1.txt", header=T)
louvain2=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/Clustering/v2_clusters_hap_IBD_louvain_iteration2.txt", header=T)
louvain3=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/Clustering/v2_clusters_hap_IBD_louvain_iteration3.txt", header=T)

#match group per individual
pca$group1=louvain1$group[match(pca$IID, louvain1$id1)]
pca$group2=louvain2$group[match(pca$IID, louvain2$id1)]
pca$group3=louvain3$group[match(pca$IID, louvain3$id1)]
#pca$group=cluster$group[match(pca$IID, cluster$IID)]


#IBD stats per pair
ibd=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/CaG_Result_2cM/all_chr.per_pair_unique.hapIBD_phaseIBD.ibd", header=T)


#include the group in the ibd per pair info
if(louvain == 2){
	print(paste0("NUMBER OF LOUVAIN ITERATION: ", louvain))
	ibd$group_id1=louvain2$group[match(ibd$id1, louvain2$id1)]
	ibd$group_id2=louvain2$group[match(ibd$id2, louvain2$id1)]
}
if(louvain == 3){
        print(paste0("NUMBER OF LOUVAIN ITERATION: ", louvain))
	ibd$group_id1=louvain3$group[match(ibd$id1, louvain3$id1)]
        ibd$group_id2=louvain3$group[match(ibd$id2, louvain3$id1)]
}



#-----------mean of distribution---------------

#empty vector for dataframe creation
average=c()
group_id1=c()
group_id2=c()
list_clusters=unique(as.numeric(as.character(ibd$group_id1)))

for (i in list_clusters){
	for (j in list_clusters){
		#print(i)
		#print(j)	
		
		#add the group ids of comparison
		group_id1=c(group_id1, i)
		group_id2=c(group_id2, j)
		
		#compute the average IBD sharing between individuals in the two clusters
		average=c(average, mean(ibd[which( ((ibd$group_id1 == i & ibd$group_id2 == j) | (ibd$group_id1 == j & ibd$group_id2 == i)) & !(is.na(ibd$size_cM_hapIBD))),]$size_cM_hapIBD))
	}
}

avg_df=data.frame(cluster1=group_id1, cluster2=group_id2, mean_IBD=average)

#sort
avg_df$cluster1= as.numeric(avg_df$cluster1)
avg_df$cluster2= as.numeric(avg_df$cluster2)
sorted_avg_df =avg_df[order(avg_df$cluster1, avg_df$cluster2),]
avg_df=sorted_avg_df

write.table(avg_df, paste0("between_cluster_mean_IBD_sharing_louvain", louvain, ".txt"), sep="\t", quote=F, row.names=F, col.names=T)


#mean of the distribution
mean(ibd[which(!(is.na(ibd$size_cM_hapIBD))),]$size_cM_hapIBD)
#16.45689

#-----------Plot-----------------

avg_df=read.table(paste0("between_cluster_mean_IBD_sharing_louvain",louvain,".txt"), header=T)



avg_df$cluster1= as.factor(avg_df$cluster1)
avg_df$cluster2= as.factor(avg_df$cluster2)
avg_df$mean_IBD= round(avg_df$mean_IBD, digits=2)

#mean IBD sahring between clusters
mean(avg_df[which(!(is.na(avg_df$mean_IBD))),]$mean_IBD)


p=ggplot(avg_df, aes(x = cluster1, y = cluster2, fill = mean_IBD)) +  geom_tile(color = "black") + scale_fill_gradientn(colors=c("blue","white","red")) #+  geom_text(aes(label = mean_IBD), color = "white", size = 1) +  coord_fixed()
ggsave(paste0("between_cluster_hapIBD_mean_length_louvain",louvain,".png"), p)

p=ggplot(avg_df, aes(x = cluster1, y = cluster2, fill = mean_IBD)) +  geom_tile(color = "black") +  geom_text(aes(label = as.numeric(mean_IBD)), color = "black", size = 2.5) +  coord_fixed() + scale_fill_gradientn(colors=c("blue","white","red"))
ggsave(paste0("between_cluster_hapIBD_mean_length_txt_louvain",louvain,".png"), width = 12, height=12, p)



print("plot generated")




































