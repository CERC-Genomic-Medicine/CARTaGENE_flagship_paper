library(ggplot2)
library(dplyr)

#--------------------------STATS within cluster------------------------------
#Cluster after FST pruning
cluster=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/Clustering/FST_pruned_clusters.txt", header=F)
names(cluster)=c("IID", "group")


#----------histogram nb sample per cluster-------

p=ggplot(cluster, aes(x=group)) + geom_histogram(fill="lightblue", color="black", binwidth=1) + xlab("Clusters") + ylab("Sample size") + theme(legend.position = "none") + scale_x_continuous(breaks=seq(1,56,1)) + theme(text = element_text(size = 20),axis.text.x = element_text(angle = -55, vjust = 0.5, hjust=1)) 
ggsave("histogram_nb_sample_per_cluster_louvain3_FST.png", width=16, height=11, p)

#------------------------------------------------

#PCA
pca=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/PCA/CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all.vcf.gz.LDpruned.no_outliers.PC", header=T)
#louvain clustering
louvain1=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/Clustering/clusters_hap_IBD_louvain_iteration1.txt", header=T)
louvain2=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/Clustering/clusters_hap_IBD_louvain_iteration2.txt", header=T)
louvain3=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/Clustering/clusters_hap_IBD_louvain_iteration3.txt", header=T)

#match group per individual
pca$group1=louvain1$group[match(pca$IID, louvain1$id1)]
pca$group2=louvain2$group[match(pca$IID, louvain2$id1)]
pca$group3=louvain3$group[match(pca$IID, louvain3$id1)]
pca$group=cluster$group[match(pca$IID, cluster$IID)]


#IBD stats per pair
ibd=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping/CaG_Result_2cM/all_chr.per_pair_unique.hapIBD_phaseIBD.ibd", header=T)


#include the group in the ibd per pair info
ibd$group_id1=cluster$group[match(ibd$id1, cluster$IID)]
ibd$group_id2=cluster$group[match(ibd$id2, cluster$IID)]

#only keep pairs within the same cluster
within_ibd=ibd[which(ibd$group_id1 == ibd$group_id2),]

within_ibd$group_id1=as.factor(within_ibd$group_id1)



#-----------mean of distribution---------------

average=c()
group_id=c()
for (i in unique(within_ibd$group_id1)){
	print(i)
	group_id=c(group_id, i)

	average=c(average, mean(within_ibd[which(within_ibd$group_id1 == i & !(is.na(within_ibd$size_cM_hapIBD))),]$size_cM_hapIBD))
}

avg_df=data.frame(group=group_id, mean_IBD=average)

write.table(avg_df, "within_cluster_mean_IBD_sharing.txt", sep="\t", quote=F, row.names=F, col.names=T)
avg_df=read.table("within_cluster_mean_IBD_sharing.txt", header=T)

df2 <- avg_df[order(avg_df$mean_IBD),]

#mean of the distribution
#mean(within_ibd[which(!(is.na(within_ibd$size_cM_hapIBD))),]$size_cM_hapIBD)
#16.45689
#without top3
#mean(within_ibd[which(within_ibd$group_id1 != 1 & within_ibd$group_id1 != 16 & within_ibd$group_id1 != 30 & !(is.na(within_ibd$size_cM_hapIBD))),]$size_cM_hapIBD)
#10.25571

#-----------Plot-----------------

#reorder the boxes
within_ibd$group_id1=factor(within_ibd$group_id1, levels=df2$group)


#compute sample size in each cluster
size_sample=data.frame(table(pca$group))
names(size_sample)=c("group_id1","size")



give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x)))
  # experiment with the multiplier to find the perfect position
}



color_custom=c("#d8d8d8", "#d0d0d0", "#c9c9c9", "#B454DA", "#bababa", "#b2b2b2", "#ababab", "#a3a3a3", "#9c9c9c", "#949494", "#8d8d8d", "#858585", "#7e7e7e", "#767676", "#6f6f6f", "#676767", "#606060", "#585858" , "#F84C3C", "#019FF4", "#4EEA70")


p<-ggplot(within_ibd, aes(x=size_cM_hapIBD, y=group_id1, fill=group_id1)) + geom_boxplot(outlier.shape=NA) +stat_summary(fun=mean, geom="point") + xlab("IBD sharing (cM)") + ylab("Clusters")+ theme(text = element_text(size = 20)) + scale_fill_manual(values=color_custom) +  theme(legend.position = "none") + xlim(c(0,100)) + theme_classic()  +  theme(legend.position = "none") + theme(text = element_text(size = 20))#+ coord_flip() 
ggsave("within_cluster_hapIBD_mean_length_classic.png", p)

p<-ggplot(within_ibd, aes(x=size_cM_hapIBD, y=group_id1, fill=group_id1)) + geom_boxplot(outlier.shape=NA) +stat_summary(fun=mean, geom="point") + xlab("IBD sharing (cM)") + ylab("Clusters") + theme(legend.position = "none") + xlim(c(0,100)) + theme(text = element_text(size = 20)) + scale_fill_manual(values=color_custom) #+ coord_flip() 
ggsave("within_cluster_hapIBD_mean_length.png", p)




print("plot generated")














































