library(ggplot2)
library(wordcloud)
library(dplyr)

#library("/home/justinp/R/x86_64-pc-linux-gnu-library/4.3/ggplot2")



#-------------------------INPUT--------------------------------------------------------pca$ETHNICITY_FR=metadata$ETHNICITY_FR[match(pca$IID, metadata$PROJECT_CODE)]


#pca=read.table("CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all.LDpruned.PC", header=T)
pca=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version/PCA/CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all.vcf.gz.LDpruned.no_outliers.PC", header=T)



#louvain clustering
louvain1=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version/Clustering/v2_clusters_hap_IBD_louvain_iteration1.txt", header=T)
louvain2=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version/Clustering/v2_clusters_hap_IBD_louvain_iteration2.txt", header=T)
louvain3=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version/Clustering/v2_clusters_hap_IBD_louvain_iteration3.txt", header=T)

#metadata country
metadata = read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version/PCA/phenotypes/QCed_phenotypes_v3.txt", header=T)

#Cluster after FST pruning
cluster=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version/Clustering/FST_pruned_clusters.txt", header=F)
names(cluster)=c("IID", "group")

# Cluster renamed
conv_cluster=read.table("FST_pruned_clusters_with_conversion.txt", header=T)

#alex data
alex = read.csv("/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version/Clustering_3cM/genotype_clusters_cob_umap_coords.csv",  header=T, na.strings=c(-1,''))
correspondance=read.csv("/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version/Clustering_3cM/Mooser_433651_imputation_codes_20221124.csv", header=T)
alex$project_code=correspondance$project_code[match(alex$file111, correspondance$imputation_code)]
cluster$umap1=alex$UMAP2D_1[match(cluster$IID, alex$project_code)]
cluster$umap2=alex$UMAP2D_2[match(cluster$IID, alex$project_code)]
cluster$cluster_alex=alex$Cluster[match(cluster$IID, alex$project_code)]



#match group per individual
pca$group1=louvain1$group[match(pca$IID, louvain1$id1)]
pca$group2=louvain2$group[match(pca$IID, louvain2$id1)]
pca$group3=louvain3$group[match(pca$IID, louvain3$id1)]
pca$group=conv_cluster$New_Cluster[match(pca$IID,conv_cluster$IID)]
pca$group_alex=alex$Cluster[match(pca$IID, alex$project_code)]

#match umap coordinates
pca$umap1=alex$UMAP2D_1[match(pca$IID, alex$project_code)]
pca$umap2=alex$UMAP2D_2[match(pca$IID, alex$project_code)]


#match country of  birth for all individuals
pca$country=metadata$COUNTRY_BIRTH[match(pca$IID, metadata$PROJECT_CODE)]
pca$country_mothers_mother=metadata$MOTHERS_MOTHER_COUNTRY_BIRTH[match(pca$IID, metadata$PROJECT_CODE)]
pca$country_mothers_father=metadata$MOTHERS_FATHER_COUNTRY_BIRTH[match(pca$IID, metadata$PROJECT_CODE)]
pca$country_fathers_mother=metadata$FATHERS_MOTHER_COUNTRY_BIRTH[match(pca$IID, metadata$PROJECT_CODE)]
pca$country_fathers_father=metadata$FATHERS_FATHER_COUNTRY_BIRTH[match(pca$IID, metadata$PROJECT_CODE)]
pca$CREGION=metadata$CREGION[match(pca$IID, metadata$PROJECT_CODE)]
pca$PROVINCE_BIRTH=metadata$PROVINCE_BIRTH[match(pca$IID, metadata$PROJECT_CODE)]
#correct for spelling errors
pca[which(pca$PROVINCE_BIRTH == "NEW_BRUNSWIC"),]$PROVINCE_BIRTH="NEW_BRUNSWICK"
pca[which(pca$PROVINCE_BIRTH == "PRINCE_EDWAR"),]$PROVINCE_BIRTH="PRINCE_EDWARD_ISLAND"
#pca$ETHNICITY_CORRECT_ME=metadata$ETHNICITY_CORRECT_ME[match(pca$IID, metadata$PROJECT_CODE)]
pca$ETHNICITY_ME=metadata$ETHNICITY_ME[match(pca$IID, metadata$PROJECT_CODE)]
pca$ETHNICITY=metadata$ETHNICITY[match(pca$IID, metadata$PROJECT_CODE)]
pca$ETHNICITY_FR=metadata$ETHNICITY_FR[match(pca$IID, metadata$PROJECT_CODE)]
pca$ETHNICITY_MOTHER=metadata$ETHNICITY_MOTHER[match(pca$IID, metadata$PROJECT_CODE)]
pca$ETHNICITY_FATHER=metadata$ETHNICITY_FATHER[match(pca$IID, metadata$PROJECT_CODE)]




#------------ACADIANS CLUSTER FLAGSHIP-------------------#

clusters=pca[which(!(is.na(pca$group))),]


cluster8 = clusters[which(clusters$group == 8),]
all_but_8 = clusters[which(clusters$group != 8),]

table(cluster8$PROVINCE_BIRTH)
table(all_but_8$PROVINCE_BIRTH)

# Overall atlantic provinces
nrow(clusters[which(clusters$PROVINCE_BIRTH == "NEW_BRUNSWICK" | clusters$PROVINCE_BIRTH == "NOVA_SCOTIA" | clusters$PROVINCE_BIRTH == "PRINCE_EDWARD_ISLAND"),])
# FC atlantic provinces 
nrow(clusters[which((clusters$PROVINCE_BIRTH == "NEW_BRUNSWICK" | clusters$PROVINCE_BIRTH == "NOVA_SCOTIA" | clusters$PROVINCE_BIRTH == "PRINCE_EDWARD_ISLAND") & (clusters$ETHNICITY_FR == "FRENCH_CANADIAN")),])
nrow(cluster8[which((cluster8$PROVINCE_BIRTH == "NEW_BRUNSWICK" | cluster8$PROVINCE_BIRTH == "NOVA_SCOTIA" | cluster8$PROVINCE_BIRTH == "PRINCE_EDWARD_ISLAND") & (cluster8$ETHNICITY_FR == "FRENCH_CANADIAN")),])


#-------------------------Histogram groups-------------------------------------------------------

p<-ggplot(pca, aes(x=group)) + geom_histogram(bins=length(unique(pca$group)),fill="white",color="black")+ geom_hline(yintercept=20, color="red")
ggsave("number_of_sample_per_group_louvain_3cM.png", p)

p<-ggplot(pca, aes(x=group2)) + geom_histogram(bins=length(unique(pca$group2)),fill="white",color="black")+ geom_hline(yintercept=20, color="red")
ggsave("number_of_sample_per_group2_louvain_3cM.png", p)

p<-ggplot(pca, aes(x=group3)) + geom_histogram(bins=length(unique(pca$group3)))+ geom_hline(yintercept=20, color="red")
ggsave("number_of_sample_per_group3_louvain_3cM.png", p)



#remove clusters with less than 10 individuals
group3_size=as.data.frame(table(pca$group3))
list_small_clusters<-group3_size[which(group3_size$Freq < 10),]$Var1
cluster_exclude=pca[which(pca$group3 %in% list_small_clusters),]$IID
write.table(data.frame(FID=cluster_exclude, IID=cluster_exclude), "/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version/Clustering/group3_small_cluster_to_exclude.txt", sep="\t", quote=F, row.names=F, col.names=F)

tmp=pca[which(!(pca$group3 %in% list_small_clusters)), ]
write.table(data.frame(FID=tmp$IID, IID=tmp$IID, Cluster=tmp$group3), "/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version/Clustering/group3_cluster_info.txt", sep="\t", quote=F, row.names=F, col.names=F)



#-------------Number of sample per cluster-------------------#
#cluster from louvain 3rd iteration with FST combining
clusters=pca[which(!(is.na(pca$group))),]

for (cluster in unique(clusters$group)){
	print(paste0("Number of samples in ", cluster, ": ",nrow(clusters[which(clusters$group == cluster),])))
}


#----------------------TAG for cluster---------------------------------

#cluster from louvain 3rd iteration with FST combining
clusters=pca[which(!(is.na(pca$group))),]

for (cluster in unique(clusters$group)){
	tmp=clusters[which(clusters$group == cluster),]
	#print(head(tmp))
	nb_sample=nrow(tmp)
	print(cluster)


	#----------Capital letter first and lowercase then
	# Custom function to capitalize the first letter of each word
	capitalize_first_letter <- function(x) {
	  sapply(x, function(y) {
	    if(is.na(y)) return(NA)
	    s <- strsplit(y, " ")[[1]]
	    paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2)),
	          sep = "", collapse = " ")
	  })
	}

	# Apply the custom function to all character columns in the dataframe
	tmp <- tmp %>%
	  mutate(across(where(is.character), capitalize_first_letter))



	#----------DF for each field-----------#

	
	dfcountry <- as.data.frame(table(tmp$country))
	dfcountry$Proportion <- dfcountry$Freq / sum(dfcountry$Freq)
	
	dfcregion <- as.data.frame(table(tmp$CREGION))
	dfcregion$Proportion <- dfcregion$Freq / sum(dfcregion$Freq)
        dfprovince <- as.data.frame(table(tmp$PROVINCE_BIRTH))
	dfprovince$Proportion <- dfprovince$Freq / sum(dfprovince$Freq)
	dfethnicityfr <- as.data.frame(table(tmp$ETHNICITY_FR))
        dfethnicityme <- as.data.frame(table(tmp$ETHNICITY_ME))
	dfethnicityme$Proportion <- dfethnicityme$Freq / sum(dfethnicityme$Freq)
	dfethnicity <- as.data.frame(table(tmp$ETHNICITY))
	dfethnicity$Proportion <- dfethnicity$Freq / sum(dfethnicity$Freq)


	#---------bar plot---------------------#

	if(F){
	p=ggplot(dfcountry, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity")+ coord_flip() + ggtitle(paste0("cluster: ",cluster, " (", nb_sample," samples)")) + xlab("Country of birth")+ ylab("Count")
	ggsave(paste0("country_birth_cluster",cluster,"_3cM.png"), p)

	p=ggplot(dfcregion, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity")+ coord_flip() + ggtitle(paste0("cluster: ",cluster, " (", nb_sample," samples)"))+ xlab("Recruitment region")+ ylab("Count")
	ggsave(paste0("cREGION_cluster",cluster,"_3cM.png"), p)

	p=ggplot(dfprovince, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity")+ coord_flip() + ggtitle(paste0("cluster: ",cluster, " (", nb_sample," samples)"))+ xlab("Province of birth")+ ylab("Count")
	ggsave(paste0("province_birth_cluster",cluster,"_3cM.png"), p)

	p=ggplot(dfethnicityfr, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity")+ coord_flip() + ggtitle(paste0("cluster: ",cluster, " (", nb_sample," samples)"))+ xlab("Ethnicity_FR") + ylab("Count")
        ggsave(paste0("ethnicity_FR_cluster",cluster,"_3cM.png"), p)

	p=ggplot(dfethnicityme, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity")+ coord_flip() + ggtitle(paste0("cluster: ",cluster, " (", nb_sample," samples)"))+ xlab("Ethnicity_me") + ylab("Count")
        ggsave(paste0("ethnicity_me_cluster",cluster,"_3cM.png"), p)	

        p=ggplot(dfethnicity, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity")+ coord_flip() + ggtitle(paste0("cluster: ",cluster, " (", nb_sample," samples)"))+ xlab("Ethnicity") + ylab("Count")
        ggsave(paste0("ethnicity_cluster",cluster,"_3cM.png"), p)

	#combine four grand-parents
	four_gp=c(tmp$country_mothers_mother, tmp$country_mothers_father, tmp$country_fathers_mother, tmp$country_fathers_father)
	df2 <- as.data.frame(table(four_gp))
	p=ggplot(df2, aes(x=reorder(four_gp, -Freq), y=Freq)) + geom_bar(stat="identity")+ coord_flip() + ggtitle(paste0("cluster: ",cluster, " (", nb_sample," samples)"))+ xlab("Grand-parents country") + ylab("Count")
        ggsave(paste0("GrandParents_country_cluster",cluster,"_3cM.png"), p)

	}

	#-------HEATMAP-------------#
	
	# country of birth
	if(cluster == 15){
		custom_color="#B454DA"
		# Add a Region column based on the country
		dfcountry$Region <- with(dfcountry, 
                         ifelse(Var1 %in% c("Barbados", "Bermuda", "Dominica", "Dominican Republic",
                                            "Grenada", "Guadeloupe", "Jamaica", "Martinique",
                                            "Saint Vincent And Grenadines", "Trinidad And Tobago"), 
                                "Caribbean",
                         ifelse(Var1 %in% c("Australia", "Guam", "Heard Mcdonald Islands"), 
                                "Oceania",
                         ifelse(Var1 %in% c("Benin", "Burkina Faso", "Burundi", "Cameroon", "Cape Verde",
                                            "Chad", "Congo Democratic Republic", "Cote Divoire", "Ghana",
                                            "Guinea", "Guinea Bissau", "Niger", "Nigeria", "Senegal",
                                            "South Africa", "Togo"), 
                                "Africa",
                         ifelse(Var1 == "Haiti", 
                                "Haiti",
                         ifelse(Var1 %in% c("France", "United Kingdom"), 
                                "Europe",
                         ifelse(Var1 %in% c("Canada", "United States"), 
                                "North America",
                         ifelse(Var1 %in% c("Cuba", "Ecuador", "Costa Rica", "Venezuela"), 
                                "Central and South America",
                         ifelse(Var1 %in% c("China", "Cambodia"), 
                                "Asia", 
                                "Other")))))))))

		# Summarize proportions by region
		dfregion <- dfcountry %>%
		  group_by(Region) %>%
  		summarise(Proportion = sum(Proportion))
		dfcountry <- as.data.frame(dfregion) 
		names(dfcountry)=c("Var1", "Proportion")
	}else{
		custom_color="#D51A10"
	}
	p = ggplot(dfcountry, aes(x = factor(cluster), y = reorder(Var1, Proportion), fill = Proportion)) +
	  geom_tile(color = "black", width = 0.9, height = 0.9) + 
	  geom_text(aes(label = sprintf("%.2f%%", Proportion * 100)), color = "black", size = 12) +  
	  scale_fill_gradient(low = "white", high = custom_color) +    
	  #scale_fill_distiller(palette = "Reds", direction = 1) +    
	  xlab("IBD cluster") +
	  ylab("Place of birth") +
	  labs(fill = "number of samples") +
	  theme(text = element_text(size = 30)) +
	  theme(legend.position = "none") +
	  theme(panel.background = element_blank(),
	        panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank())
  	ggsave(paste0("heatmap_country_of_birth_",cluster,"_2cM.png"), p, dpi = 300)



	# Province of birth
	if(cluster == 8){
                custom_color="#F84C3C"
        }else{
                custom_color="#D51A10"
        }
        p = ggplot(dfprovince, aes(x = factor(cluster), y = reorder(Var1, Proportion), fill = Proportion)) +
          geom_tile(color = "black", width = 0.9, height = 0.9) +
          geom_text(aes(label = sprintf("%.2f%%", Proportion * 100)), color = "black", size = 12) +
          scale_fill_gradient(low = "white", high = custom_color) + 
	  #scale_fill_distiller(palette = "Reds", direction = 1) +
          xlab("IBD cluster") +
          ylab("Province of birth") +
          labs(fill = "number of samples") +
          theme(text = element_text(size = 30)) +
          theme(legend.position = "none") +
          theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())

	ggsave(paste0("heatmap_province_of_birth_",cluster,"_2cM.png"), p, dpi = 300)


	# CREGION
	if(cluster == 2){
                custom_color="#4EEA70"
        }else{
                custom_color="#D51A10"
        }
        p = ggplot(dfcregion, aes(x = factor(cluster), y = reorder(Var1, Proportion), fill = Proportion)) +
          geom_tile(color = "black", width = 0.9, height = 0.9) +
          geom_text(aes(label = sprintf("%.2f%%", Proportion * 100)), color = "black", size = 12) +
          scale_fill_gradient(low = "white", high = custom_color) + 
	  #scale_fill_distiller(palette = "Reds", direction = 1) +
          xlab("IBD cluster") +
          ylab("Recruitment center") +
          labs(fill = "number of samples") +
          theme(text = element_text(size = 30)) +
          theme(legend.position = "none") +
          theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())

	ggsave(paste0("heatmap_cregion_",cluster,"_2cM.png"), p, dpi = 300)


	# ethnicity_me
	if(cluster == 13){
                custom_color="#019FF4"

		# Replace "Whitejewish" by "White Jewish"
		dfethnicityme[which(dfethnicityme$Var1 == "Whitejewish"),]$Var1 = "White Jewish"

		# Combine categories with less than 1% into 'Other'
		dfethnicityme$Var1[dfethnicityme$Proportion < 0.01 & dfethnicityme$Var1 != "Other"] <- "Other"

		# Recalculate frequencies and proportions
		dfethnicityme_combined <- aggregate(cbind(Freq, Proportion) ~ Var1, data = dfethnicityme, sum)

		# Display the updated dataframe
		dfethnicityme = dfethnicityme_combined
        }else{
                custom_color="#D51A10"
        }
        p = ggplot(dfethnicityme, aes(x = factor(cluster), y = reorder(Var1, Proportion), fill = Proportion)) +
          geom_tile(color = "black", width = 0.9, height = 0.9) +
          geom_text(aes(label = sprintf("%.2f%%", Proportion * 100)), color = "black", size = 12) + 
	  scale_fill_gradient(low = "white", high = custom_color) + 
	  #scale_fill_distiller(palette = "Reds", direction = 1) +
          xlab("IBD cluster") +
          ylab("Self-declared ethnicity") +
          labs(fill = "number of samples") +
          theme(text = element_text(size = 30)) +
          theme(legend.position = "none") +
          theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
	
	ggsave(paste0("heatmap_ethnicity_",cluster,"_2cM.png"), p, dpi = 300)




	#-------Word cloud-----------#

	df_freq=data.frame(table(tmp$country))
	df_freq$Var1=tolower(df_freq$Var1)
	png(paste0("word_cloud_country", cluster, "_3cM.png"), width=12,height=8, units='in', res=500)
	p=wordcloud(words = df_freq$Var1, freq = df_freq$Freq, min.freq = 1, max.words=200, random.order=FALSE, rot.per=0.35)#, family = "mono")
	print(p)
	dev.off()

	df_freq=data.frame(table(tmp$CREGION))
        df_freq$Var1=tolower(df_freq$Var1)
        png(paste0("word_cloud_cregion", cluster, "_3cM.png"), width=12,height=8, units='in', res=500)
        p=wordcloud(words = df_freq$Var1, freq = df_freq$Freq, min.freq = 1, max.words=200, random.order=FALSE, rot.per=0.35)#, family = "mono")
        print(p)
        dev.off()

	df_freq=data.frame(table(tmp$ETHNICITY_CORRECT_ME))
	df_freq$Var1=tolower(df_freq$Var1)
        png(paste0("word_cloud_ethnicity", cluster, "_3cM.png"), width=12,height=8, units='in', res=500)
        p=wordcloud(words = df_freq$Var1, freq = df_freq$Freq, min.freq = 1, max.words=200, random.order=FALSE, rot.per=0.35)#, family = "mono")
        print(p)
        dev.off()

	df_freq=data.frame(table(tmp$PROVINCE_BIRTH))
	df_freq$Var1=tolower(df_freq$Var1)
        png(paste0("word_cloud_province", cluster, "_3cM.png"), width=12,height=8, units='in', res=500)
        p=wordcloud(words = df_freq$Var1, freq = df_freq$Freq, min.freq = 1, max.words=200, random.order=FALSE, rot.per=0.35 )#, family = "mono")
        print(p)
        dev.off()
}


#-------specific wordclouds---------

#ACADIANS
cluster=8
tmp=clusters[which(clusters$group == cluster),]

df_freq=data.frame(table(tmp$PROVINCE_BIRTH))
df_freq$Var1=tolower(df_freq$Var1)
test=c("black", "#F84C3C", "black", "black", "black", "black", "black", "black", "black", "black", "black")

#figure
png(paste0("word_cloud_province", cluster, "_3cM.png"), width=12,height=8, units='in', res=500)
p=wordcloud(words = df_freq$Var1, freq = df_freq$Freq, min.freq = 1, max.words=200, random.order=FALSE, rot.per=0.35, colors=test )#, family = "mono")
print(p)
dev.off()


#JEWISH
cluster=13
tmp=clusters[which(clusters$group == cluster),]
df_freq=data.frame(table(tmp$ETHNICITY_ME))
df_freq$Var1=tolower(df_freq$Var1)
test=c( "black", "black", "#019FF4", "black", "black", "black", "black", "black", "black", "black", "black")


png(paste0("word_cloud_ethnicity", cluster, "_3cM.png"), width=12,height=8, units='in', res=500)
p=wordcloud(words = df_freq$Var1, freq = df_freq$Freq, min.freq = 1, max.words=200, random.order=FALSE, random.color = F, rot.per=0.35, colors="black" )#, family = "mono")
print(p)
dev.off()





#-------------COMPARISON ALEX CLUSTER VS IBD CLUSTER-------------------------------



#empty vector for dataframe creation
count=c()
#total=c()
total_IBD=c()
total_UMAP=c()
group_id1=c()
group_id2=c()
list_clusters_tmp=unique(as.numeric(as.character(pca$group)))
list_clusters<-list_clusters_tmp[!is.na(list_clusters_tmp)]

for (i in list_clusters){
        for (j in unique(as.numeric(as.character(pca$group_alex)))){
                print(i)
                print(j)

                #add the group ids of comparison
                group_id1=c(group_id1, i)
                group_id2=c(group_id2, j)

                #count of individuals intersecting in these clusters
                tmp_count=nrow(pca[which(pca$group==i & pca$group_alex == j),])
		count=c(count, tmp_count)
		
		#count the total number in FST group and alex cluster
		total_IBD=c(total_IBD, nrow(pca[which(pca$group==i),]))
		total_UMAP=c(total_UMAP, nrow(pca[which(pca$group_alex==j),]))
		#total=c(total, (nrow(pca[which(pca$group==i),]) + nrow(pca[which(pca$group_alex==j),])))

        }
}

avg_df=data.frame(cluster_IBD=group_id1, cluster_UMAP=group_id2, nb_sample=count, total_IBD=total_IBD, total_UMAP=total_UMAP)

avg_df$cluster_IBD=as.character(avg_df$cluster_IBD)
avg_df$cluster_IBD=factor(avg_df$cluster_IBD, levels=as.character(sort(as.numeric(unique(avg_df$cluster_IBD)))))
avg_df$cluster_UMAP=as.character(avg_df$cluster_UMAP)
avg_df$cluster_UMAP=factor(avg_df$cluster_UMAP, levels=as.character(sort(as.numeric(unique(avg_df$cluster_UMAP)))))
avg_df$relative_count_IBD=round((avg_df$nb_sample/avg_df$total_IBD), digits=2)
avg_df$relative_count_UMAP=round((avg_df$nb_sample/avg_df$total_UMAP), digits=2)



p=ggplot(avg_df, aes(x=cluster_IBD , y = cluster_UMAP, fill = relative_count_IBD)) +  geom_tile(color = "black") +  geom_text(aes(label = relative_count_IBD), color = "black", size = 7)  + scale_fill_gradientn(colors=c("white", "#7F7F7F")) +xlab("IBD cluster") + ylab("HDBSCAN cluster") + labs(fill = "number of samples") + theme(text = element_text(size=25)) + theme(legend.title.align=0.5) +  theme(legend.position="none") 
ggsave("comparison_umap_vs_ibd_cluster_matrix_IBD_relative.png", height=10, width=20, p)

p=ggplot(avg_df, aes(x=cluster_IBD , y = cluster_UMAP, fill = relative_count_UMAP)) +  geom_tile(color = "black") +  geom_text(aes(label = relative_count_UMAP), color = "black", size = 5) +  coord_fixed() + scale_fill_gradientn(colors=c("white", "#BC0010")) +xlab("IBD cluster") + ylab("HDBSCAN cluster") + labs(fill = "number of samples") + theme(text = element_text(size=20)) + theme(legend.title.align=0.5) +  theme(legend.position="none") 
ggsave("comparison_umap_vs_ibd_cluster_matrix_UMAP_relative.png", height=10, width=10, p)































#-------------Color per group on PC1-PC2-------------------------------

clusters=pca[which(!(is.na(pca$group))),]

for (cluster in unique(clusters$group)){
	
	tmp=clusters
	tmp$cluster="other"
	tmp[which(tmp$group == cluster),]$cluster=cluster
	tmp$cluster=factor(tmp$cluster, levels=c("other", cluster))
	pop_colors=c("#354446","#F8766D")

	print(table(tmp$cluster))
	
	if(cluster==2){
                print("saguenay")
                pop_colors=c("#354446","#4EEA70")
        }
        if(cluster==13){
                print("jewish")
                pop_colors=c("#354446","#019FF4")
        }
        if(cluster==8){
                print("acadian")
                pop_colors=c("#354446","#F84C3C")
        }
        if(cluster==15){ 
                print("haitian")
                pop_colors=c("#354446","#B454DA")
        }

	names(pop_colors)=levels(tmp$cluster)
        part1=tmp[which(tmp$group == cluster),]
        part2=tmp[which(tmp$group != cluster),]
        tmp=rbind(part2, part1)


	p1<-ggplot(tmp, aes(x=PC1, y=PC2, color=cluster)) + geom_point(alpha=0.8, size = 4) +theme_classic()+scale_colour_manual(values=pop_colors[tmp$cluster]) +theme(text = element_text(size=30))+
	 scale_fill_discrete(name = "IBD Cluster")
	ggsave(paste0("CaG_genotyping_PC1_PC2_cluster",cluster,"_2cM.png"),p1, dpi=300)


}







#without FST merging
clusters=pca[which(!(is.na(pca$group3))),]

for (cluster in unique(clusters$group3)){

        tmp=clusters
        tmp$cluster="other"
        tmp[which(tmp$group3 == cluster),]$cluster=cluster
        tmp$cluster=factor(tmp$cluster, levels=c("other", cluster))
        pop_colors=c("#354446","#F8766D")
        names(pop_colors)=levels(tmp$cluster)
        part1=tmp[which(tmp$group3 == cluster),]
        part2=tmp[which(tmp$group3 != cluster),]
        tmp=rbind(part2, part1)

        print(table(tmp$cluster))

        p1<-ggplot(tmp, aes(x=PC1, y=PC2, color=cluster)) + geom_point(alpha=0.8) +theme_classic()+scale_colour_manual(values=pop_colors[tmp$cluster])
        ggsave(paste0("group3_CaG_genotyping_PC1_PC2_cluster",cluster,"_3cM.png"),p1)


}

#color morrocan
clusters$cluster="other"
clusters[which(clusters$country == "MOROCCO"),]$cluster="morocco"
clusters$cluster=factor(tmp$cluster, levels=c("other", "morocco"))
pop_colors=c("#354446","#F8766D")
names(pop_colors)=levels(clusters$cluster)
part1=clusters[which(clusters$cluster == "morocco"),]
part2=clusters[which(clusters$cluster != "morocco"),]
tmp=rbind(part2, part1)
p1<-ggplot(tmp, aes(x=PC1, y=PC2, color=cluster)) + geom_point(alpha=0.8) +theme_classic()+scale_colour_manual(values=pop_colors[tmp$cluster])
ggsave(paste0("CaG_genotyping_PC1_PC2_MORROCAN_3cM.png"),p1)

#------UMAP Alex----------

#FST clusters
clusters=pca[which(!(is.na(pca$group))),]

for (cluster in c(1,16,30,24)){
#for (cluster in unique(clusters$group)){

        tmp=clusters
        tmp$cluster="other"
        tmp[which(tmp$group == cluster),]$cluster=cluster
        tmp$cluster=factor(tmp$cluster, levels=c("other", cluster))
        print(cluster)
	pop_colors=c("#354446","#F8766D")
	if(cluster==1){
		print("saguenay")
                pop_colors=c("#354446","#4EEA70")
        }
        if(cluster==16){
		print("jewish")
                pop_colors=c("#354446","#019FF4")
        }
        if(cluster==30){
		print("acadian")
                pop_colors=c("#354446","#F84C3C")
        }
	if(cluster==24){
                print("haitian")
                pop_colors=c("#354446","#B454DA")
        }
	names(pop_colors)=levels(tmp$cluster)
        part1=tmp[which(tmp$group == cluster),]
        part2=tmp[which(tmp$group != cluster),]
        tmp=rbind(part2, part1)

        print(table(tmp$cluster))

        p1<-ggplot(tmp, aes(x=umap1, y=umap2, color=cluster)) + geom_point(alpha=0.8) +theme_classic()+scale_colour_manual(values=pop_colors[tmp$cluster]) + xlab("UMAP 1") +ylab("UMAP 2") +theme(text = element_text(size=20))
        ggsave(paste0("CaG_genotyping_Alex_UMAP_cluster",cluster,"_3cM.png"),p1)

}

#----ADMIXED highlight-----#
admx_cluster1=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version_final/ADMIXTURE/admixed_cluster1.txt", header=T)

#FST clusters
clusters=pca[which(!(is.na(pca$group))),]

for (cluster in c(1,16,30,24)){
#for (cluster in unique(clusters$group)){

        tmp=clusters
        tmp$cluster="other"
        tmp[which(tmp$group == cluster),]$cluster=cluster
        print(cluster)
        pop_colors=c("#354446","#F8766D")
        if(cluster==1){
                print("saguenay")
		tmp[which(tmp$IID %in% admx_cluster1$V4),]$cluster="cluster1_admixed"
        	tmp$cluster=factor(tmp$cluster, levels=c("other", cluster, "cluster1_admixed"))
                print(tmp[which(tmp$cluster == "cluster1_admixed"),])
		pop_colors=c("#354446","#4EEA70", "#FFE346")
        }
        if(cluster==16){
                print("jewish")
                pop_colors=c("#354446","#019FF4", "#FFE346")
        }
        if(cluster==30){
                print("acadian")
                pop_colors=c("#354446","#F84C3C", "#FFE346")
        }
        if(cluster==24){
                print("haitian")
                pop_colors=c("#354446","#B454DA", "#FFE346")
        }
        names(pop_colors)=levels(tmp$cluster)
        part1=tmp[which(tmp$group == cluster),]
        part2=tmp[which(tmp$group != cluster),]
        tmp=rbind(part2, part1)

        print(table(tmp$cluster))

        p1<-ggplot(tmp, aes(x=umap1, y=umap2, color=cluster)) + geom_point(alpha=0.8) +theme_classic()+scale_colour_manual(values=pop_colors[tmp$cluster]) + xlab("UMAP 1") +ylab("UMAP 2") +theme(text = element_text(size=20))
        ggsave(paste0("CaG_genotyping_Alex_UMAP_ADMIXED_cluster",cluster,"_2cM.png"),p1)

}




















#Louvain3 clusters
clusters=pca[which(!(is.na(pca$group3))),]

for (cluster in unique(clusters$group3)){

        tmp=clusters
        tmp$cluster="other"
        tmp[which(tmp$group3 == cluster),]$cluster=cluster
        tmp$cluster=factor(tmp$cluster, levels=c("other", cluster))
        pop_colors=c("#354446","#F8766D")
        names(pop_colors)=levels(tmp$cluster)
        part1=tmp[which(tmp$group3 == cluster),]
        part2=tmp[which(tmp$group3 != cluster),]
        tmp=rbind(part2, part1)

        print(table(tmp$cluster))

        p1<-ggplot(tmp, aes(x=umap1, y=umap2, color=cluster)) + geom_point(alpha=0.8) +theme_classic()+scale_colour_manual(values=pop_colors[tmp$cluster])
        ggsave(paste0("CaG_genotyping_Alex_UMAP_louvain3_cluster",cluster,"_3cM.png"),p1)

}

#Louvain2 clusters
clusters=pca[which(!(is.na(pca$group2))),]

for (cluster in unique(clusters$group2)){

        tmp=clusters
        tmp$cluster="other"
        tmp[which(tmp$group2 == cluster),]$cluster=cluster
        tmp$cluster=factor(tmp$cluster, levels=c("other", cluster))
        pop_colors=c("#354446","#F8766D")
        names(pop_colors)=levels(tmp$cluster)
        part1=tmp[which(tmp$group2 == cluster),]
        part2=tmp[which(tmp$group2 != cluster),]
        tmp=rbind(part2, part1)

        print(table(tmp$cluster))

        p1<-ggplot(tmp, aes(x=umap1, y=umap2, color=cluster)) + geom_point(alpha=0.8) +theme_classic()+scale_colour_manual(values=pop_colors[tmp$cluster])
        ggsave(paste0("CaG_genotyping_Alex_UMAP_louvain2_cluster",cluster,"_3cM.png"),p1)

}





#-------ADMIXED in cluster 24 - Haiti----------#

#get outliers 24
outliers1_cluster24=pca[which(pca$group == 24 & (pca$umap1 >8)),]
outliers2_cluster24=pca[which(pca$group == 24 & (pca$umap2 <9)),]
write.table(outliers_cluster24, "outliers_cluster24_FST.txt", quote=F, row.names=F, sep="\t")
#cut -f1,30-34,37,38 outliers_cluster24_FST.txt

#cluster21
cluster21=pca[which(pca$group == 21),]

#branch under latinos (26)
branch_weird=pca[which(pca$umap1 <=5.25 & pca$umap2<=7 & pca$umap2>5),]
#probably admixed indian descent


#cluster down middle
cluster_bottom=pca[which(pca$umap1 <=10 & pca$umap2<=-9.5 & pca$umap1>6.5),]



#cluster1 asian
cluster_bottom=pca[which(pca$umap1 <=2.5 & pca$umap2>=-1  &  pca$umap2<=1 $ pca$group == 1),]






























#------UMAP----------
library(umap)


PC1_to_PC10<-pca[,3:7]
#PC1_to_PC10<-pca[,3:12]
umap_cag<-umap(PC1_to_PC10)

data_umap=data.frame(umap1=umap_cag$layout[,1], umap2=umap_cag$layout[,2])
data_umap$country=pca$country
data_umap$group=pca$group
data_umap$group3=pca$group3
data_umap$IID=pca$IID


clusters=data_umap[which(!(is.na(data_umap$group))),]

for (cluster in unique(clusters$group)){

        tmp=clusters
        tmp$cluster="other"
        tmp[which(tmp$group == cluster),]$cluster=cluster
        tmp$cluster=factor(tmp$cluster, levels=c("other", cluster))
        pop_colors=c("#354446","#F8766D")
        names(pop_colors)=levels(tmp$cluster)
        part1=tmp[which(tmp$group == cluster),]
        part2=tmp[which(tmp$group != cluster),]
        tmp=rbind(part2, part1)

        print(table(tmp$cluster))

        p1<-ggplot(tmp, aes(x=umap1, y=umap2, color=cluster)) + geom_point(alpha=0.8) +theme_classic()+scale_colour_manual(values=pop_colors[tmp$cluster])
        #ggsave(paste0("CaG_umap_PC1_PC10_cluster",cluster,"_3cM.png"),p1)
        ggsave(paste0("CaG_umap_PC1_PC5_cluster",cluster,"_3cM.png"),p1)



}

#without FST merging
clusters=data_umap[which(!(is.na(data_umap$group3))),]

for (cluster in unique(clusters$group3)){

        tmp=clusters
        tmp$cluster="other"
        tmp[which(tmp$group3 == cluster),]$cluster=cluster
        tmp$cluster=factor(tmp$cluster, levels=c("other", cluster))
        pop_colors=c("#354446","#F8766D")
        names(pop_colors)=levels(tmp$cluster)
        part1=tmp[which(tmp$group3 == cluster),]
        part2=tmp[which(tmp$group3 != cluster),]
        tmp=rbind(part2, part1)

        print(table(tmp$cluster))

        p1<-ggplot(tmp, aes(x=umap1, y=umap2, color=cluster)) + geom_point(alpha=0.8) +theme_classic()+scale_colour_manual(values=pop_colors[tmp$cluster])
        ggsave(paste0("group3_CaG_umap_PC1_PC10_cluster",cluster,"_3cM.png"),p1)


}




#--------------------------STATS within cluster------------------------------


#IBD stats per pair
ibd=read.table("/lustre07/scratch/justinp/NextFlow/CaG_genotyping_only_version/CaG_Result_3cM/all_chr.per_pair_unique.hapIBD_phaseIBD.ibd", header=T)

#get segment for all the individuals in the clusters even the one without a shared segment
test=expand.grid(unique(cluster$IID),unique(cluster$IID))
names(test)=c("id1", "id2")
new_df=test[which(test$id1 != test$id2),]


#include the group in the ibd per pair info
new_df$group_id1=cluster$group[match(new_df$id1, cluster$IID)]
new_df$group_id2=cluster$group[match(new_df$id2, cluster$IID)]

#only keep pairs within the same cluster
within_ibd=[which(new_df$group_id1 == new_df$group_id2),]

#add info on IBD segments
within_ibd$group_id=paste0(within_ibd$id1,"_",within_ibd$id2)
ibd$group_id=paste0(ibd$id1,"_",ibd$id2)
within_ibd$size_cM_hapIBD=ibd$size_cM_hapIBD[match(within_ibd$group_id, ibd$group_id)]
within_ibd[which(is.na(within_ibd$size_cM_hapIBD)),]$size_cM_hapIBD=0
within_ibd$size_cM_phaseIBD=ibd$size_cM_phaseIBD[match(within_ibd$group_id, ibd$group_id)]
within_ibd[which(is.na(within_ibd$size_cM_phaseIBD)),]$size_cM_phaseIBD=0


#IBD within cluster
avg_ibd=as.data.frame(aggregate(within_ibd$size_cM_hapIBD, list(within_ibd$group), FUN=mean))
names(avg_ibd)=c("cluster", "mean_size_cM_hapIBD")












#-----------------------Subset the europeans------------------------


eur = pca[which((pca$PC2 >= -0.023 & pca$PC2 < 0.025) & (pca$PC1 <= 0.04)),]


#test on pca
pca$test=0
pca[which(pca$IID %in% eur$IID),]$test=1
p2<-ggplot(pca, aes(x=PC1, y=PC2, color=as.factor(test))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC1_PC2_test_EUR_3cM.png",p2)


eur_data=data.frame(IID=eur$IID, FID=eur$FID)
write.table(eur_data, "EUR_CAG_genotyping.txt", row.names=F, col.names=F, quote=F)



#------------PCA results-----------


pca = read.table("CARTaGENE_hg38_shared.phased.maf.hwe.missingness.chr_all.vcf.gz.LDpruned.EUR.PC", header=T) 
#match country of  birth for all individuals
pca$country=metadata$COUNTRY_BIRTH[match(pca$IID, metadata$PROJECT_CODE)]

#----outliers----
outliers_pca<-fc.pca.ID.outliers.pcs(pca, npc=10, nsd=10)
#write.table(outliers_pca, "outliers_CaG_genotyping_pca_EUR.txt", row.name=F, col.name=F, quote=F)


p1<-ggplot(pca, aes(x=PC1, y=PC2, color=as.factor(country))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC1_PC2_EUR_3cM.png",p1)

p2<-ggplot(pca, aes(x=PC3, y=PC4, color=as.factor(country))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC3_PC4_EUR_3cM.png",p2)

p3<-ggplot(pca, aes(x=PC5, y=PC6, color=as.factor(country))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC5_PC6_EUR_3cM.png",p3)

p4<-ggplot(pca, aes(x=PC7, y=PC8, color=as.factor(country))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC7_PC8_EUR_3cM.png",p4)







#------------------Color by cluster of IBD with louvain----------------------------




p1<-ggplot(pca, aes(x=PC1, y=PC2, color=as.factor(group1))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC1_PC2_IBD1_3cM.png",p1)
p1<-ggplot(pca, aes(x=PC1, y=PC2, color=as.factor(group2))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC1_PC2_IBD2_3cM.png",p1)
p1<-ggplot(pca, aes(x=PC1, y=PC2, color=as.factor(group3))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC1_PC2_IBD3_3cM.png",p1)

p2<-ggplot(pca, aes(x=PC3, y=PC4, color=as.factor(group1))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC3_PC4_IBD1_3cM.png",p2)
p2<-ggplot(pca, aes(x=PC3, y=PC4, color=as.factor(group2))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC3_PC4_IBD2_3cM.png",p2)
p2<-ggplot(pca, aes(x=PC3, y=PC4, color=as.factor(group3))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC3_PC4_IBD3_3cM.png",p2)

p2<-ggplot(pca, aes(x=PC5, y=PC6, color=as.factor(group1))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC5_PC6_IBD1_3cM.png",p2)
p2<-ggplot(pca, aes(x=PC5, y=PC6, color=as.factor(group2))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC5_PC6_IBD2_3cM.png",p2)
p2<-ggplot(pca, aes(x=PC5, y=PC6, color=as.factor(group3))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC5_PC6_IBD3_3cM.png",p2)

p2<-ggplot(pca, aes(x=PC7, y=PC8, color=as.factor(group1))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC7_PC8_IBD1_3cM.png",p2)
p2<-ggplot(pca, aes(x=PC7, y=PC8, color=as.factor(group2))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC7_PC8_IBD2_3cM.png",p2)
p2<-ggplot(pca, aes(x=PC7, y=PC8, color=as.factor(group3))) + geom_point(alpha=0.8)
ggsave("CaG_genotyping_PC7_PC8_IBD3_3cM.png",p2)



}
