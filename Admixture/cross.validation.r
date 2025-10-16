library(stringr)
library(ggplot2)
library(dplyr)

# grep -w "CV error" *.out | cut -f2,3 -d":" > cross_validation.txt

cv <- read.table("cross_validation.txt")

cv$K <-gsub("[\\(K=\\)]", "", regmatches(cv$V3, gregexpr("\\(.*?\\)", cv$V3)))
CV <- cv %>% mutate(K=as.numeric(K)) %>% select(V4,K)

#CV <- CV[seq(1,18,3),]

colnames(CV) <- c("cross.validation","K")

graph_title="Cross-Validation plot"
x_title="K"
y_title="Cross-validation error"
graph_1 <- ggplot(CV,aes(x=K,y=cross.validation)) + geom_point()

graph_1+geom_line()+
  labs(title=graph_title)+
  labs(x=x_title)+
  labs(y=y_title)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour = "black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour ="black", fill=NA, linewidth =3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))


ggsave("Admixture_cross-validation.png",width=7,height=5,dpi=600)
dev.off()
