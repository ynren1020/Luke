#######################2019-09-03########################
##CHUK FPKM in EGFR mutant or wt patients distribution###
##TCGA cohort LUAD and LUSC #############################
#########################################################

library(tidyverse)
library(ggpubr)
input1<-"LUAD.FPKM.CHUK.GeneID.txt"
input2<-"LUAD.mut.txt"
input3<-"LUSC.FPKM.CHUK.GeneID.txt"
input4<-"LUSC.mut.txt"

luad.fpkm<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
luad.mut<-read.delim(input2,header = FALSE,stringsAsFactors = FALSE,check.names = FALSE)
lusc.fpkm<-read.delim(input3,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
lusc.mut<-read.delim(input4,header = FALSE,stringsAsFactors = FALSE,check.names = FALSE)


luad.fpkmL<-gather(luad.fpkm,"patient","FPKM",-GeneID)
luad.fpkmL$id<-str_sub(luad.fpkmL$patient,start = 1,end=12)
luad.mut$id<-str_sub(luad.mut$V4,start = 1,end=12)

for (i in 1:nrow(luad.fpkmL)){
  luad.fpkmL$group[i]<-ifelse(luad.fpkmL$id[i]%in%luad.mut$id,"EGFR mutant", "WT")
}


chuk<-function(input1,input2){
  fpkm<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
  mut<-read.delim(input2,header = FALSE,stringsAsFactors = FALSE,check.names = FALSE)
  fpkmL<-gather(fpkm,"patient","FPKM",-GeneID)
  fpkmL$id<-str_sub(fpkmL$patient,start = 1,end=12)
  mut$id<-str_sub(mut$V4,start = 1,end=12)
  
  for (i in 1:nrow(fpkmL)){
    fpkmL$group[i]<-ifelse(fpkmL$id[i]%in%mut$id,"EGFR mutant", "WT")
  }
  return(fpkmL)
}

lusc.fpkmL<-chuk(input3,input4)

##add cohort##
lusc.fpkmL$cohort<-"LUSC"
luad.fpkmL$cohort<-"LUAD"

##rbind for boxplot ##
luad.lusc.join<-bind_rows(luad.fpkmL,lusc.fpkmL)
luad.lusc.join$log2FPKM<-log2(luad.lusc.join$FPKM)
##boxplot##
p<-ggboxplot(luad.lusc.join,"group","log2FPKM",
             fill = "group",
             facet.by = "cohort",
             xlab = FALSE,
             ylab = "CHUK expression (log2FPKM)")
#ggpar(p,legend = "none",yscale = "log2")
p1<-p+stat_compare_means(method = "t.test",label.y = 4.3,label.x = 1.2)
ggpar(p1,legend = "none")
ggsave("LUAD.LUSC.CHUK.FPKM.boxplot.pdf",width = 5,height = 5,dpi = 300)




