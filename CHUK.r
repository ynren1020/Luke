#######################2019-09-03########################
##CHUK FPKM in EGFR mutant or wt patients distribution###
##TCGA cohort LUAD and LUSC #############################
#########################################################

library(tidyverse)
library(ggpubr)
library(gridExtra)
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


chuk<-function(input1,input2,a,b){
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
luad.fpkmL<-chuk(input1,input2)

##add cohort##
lusc.fpkmL$cohort<-"LUSC" #WT 457 mutant 12
luad.fpkmL$cohort<-"LUAD" #WT 441 mutant 65

lusc.fpkmL$log2FPKM<-log2(lusc.fpkmL$FPKM)
luad.fpkmL$log2FPKM<-log2(luad.fpkmL$FPKM)
##rbind for boxplot ##
#luad.lusc.join<-bind_rows(luad.fpkmL,lusc.fpkmL)
#luad.lusc.join$log2FPKM<-log2(luad.lusc.join$FPKM)
##boxplot##
p1<-ggboxplot(luad.fpkmL,"group","log2FPKM",
             fill = "group",
             facet.by = "cohort",
             xlab = FALSE,
             ylab = "CHUK expression (log2FPKM)")
#ggpar(p,legend = "none",yscale = "log2")
p2<-p1+stat_compare_means(method = "t.test",label.x = 1.4,label.y = 4.3)+
  scale_x_discrete(labels=c("EGFR mutant" = "EGFR mutant\n(N=65)", "WT" = "WT\n(N=441)"))+rotate_x_text(45)
p22<-ggpar(p2,legend = "none",ylim = c(0.5,4.7),font.y = c(14,"bold"),font.tickslab = c(12,"bold"))

  
#LUSC##
mytheme <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)

p3<-ggboxplot(lusc.fpkmL,"group","log2FPKM",
                fill = "group",
                facet.by = "cohort",
                xlab = FALSE,
                ylab = FALSE)
  #ggpar(p,legend = "none",yscale = "log2")
p4<-p3+stat_compare_means(method = "t.test",label.x = 1.4,label.y = 4.3)+
    scale_x_discrete(labels=c("EGFR mutant" = "EGFR mutant\n(N=12)", "WT" = "WT\n(N=457)"))+rotate_x_text(45)
p44<-ggpar(p4,legend = "none",ylim=c(0.5,4.7),font.xtickslab =  c(12,"bold"))+mytheme

pdf("LUAD.LUSC.CHUK.FPKM.boxplot_withsamplesize.pdf",width = 5,height = 5)
grid.arrange(p22,p44,ncol=2)
dev.off()



##the way to save grid.arrange output##
pdf("filename.pdf", width = 8, height = 12) # Open a new pdf file
grid.arrange(plot1, plot2, plot3, nrow=3) # Write the grid.arrange in the file
dev.off() # Close the file

