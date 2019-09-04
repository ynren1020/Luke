###################2019-09-04###################
##salmon output of LUAD and LUSC for boxplot####
##DARPP-32(long),t-DARPP(short)2 isoforms ######
################################################

library(tidyverse)
library(ggpubr)

isof<-read.delim("luad.lusc.isoform.txt",header = FALSE,stringsAsFactors = FALSE)
names(isof)<-c("Name","Length","Effectivelength","TPM","NumReads","id","gene")
meta<-read.delim("TCGA.metainfo.txt",header = TRUE,stringsAsFactors = FALSE)

##match id to get cohort and patient info##
isof.meta<-left_join(isof,meta,by=c("id"="FILE_ID"))

##EGFR mutant info##
luad.mut<-read.delim("LUAD.mut.txt",header = FALSE,stringsAsFactors = FALSE,check.names = FALSE)
lusc.mut<-read.delim("LUSC.mut.txt",header = FALSE,stringsAsFactors = FALSE,check.names = FALSE)

mut<-bind_rows(luad.mut,lusc.mut)
mut$patient<-str_sub(mut$V4,1,12)

isof.meta$group<-ifelse(isof.meta$PATIENT_ID%in%mut$patient,"EGFR mutant","WT")
isof.meta$cohort<-str_sub(isof.meta$PROJECT,6,9)
write.table(isof.meta,"luad.lusc.isoform.EGFRmut.quant.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")

##subset cohort##
isof.meta.sub<-isof.meta[,c(4,7,9,12,13)]

isof.meta.sub.luad<-isof.meta.sub[isof.meta.sub$cohort=="LUAD",]
isof.meta.sub.lusc<-isof.meta.sub[isof.meta.sub$cohort=="LUSC",]

##subset long and short isoform##
##luad##
isof.meta.sub.luad.l<-isof.meta.sub.luad[isof.meta.sub.luad$gene=="DARPP-32",]
isof.meta.sub.luad.s<-isof.meta.sub.luad[isof.meta.sub.luad$gene=="t-DARPP",]

isof.meta.sub.luad.l<-isof.meta.sub.luad.l%>%rename("DARPP32"="TPM")
isof.meta.sub.luad.s<-isof.meta.sub.luad.s%>%rename("tDARPP"="TPM")
                                                    
isof.meta.sub.luad.l<-isof.meta.sub.luad.l[,-2]
isof.meta.sub.luad.s<-isof.meta.sub.luad.s[,-2]

isof.meta.sub.luad.join<-full_join(isof.meta.sub.luad.l,isof.meta.sub.luad.s,by=c("PATIENT_ID","group","cohort"))

isof.meta.sub.luad.join$ratio<-isof.meta.sub.luad.join$tDARPP/isof.meta.sub.luad.join$DARPP32
##WT 679 mutant 116##

##lusc##
isof.meta.sub.lusc.l<-isof.meta.sub.lusc[isof.meta.sub.lusc$gene=="DARPP-32",]
isof.meta.sub.lusc.s<-isof.meta.sub.lusc[isof.meta.sub.lusc$gene=="t-DARPP",]

isof.meta.sub.lusc.l<-isof.meta.sub.lusc.l%>%rename("DARPP32"="TPM")
isof.meta.sub.lusc.s<-isof.meta.sub.lusc.s%>%rename("tDARPP"="TPM")

isof.meta.sub.lusc.l<-isof.meta.sub.lusc.l[,-2]
isof.meta.sub.lusc.s<-isof.meta.sub.lusc.s[,-2]

isof.meta.sub.lusc.join<-full_join(isof.meta.sub.lusc.l,isof.meta.sub.lusc.s,by=c("PATIENT_ID","group","cohort"))

isof.meta.sub.lusc.join$ratio<-isof.meta.sub.lusc.join$tDARPP/isof.meta.sub.lusc.join$DARPP32
##WT 636 mutant 15##
##logscale##
isof.meta.sub.luad.join$logDARPP32<-log10(isof.meta.sub.luad.join$DARPP32+0.5)
isof.meta.sub.luad.join$logtDARPP<-log10(isof.meta.sub.luad.join$tDARPP+0.5)
isof.meta.sub.luad.join$ratio2<-(isof.meta.sub.luad.join$tDARPP+0.5)/(isof.meta.sub.luad.join$DARPP32+0.5)
#isof.meta.sub.luad.join$logratio<-isof.meta.sub.luad.join$logtDARPP/isof.meta.sub.luad.join$logDARPP32
isof.meta.sub.luad.join$logratio<-log2(isof.meta.sub.luad.join$ratio2)
isof.meta.sub.luad.join$ratioR<-1/isof.meta.sub.luad.join$ratio

isof.meta.sub.lusc.join$logDARPP32<-log10(isof.meta.sub.lusc.join$DARPP32+0.5)
isof.meta.sub.lusc.join$logtDARPP<-log10(isof.meta.sub.lusc.join$tDARPP+0.5)
isof.meta.sub.lusc.join$ratio2<-(isof.meta.sub.lusc.join$tDARPP+0.5)/(isof.meta.sub.lusc.join$DARPP32+0.5)
#isof.meta.sub.lusc.join$logratio<-isof.meta.sub.lusc.join$logtDARPP/isof.meta.sub.lusc.join$logDARPP32
isof.meta.sub.lusc.join$logratio<-log2(isof.meta.sub.lusc.join$ratio2)
isof.meta.sub.lusc.join$ratioR<-1/isof.meta.sub.lusc.join$ratio

#rename group to LUAD or LUSC##
isof.meta.sub.luad.join<-isof.meta.sub.luad.join%>%rename("group"="LUAD")
isof.meta.sub.lusc.join<-isof.meta.sub.lusc.join%>%rename("group"="LUSC")


##boxplot##
mythemey <- theme(
  #axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)
mythemexy <- theme(
  #axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.x = element_blank()
  #axis.ticks.x = element_blank()
)

mythemex <- theme(
  #axis.title.y = element_blank(),
  #axis.text.y = element_blank(),
  #axis.ticks.y = element_blank(),
  axis.text.x = element_blank()
  #axis.ticks.x = element_blank()
)

##long luad##
p1<-ggboxplot(isof.meta.sub.luad.join,"group","DARPP32",
              fill = "group",
              #facet.by = "cohort",
              xlab = FALSE,
              ylab = "DARPP-32 expression (TPM)")
#ggpar(p,legend = "none",yscale = "log2")
#
p2<-p1+stat_compare_means(method = "t.test",label.x = 1.3)
 # scale_x_discrete(labels=c("EGFR mutant" = "EGFR mutant\n(N=116)", "WT" = "WT\n(N=679)"))+rotate_x_text(45)
p22<-ggpar(p2,legend = "none",font.y = c(12,"bold"),font.ytickslab = c(12,"bold"))+mythemex
#ylim = c(0.5,4.7),

##long lusc##
p1s<-ggboxplot(isof.meta.sub.lusc.join,"group","DARPP32",
              fill = "group",
              #facet.by = "cohort",
              width = 0.6,
              xlab = FALSE,
              ylab = FALSE)
#ggpar(p,legend = "none",yscale = "log2")
p2s<-p1s+stat_compare_means(method = "t.test",label.x = 1.3) #method = "t.test"
#  scale_x_discrete(labels=c("EGFR mutant" = "EGFR mutant\n(N=15)", "WT" = "WT\n(N=636)"))+rotate_x_text(45)
p22s<-ggpar(p2s,legend = "none")+mythemexy
#ylim = c(0.5,4.7),

##short##
##luad##
p1t<-ggboxplot(isof.meta.sub.luad.join,"group","tDARPP",
              fill = "group",
              #facet.by = "cohort",
              xlab = FALSE,
              ylab = "t-DARPP expression (TPM)")
#ggpar(p,legend = "none",yscale = "log2")
p2t<-p1t+stat_compare_means(method = "t.test",label.x = 1.3) #method = "t.test"
#  scale_x_discrete(labels=c("EGFR mutant" = "EGFR mutant\n(N=116)", "WT" = "WT\n(N=679)"))+rotate_x_text(45)
p22t<-ggpar(p2t,legend = "none",font.y = c(12,"bold"),font.ytickslab = c(12,"bold"))+mythemex
#ylim = c(0.5,4.7),

##long lusc##
p1st<-ggboxplot(isof.meta.sub.lusc.join,"group","tDARPP",
               fill = "group",
               #facet.by = "cohort",
               width = 0.6,
               xlab = FALSE,
               ylab = FALSE)
#ggpar(p,legend = "none",yscale = "log2")
p2st<-p1st+stat_compare_means(method = "t.test",label.x = 1.3) #method = "t.test"
 # scale_x_discrete(labels=c("EGFR mutant" = "EGFR mutant\n(N=15)", "WT" = "WT\n(N=636)"))+rotate_x_text(45)
p22st<-ggpar(p2st,legend = "none")+mythemexy
#ylim = c(0.5,4.7),



##ratio##
##luad##
p1tR<-ggboxplot(isof.meta.sub.luad.join,"group","ratio",
               fill = "group",
               #facet.by = "cohort",
               xlab = "LUAD",
               ylab = "t-DARPP/DARPP-32")
#ggpar(p,legend = "none",yscale = "log2")
#method = "t.test"
p2tR<-p1tR+stat_compare_means(method = "t.test",label.x = 1.3,label.y=140)+
  scale_x_discrete(labels=c("EGFR mutant" = "EGFR mutant\n(N=116)", "WT" = "WT\n(N=679)"))+rotate_x_text(45)
p22tR<-ggpar(p2tR,legend = "none",font.y = c(12,"bold"),ylim = c(0.5,150),font.x = c(12,"bold"),font.tickslab = c(12,"bold"))
#

##lusc##
p1stR<-ggboxplot(isof.meta.sub.lusc.join,"group","ratio",
                fill = "group",
                #facet.by = "cohort",
                width = 0.6,
                xlab = "LUSC",
                ylab = FALSE)
#ggpar(p,legend = "none",yscale = "log2")
#method = "t.test"
p2stR<-p1stR+stat_compare_means(method = "t.test",label.x = 1.3,label.y = 140)+
  scale_x_discrete(labels=c("EGFR mutant" = "EGFR mutant\n(N=15)", "WT" = "WT\n(N=636)"))+rotate_x_text(45)
p22stR<-ggpar(p2stR,legend = "none",ylim = c(0.5,150),font.x = c(12,"bold"),font.xtickslab = c(12,"bold"))+mythemey
#

pdf("LUAD.LUSC.PPP1R1B_isoform.TPM.boxplot_withsamplesize.pdf",width = 8,height = 8)
grid.arrange(p22,p22s,p22t,p22st,p22tR,p22stR,ncol=2)
dev.off()














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







