###########################2019-09-04###########################
##modify quant.sf file from salmon with file_id as one column###
################################################################

#library(tidyverse)

args<-commandArgs(TRUE)

#dat<-read.delim("ffef03e7-ce56-408e-a6d6-2357f597da28_quant.sf",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
#dat$id<-strsplit("ffef03e7-ce56-408e-a6d6-2357f597da28_quant.sf","_")[[1]][1]
#dat$gene<-ifelse(dat$Name==0,"DARPP-32","t-DARPP")

#write.table(dat,"ffef03e7-ce56-408e-a6d6-2357f597da28_quant.sf.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")

##for Rscript##
dat<-read.delim(args[1],header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
dat$id<-strsplit(args[1],"_")[[1]][1]
dat$gene<-ifelse(dat$Name==0,"DARPP-32","t-DARPP")
write.table(dat,paste0(args[1],".txt"),quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")