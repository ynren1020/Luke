###############2019-08-02################
##from CI to SE #########################
#########################################
input<-"shDP32_03.degs.txt"
input<-"shDP32_04.degs.txt"


secal<-function(input){
  dat<-read.delim(input,header = TRUE,stringsAsFactors = FALSE)
  dat$SE<-(dat$CI.R-dat$CI.L)/(1.96*2)
  dat$gene_id<-row.names(dat)
  dat<-as.data.frame(dat)
  write.table(dat,file = paste0(strsplit(input,"[.]")[[1]][1],".meta.txt"),quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
}

secal(input="shDP32_03.degs.txt")
secal(input="shDP32_04.degs.txt")

