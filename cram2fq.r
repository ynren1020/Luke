###################2019-09-03######################
##TCGA cohort: LUAD LUSC cram file to fastq file###
##for transcript quantification by salmon or RSEM##
###################################################
library(tidyverse)

input<-"tcga_cmd_luad.cram2fq.txt"
input<-"tcga_cmd_lusc.cram2fq.txt"

luad<-read.delim(input,header = FALSE,stringsAsFactors = FALSE,sep = " ")
#samtools fastq -@4 -1 sample1_R1.fastq -2 sample1_R2.fastq /data/ryang/projects/tcga/lung_luad/data/247f245e-7878-4d0f-b550-a4ea6be16fb4/56466297-cbb3-4015-b7e5-5dc386bafaee_gdc_realn_rehead.cram
luad.sub<-luad[,c(1,3:5)]

for (i in 1:nrow(luad.sub)){
  luad.sub$sample[i]<-strsplit(luad.sub$V5[i],"/")[[1]][8]
  
}
luad.sub$V2<-"fastq"
luad.sub$out1<-paste0("-1 ",luad.sub$sample,"_R1.fastq")
luad.sub$out2<-paste0("-2 ",luad.sub$sample,"_R2.fastq")

luad.sub.out<-luad.sub[,c(1,6,2,3,7,8,4,5)]

##salmon quant##
#&& salmon quant -p 8 -i transcripts_index -l A -1 reads1.fq -2 reads2.fq --validateMappings -o transcripts_quant

luad.sub.out$salmon<-"&& salmon quant -p 8 -i transcripts_index -l A"
luad.sub.out$output<-paste0("--validateMappings -o ",luad.sub.out$sample,".transcripts_quant")

luad.sub.out2<-luad.sub.out[,c(1:7,9,5,6,10)]
write.table(luad.sub.out2,"LUAD.cram2fq.salmon.cmd.txt",sep = " ",quote = FALSE,col.names = FALSE,row.names = FALSE)



