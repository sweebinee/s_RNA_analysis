##R
setwd("/media/subin/fec664d7-6f05-4f07-86cb-b7b76cd331b2/ChimerDB4/ChimerKB")

KB <- read.delim("ChimerDB3.0_ChimerKB.tsv",header =T, stringsAsFactors = F, sep = '\t')

Pub4 <- read.delim("190620_ChimerPub_symbol.tsv",header =T, stringsAsFactors = F, sep = '\t')
Pub3 <- read.delim("ChimerDB3.0_ChimerPub.tsv",header =T, stringsAsFactors = F, sep = '\t')
Pub4 <- Pub4[,c('H_gene','T_gene','Score','Validation','Translocation')]
Pub3 <- Pub3[,c('H_gene','T_gene','Score','Validation','Translocation')]
Pub4$ver <- rep(4,nrow(Pub4))
Pub3$ver <- rep(4,nrow(Pub3))
Pub <- rbind(Pub3,Pub4)

Seq <- read.delim("/media/subin/fec664d7-6f05-4f07-86cb-b7b76cd331b2/ChimerDB4/ChimerSeq/ChimerSeq.txt",header =F, stringsAsFactors = F)
colnames(Seq) <- c('id','ChimerDB_Type','Source','webSource','Fusion_pair','H_gene','H_chr','H_position','H_strand','T_gene','T_chr','T_position','T_strand','Breakpoint','Genome_Build_Version','Cancertype','BarcodeID','Seed_read','Spanning_read','Junction_read','Frame','Chr_info','H_locus','H_Kinase','H_Oncogene','H_Tumor_suppressor','H_Receptor','H_Transcription_factor','T_locus','T_Kinase','T_OncogeneT_Tumor_suppressor','T_Receptor','T_Transcription_factor','ChimerKB','ChimerPub')

#Pub99 <- Pub[Pub$Score>=0.9,]

KB_F <- KB$Fusion_pair
Pub_F <- paste(Pub$H_gene,Pub$T_gene,sep='_')
Seq_F <- paste(Seq$H_gene,Seq$T_gene,sep='_')

library(venn)
Input = list(
	    Pub_3.0_4.0 = Pub_F,
        KB_3.0 = KB_F,
        Seq_4.0 = Seq_F
        )
venn(Input,ilab=TRUE, zcolor = "style",opacity = 0.3, size = 15, cexil = 1.5, cexsn = 1.5)

Pub_only = setdiff(setdiff(Pub_F,KB_F),Seq_F)
Pub_only_H = vector()
Pub_only_T = vector()
for(i in c(1:length(Pub_only))){
	Pub_only_H <- append(Pub_only_H,strsplit(Pub_only,'_')[[i]][1])
	Pub_only_T <- append(Pub_only_T,strsplit(Pub_only,'_')[[i]][2])
}
Pub_only_F = data.frame(Pub_only_H, Pub_only_T)
colnames(Pub_only_F) <- c('H_gene','T_gene')

Pub_only.df <- merge(Pub,Pub_only_F,by= c('H_gene','T_gene'), all=FALSE)

#test <- paste(Pub_only.df$H_gene,Pub_only.df$T_gene,sep='_')
#> length(unique(test))
#[1] 4005

###############Pub_only validation step있는가???
Pub_only_val <- Pub_only.df[Pub_only.df$Validation!="",]

#test <- paste(Pub_only_val$H_gene,Pub_only_val$T_gene,sep='_')
#length(unique(test))

###############Pub_only translocation 정보 있는가???
Pub_only_val_trans <- Pub_only_val[Pub_only_val$Translocation!="",]

#test <- paste(Pub_only_val_trans$H_gene,Pub_only_val_trans$T_gene,sep='_')
#length(unique(test))

####################################################################3
####Pub n Seq 108 fusions
PubSeq = setdiff(intersect(Pub_F,Seq_F),KB_F)
PubSeq_H = vector()
PubSeq_T = vector()
for(i in c(1:length(PubSeq))){
	PubSeq_H <- append(PubSeq_H,strsplit(PubSeq,'_')[[i]][1])
	PubSeq_T <- append(PubSeq_T,strsplit(PubSeq,'_')[[i]][2])
}
PubSeq_F = data.frame(PubSeq_H, PubSeq_T)
colnames(PubSeq_F) <- c('H_gene','T_gene')

PubSeq.df <- merge(Pub,PubSeq_F,by= c('H_gene','T_gene'), all=FALSE)

PubSeq_val <- PubSeq.df[PubSeq.df$Validation!="",]
#test <- paste(PubSeq_val$H_gene,PubSeq_val$T_gene,sep='_')
#length(unique(test))


test <- Pub[Pub$Validation==""&Pub$Translocation=="",]
test <- paste(test$H_gene,test$T_gene,sep='_')