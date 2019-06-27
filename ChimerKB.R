##R
setwd("/media/subin/fec664d7-6f05-4f07-86cb-b7b76cd331b2/ChimerDB4/ChimerKB")

KB <- read.delim("ChimerDB3_new_2.tsv",header =T, stringsAsFactors = F, sep = '\t')

Pub4_S <- read.delim("190621_ChimerPub_symbol.tsv",header =T, stringsAsFactors = F, sep = '\t',quote = "")
Pub4_NS <- read.delim("190625_ChimerPub_nonsymbol_.tsv",header =T, stringsAsFactors = F, sep = '\t')
Pub3 <- read.delim("ChimerDB3.0_ChimerPub.tsv",header =T, stringsAsFactors = F, sep = '\t')

Pub4_S <- Pub4_S[,c('H_gene','T_gene','Score','Validation','Translocation')]
Pub4_NS <- Pub4_NS[,c('H_gene','T_gene','Score','Validation','Translocation')]
Pub3 <- Pub3[,c('H_gene','T_gene','Score','Validation','Translocation')]

Pub4_S$ver <- rep(4,nrow(Pub4_S))
Pub4_NS$ver <- rep(4,nrow(Pub4_NS))
Pub3$ver <- rep(4,nrow(Pub3))

Pub <- rbind(Pub3,Pub4_S,Pub4_NS)

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

####################################################################
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

####################################################################
##Chimer DB 4.0
####################################################################

Pub4_S <- read.delim("190621_ChimerPub_symbol.tsv",header =T, stringsAsFactors = F, sep = '\t',quote = "")
Pub4_NS <- read.delim("190625_ChimerPub_nonsymbol_.tsv",header =T, stringsAsFactors = F, sep = '\t')
Pub3 <- read.delim("ChimerDB3.0_ChimerPub.tsv",header =T, stringsAsFactors = F, sep = '\t')
colnames(Pub3)[16] <- 'ChimerSeq_supported'

Pub3$Source <- rep("ChimerPub 3.0",nrow(Pub3))
Pub3$webSource <- rep("ChimerPub 3.0",nrow(Pub3))
Pub4_S$Source <- rep("ChimerPub 4.0 Symbol",nrow(Pub4_S))
Pub4_S$webSource <- rep("ChimerPub 4.0 Symbol",nrow(Pub4_S))
Pub4_NS$Source <- rep("ChimerPub 4.0 NonSymbol",nrow(Pub4_NS))
Pub4_NS$webSource <- rep("ChimerPub 4.0 NonSymbol",nrow(Pub4_NS))

Pub4_S$Fusion_pair <- paste(Pub4_S$H_gene,Pub4_S$T_gene,sep='-')
Pub4_NS$Fusion_pair <- paste(Pub4_NS$H_gene,Pub4_NS$T_gene,sep='-')

Pub3$ChimerPub_supported <- rep("Pub",nrow(Pub3))
Pub4_S$ChimerPub_supported <- rep("Pub",nrow(Pub4_S))
Pub4_NS$ChimerPub_supported <- rep("Pub",nrow(Pub4_NS))

Pub4 <- merge(Pub4_S,Pub4_NS,by=c('Source','webSource','Fusion_pair','H_gene','Translocation','T_gene','PMID','Disease','Validation','ChimerPub_supported'),all=TRUE)
Pub4 <- rename(Pub4,'Transcription_Factor'='Transcription_Faction')
Pub <- merge(Pub3, Pub4, by=c('Source','webSource','Fusion_pair','H_gene','Translocation','T_gene','PMID','Disease','Validation','Kinase','Oncogene','Tumor_suppressor','Receptor','Transcription_Factor','ChimerPub_supported'),all=TRUE)

#Pub_only_val_trans : 563 fusions
Pub_only_val_trans_F <- paste(Pub_only_val_trans$H_gene,Pub_only_val_trans$T_gene,sep='-')
Pub_only_KB <- Pub[Pub$Fusion_pair%in%Pub_only_val_trans_F,]
#PubSeq_val : 57 fusions
PubSeq_val_F <- paste(PubSeq_val$H_gene,PubSeq_val$T_gene,sep='-')
PubSeq_KB <- Pub[Pub$Fusion_pair%in%PubSeq_val_F,]
PubSeq_KB$ChimerSeq_supported <- rep("Seq",nrow(PubSeq_KB))

PubSeq_KB$H_chr <- unique(Seq[PubSeq_KB$Fusion_pair==Seq$Fusion_pair,'H_chr'])
PubSeq_KB$T_chr <- unique(Seq[PubSeq_KB$Fusion_pair==Seq$Fusion_pair,'T_chr'])

add <- merge(Pub_only_KB,PubSeq_KB,by=c('Source','webSource','Fusion_pair','H_gene','Translocation','T_gene','PMID','Disease','Validation','Kinase','Oncogene','Tumor_suppressor','Receptor','Transcription_Factor','ChimerPub_supported','ChimerSeq_supported'),all=TRUE)
add <- add[,c(1:16,53,54)]
add$id <- rep("",nrow(add))
add$exon_info <- rep("",nrow(add))
add$H_position <- rep("",nrow(add))
add$H_strand <- rep("",nrow(add))
add$T_position <- rep("",nrow(add))
add$T_strand <- rep("",nrow(add))
add$Breakpoint_Type <- rep("",nrow(add))
add$Genome_Build_Version <- rep("",nrow(add))

newKB <- add[,c(19,1,2,3,4,5,17,20,21,22,6,18,23,24,25,26,7,8,9,10,11,12,13,14,15,16)]

write.table(newKB, "ChimerKB4.0_additional.txt", sep = '\t')

####################################################
####################################################
KB <- read.delim("ChimerDB3_new_2.tsv",header =T, stringsAsFactors = F, sep = '\t')
KB$Fusion_pair <- paste(KB$H_gene,KB$T_gene,sep='-')
colnames(KB)[6] <- 'Translocation'
colnames(newKB)[8] <- 'exon.info'
colnames(newKB)[25] <- 'ChimerPub.supported'
colnames(newKB)[26] <- 'ChimerSeq.supported'

KB4  <- rbind(KB,newKB) 

KB_F <- KB4$Fusion_pair
Pub_F <- paste(Pub$H_gene,Pub$T_gene,sep='-')
Seq_F <- paste(Seq$H_gene,Seq$T_gene,sep='-')

library(venn)
Input = list(
	    Pub_3.0_4.0 = Pub_F,
        KB_4.0 = KB_F,
        Seq_4.0 = Seq_F
        )
venn(Input,ilab=TRUE, zcolor = "style",opacity = 0.3, size = 15, cexil = 1.5, cexsn = 1.5)

ChimerDB4_F <- union(union(KB_F,Pub_F),Seq_F) 
####################################################
####################################################

setwd('/media/subin/fec664d7-6f05-4f07-86cb-b7b76cd331b2/ChimerDB4/Final_data')

SF <- read.delim("ChimerSeq_STAR_all2_rename_final.csv",header =F, stringsAsFactors = F, sep = ',')
FS <- read.delim("ChimerSeq_FS_all2_rename_final.csv",header =F, stringsAsFactors = F, sep = ',')
cell <- read.delim("ChimerSeq_cell_report_all2_rename_final.csv",header =F, stringsAsFactors = F, sep = ',')
PRADA <- read.delim("ChimerSeq_PRADA_all2_rename_final.csv",header =F, stringsAsFactors = F, sep = ',')
TF <- read.delim("ChimerSeq_Tophat_all2_rename_final.csv",header =F, stringsAsFactors = F, sep = ',')
CHI <- read.delim("ChimerSeq_Chitars_chimer2_all3_final.csv",header =F, stringsAsFactors = F, sep = ',')

Seq <- rbind(rbind(rbind(rbind(rbind(SF,FS),cell),PRADA),TF),CHI)

> length(unique(Seq[Seq$V36=='ChimerPub','V5']))
[1] 264
> length(unique(Seq[Seq$V35=='ChimerKB','V5']))
[1] 235
> length(unique(Seq[Seq$V35==''&Seq$V36=='','V5']))
[1] 66087

> length(unique(Seq[Seq$V35==''&Seq$V36==''&Seq$V4=='TCGA_RNA-Seq','V5']))
[1] 49911
> length(unique(Seq[Seq$V35==''&Seq$V36==''&Seq$V4=='ChiTaRs','V5']))
[1] 16131
> length(unique(Seq[Seq$V35==''&Seq$V36==''&Seq$V4=='ChimerDB 2.0','V5']))
[1] 142

###############
P1 <- read.delim("Final_Pub_1st.csv",header =F, stringsAsFactors = F, sep = ',',quote = "")
P2 <- read.delim("Final_Pub_2nd.csv",header =F, stringsAsFactors = F, sep = ',',quote = "")
P3 <- read.delim("Final_Pub_3rd.csv",header =F, stringsAsFactors = F, sep = ',',quote = "")

Pub <- rbind(rbind(P1,P2),P3)