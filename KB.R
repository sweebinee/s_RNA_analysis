KB <- read.delim('ChimerKB - ChimerKB.tsv', header = F)
KB <- KB[-1,]

KB3 <- KB[KB$V2=="ChimerKB3.1",]

KB_new <- KB[KB$V3%in%c('ChimerPub 4.0 Symbol','ChimerPub 4.0 NonSymbol','ChimerPub 3.0'),]


length(unique(KB[KB$V4%in%c('ChimerPub 4.0 Symbol','ChimerPub 4.0 NonSymbol','ChimerPub 3.0'),]$V5))


KB_PUB <- KB[KB$V4%in%c('ChimerPub 4.0 Symbol','ChimerPub 4.0 NonSymbol','ChimerPub 3.0'),]

KB_PUB[KB_PUB$V30=='Pub'&KB_PUB$V31=='',]




library(venn)
Input = list(
        KB_PUB_F = KB_PUB_F,
        KB_PUB_sup_F = KB_PUB_sup_F,
        KB_PUB_Seq_sup_F = KB_PUB_Seq_sup_F
        )
venn(Input,ilab=TRUE, zcolor = "style",opacity = 0.3, size = 15, cexil = 0.85, cexsn = 1)


Pub<-read.delim('Pub.tsv', header = F, quote="")

unique(Pub[Pub$V3!='',]$V2)

unique(Pub[Pub$V15=='ChimerKB',]$V2)

###
Seq<-read.delim('Seq.tsv', header = F)
Seq<-Seq[-1,]

>unique(Seq$V3)
[1]                    FusionScan         STARFusion         PRADA(TumorFusion)
[5] Tophat-Fusion      Gao et al          ChimerDB 2.0       ChiTaRs  

length(unique(Seq[Seq$V3=="ChiTaRs",]$V5))
length(unique(Seq[Seq$V4=="TCGA_RNA-Seq",]$V5))

colnames(Seq) <- c('id','ChimerDB_Type','Source','webSource','Fusion_pair','H_gene','H_chr','H_position','H_strand','T_gene','T_chr','T_position','T_strand','Breakpoint','Genome_Build_Version','Cancertype','BarcodeID','Seed_read','Spanning_read','Junction_read','Frame','Chr_info','H_locus','H_Kinase','H_Oncogene','H_Tumor_suppressor','H_Receptor','H_Transcription_factor','T_locus','T_Kinase','T_Oncogene','T_Tumor_suppressor','T_Receptor','T_Transcription_factor','ChimerKB','ChimerPub')

FS <- unique(Seq[Seq$Source=='FusionScan','Fusion_pair'])
Gao <- unique(Seq[Seq$Source=='Gao et al','Fusion_pair'])
PRADA <- unique(Seq[Seq$Source=='PRADA(TumorFusion)','Fusion_pair'])
SF <- unique(Seq[Seq$Source=='STARFusion','Fusion_pair'])
TF <- unique(Seq[Seq$Source=='Tophat-Fusion','Fusion_pair'])
DB2 <- unique(Seq[Seq$Source=='ChimerDB 2.0','Fusion_pair'])
CT <- unique(Seq[Seq$Source=='ChiTaRs','Fusion_pair'])

library(venn)
Input = list(
        FusionScan = FS,
        Gao = Gao,
        PRADA = PRADA,
        STARFusion = SF,
        TophatFusion = TF
        )
venn(Input,ilab=TRUE, zcolor = "style",opacity = 0.3, size = 15, cexil = 0.85, cexsn = 1)


Seq_normal <- read.delim("Normal.tsv", header=F)

length(unique(Seq[Seq$V4=="ChiTaRs"&Seq$V35==""&Seq$V36=='',]$V5))

cancerType = c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')

for(i in cancerType){
	print(i)
	print(length(unique(Seq[Seq$V3=='Tophat-Fusion'&Seq$V16==i,'V17'])))
}


for(i in cancerType){
	print(i)
	print(length(unique(Seq[Seq$V4=='TCGA_RNA-Seq'&Seq$V16==i,'V17'])))
}

## table # of fusions
library(plotly)
library(dplyr)
library(MASS)
TCGA <- Seq[Seq$webSource=="TCGA_RNA-Seq",]


fusion <- as.data.frame(subset(TCGA,!duplicated(Fusion_pair)) %>% count(Cancertype)) 
#F_data <- F_data[order(F_data$fusion,decreasing=TRUE),]
F_data <- fusion[order(fusion$n,decreasing=TRUE),]
F_data$Cancertype <- factor(F_data$Cancertype, levels = F_data[["Cancertype"]])

p <- plot_ly(F_data, x = F_data$Cancertype, y = F_data$n, type = 'bar', name = '# of fusions',marker = list(color = 'rgb(56,118,29)')) %>%
  layout(xaxis = list(title = 'CancerTypes'),
    yaxis = list(title = '# of fusions'))


###
KB_F <- unique(KB$V5)
Pub_F <- unique(Pub$V4)
Seq_F <- unique(Seq$V5)

Input = list(
        KB = KB_F,
        Pub = Pub_F,
        Seq = Seq_F
        )
venn(Input,ilab=TRUE, zcolor = "style",opacity = 0.3, size = 15, cexil = 0.85, cexsn = 1)





## 
#HRS statistics
length(unique(HRS[HRS$ChimerPub=='ChimerPub'|HRS$ChimerKB=='ChimerKB',]$Fusion_pair))
length(unique(HRS[HRS$ChimerPub==''&HRS$ChimerKB=='',]$Fusion_pair))


length(unique(HRS[HRS$T_Oncogene=='Oncogene'|HRS$H_Oncogene=='Oncogene',]$Fusion_pair))
length(unique(HRS[HRS$T_Kinase!=''|HRS$H_Kinase!='',]$Fusion_pair))
length(unique(HRS[HRS$T_Tumor_suppressor!=''|HRS$H_Tumor_suppressor!='',]$Fusion_pair))
length(unique(HRS[HRS$T_Receptor!=''|HRS$H_Receptor!='',]$Fusion_pair))

F1 <- unique(HRS[HRS$T_Oncogene=='Oncogene'|HRS$H_Oncogene=='Oncogene',]$Fusion_pair)
F2 <- unique(HRS[HRS$T_Kinase!=''|HRS$H_Kinase!='',]$Fusion_pair)
F3 <- unique(HRS[HRS$T_Tumor_suppressor!=''|HRS$H_Tumor_suppressor!='',]$Fusion_pair)
F4 <- unique(HRS[HRS$T_Receptor!=''|HRS$H_Receptor!='',]$Fusion_pair)
Func <- union(union(union(F1,F2),F3),F4)
setdiff(unique(HRS$Fusion_pair),Func)


#ranke table
rank <- data.frame(table(HRS$Fusion_pair))
rank<-rank[order(rank$Freq,decreasing=TRUE),]

HRS_rank <- as.data.frame(HRS[,c('Source','Fusion_pair')])
HRS_rank_unique <- unique(HRS_rank)
HRS_rank_unique_count <- data.frame(table(HRS_rank_unique$Fusion_pair))
HRS_rank_unique_count[order(HRS_rank_unique_count$Freq,decreasing=TRUE),]
colnames(HRS_rank_unique_count)[2]<-'count'

rank_table <- merge(rank,HRS_rank_unique_count,by='Var1')

rank_table <- rank_table[rank_table$Freq!=0&rank_table$count!=0,]
rank_table[order(rank_table$count,rank_table$Freq,decreasing=TRUE),]

HRS[HRS$Fusion_pair=='ERBB2-PPP1R1B',]
