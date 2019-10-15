cat ChimerSeq_* > 1.txt
sed '/Fusion_pair/d' 1.txt > ChimerSeq.txt

###############################################################################################
## Finding Highly reliable seq
### venndiagram of unique fusion pairs from different programs
##R

#setwd("/media/subin/fec664d7-6f05-4f07-86cb-b7b76cd331b2/ChimerDB4/ChimerSeq")
#setwd("/home/subin/Desktop/ChimerDB")
#setwd("/media/subin/fec664d7-6f05-4f07-86cb-b7b76cd331b2/ChimerDB4")
#setwd('/Users/subin/Desktop/ChimerDB04')
#setwd('/Users/subin/Downloads')


Seq <- read.delim("ChimerSeq4_final_result.csv",header =T, stringsAsFactors = F,sep=',')
#colnames(Seq) <- c('id','ChimerDB_Type','Source','webSource','Fusion_pair','H_gene','H_chr','H_position','H_strand','T_gene','T_chr','T_position','T_strand','Breakpoint','Genome_Build_Version','Cancertype','BarcodeID','Seed_read','Spanning_read','Junction_read','Frame','Chr_info','H_locus','H_Kinase','H_Oncogene','H_Tumor_suppressor','H_Receptor','H_Transcription_factor','T_locus','T_Kinase','T_Oncogene','T_Tumor_suppressor','T_Receptor','T_Transcription_factor','ChimerKB','ChimerPub')
#Seq <- Seq[,-36]

FS <- unique(Seq[Seq$Source=='FusionScan'&Seq$Cancertype!='Non-Cancer','Fusion_pair'])
Gao <- unique(Seq[Seq$Source=='Gao et al','Fusion_pair'])
PRADA <- unique(Seq[Seq$Source=='PRADA(TumorFusion)','Fusion_pair'])
SF <- unique(Seq[Seq$Source=='STARFusion'&Seq$Cancertype!='Non-Cancer','Fusion_pair'])
TF <- unique(Seq[Seq$Source=='Tophat-Fusion','Fusion_pair'])
DB2 <- unique(Seq[Seq$Source=='ChimerDB 2.0','Fusion_pair'])
CT <- unique(Seq[Seq$Source=='ChiTaRs','Fusion_pair'])
NC<- unique(Seq[Seq$Cancertype=='Non-Cancer','Fusion_pair'])

total_F <- union(union(union(union(union(union(FS, Gao),PRADA),SF),TF),DB2),CT)
Seq_FP <- total_F
length(unique(Seq$Fusion_pair))

KB_s <- unique(Seq[Seq$ChimerKB=='ChimerKB','Fusion_pair'])
Pub_s <- unique(Seq[Seq$ChimerPub=='ChimerPub','Fusion_pair'])

TCGA_novel <- unique(Seq[Seq$Fusion_pair %in% union_fusion & Seq$ChimerKB=='','Fusion_pair'])
TCGA_known <- unique(Seq[Seq$Fusion_pair %in% union_fusion & Seq$ChimerKB=='ChimerKB','Fusion_pair'])

CT_novel <- unique(Seq[Seq$Fusion_pair %in% CT & Seq$ChimerKB=='','Fusion_pair'])


library(venn)
Input = list(
        FusionScan = FS,
        TCGA_FAWG = Gao,
        TumorFusions = PRADA,
        STARFusion = SF,
        TophatFusions = TF
        )
venn(Input,ilab=TRUE, zcolor = "style",opacity = 0.35, size = 15, cexil = 1.1, cexsn = 0.5,device=cairo_ps)

union_fusion = union(union(union(union(FS, Gao),PRADA),SF),TF)
FS_only = setdiff(setdiff(setdiff(setdiff(FS, Gao),PRADA),SF),TF)
SF_only = setdiff(setdiff(setdiff(setdiff(SF, Gao),PRADA),FS),TF)
PRADA_only = setdiff(setdiff(setdiff(setdiff(PRADA, Gao),SF),FS),TF)
Gao_only = setdiff(setdiff(setdiff(setdiff(Gao, PRADA),SF),FS),TF)
TF_only = setdiff(setdiff(setdiff(setdiff(TF, PRADA),SF),FS),Gao)

HRS_fusion = setdiff(setdiff(setdiff(setdiff(setdiff(union_fusion,FS_only),SF_only),PRADA_only),Gao_only),TF_only)

HRS <- Seq[Seq$Fusion_pair %in% HRS_fusion,]

write.table(HRS, "Highly_Reliable_Seq.txt", sep = '\t')

###
#기여도 테스트
table <- data.frame(HRS_fusion)
table$SF <- 0
table$FAWG <- 0
table$PRADA <- 0
table$FS <- 0
table$TF <- 0
table$sum <- 0

for(i in c(1:nrow(table))){
  if(table$HRS_fusion[i] %in% SF){
    table$SF[i] <- 1
  }
  if(table$HRS_fusion[i] %in% Gao){
    table$FAWG[i] <- 1
  }
  if(table$HRS_fusion[i] %in% PRADA){
    table$PRADA[i] <- 1
  }
  if(table$HRS_fusion[i] %in% FS){
    table$FS[i] <- 1
  }
  if(table$HRS_fusion[i] %in% TF){
    table$TF[i] <- 1
  }
  table$sum[i] <- table$SF[i]+table$FAWG[i]+table$PRADA[i]+table$FS[i]+table$TF[1]
}

SF_HRS <- length(table[table$SF==1,'HRS_fusion'])
FAWG_HRS <- length(table[table$FAWG==1,'HRS_fusion'])
PRADA_HRS <- length(table[table$PRADA==1,'HRS_fusion'])
FS_HRS <- length(table[table$FS==1,'HRS_fusion'])
TF_HRS <- length(table[table$TF==1,'HRS_fusion'])

eFAWG <- table[table$FAWG==1&table$sum==2,]
ePRADA <- table[table$PRADA==1&table$sum==2,]
eFS <- table[table$FS==1&table$sum==2,]
eTF <- table[table$TF==1&table$sum==2,]

eSF <- table[table$SF==1&table$sum==2,]
eFAWG <- table[table$FAWG==1&table$sum==2,]
ePRADA <- table[table$PRADA==1&table$sum==2,]
eFS <- table[table$FS==1&table$sum==2,]
eTF <- table[table$TF==1&table$sum==2,]

length(eSF$HRS_fusion)
length(eFAWG$HRS_fusion)
length(ePRADA$HRS_fusion)
length(eFS$HRS_fusion)
length(eTF$HRS_fusion)



###############################################################################################
## KB, Pub, Seq plot
KB <- read.delim("Final_ChimerKB_for_db - Final_ChimerKB_for_db.tsv",header =T, stringsAsFactors = F)
KB_FP <- unique(KB$Fusion_pair)

Pub <- read.delim("ChimerPub_Final2.csv",header =T, stringsAsFactors = F, sep=',')
Pub <- Pub[,c(1:17)]
colnames(Pub) <- c('num','Fusion_pair','Translocation','5gene','3gene','PMID','score','disase','validation','kinase','oncogene','TumorSuppressor','','TF','ChimerKB','ChimerSeq','ChimerSeqPlus')
Pub_FP <- unique(Pub$Fusion_pair)

ChimerDB_FP <- union(union(KB_FP,Pub_FP),Seq_FP)

Input = list(
        KB = KB_FP,
        Pub = Pub_FP,
        Seq = Seq_FP
        )
venn(Input,ilab=TRUE, zcolor = "style",opacity = 0.35, size = 15, cexil = 1.1, cexsn = 0.5,device=cairo_ps)



###############################################################################################
## bar plot
## # of samples
library(plotly)
library(dplyr)
library(MASS)

cancerType = unique(Seq$Cancertype)
FS <- as.data.frame(subset(Seq[Seq$Source=='FusionScan',],!duplicated(BarcodeID))[,c('BarcodeID','Cancertype')] %>% count(Cancertype))
SF <- as.data.frame(subset(Seq[Seq$Source=='STARFusion',],!duplicated(BarcodeID))[,c('BarcodeID','Cancertype')] %>% count(Cancertype))
PRADA <- as.data.frame(subset(Seq[Seq$Source=='PRADA(TumorFusion)',],!duplicated(BarcodeID))[,c('BarcodeID','Cancertype')] %>% count(Cancertype))
Gao <- as.data.frame(subset(Seq[Seq$Source=='Gao et al',],!duplicated(BarcodeID))[,c('BarcodeID','Cancertype')] %>% count(Cancertype))
TF <- as.data.frame(subset(Seq[Seq$Source=='Tophat-Fusion',],!duplicated(BarcodeID))[,c('BarcodeID','Cancertype')] %>% count(Cancertype))
data <- merge(merge(merge(merge(FS,SF,by='Cancertype',all=TRUE),PRADA,by='Cancertype', all=TRUE),Gao,by='Cancertype', all=TRUE),TF,by='Cancertype', all=TRUE)
colnames(data) <- c("Cancertype","FusionScan","STARFusion","PRADA","Gao","TophatFusion")
data[is.na(data)] <- 0
data[data$Cancertype=="LAML",'STARFusion'] <- length(unique(substr(Seq[Seq$Source=='STARFusion'&Seq$Cancertype=='LAML','BarcodeID'],1,12)))

data$samples <- c(79,408,1026,304,36,458,48,184,167,521,66,533,290,179,516,318,517,501,87,430,179,497,167,259,446,419,505,120,557,57,80,150,178)

p <- plot_ly(data, x = ~Cancertype, y = ~samples, type = 'bar', name = '# of TCGA samples',marker = list(color = 'rgb(204,204,204)')) %>%
  add_trace(y = ~FusionScan, name = 'FusionScan',marker = list(color = 'rgb(255,102,102)')) %>%
  add_trace(y = ~STARFusion, name = 'STAR-Fusion',marker = list(color = 'rgb(255,178,102)')) %>%
  add_trace(y = ~PRADA, name = 'PRADA',marker = list(color = 'rgb(255,255,102)')) %>%
  add_trace(y = ~Gao, name = 'Gao',marker = list(color = 'rgb(178,255,102)')) %>%
  add_trace(y = ~TophatFusion, name = 'Tophat-Fusion',marker = list(color = 'rgb(102,178,255)')) %>%
  layout(yaxis = list(title = 'Count'), barmode = 'group')
##
## # of fusions
setwd('/media/subin/fec664d7-6f05-4f07-86cb-b7b76cd331b2/ChimerDB4/Final_data')

Seq <- read.delim("ChimerSeq_TCGA.txt",header =F, stringsAsFactors = F,sep=',')
colnames(Seq) <- c('id','ChimerDB_Type','Source','webSource','Fusion_pair','H_gene','H_chr','H_position','H_strand','T_gene','T_chr','T_position','T_strand','Breakpoint','Genome_Build_Version','Cancertype','BarcodeID','Seed_read','Spanning_read','Junction_read','Frame','Chr_info','H_locus','H_Kinase','H_Oncogene','H_Tumor_suppressor','H_Receptor','H_Transcription_factor','T_locus','T_Kinase','T_OncogeneT_Tumor_suppressor','T_Receptor','T_Transcription_factor','ChimerKB','ChimerPub')
Seq <- Seq[,-36]

F_data <- data[,c('Cancertype','samples')]
F_data$fusion <- as.data.frame(subset(Seq,!duplicated(Fusion_pair)) %>% count(Cancertype))[,2] 
#F_data <- F_data[order(F_data$fusion,decreasing=TRUE),]
F_data <- F_data[order(F_data$samples,decreasing=TRUE),]
F_data$Cancertype <- factor(F_data$Cancertype, levels = F_data[["Cancertype"]])

p <- plot_ly(F_data, x = F_data$Cancertype, y = F_data$fusion, type = 'bar', name = '# of fusions',marker = list(color = 'rgb(204,0,0)')) %>%
  layout(xaxis = list(title = 'CancerTypes'),
    yaxis = list(title = '# of fusions'))


p <- plot_ly(F_data, x = F_data$Cancertype, y = F_data$samples, type = 'bar', name = '# of TCGA samples',marker = list(color = 'rgb(204,204,204)')) %>%
  layout(xaxis = list(title = 'CancerTypes'),
    yaxis = list(title = '# of TCGA samples'))



p <- plot_ly(F_data, x = F_data$Cancertype, y = F_data$fusion, type = 'bar', name = '# of fusions',marker = list(color = 'rgb(204,0,0)')) %>%
  add_trace(y = F_data$samples, name = '# of TCGA samples',marker = list(color = 'rgb(204,204,204)'), yaxis='y2') %>%
  layout(yaxis2 = list(title = 'samples',overlaying = 'y',side='right'),
    xaxis = list(title = 'CancerTypes'),
    yaxis = list(title = 'fusions'),
    margin = list(b=0),
    barmode = '')


import plotly
import plotly.graph_objs as go
cancerType = ['BRCA','STAD','SARC','OV','LUSC','PRAD','LUAD','BLCA','SKCM','ESCA','HNSC','UCEC','GBM','LIHC','LGG','COAD','CESC','KIRC','UCS','ACC','THCA','PAAD','READ','KIRP','PCPG','MESO','TGCT','LAML','CHOL','THYM','DLBC','KICH','UVM']
trace = [go.Bar(x = cancerType,y = [8643,6986,3037,2916,2867,2637,2540,2173,2070,1873,1821,1683,1424,1351,1342,1059,849,630,554,487,430,376,353,315,245,227,209,186,111,100,92,74,58],name = '# of fusions'),
  go.Bar(x = ['BRCA'], y = [0], name = 'y_dummy', showlegend=False),
  go.Bar(x = ['BRCA'], y = [0], name = 'y2_dummy', yaxis = 'y2', showlegend=False),
  go.bar(x = cancerType, y = [1026,505,446,430,501,167,517,408,419,184,521,80,167,318,516,458,304,533,150,79,557,179,259,290,497,87,120,179,36,57,48,66,178], yaxis = 'y2', name = '# of TCGA samples')]
layout = go.layout(barmode = 'group',
  yaxis = dict(title = '# of fusions', overlaying='y2'),
  yaxis2 = dict(title = '# of TCGA samples', side = 'right'))
fig = go.Figure(data = traces, layout = layout)
plotly.offline.iplot(fig)
###############################################################################################



