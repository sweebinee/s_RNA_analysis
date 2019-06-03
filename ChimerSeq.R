
###############################################################################################
## Finding Highly reliable seq
### venndiagram of unique fusion pairs from different programs
##R
setwd("/media/subin/fec664d7-6f05-4f07-86cb-b7b76cd331b2/ChimerDB4")

biocLite("VennDiagram")
library(VennDiagram)

Seq <- read.delim("ChimerSeq.txt",header =F, stringsAsFactors = F)
colnames(Seq) <- c('id','Source','Fusion_pair','H_gene','H_chr','H_position','H_strand','T_gene','T_chr','T_position','T_strand','Breakpoint','Genome_Build_Version','Cancertype','BarcodeID','Seed_read','Spanning_read','Junction_read','Frame','Chr_info','H_locus','H_Kinase','H_Oncogene','H_Tumor_suppressor','H_Receptor','H_Transcription_factor','T_locus','T_Kinase','T_Oncogene','T_Tumor_suppressor','T_Receptor','T_Transcription_factor','ChimerKB','ChimerPub')
> table(Seq$Source)
        FusionScan          Gao et al PRADA(TumorFusion)         STARFusion	Tophat-Fusion 
             15129              25664              20731              56079          4482

FS <- Seq[Seq$Source=='FusionScan','Fusion_pair']
Gao <- Seq[Seq$Source=='Gao et al','Fusion_pair']
PRADA <- Seq[Seq$Source=='PRADA(TumorFusion)','Fusion_pair']
SF <- Seq[Seq$Source=='STARFusion','Fusion_pair']
TF <- gsub('_','-',Seq[Seq$Source=='Tophat-Fusion','Fusion_pair'])

library(venn)
Input = list(
        FusionScan = FS,
        Gao = Gao,
        PRADA = PRADA,
        STARFusion = SF,
        TophatFusion = TF
        )
venn(Input,ilab=TRUE, zcolor = "style",opacity = 0.3, size = 15, cexil = 0.85, cexsn = 1)

###############################################################################################
## bar plot
library(plotly)
library(dplyr)
library(MASS)

cancerType = unique(Seq$Cancertype)
FS <- as.data.frame(Seq[Seq$Source=='FusionScan',] %>% count(Cancertype))
SF <- as.data.frame(Seq[Seq$Source=='STARFusion',] %>% count(Cancertype))
PRADA <- as.data.frame(Seq[Seq$Source=='PRADA(TumorFusion)',] %>% count(Cancertype))
Gao <- as.data.frame(Seq[Seq$Source=='Gao et al',] %>% count(Cancertype))
TF <- as.data.frame(Seq[Seq$Source=='Tophat-Fusion',] %>% count(Cancertype))
data <- merge(merge(merge(merge(FS,SF,by='Cancertype',all=TRUE),PRADA,by='Cancertype', all=TRUE),Gao,by='Cancertype', all=TRUE),TF,by='Cancertype', all=TRUE)
colnames(data) <- c("Cancertype","FusionScan","STARFusion","PRADA","Gao","TophatFusion")
is.na(data) <- 0

p <- plot_ly(data, x = ~Animals, y = ~SF_Zoo, type = 'bar', name = 'SF Zoo') %>%
  add_trace(y = ~LA_Zoo, name = 'LA Zoo') %>%
  layout(yaxis = list(title = 'Count'), barmode = 'group')