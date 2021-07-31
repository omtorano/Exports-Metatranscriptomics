BiocManager::install('tximport')
library(tximport)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(stringr)


samples<-list.files(path="//pine/scr/o/m/omtorano/alignment/trinity/", full.names=T)
files<-file.path(samples,"quant.sf")
names(files)<-str_replace(samples, "/pine/scr/o/m/omtorano/alignment/trinity/${input1}_quants_trinity","")%>%str_replace(".salmon","")
tx2gene <- read.delim("/proj/marchlab/projects/PUPCYCLE_2019/mega_annotations/annotations/kegg_annotations/mega_kegg.tsv")
txi.keggannot.transcriptlevel<-tximport(files,type="salmon",tx2gene=tx2gene[,c('query', 'KO')], txOut = TRUE)
txi.keggannot.transcriptlevel.1.rds <- saveRDS(txi.keggannot.transcriptlevel, "txi.keggannot.transcriptlevel.1.rds")
kegg_tl<-readRDS('txi.keggannot.transcriptlevel.1.rds')
rownames(kegg_tl) <- kegg_tl$X
colnames(kegg_tl)[colnames(kegg_tl) == 'X'] <-'query'
#phylodb annotations
tx2gene1<-read.delim("/proj/marchlab/projects/PUPCYCLE_2019/mega_annotations/annotations/phylodb_annotations/megaphylodb.tsv")
#import
txi.phylodb.transcriptlevel<-tximport(files, type="salmon", tx2gene =tx2gene1[,c("TrinityID", "Organism")],txOut=TRUE)
#save r data
txi.phylodb.transcriptlevel.1.rds<-saveRDS(txi.phylodb.transcriptlevel, "txi.phylodb.transcriptlevel.1.csv")

phylo<-read.csv('sub.phylo.csv', header=T)
kegg<-read.csv('sub.phylo.csv', header=T)
rownames(kegg)<-kegg$X
colnames(kegg)[colnames(kegg)=='X']<-'TrinityID'

p<-read.delim('/proj/marchlab/projects/PUPCYCLE_2019/mega_annotations/annotations/phylodb_annotations/megaphylodb.tsv')
k<-read.delim('/proj/marchlab/projects/PUPCYCLE_2019/mega_annotations/annotations/kegg_annotations/mega_kegg.tsv')
colnames(keggannot)[colnames(keggannot)=='query']<-'TrinityID'
tpmphy<-merge(kegg, p, by ='TrinityID')
tpmphy<-merge(tpmphy,k, by ='TrinityID')