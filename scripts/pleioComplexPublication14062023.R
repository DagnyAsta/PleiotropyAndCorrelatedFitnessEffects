#This is an edited script from:
#analysis for portugal D.mel populations
##written by Chaimae and edited by SKHsu
##20190801_V3
#Edidet by Dagny Asta Runarsdottir in order to compare d.mel and d.sim transcriptome response
#For manuscript: 


rm(list=ls())
library(edgeR)
library(topGO)
library(emmeans)
library(stringr)
library(ggplot2)
library(Rmisc)
library(poolr)
library(ggsignif)
library(viridis)
library(dplyr)
library(ggbreak)
setwd("/Volumes/Data/Dagny/Dropbox (PopGen)/DagnyPhD/Publications/PleiotropyInDifferentEnvironments/scripts/")

#public data from Zhao et.al. 2015 download: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP047141&o=acc_s%3Aa


####InputRNASeqCountTables####
##simulans
simcounts = read.csv("../data/laboratory/count_table_Dsim_Pt_redo_SK.csv", header = T, stringsAsFactors=F,row.names = 1)
simcounts <- simcounts[,c(4:6,8,10:15)] #exlude base
head(simcounts)
nrow(simcounts)
simcounts_use=simcounts[apply(cpm(simcounts),1,function(x) !sum(x<1)>=1),]
nrow(simcounts_use)
simevo=str_sub(colnames(simcounts_use),1,-3)#translate labels to evolutionary states (B:ancestral; H: evolved)
simrep <- str_sub(colnames(simcounts_use),-1)


#melanogaster
melcounts = read.csv("../data/laboratory/Portugal_Dmel_readcounts_r5.49.csv", header = T, stringsAsFactors=F,row.names = 1)
head(melcounts)
nrow(melcounts)
melcounts <- melcounts[,c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14)] #take only cold and hot
selectSamples <- colnames(melcounts)
melcountsmRNA = read.csv("../data/laboratory/readcountsmRNA_r5.49.csv", header = T, stringsAsFactors=F,row.names = 1) #recount by mRNA
colnames(melcountsmRNA)
melcountsmRNA <- melcountsmRNA[,selectSamples]
head(melcountsmRNA)
nrow(melcountsmRNA)
melcounts_use=melcountsmRNA[apply(cpm(melcountsmRNA),1,function(x) !sum(x<1)>=1),]
nrow(melcounts_use)
melevo=rep(c("cold","hot"),5)


#normalise the genes extracted
simy2=DGEList(counts=simcounts_use,group = simevo)
simy2=calcNormFactors(simy2)

mely2=DGEList(counts=melcounts_use,group = melevo)
mely2=calcNormFactors(mely2)

nrow(mely2)
nrow(simy2)

####VarianceDecoposition####
plotMDS(simy2)
simpca=(prcomp(t(log(cpm(simy2)))))
simpcaPlot=as.data.frame(simpca$x)
simpcaPlot$evolution=simevo
ve=simpca$sdev^2/sum(simpca$sdev^2)
tiff("../output/pcaSim.tiff", height = 10, width = 10, res = 300, units = "cm")
ggplot(simpcaPlot, aes(x = PC1, y = PC2, color = evolution)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title = "A") +
  theme(legend.title = element_blank()) +
  scale_color_manual(labels = c("Cold Evolved", "Hot Evolved"), values = viridis(n=4)[c(1,3)])+
  xlab(paste0("PC1 (",round(ve[1]*100,2),"%)")) +
  ylab(paste0("PC2 (",round(ve[2]*100,2),"%)")) +
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),  
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) 
dev.off()
  
plotMDS(mely2)
melpca=prcomp(t(log(cpm(mely2))))
melpcaPlot=as.data.frame(melpca$x)
melpcaPlot$evolution=melevo
ve=melpca$sdev^2/sum(melpca$sdev^2)
tiff("../output/pcaMel.tiff", height = 10, width = 10, res = 300, units = "cm")
ggplot(melpcaPlot, aes(x = PC1, y = PC2, color = evolution)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title = "B") +
  theme(legend.title = element_blank()) +
  scale_color_manual(labels = c("Cold Evolved", "Hot Evolved"), values =viridis(n=4)[c(2,4)])+
  xlab(paste0("PC1 (",round(ve[1]*100,2),"%)")) +
  ylab(paste0("PC2 (",round(ve[2]*100,2),"%)")) +
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),  
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) 
dev.off()



####ExperimentalData####
##simulans
simModelDesign=model.matrix(~0+simevo) #is the intercept set to 0 in able to make constrast analysis possable?
simDGE2=estimateDisp(simy2,design = simModelDesign,robust = T)
simGLM=glmFit(simDGE2,design = simModelDesign)
simmycontrast=makeContrasts("HC"=simevohot-simevocold, #don't quiet get the lab contrast, the average of the cold and hot compared to base, but could it be written as follows: Base-(Hot+Cold)/2?
                            levels = simModelDesign)
simLRT_res_HC=glmLRT(simGLM,contrast = simmycontrast[,"HC"])
simres_table_HC=simLRT_res_HC$table
simres_table_HC$padj=p.adjust(simres_table_HC$PValue,method = "BH")
#write.csv(res_table_interaction,"./DE_test_TA.csv",quote = F)
#write.csv(res_table_LB,"./DE_test_LB.csv",quote = F)

##melanogaster
melModelDesign=model.matrix(~0+melevo)
melDGE2=estimateDisp(mely2,design = melModelDesign,robust = T)
melGLM=glmFit(melDGE2,design = melModelDesign)
melmycontrast=makeContrasts("HC"=melevohot-melevocold,
                            levels = melModelDesign)
melLRT_res_HC=glmLRT(melGLM,contrast = melmycontrast[,"HC"])
melres_table_HC=melLRT_res_HC$table
melres_table_HC$padj=p.adjust(melres_table_HC$PValue,method = "BH")
#write.csv(res_table_interaction,"./DE_test_TA.csv",quote = F)
#write.csv(res_table_LB,"./DE_test_LB.csv",quote = F)

simTempNr <-  row.names(simres_table_HC[simres_table_HC$padj<0.05, ])
simTempNr001 <-  row.names(simres_table_HC[simres_table_HC$padj<0.01, ])
simTempNr01 <-  row.names(simres_table_HC[simres_table_HC$padj<0.1, ])
simTempNr15 <-  row.names(simres_table_HC[simres_table_HC$padj<0.15, ])
simTempNr20 <-  row.names(simres_table_HC[simres_table_HC$padj<0.2, ])

simTempNrUP <- row.names(simres_table_HC[simres_table_HC$padj<0.05&simres_table_HC$logFC>0, ])
simTempNrUP01 <- row.names(simres_table_HC[simres_table_HC$padj<0.1&simres_table_HC$logFC>0, ])
simTempNrUP015 <- row.names(simres_table_HC[simres_table_HC$padj<0.15&simres_table_HC$logFC>0, ])
simTempNrUP02 <- row.names(simres_table_HC[simres_table_HC$padj<0.2&simres_table_HC$logFC>0, ])

simTempNrDOWN <- row.names(simres_table_HC[simres_table_HC$padj<0.05&simres_table_HC$logFC<0, ])
simTempNrDOWN01 <- row.names(simres_table_HC[simres_table_HC$padj<0.1&simres_table_HC$logFC<0, ])
simTempNrDOWN015 <- row.names(simres_table_HC[simres_table_HC$padj<0.15&simres_table_HC$logFC<0, ])
simTempNrDOWN02 <- row.names(simres_table_HC[simres_table_HC$padj<0.2&simres_table_HC$logFC<0, ])

length(simTempNr)
length(simTempNrUP)
length(simTempNrDOWN)

melTempNr <- row.names(melres_table_HC[melres_table_HC$padj<0.05, ])
melTempNr001 <- row.names(melres_table_HC[melres_table_HC$padj<0.01, ])
melTempNr01 <- row.names(melres_table_HC[melres_table_HC$padj<0.1, ])
melTempNr015 <- row.names(melres_table_HC[melres_table_HC$padj<0.15, ])
melTempNr02 <- row.names(melres_table_HC[melres_table_HC$padj<0.2, ])

melTempNrUP <- row.names(melres_table_HC[melres_table_HC$padj<0.05&melres_table_HC$logFC>0, ])
melTempNrUP01 <- row.names(melres_table_HC[melres_table_HC$padj<0.1&melres_table_HC$logFC>0, ])
melTempNrUP015 <- row.names(melres_table_HC[melres_table_HC$padj<0.15&melres_table_HC$logFC>0, ])
melTempNrUP02 <- row.names(melres_table_HC[melres_table_HC$padj<0.2&melres_table_HC$logFC>0, ])

melTempNrDOWN <- row.names(melres_table_HC[melres_table_HC$padj<0.05&melres_table_HC$logFC<0, ])
melTempNrDOWN01 <- row.names(melres_table_HC[melres_table_HC$padj<0.1&melres_table_HC$logFC<0, ])
melTempNrDOWN015 <- row.names(melres_table_HC[melres_table_HC$padj<0.15&melres_table_HC$logFC<0, ])
melTempNrDOWN02 <- row.names(melres_table_HC[melres_table_HC$padj<0.2&melres_table_HC$logFC<0, ])

length(melTempNr)
length(melTempNrUP)
length(melTempNrDOWN)

#sim DE
write.table(simTempNr, file = "../output/DEgeneLists/simExp.txt", quote = F, sep = "\n", row.names = F, col.names = F)
write.table(simTempNrUP, file = "../output/DEgeneLists/simExpUp.txt", quote = F, sep = "\n", row.names = F, col.names = F)
write.table(simTempNrDOWN, file = "../output/DEgeneLists/simExpDown.txt", quote = F, sep = "\n", row.names = F, col.names = F)

#mel DE
write.table(melTempNr, file = "../output/DEgeneLists/melExp.txt", quote = F, sep = "\n", row.names = F, col.names = F)
write.table(melTempNrUP, file = "../output/DEgeneLists/melExpUp.txt", quote = F, sep = "\n", row.names = F, col.names = F)
write.table(melTempNrDOWN, file = "../output/DEgeneLists/melExpDown.txt", quote = F, sep = "\n", row.names = F, col.names = F)

####ClinalData####
##reading and joining the tables for 
#melanogaster
melClinal <- c("SRR1576460", "SRR1576461", "SRR1576462", "SRR1576514", "SRR1576515", "SRR1576516")
melClinalCounts <- read.csv("../data/natureZhao2015/countReads/qcsorted_mapped_trimmed_mel_SRR1576460..tsv", header = T, sep = "\t")
nrow(melClinalCounts)
for(i in melClinal[2:6]) {
  temp <- read.csv(paste("../data/natureZhao2015/countReads/qcsorted_mapped_trimmed_mel_", i, "..tsv", sep = ""), header = T, sep = "\t")
  if(setequal(row.names(temp), row.names(melClinalCounts))) {
    melClinalCounts <- cbind(melClinalCounts, temp)
  }
}

colnames(melClinalCounts) <- melClinal
head(melClinalCounts)
str(melClinalCounts)

melClinalCounts_use=melClinalCounts[apply(cpm(melClinalCounts),1,function(x) !sum(x<1)>=1),]
nrow(melClinalCounts_use)

#simulans
simClinal <- c("SRR1576520", "SRR1576522", "SRR1576523", "SRR1576527", "SRR1576528", "SRR1576529")
simClinalCounts <- read.csv("../data/natureZhao2015/countReads/qcsorted_mapped_trimmed_sim_SRR1576520..tsv", header = T, sep = "\t")
nrow(simClinalCounts)
for(i in simClinal[2:6]) {
  temp <- read.csv(paste("../data/natureZhao2015/countReads/qcsorted_mapped_trimmed_sim_", i, "..tsv", sep = ""), header = T, sep = "\t")
  if(setequal(row.names(temp), row.names(simClinalCounts))) {
    simClinalCounts <- cbind(simClinalCounts, temp)
  }
}

colnames(simClinalCounts) <- simClinal
head(simClinalCounts)
str(simClinalCounts)

simClinalCounts_use=simClinalCounts[apply(cpm(simClinalCounts),1,function(x) !sum(x<1)>=1),]
nrow(simClinalCounts_use)

#normalise the genes extracted
melSimSampleTable <- read.table("../data/natureZhao2015/filereport_read_run_PRJNA260940_tsv-3.txt", sep = "\t", header = T)
head(melSimSampleTable)
melSimSampleTable <- melSimSampleTable[,c("run_accession", "sample_alias")]
sampleInfo <- as.data.frame(t(as.data.frame(strsplit(melSimSampleTable$sample_alias, split = "-"))))
melSimSampleTable$species <- sampleInfo$V1
melSimSampleTable$population <- sampleInfo$V2
melSimSampleTable$reartemp <- sampleInfo$V3
melSimSampleTable$replicate <- sampleInfo$V4
melSimSampleTable <- melSimSampleTable[melSimSampleTable$reartemp == "29C", ]
head(melSimSampleTable)

table(melSimSampleTable[melSimSampleTable$species == "Dsim", ]$run_accession == colnames(simClinalCounts_use)) #good?
simClinalSamples <- melSimSampleTable[melSimSampleTable$species == "Dsim", ]$population
clinalsimy2=DGEList(counts=simClinalCounts_use,group = simClinalSamples)
clinalsimy2=calcNormFactors(clinalsimy2)
clinalsimModelDesign=model.matrix(~0+simClinalSamples)
clinalsimDGE2=estimateDisp(clinalsimy2,design = clinalsimModelDesign,robust = T)
clinalsimGLM=glmFit(clinalsimDGE2,design = clinalsimModelDesign)
clinalsimmycontrast=makeContrasts("PM"=simClinalSamplesPanama-simClinalSamplesMaine,
                                  levels = clinalsimModelDesign)
clinalsimLRT_res_PM=glmLRT(clinalsimGLM,contrast = clinalsimmycontrast[,"PM"])
clinalsimres_table_PM=clinalsimLRT_res_PM$table
clinalsimres_table_PM$padj=p.adjust(clinalsimres_table_PM$PValue,method = "BH")
clinalsimres_table_PM <- clinalsimres_table_PM[order(row.names(clinalsimres_table_PM)), ]
head(clinalsimres_table_PM)

clinalSimSigGenes <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.05, ])
clinalSimSigGenes001 <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.01, ])
clinalSimSigGenes01 <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.1, ])
clinalSimSigGenes015 <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.15, ])
clinalSimSigGenes02 <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.2, ])

clinalSimSigGenesPanamaUp <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.05 & clinalsimres_table_PM$logFC>0, ])
clinalSimSigGenesPanamaUp01 <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.1 & clinalsimres_table_PM$logFC>0, ])
clinalSimSigGenesPanamaUp015 <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.15 & clinalsimres_table_PM$logFC>0, ])
clinalSimSigGenesPanamaUp02 <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.2 & clinalsimres_table_PM$logFC>0, ])

clinalSimSigGenesPanamaDown <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.05 & clinalsimres_table_PM$logFC<0, ])
clinalSimSigGenesPanamaDown01 <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.1 & clinalsimres_table_PM$logFC<0, ])
clinalSimSigGenesPanamaDown015 <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.15 & clinalsimres_table_PM$logFC<0, ])
clinalSimSigGenesPanamaDown02 <- row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj<0.2 & clinalsimres_table_PM$logFC<0, ])

length(clinalSimSigGenes)
length(clinalSimSigGenesPanamaUp)
length(clinalSimSigGenesPanamaDown)

table(melSimSampleTable[melSimSampleTable$species == "Dmel", ]$run_accession == colnames(melClinalCounts_use)) #good?
melClinalSamples <- melSimSampleTable[melSimSampleTable$species == "Dmel", ]$population
clinalmely2=DGEList(counts=melClinalCounts_use,group = melClinalSamples)
clinalmely2=calcNormFactors(clinalmely2)
clinalmelModelDesign=model.matrix(~0+melClinalSamples)
clinalmelDGE2=estimateDisp(clinalmely2,design = clinalmelModelDesign,robust = T)
clinalmelGLM=glmFit(clinalmelDGE2,design = clinalmelModelDesign)
clinalmelmycontrast=makeContrasts("PM"=melClinalSamplesPanama-melClinalSamplesMaine,
                                  levels = clinalmelModelDesign)
clinalmelLRT_res_PM=glmLRT(clinalmelGLM,contrast = clinalmelmycontrast[,"PM"])
clinalmelres_table_PM=clinalmelLRT_res_PM$table
clinalmelres_table_PM$padj=p.adjust(clinalmelres_table_PM$PValue,method = "BH")
clinalmelres_table_PM <- clinalmelres_table_PM[order(row.names(clinalmelres_table_PM)), ]
head(clinalmelres_table_PM)

clinalMelSigGenes <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.05, ])
clinalMelSigGenes001 <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.01, ])
clinalMelSigGenes01 <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.1, ])
clinalMelSigGenes015 <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.15, ])
clinalMelSigGenes02 <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.2, ])

clinalMelSigGenesPanamaUp <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.05 & clinalmelres_table_PM$logFC>0, ])
clinalMelSigGenesPanamaUp01 <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.1 & clinalmelres_table_PM$logFC>0, ])
clinalMelSigGenesPanamaUp015 <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.15 & clinalmelres_table_PM$logFC>0, ])
clinalMelSigGenesPanamaUp02 <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.2 & clinalmelres_table_PM$logFC>0, ])

clinalMelSigGenesPanamaDown <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.05 & clinalmelres_table_PM$logFC<0, ])
clinalMelSigGenesPanamaDown01 <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.1 & clinalmelres_table_PM$logFC<0, ])
clinalMelSigGenesPanamaDown015 <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.15 & clinalmelres_table_PM$logFC<0, ])
clinalMelSigGenesPanamaDown02 <- row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj<0.2 & clinalmelres_table_PM$logFC<0, ])

length(clinalMelSigGenes)
length(clinalMelSigGenesPanamaUp)
length(clinalMelSigGenesPanamaDown)

nrow(clinalsimy2)
nrow(clinalmely2)

####GeneOntology####
##enrichment tests
#make a list of the gene sets I want to make an enrichment with, for the experimental part

#simulans exp.
DEGeneListExperiSim <-
  list(
    simHC = row.names(simy2)[simres_table_HC$padj < 0.05],
    simHCUp = row.names(simy2)[simres_table_HC$padj <
                                 0.05 & simres_table_HC$logFC > 0],
    simHCDown = row.names(simy2)[simres_table_HC$padj <
                                   0.05 & simres_table_HC$logFC < 0]
  )
background=rownames(simy2)
ExperiGOResSim <- list()
for (i in 1:length(DEGeneListExperiSim)){
  tmp5=factor(as.integer(rownames(simy2)%in%DEGeneListExperiSim[[i]]))
  names(tmp5)=rownames(simy2)#genelist
  tgd5=new( "topGOdata", ontology="BP", allGenes = tmp5, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")#data preparation#
  genesInGO = genesInTerm(tgd5)
  genesInGOPaste = c()
  for(j in 1:length(genesInGO)){
    genesInGOPaste = c(genesInGOPaste, paste(genesInGO[[j]], collapse = ","))
  }
  resTopGO.classic=runTest(tgd5, algorithm = "classic", statistic = "Fisher")#enrichment test#
  resTopGO.weight01=runTest(tgd5, algorithm = "weight01", statistic = "Fisher")
  tmp_res5=GenTable(tgd5,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)#analysis of results#
  tmp_res5=tmp_res5[order(tmp_res5$GO.ID),]
  if((sum((tmp_res5$GO.ID == names(genesInGO)) == TRUE)) == length(names(genesInGO))) {
    tmp_res5$Genes=genesInGOPaste
  }
  ExperiGOResSim[[i]]=tmp_res5
  ExperiGOResSim[[i]]$Fisher.classic <- as.numeric(ExperiGOResSim[[i]]$Fisher.classic)
  ExperiGOResSim[[i]]$Fisher.weight01 <- as.numeric(ExperiGOResSim[[i]]$Fisher.weight01)
}
names(ExperiGOResSim) <- names(DEGeneListExperiSim)

#melanogaster exp.
DEGeneListExperiMel <-
  list(
    melHC = row.names(mely2)[melres_table_HC$padj < 0.05],
    melHCUp = row.names(mely2)[melres_table_HC$padj <
                                 0.05 & melres_table_HC$logFC > 0],
    melHCDown = row.names(mely2)[melres_table_HC$padj <
                                   0.05 & melres_table_HC$logFC < 0]
  )
background=rownames(mely2)
ExperiGOResMel <- list()
for (i in 1:length(DEGeneListExperiMel)){
  tmp5=factor(as.integer(rownames(mely2)%in%DEGeneListExperiMel[[i]]))
  names(tmp5)=rownames(mely2)#genelist
  tgd5=new( "topGOdata", ontology="BP", allGenes = tmp5, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")#data preparation#
  genesInGO = genesInTerm(tgd5)
  genesInGOPaste = c()
  for(j in 1:length(genesInGO)){
    genesInGOPaste = c(genesInGOPaste, paste(genesInGO[[j]], collapse = ","))
  }
  resTopGO.classic=runTest(tgd5, algorithm = "classic", statistic = "Fisher")#enrichment test#
  resTopGO.weight01=runTest(tgd5, algorithm = "weight01", statistic = "Fisher")
  tmp_res5=GenTable(tgd5,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)#analysis of results#
  tmp_res5=tmp_res5[order(tmp_res5$GO.ID),]
  if((sum((tmp_res5$GO.ID == names(genesInGO)) == TRUE)) == length(names(genesInGO))) {
    tmp_res5$Genes=genesInGOPaste
  }
  ExperiGOResMel[[i]]=tmp_res5
  ExperiGOResMel[[i]]$Fisher.classic <- as.numeric(ExperiGOResMel[[i]]$Fisher.classic)
  ExperiGOResMel[[i]]$Fisher.weight01 <- as.numeric(ExperiGOResMel[[i]]$Fisher.weight01)
}
names(ExperiGOResMel) <- names(DEGeneListExperiMel)

#simulans cli.
DEGeneListClineSim <-
  list(
    simPM = row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj < 0.05,]),
    simPMUp = row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj <
                                                0.05 & clinalsimres_table_PM$logFC > 0,]),
    simPMDown = row.names(clinalsimres_table_PM[clinalsimres_table_PM$padj <
                                                  0.05 & clinalsimres_table_PM$logFC < 0, ])
  )
background=rownames(clinalsimres_table_PM) #the background
ClineGOResSim=list()
for (i in 1:length(DEGeneListClineSim)){
  tmp5=factor(as.integer(rownames(clinalsimres_table_PM)%in%DEGeneListClineSim[[i]]))
  names(tmp5)=rownames(clinalsimres_table_PM)#genelist
  tgd5=new("topGOdata", ontology="BP", allGenes = tmp5, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")#data preparation#
  genesInGO = genesInTerm(tgd5)
  genesInGOPaste = c()
  for(j in 1:length(genesInGO)){
    genesInGOPaste = c(genesInGOPaste, paste(genesInGO[[j]], collapse = ","))
  }
  resTopGO.classic=runTest(tgd5, algorithm = "classic", statistic = "Fisher")#enrichment test#
  resTopGO.weight01=runTest(tgd5, algorithm = "weight01", statistic = "Fisher")
  tmp_res5=GenTable(tgd5,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)#analysis of results#
  tmp_res5=tmp_res5[order(tmp_res5$GO.ID),]
  if((sum((tmp_res5$GO.ID == names(genesInGO)) == TRUE)) == length(names(genesInGO))) {
    tmp_res5$Genes=genesInGOPaste
  }
  ClineGOResSim[[i]]=tmp_res5
  ClineGOResSim[[i]]$Fisher.classic <- as.numeric(ClineGOResSim[[i]]$Fisher.classic)
  ClineGOResSim[[i]]$Fisher.weight01 <- as.numeric(ClineGOResSim[[i]]$Fisher.weight01)
}
names(ClineGOResSim) <- names(DEGeneListClineSim)

#melanogaster cli.
DEGeneListClineMel <-
  list(
    melPM = row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj < 0.05,]),
    melPMUp = row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj <
                                                0.05 & clinalmelres_table_PM$logFC > 0, ]),
    melPMDown = row.names(clinalmelres_table_PM[clinalmelres_table_PM$padj <
                                                  0.05 & clinalmelres_table_PM$logFC < 0, ])
  )
background=rownames(clinalmelres_table_PM) #the background
ClineGOResMel=list()
for (i in 1:length(DEGeneListClineMel)){
  tmp5=factor(as.integer(rownames(clinalmelres_table_PM)%in%DEGeneListClineMel[[i]]))
  names(tmp5)=rownames(clinalmelres_table_PM)#genelist
  tgd5=new("topGOdata", ontology="BP", allGenes = tmp5, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")#data preparation#
  genesInGO = genesInTerm(tgd5)
  genesInGOPaste = c()
  for(j in 1:length(genesInGO)){
    genesInGOPaste = c(genesInGOPaste, paste(genesInGO[[j]], collapse = ","))
  }
  resTopGO.classic=runTest(tgd5, algorithm = "classic", statistic = "Fisher")#enrichment test#
  resTopGO.weight01=runTest(tgd5, algorithm = "weight01", statistic = "Fisher")
  tmp_res5=GenTable(tgd5,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)#analysis of results#
  tmp_res5=tmp_res5[order(tmp_res5$GO.ID),]
  if((sum((tmp_res5$GO.ID == names(genesInGO)) == TRUE)) == length(names(genesInGO))) {
    tmp_res5$Genes=genesInGOPaste
  }
  ClineGOResMel[[i]]=tmp_res5
  ClineGOResMel[[i]]$Fisher.classic <- as.numeric(ClineGOResMel[[i]]$Fisher.classic)
  ClineGOResMel[[i]]$Fisher.weight01 <- as.numeric(ClineGOResMel[[i]]$Fisher.weight01)
}
names(ClineGOResMel) <- names(DEGeneListClineMel)

#what are the significantly enriched terms
ExperiGOResSim$simHC$Fisher.weight01 <- as.numeric(ExperiGOResSim$simHC$Fisher.weight01)
ExperiGOResSim$simHC[ExperiGOResSim$simHC$Fisher.weight01 < 0.05, ]$Term
ExperiGOResSim$simHCUp[ExperiGOResSim$simHCUp$Fisher.weight01 < 0.05, ]$Term
ExperiGOResSim$simHCDown[ExperiGOResSim$simHCDown$Fisher.weight01 < 0.05, ]$Term

ExperiGOResSim$simHC <- ExperiGOResSim$simHC[order(ExperiGOResSim$simHC$Fisher.weight01), ]
write.table(ExperiGOResSim$simHC[,1:8], quote = F, file = "../output/GOtables/ExperiSim.csv", sep = ";", row.names = F, dec = ",")

ExperiGOResMel$melHC$Fisher.weight01 <- as.numeric(ExperiGOResMel$melHC$Fisher.weight01)
ExperiGOResMel$melHC[ExperiGOResMel$melHC$Fisher.weight01 < 0.05, ]$Term
ExperiGOResMel$melHCUp[ExperiGOResMel$melHCUp$Fisher.weight01 < 0.05, ]$Term
ExperiGOResMel$melHCDown[ExperiGOResMel$melHCDown$Fisher.weight01 < 0.05, ]$Term

ExperiGOResMel$melHC <- ExperiGOResMel$melHC[order(ExperiGOResMel$melHC$Fisher.weight01), ]
write.table(ExperiGOResMel$melHC[,1:8], quote = F, file = "../output/GOtables/ExperiMel.csv", sep = ";", row.names = F, dec = ",")

ClineGOResSim$simPM$Fisher.weight01 <- as.numeric(ClineGOResSim$simPM$Fisher.weight01)
ClineGOResSim$simPM[ClineGOResSim$simPM$Fisher.weight01 < 0.05, ]$Term
ClineGOResSim$simPMUp[ClineGOResSim$simPMUp$Fisher.weight01 < 0.05, ]$Term
ClineGOResSim$simPMDown[ClineGOResSim$simPMDown$Fisher.weight01 < 0.05, ]$Term

ClineGOResSim$simPM <- ClineGOResSim$simPM[order(ClineGOResSim$simPM$Fisher.weight01), ]
write.table(ClineGOResSim$simPM[,1:8], quote = F, file = "../output/GOtables/ClineSim.csv", sep = ";", row.names = F, dec = ",")

ClineGOResMel$melPM$Fisher.weight01 <- as.numeric(ClineGOResMel$melPM$Fisher.weight01)
ClineGOResMel$melPM[ClineGOResMel$melPM$Fisher.weight01 < 0.05, ]$Term
ClineGOResMel$melPMUp[ClineGOResMel$melPMUp$Fisher.weight01 < 0.05, ]$Term
ClineGOResMel$melPMDown[ClineGOResMel$melPMDown$Fisher.weight01 < 0.05, ]$Term

ClineGOResMel$melPM <- ClineGOResMel$melPM[order(ClineGOResMel$melPM$Fisher.weight01), ]
write.table(ClineGOResMel$melPM[,1:8], quote = F, file = "../output/GOtables/ClineMel.csv", sep = ";", row.names = F, dec = ",")


#write some tables


####TissueSpecificity####
tau_bg = read.table("../data/pleiotropyMeasures/tau_flyatlas2_bg.txt",sep = "\t",header = T,stringsAsFactors = F) # tau estimated from flyatlas2 (Scott's analysis 2018)
#png(filename = "tau_tissue_boxplot.png",height = 12,width = 22,units = "cm",res = 300,pointsize = 10)
par(mar=c(5,5,1,2))
boxplot(tau_bg$tau_male_4~tau_bg$max_tissue,las=2,
        ylab= "tissue specificity(tau)", axes=F,
        xlab=" ",
        names=rep('',6),
        las=2 )
axis(2,at=seq(0,1,0.2), labels=c('0','0.2','0.4','0.6','0.8',"1.0"),las=1)
labels=c("Accessory glands","Brain","Carcass","Crop","Eye",
         "Head","Hindgut","Malpighian Tubules","Midgut","Rectal pad",
         "Salivary gland","Testis","Thoracicoabdominal ganglion")
text(x=c(1:13)+.1, y=0.09, cex=.7, srt = 30, adj = 1.1,
     labels = labels, xpd = TRUE)
box()
#dev.off()

####NetworkConnectivity####
net.con=read.table("../data/pleiotropyMeasures/flynet_supervised_0.6.txt",header=F,stringsAsFactors = F) #reference Marbach et al., 2012 genome research
#out-degree
ind=unique(net.con$V1)
out_degree=as.data.frame(matrix(NA,length(ind),2))
out_degree[,1]=as.character(unique(net.con$V1))
for (i in 1:length(ind)) {
  out_degree[i,2]=sum(net.con$V1==ind[i])
}
colnames(out_degree)=c("FBgn","out_de")
hist(as.numeric(paste(out_degree[,2])))

#in-degree
ind=unique(net.con$V2)
in_degree=as.data.frame(matrix(NA,length(ind),2))
in_degree[,1]=as.character(unique(net.con$V2))
for (i in 1:length(ind)) {
  in_degree[i,2]=sum(net.con$V2==ind[i])
}
colnames(in_degree)=c("FBgn","in_de")
hist(as.numeric(paste(in_degree[,2])))

conn_use=merge(in_degree,out_degree,by = "FBgn",all = T)
conn_use[which(is.na(conn_use$out_de)),3]=0
conn_use[which(is.na(conn_use$in_de)),2]=0 #there are two values NA after merging - some genes which only have an "in" connection.
conn_use$conn=conn_use$in_de+conn_use$out_de

PleioIndexTau <- tau_bg[,c(1,3)]
PleioIndexTau$tau_male_4 <- 1-PleioIndexTau$tau_male_4
colnames(PleioIndexTau) <-  c("FBgn","tau")
row.names(PleioIndexTau) <- PleioIndexTau$FBgn
PleioIndexTau <- PleioIndexTau[order(PleioIndexTau$FBgn), ]

PleioIndexConn <- conn_use[,c(1,4)]
colnames(PleioIndexConn) <-  c("FBgn","connectivity")
row.names(PleioIndexConn) <- PleioIndexConn$FBgn
PleioIndexConn <- PleioIndexConn[order(PleioIndexConn$FBgn), ]

PleioIndex=merge(conn_use[,c(1,4)],tau_bg[,c(1,3)],by="FBgn")
PleioIndex$tau_male_4=1-PleioIndex$tau_male_4
colnames(PleioIndex) <-  c("FBgn","connectivity", "tau")
row.names(PleioIndex) <- PleioIndex$FBgn
PleioIndex <- PleioIndex[order(PleioIndex$FBgn), ]

#are the estimates correlated? 
cor.test(PleioIndex$connectivity, PleioIndex$tau, method = "spearman")

####PleiotropicCost####

tau <- ggplot(data = PleioIndexTau, aes(x=tau)) +
  geom_histogram(bins = 20, color=viridis(n=10)[3], fill="white", size = 1.5) +
  theme_classic() +
  xlab("1-tau") +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20)) 
tau

conn <- ggplot(data = PleioIndexConn, aes(x=log10(connectivity))) +
  geom_histogram(bins = 20, color=viridis(n=10)[6], fill="white", size = 1.5) +
  theme_classic() +
  xlab("log10(connectivity)") +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20)) 
conn

tiff("../output/distOfTauAndConnectivity.tiff", height = 30, width = 30, res = 300, units = "cm")
multiplot(tau, conn, cols = 1)
dev.off()



####does pleiotropy corrilate with logFC

#does pleiotropy correlate with logFC
####Tau####
pleioIndexs <- PleioIndexTau

simExDataIndexes <- pleioIndexs[row.names(pleioIndexs) %in% row.names(simres_table_HC), ]
simExDataPleio <-  simres_table_HC[row.names(simres_table_HC) %in% simExDataIndexes$FBgn, ]
table(row.names(simExDataPleio) == row.names(simExDataIndexes))
simExDataPleio$tau <- simExDataIndexes$tau

melExDataIndexes <- pleioIndexs[row.names(pleioIndexs) %in% row.names(melres_table_HC), ]
melExDataPleio <-  melres_table_HC[row.names(melres_table_HC) %in% melExDataIndexes$FBgn, ]
table(row.names(melExDataPleio) == row.names(melExDataIndexes))
melExDataPleio$tau <- melExDataIndexes$tau

clinalsimres_table_PM <- clinalsimres_table_PM[order(row.names(clinalsimres_table_PM)), ]
simClDataIndexes <- pleioIndexs[row.names(pleioIndexs) %in% row.names(clinalsimres_table_PM), ]
simClDataPleio <-  clinalsimres_table_PM[row.names(clinalsimres_table_PM) %in% simClDataIndexes$FBgn, ]
table(row.names(simClDataPleio) == row.names(simClDataIndexes))
simClDataPleio$tau <- simClDataIndexes$tau

clinalmelres_table_PM <- clinalmelres_table_PM[order(row.names(clinalmelres_table_PM)), ]
melClDataIndexes <- pleioIndexs[row.names(pleioIndexs) %in% row.names(clinalmelres_table_PM), ]
melClDataPleio <-  clinalmelres_table_PM[row.names(clinalmelres_table_PM) %in% melClDataIndexes$FBgn, ]
table(row.names(melClDataPleio) == row.names(melClDataIndexes))
melClDataPleio$tau <- melClDataIndexes$tau

simClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalSimSigGenes, ] #clinalSimSigGenes
nrow(simClSigGenesPleio)
#simClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalSimSigGenes01, ] #clinalSimSigGenes
#simClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalSimSigGenes015, ] #clinalSimSigGenes
#simClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalSimSigGenes02, ] #clinalSimSigGenes
simClSigGenesPleio$group <- rep("SimClinal", nrow(simClSigGenesPleio))
simClSigGenesPleio$species <- rep("Sim", nrow(simClSigGenesPleio))


melClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalMelSigGenes, ] #clinalSimSigGenes
nrow(melClSigGenesPleio)
#melClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalMelSigGenes01, ] #clinalSimSigGenes
#melClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalMelSigGenes015, ] #clinalSimSigGenes
#melClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalMelSigGenes02, ] #clinalSimSigGenes
melClSigGenesPleio$group <- rep("MelClinal", nrow(melClSigGenesPleio))
melClSigGenesPleio$species <- rep("Mel", nrow(melClSigGenesPleio))


simExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% simTempNr, ] #clinalSimSigGenes
nrow(simExSigGenesPleio)
#simExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% simTempNr01, ] #clinalSimSigGenes
#simExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% simTempNr15, ] #clinalSimSigGenes
#simExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% simTempNr20, ] #clinalSimSigGenes
simExSigGenesPleio$group <- rep("SimExperi", nrow(simExSigGenesPleio))
simExSigGenesPleio$species <- rep("Sim", nrow(simExSigGenesPleio))


melExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% melTempNr, ] #clinalSimSigGenes
nrow(melExSigGenesPleio)
#melExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% melTempNr01, ] #clinalSimSigGenes
#melExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% melTempNr015, ] #clinalSimSigGenes
#melExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% melTempNr02, ] #clinalSimSigGenes
melExSigGenesPleio$group <- rep("MelExperi", nrow(melExSigGenesPleio))
melExSigGenesPleio$species <- rep("Mel", nrow(melExSigGenesPleio))

pleiotropyTauWilcox <- c()
pleiotropyTauWilcox <- c(pleiotropyTauWilcox, wilcox.test(x = simClSigGenesPleio$tau, y=simExSigGenesPleio$tau,  alternative = c("two.sided"))$p.value)
pleiotropyTauWilcox<- c(pleiotropyTauWilcox, wilcox.test(x = melClSigGenesPleio$tau, y=melExSigGenesPleio$tau,  alternative = c("two.sided"))$p.value)
names(pleiotropyTauWilcox) <- c("sim", "mel")
pleiotropyTauWilcox

pleioDataSig <- rbind(simClSigGenesPleio, melClSigGenesPleio, simExSigGenesPleio, melExSigGenesPleio)

pleioDataSig$group <- as.factor(pleioDataSig$group)


tiff("../output/PleioInClineExperiBlueBrightGreen.tiff", height = 15, width = 15, res = 300, units = "cm")
Figure2 <- ggplot(data = pleioDataSig, aes(x = group, y = tau, colour = group, group = group)) +
  geom_boxplot(size = 1.5) + 
  geom_jitter(width = 0.3,  shape=1, size = 0.5) +
  theme_classic() +
  scale_color_manual(values=c(
    viridis(n=40, option = "turbo")[3],
    viridis(n=40, option = "turbo")[37],
    viridis(n=40,option = "turbo")[8],
    viridis(n=40, option = "turbo")[30])) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),  
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.position='none') +
  geom_signif(comparisons = list(c("MelClinal", "MelExperi")),
              annotations = "3.4e-05",  color = "black") +
  geom_signif(comparisons = list(c("SimClinal", "SimExperi")),
              annotations = "6.4e-04",color = "black")+
  xlab("") +
  ylab("Pleiotropy (1-tau)") 
dev.off()




pleiotropyTauWilcox

##Corrilation
simExCorrFCTau <- cor.test(simExDataPleio$tau, abs(simExDataPleio$logFC), method = "spearman", exact = F)
melExCorrFCTau <- cor.test(melExDataPleio$tau, abs(melExDataPleio$logFC), method = "spearman",  exact = F)
simCLCorrFCTau <- cor.test(simClDataPleio$tau, abs(simClDataPleio$logFC), method = "spearman",  exact = F)
melCLCorrFCTau <- cor.test(melClDataPleio$tau, abs(melClDataPleio$logFC), method = "spearman",  exact = F)

tauLogFCTauadj <- data.frame(
  rho = c(
    simExCorrFCTau$estimate,
    melExCorrFCTau$estimate,
    simCLCorrFCTau$estimate,
    melCLCorrFCTau$estimate
  ),
  pval = c(
    simExCorrFCTau$p.value,
    melExCorrFCTau$p.value,
    simCLCorrFCTau$p.value,
    melCLCorrFCTau$p.value
  ),
  padj = p.adjust(
    c(
      simExCorrFCTau$p.value,
      melExCorrFCTau$p.value,
      simCLCorrFCTau$p.value,
      melCLCorrFCTau$p.value
    ),
    method = "bonferroni"
  ),
  group = c(
    "simExCorrFCTau",
    "melExCorrFCTau",
    "simCLCorrFCTau",
    "melCLCorrFCTau"
  )
)

tauLogFCTauadj


###plotting
simClDataPleioBins <- simClDataPleio %>% mutate(Tau_bin = ntile(tau, n=10)) 
simClDataPleioBins$Tau_bin <- as.factor(simClDataPleioBins$Tau_bin)

melClDataPleioBins <- melClDataPleio %>% mutate(Tau_bin = ntile(tau, n=10)) 
melClDataPleioBins$Tau_bin <- as.factor(melClDataPleioBins$Tau_bin)

simExDataPleioBins <- simExDataPleio %>% mutate(Tau_bin = ntile(tau, n=10)) 
simExDataPleioBins$Tau_bin <- as.factor(simExDataPleioBins$Tau_bin)

melExDataPleioBins <- melExDataPleio %>% mutate(Tau_bin = ntile(tau, n=10)) 
melExDataPleioBins$Tau_bin <- as.factor(melExDataPleioBins$Tau_bin)


###
simClDataPleioBins <- simClDataPleio %>% mutate(Tau_bin = cut(tau, breaks=c(0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0, 1.1))) #breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)))
simClDataPleioBins$Tau_bin <- as.factor(simClDataPleioBins$Tau_bin)

melClDataPleioBins <- melClDataPleio %>% mutate(Tau_bin = cut(tau, breaks=c(0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))) #breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)))
melClDataPleioBins$Tau_bin <- as.factor(melClDataPleioBins$Tau_bin)

simExDataPleioBins <- simExDataPleio %>% mutate(Tau_bin = cut(tau, breaks=c(0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))) #breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)))
simExDataPleioBins$Tau_bin <- as.factor(simExDataPleioBins$Tau_bin)

melExDataPleioBins <- melExDataPleio %>% mutate(Tau_bin = cut(tau, breaks=c(0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))) #breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)))
melExDataPleioBins$Tau_bin <- as.factor(melExDataPleioBins$Tau_bin)

###
simClDataPleioBins <- simClDataPleio %>% mutate(Tau_bin = cut(tau, breaks=10)) 
simClDataPleioBins$Tau_bin <- as.factor(simClDataPleioBins$Tau_bin)

melClDataPleioBins <- melClDataPleio %>% mutate(Tau_bin = cut(tau, breaks=10))
melClDataPleioBins$Tau_bin <- as.factor(melClDataPleioBins$Tau_bin)

simExDataPleioBins <- simExDataPleio %>% mutate(Tau_bin = cut(tau, breaks=10)) 
simExDataPleioBins$Tau_bin <- as.factor(simExDataPleioBins$Tau_bin)

melExDataPleioBins <- melExDataPleio %>% mutate(Tau_bin = cut(tau, breaks=10))
melExDataPleioBins$Tau_bin <- as.factor(melExDataPleioBins$Tau_bin)

###
simClDataPleioBins <- simClDataPleio %>%  mutate(Tau_bin=cut_width(tau, width=0.1, boundary=0))
simClDataPleioBins$Tau_bin <- as.factor(simClDataPleioBins$Tau_bin)

melClDataPleioBins <- melClDataPleio %>% mutate(Tau_bin=cut_width(tau, width=0.1, boundary=0))
melClDataPleioBins$Tau_bin <- as.factor(melClDataPleioBins$Tau_bin)

simExDataPleioBins <- simExDataPleio %>% mutate(Tau_bin=cut_width(tau, width=0.1, boundary=0))
simExDataPleioBins$Tau_bin <- as.factor(simExDataPleioBins$Tau_bin)

melExDataPleioBins <- melExDataPleio %>% mutate(Tau_bin=cut_width(tau, width=0.1, boundary=0))
melExDataPleioBins$Tau_bin <- as.factor(melExDataPleioBins$Tau_bin)


ai <- ggplot(simClDataPleioBins, aes(y= abs(logFC), x = Tau_bin, group = Tau_bin)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[8], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("simulans cline", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "simCLCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("A",subtitle = "Nature (D. simulans)") +
  theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (1-tau)")
ai

bi <- ggplot(melClDataPleioBins, aes(y= abs(logFC), x = Tau_bin)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[3], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("melanogaster cline", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "melCLCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("B",subtitle = "Nature (D. melanogaster)") +
  theme(axis.text.x = element_text(size=15, angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (1-tau)")
bi

ci <- ggplot(simExDataPleioBins, aes(y= abs(logFC), x = Tau_bin)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[30], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("simulans exp.", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "simExCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("C", subtitle = "Laboratory (D. simulans)") +
  theme(axis.text.x = element_text(size=15, angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (1-tau)")
ci

di <- ggplot(melExDataPleioBins, aes(y= abs(logFC), x = Tau_bin)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[37], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("melanogaster exp.", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "melExCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("D", subtitle = "Laboratory (D. melanogaster)") +
  theme(axis.text.x = element_text(size=15, angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (1-tau)")
di

tiff("../output/PleioByLogFCBinnedByBreakes.tiff", height = 30, width = 30, res = 300, units = "cm")
multiplot(ai, ci, bi, di, cols = 2)
dev.off()

tiff("../output/PleioByLogFCBinned.tiff", height = 30, width = 30, res = 300, units = "cm")
multiplot(ai, ci, bi, di, cols = 2)
dev.off()


#creating a variable to plot the mean/median of the bin for the boxplot
simClDataPleioBins <- simClDataPleio %>% mutate(Tau_bin = ntile(tau, n=10)) 
simClDataPleioBins$Tau_bin <- as.factor(simClDataPleioBins$Tau_bin)
repTimes <- c(table(simClDataPleioBins$Tau_bin))
str(repTimes)

tauMean <- as.data.frame(simClDataPleioBins %>%
                           group_by(Tau_bin) %>%
                           summarize(MeanTau = mean(tau, digits = 8)))$MeanTau
tauMean <- round(tauMean, digits = 4)
tauMedian <- rep(tauMean, times = repTimes)
simClDataPleioBins <- simClDataPleioBins[order(simClDataPleioBins$Tau_bin), ]
simClDataPleioBins$meanTau <- as.factor(tauMedian)
str(simClDataPleioBins)
ai <- ggplot(simClDataPleioBins, aes(y= abs(logFC), x = as.factor(meanTau), group = as.factor(meanTau))) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[8], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("simulans cline", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "simCLCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("A",subtitle = "Nature (D. simulans)") +
  theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (1-tau)")
ai


melClDataPleioBins <- melClDataPleio %>% mutate(Tau_bin = ntile(tau, n=10)) 
melClDataPleioBins$Tau_bin <- as.factor(melClDataPleioBins$Tau_bin)
repTimes <- c(table(melClDataPleioBins$Tau_bin))
str(repTimes)
tauMean <- as.data.frame(melClDataPleioBins %>%
                           group_by(Tau_bin) %>%
                           summarize(MeanTau = mean(tau, digits = 8)))$MeanTau
tauMean <- round(tauMean, digits = 4)
tauMedian <- rep(tauMean, times = repTimes)
melClDataPleioBins <- melClDataPleioBins[order(melClDataPleioBins$Tau_bin), ]
melClDataPleioBins$meanTau <- as.factor(tauMedian)
str(melClDataPleioBins)
bi <- ggplot(melClDataPleioBins, aes(y= abs(logFC), x = as.factor(meanTau), group = meanTau)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[3], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("melanogaster cline", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "melCLCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("B",subtitle = "Nature (D. melanogaster)") +
  theme(axis.text.x = element_text(size=15, angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (1-tau)")
bi

simExDataPleioBins <- simExDataPleio %>% mutate(Tau_bin = ntile(tau, n=10)) 
simExDataPleioBins$Tau_bin <- as.factor(simExDataPleioBins$Tau_bin)
repTimes <- c(table(simExDataPleioBins$Tau_bin))
str(repTimes)
tauMean <- as.data.frame(simExDataPleioBins %>%
                           group_by(Tau_bin) %>%
                           summarize(MeanTau = mean(tau, digits = 8)))$MeanTau
tauMean <- round(tauMean, digits = 4)
tauMedian <- rep(tauMean, times = repTimes)
simExDataPleioBins <- simExDataPleioBins[order(simExDataPleioBins$Tau_bin), ]
simExDataPleioBins$meanTau <- as.factor(tauMedian)
str(simExDataPleioBins)
ci <- ggplot(simExDataPleioBins, aes(y= abs(logFC), x = as.factor(meanTau), group = meanTau)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[30], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("simulans exp.", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "simExCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("C", subtitle = "Laboratory (D. simulans)") +
  theme(axis.text.x = element_text(size=15, angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (1-tau)")
ci


melExDataPleioBins <- melExDataPleio %>% mutate(Tau_bin = ntile(tau, n=10)) 
melExDataPleioBins$Tau_bin <- as.factor(melExDataPleioBins$Tau_bin)
repTimes <- c(table(melExDataPleioBins$Tau_bin))
str(repTimes)
tauMean <- as.data.frame(melExDataPleioBins %>%
                           group_by(Tau_bin) %>%
                           summarize(MeanTau = mean(tau, digits = 8)))$MeanTau
tauMean <- round(tauMean, digits = 4)
tauMedian <- rep(tauMean, times = repTimes)
melExDataPleioBins <- melExDataPleioBins[order(melExDataPleioBins$Tau_bin), ]
melExDataPleioBins$meanTau <- as.factor(tauMedian)
str(melExDataPleioBins)
di <- ggplot(melExDataPleioBins, aes(y= abs(logFC), x = as.factor(meanTau), group = meanTau)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[37], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("melanogaster exp.", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "melExCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("D", subtitle = "Laboratory (D. melanogaster)") +
  theme(axis.text.x = element_text(size=15, angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (1-tau)")
di



tiff("../output/PleioByLogFCBinnedByBinMeanForEachBin.tiff", height = 30, width = 30, res = 300, units = "cm")
multiplot(ai, ci, bi, di, cols = 2)
dev.off()



#density distribution of the DE genes relative to background of pleiotropy
pleioIndexs$group <- rep("Background", nrow(pleioIndexs))
masterPleio <- rbind(pleioDataSig [, 1:3], pleioIndexs)
masterPleio$group <- factor(masterPleio$group, 
                            levels = c("MelClinal", "MelExperi", "SimClinal", "SimExperi", "Background"))

#try to compare the distribution with densityplots
i <- ggplot(data = masterPleio, aes(x = tau, colour = group, group = group, linetype=group)) +
  geom_density(size = 1.5) + 
  theme_classic() +
  scale_linetype_manual(values=c("solid","solid", "solid","solid", "dashed"))+
  scale_color_manual(values=c(
    viridis(n=40, option = "turbo")[3],
    viridis(n=40, option = "turbo")[37],
    viridis(n=40,option = "turbo")[8],
    viridis(n=40, option = "turbo")[30], "black")) +
  ggtitle("A") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=15),
        title = element_text(size=20)) +
  xlab("Pleiotropy (1-tau)")
i


####Connectivity####
pleioIndexs <- PleioIndexConn

simExDataIndexes <- pleioIndexs[row.names(pleioIndexs) %in% row.names(simres_table_HC), ]
simExDataPleio <-  simres_table_HC[row.names(simres_table_HC) %in% simExDataIndexes$FBgn, ]
table(row.names(simExDataPleio) == row.names(simExDataIndexes))
simExDataPleio$tau <- simExDataIndexes$tau
simExDataPleio$connectivity <- simExDataIndexes$connectivity

melExDataIndexes <- pleioIndexs[row.names(pleioIndexs) %in% row.names(melres_table_HC), ]
melExDataPleio <-  melres_table_HC[row.names(melres_table_HC) %in% melExDataIndexes$FBgn, ]
table(row.names(melExDataPleio) == row.names(melExDataIndexes))
melExDataPleio$tau <- melExDataIndexes$tau
melExDataPleio$connectivity <- melExDataIndexes$connectivity

clinalsimres_table_PM <- clinalsimres_table_PM[order(row.names(clinalsimres_table_PM)), ]
simClDataIndexes <- pleioIndexs[row.names(pleioIndexs) %in% row.names(clinalsimres_table_PM), ]
simClDataPleio <-  clinalsimres_table_PM[row.names(clinalsimres_table_PM) %in% simClDataIndexes$FBgn, ]
table(row.names(simClDataPleio) == row.names(simClDataIndexes))
simClDataPleio$tau <- simClDataIndexes$tau
simClDataPleio$connectivity <- simClDataIndexes$connectivity

clinalmelres_table_PM <- clinalmelres_table_PM[order(row.names(clinalmelres_table_PM)), ]
melClDataIndexes <- pleioIndexs[row.names(pleioIndexs) %in% row.names(clinalmelres_table_PM), ]
melClDataPleio <-  clinalmelres_table_PM[row.names(clinalmelres_table_PM) %in% melClDataIndexes$FBgn, ]
table(row.names(melClDataPleio) == row.names(melClDataIndexes))
melClDataPleio$tau <- melClDataIndexes$tau
melClDataPleio$connectivity <- melClDataIndexes$connectivity


simClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalSimSigGenes, ] #clinalSimSigGenes
nrow(simClSigGenesPleio)
#simClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalSimSigGenes01, ] #clinalSimSigGenes
#simClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalSimSigGenes015, ] #clinalSimSigGenes
#simClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalSimSigGenes02, ] #clinalSimSigGenes
simClSigGenesPleio$group <- rep("SimClinal", nrow(simClSigGenesPleio))
simClSigGenesPleio$species <- rep("Sim", nrow(simClSigGenesPleio))

melClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalMelSigGenes, ] #clinalSimSigGenes
nrow(melClSigGenesPleio)
#melClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalMelSigGenes01, ] #clinalSimSigGenes
#melClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalMelSigGenes015, ] #clinalSimSigGenes
#melClSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% clinalMelSigGenes02, ] #clinalSimSigGenes
melClSigGenesPleio$group <- rep("MelClinal", nrow(melClSigGenesPleio))
melClSigGenesPleio$species <- rep("Mel", nrow(melClSigGenesPleio))

simExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% simTempNr, ] #clinalSimSigGenes
nrow(simExSigGenesPleio)
#simExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% simTempNr01, ] #clinalSimSigGenes
#simExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% simTempNr15, ] #clinalSimSigGenes
#simExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% simTempNr20, ] #clinalSimSigGenes
simExSigGenesPleio$group <- rep("SimExperi", nrow(simExSigGenesPleio))
simExSigGenesPleio$species <- rep("Sim", nrow(simExSigGenesPleio))

melExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% melTempNr, ] #clinalSimSigGenes
nrow(melExSigGenesPleio)
#melExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% melTempNr01, ] #clinalSimSigGenes
#melExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% melTempNr015, ] #clinalSimSigGenes
#melExSigGenesPleio <- pleioIndexs[pleioIndexs$FBgn %in% melTempNr02, ] #clinalSimSigGenes
melExSigGenesPleio$group <- rep("MelExperi", nrow(melExSigGenesPleio))
melExSigGenesPleio$species <- rep("Mel", nrow(melExSigGenesPleio))

pleiotropyConnWilcox <- c()
pleiotropyConnWilcox <- c(pleiotropyConnWilcox, wilcox.test(x = simClSigGenesPleio$connectivity, y=simExSigGenesPleio$connectivity,  alternative = c("two.sided"))$p.value)
pleiotropyConnWilcox<- c(pleiotropyConnWilcox, wilcox.test(x = melClSigGenesPleio$connectivity, y=melExSigGenesPleio$connectivity,  alternative = c("two.sided"))$p.value)
names(pleiotropyConnWilcox) <- c("sim", "mel")
pleiotropyConnWilcox

pleioDataSig <- rbind(simClSigGenesPleio, melClSigGenesPleio, simExSigGenesPleio, melExSigGenesPleio)

pleioDataSig$group <- as.factor(pleioDataSig$group)

tiff("../output/PleioByLogFCBinnedConnectiviti.tiff", height = 15, width = 15, res = 300, units = "cm")
ggplot(data = pleioDataSig, aes(x = group, y = log10(connectivity), colour = group, group = group)) +
  geom_boxplot(size = 1) + 
  geom_jitter(width = 0.3,  shape=1, size = 0.5) +
  theme_classic() +
  scale_color_manual(values=c(
    viridis(n=40, option = "turbo")[3],
    viridis(n=40, option = "turbo")[37],
    viridis(n=40,option = "turbo")[8],
    viridis(n=40, option = "turbo")[30], "black")) +  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),  
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.position='none') +
  geom_signif(comparisons = list(c("MelClinal", "MelExperi")),
              annotations = "0.0074",  color = "black") +
  geom_signif(comparisons = list(c("SimClinal", "SimExperi")),
              annotations = "0.15",color = "black")+
  xlab("") +
  ylab("Pleiotropy (log10(connectivity))") 
dev.off()
pleiotropyConnWilcox



###combined
#combined p-value MELANOGASTER
stouffer(c(pleiotropyTauWilcox["sim"], pleiotropyConnWilcox["sim"]), side = 1) #check with Marlies

#combined p-value MELANOGASTER
stouffer(c(pleiotropyTauWilcox["mel"], pleiotropyConnWilcox["mel"]), side = 1)


###Corrilations
simExCorrFCConn <- cor.test(simExDataPleio$connectivity, abs(simExDataPleio$logFC), method = "spearman", exact = F)
melExCorrFCConn <- cor.test(melExDataPleio$connectivity, abs(melExDataPleio$logFC), method = "spearman",  exact = F)
simCLCorrFCConn <- cor.test(simClDataPleio$connectivity, abs(simClDataPleio$logFC), method = "spearman",  exact = F)
melCLCorrFCConn <- cor.test(melClDataPleio$connectivity, abs(melClDataPleio$logFC), method = "spearman",  exact = F)


ConnLogFCPadj <- data.frame(
  rho = c(
    simExCorrFCConn$estimate,
    melExCorrFCConn$estimate,
    simCLCorrFCConn$estimate,
    melCLCorrFCConn$estimate
  ),
  pval = c(
    simExCorrFCConn$p.value,
    melExCorrFCConn$p.value,
    simCLCorrFCConn$p.value,
    melCLCorrFCConn$p.value
  ),
  padj = p.adjust(
    c(
      simExCorrFCConn$p.value,
      melExCorrFCConn$p.value,
      simCLCorrFCConn$p.value,
      melCLCorrFCConn$p.value
    ),
    method = "BH"
  ),
  group = c(
    "simExCorrFCConn",
    "melExCorrFCConn",
    "simCLCorrFCConn",
    "melCLCorrFCConn"
  )
)

ConnLogFCPadj

###plotting
simClDataPleioBins <- simClDataPleio %>% mutate(Connectivity_bin = ntile(connectivity, n=10))
simClDataPleioBins$Connectivity_bin <- as.factor(simClDataPleioBins$Connectivity_bin)

melClDataPleioBins <- melClDataPleio %>% mutate(Connectivity_bin = ntile(connectivity, n=10))
melClDataPleioBins$Connectivity_bin <- as.factor(melClDataPleioBins$Connectivity_bin)

simExDataPleioBins <- simExDataPleio %>% mutate(Connectivity_bin = ntile(connectivity, n=10))
simExDataPleioBins$Connectivity_bin <- as.factor(simExDataPleioBins$Connectivity_bin)

melExDataPleioBins <- melExDataPleio %>% mutate(Connectivity_bin = ntile(connectivity, n=10))
melExDataPleioBins$Connectivity_bin <- as.factor(melExDataPleioBins$Connectivity_bin)


aj <- ggplot(simClDataPleioBins, aes(y= abs(logFC), x = Connectivity_bin, group = Connectivity_bin)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[8], size = 1.5, outlier.shape = NA)+
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("simulans cline", paste("rho", round(ConnLogFCPadj$rho[which(ConnLogFCPadj$group == "simCLCorrFCConn")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("A") +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        title = element_text(size=20)) +
  ylim(0,0.75) +
  xlab("Pleiotropy (Connectivity)")
aj
bj <- ggplot(melClDataPleioBins, aes(y= abs(logFC), x = Connectivity_bin)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[3], size = 1.5, outlier.shape = NA)+
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("melanogaster cline", paste("rho", round(ConnLogFCPadj$rho[which(ConnLogFCPadj$group == "melCLCorrFCConn")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("B") +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        title = element_text(size=20)) +
  ylim(0,0.75) +
  xlab("Pleiotropy (Connectivity)")
bj
cj <- ggplot(simExDataPleioBins, aes(y= abs(logFC), x = Connectivity_bin)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[30], size = 1.5, outlier.shape = NA)+
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("simulans exp.", paste("rho", round(ConnLogFCPadj$rho[which(ConnLogFCPadj$group == "simExCorrFCConn")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("C") +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        title = element_text(size=20)) +
  ylim(0,0.75) +
  xlab("Pleiotropy (Connectivity)")
cj
dj <- ggplot(melExDataPleioBins, aes(y= abs(logFC), x = Connectivity_bin)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[37], size = 1.5, outlier.shape = NA)+
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("melanogaster exp.", paste("rho", round(ConnLogFCPadj$rho[which(ConnLogFCPadj$group == "melExCorrFCConn")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("D") +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        title = element_text(size=20)) +
  ylim(0,0.75) +
  xlab("Pleiotropy (Connectivity)")
dj

tiff("../output/PleioByLogFCBinnedConnectivity.tiff", height = 30, width = 30, res = 300, units = "cm")
multiplot(aj, cj, bj, dj, cols = 2)
dev.off()



#creating a variable to plot the mean/median of the bin for the boxplot
repTimes <- c(table(simClDataPleioBins$Connectivity_bin))
str(repTimes)

conMean <- as.data.frame(simClDataPleioBins %>%
                           group_by(Connectivity_bin) %>%
                           summarize(MeanCon = mean(as.numeric(connectivity))))$MeanCon
conMean <- round(conMean, digits = 4)
conMean <- rep(conMean, times = repTimes)
simClDataPleioBins <- simClDataPleioBins[order(simClDataPleioBins$Connectivity_bin), ]
simClDataPleioBins$meanCon <- as.factor(conMean)
str(simClDataPleioBins)
aj <- ggplot(simClDataPleioBins, aes(y= abs(logFC), x = as.factor(meanCon), group = as.factor(meanCon))) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[8], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("simulans cline", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "simCLCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("A",subtitle = "Nature (D. simulans)") +
  theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (Connectivity)")
aj


repTimes <- c(table(melClDataPleioBins$Connectivity_bin))
str(repTimes)
conMean <- as.data.frame(melClDataPleioBins %>%
                           group_by(Connectivity_bin) %>%
                           summarize(MeanCon = mean(as.numeric(connectivity))))$MeanCon
conMean <- round(conMean, digits = 4)
conMean <- rep(conMean, times = repTimes)
melClDataPleioBins <- melClDataPleioBins[order(melClDataPleioBins$Connectivity_bin), ]
melClDataPleioBins$meanCon <- as.factor(conMean)
str(melClDataPleioBins)
bj <- ggplot(melClDataPleioBins, aes(y= abs(logFC), x = as.factor(meanCon), group = meanCon)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[3], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("melanogaster cline", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "melCLCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("B",subtitle = "Nature (D. melanogaster)") +
  theme(axis.text.x = element_text(size=15, angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (Connectivity)")
bj


repTimes <- c(table(simExDataPleioBins$Connectivity_bin))
str(repTimes)
conMean <- as.data.frame(simExDataPleioBins %>%
                           group_by(Connectivity_bin) %>%
                           summarize(MeanCon = mean(as.numeric(connectivity))))$MeanCon
conMean <- round(conMean, digits = 4)
conMean <- rep(conMean, times = repTimes)
simExDataPleioBins <- simExDataPleioBins[order(simExDataPleioBins$Connectivity_bin), ]
simExDataPleioBins$meanCon <- as.factor(conMean)
str(simExDataPleioBins)
cj <- ggplot(simExDataPleioBins, aes(y= abs(logFC), x = as.factor(meanCon), group = meanCon)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[30], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("simulans exp.", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "simExCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("C", subtitle = "Laboratory (D. simulans)") +
  theme(axis.text.x = element_text(size=15, angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (Connectivity)")
cj


repTimes <- c(table(melExDataPleioBins$Connectivity_bin))
str(repTimes)
conMean <- as.data.frame(melExDataPleioBins %>%
                           group_by(Connectivity_bin) %>%
                           summarize(MeanCon = mean(as.numeric(connectivity))))$MeanCon
conMean <- round(conMean, digits = 4)
conMean <- rep(conMean, times = repTimes)
melExDataPleioBins <- melExDataPleioBins[order(melExDataPleioBins$Connectivity_bin), ]
melExDataPleioBins$meanCon <- as.factor(conMean)
str(melExDataPleioBins)
dj <- ggplot(melExDataPleioBins, aes(y= abs(logFC), x = as.factor(meanCon), group = meanCon)) +
  geom_boxplot(colour = viridis(n=40, option = "turbo")[37], size = 1.5, outlier.shape = NA)+
  #geom_jitter(width = 0.15,  shape=1, size = 0.01, alpha = 0.3) +
  theme_classic() +
  #xlim(0,2) +
  ylab("absolute logFC") +
  #ggtitle(paste("melanogaster exp.", paste("rho", round(tauLogFCTauadj$rho[which(tauLogFCTauadj$group == "melExCorrFCTau")],digits = 3), sep = ": "), sep = "\n")) +
  ggtitle("D", subtitle = "Laboratory (D. melanogaster)") +
  theme(axis.text.x = element_text(size=15, angle = 45, hjust = 1),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=15, hjust = 1,face = "italic")) +
  ylim(0,0.75) +
  xlab("Pleiotropy (Connectivity)")
dj

tiff("../output/PleioByLogFCBinnedConnectivityMeanofEachBin.tiff", res = 600, height = 30, width = 30, units = "cm")
multiplot(aj, cj, bj, dj, cols = 2)
dev.off()



#density distribution of the DE genes relative to background of pleiotropy
pleioIndexs$group <- rep("Background", nrow(pleioIndexs))
masterPleio <- rbind(pleioDataSig [, 1:3], pleioIndexs)
masterPleio$group <- factor(masterPleio$group, 
                            levels = c("MelClinal", "MelExperi", "SimClinal", "SimExperi", "Background"))

#try to compare the distribution with densityplots
ii <- ggplot(data = masterPleio, aes(x = log10(connectivity), colour = group, group = group, linetype=group)) +
  geom_density(size = 1.5) + 
  theme_classic() +
  scale_linetype_manual(values=c("solid","solid", "solid","solid", "dashed"))+
  scale_color_manual(values=c(
    viridis(n=40, option = "turbo")[3],
    viridis(n=40, option = "turbo")[37],
    viridis(n=40,option = "turbo")[8],
    viridis(n=40, option = "turbo")[30], "black")) +
  ggtitle("B") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),  
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=15),
        title = element_text(size=20)) +
  xlab("Pleiotropy (log10(connectivity))")
ii




tiff("../output/PleioDensitySigGenesAllPairs.tiff", height = 15, width = 30, res = 600, dpi = 600, units = "cm")
multiplot(i, ii, cols = 2)
dev.off()



                                                                                                                       