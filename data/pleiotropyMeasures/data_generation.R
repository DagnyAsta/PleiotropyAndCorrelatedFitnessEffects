# Rscript of RNASeq data analysis and pleiotropy measurement processing
# author: Wei-Yun Lai
# date created: 2021.04.
# date last edited: 2022.08.19
# reviewed by Sheng-Kai Hsu 2022.08.21
# suggestion: Avoid using the same object name

setwd("~/Dropbox (PopGen)/backup/Wei-Yun/project2/pleiotropy/")
library(edgeR)
rm(list=ls())

####tissue specificity####
tau_bg = read.table("./tau_flyatlas2_bg.txt",sep = "\t",header = T,stringsAsFactors = F) # tau estimated from flyatlas2 (Scott's analysis 2018)
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

####network connectivity####
net.con=read.table("./flynet_supervised_0.6.txt",header=F,stringsAsFactors = F) #reference Marbach et al., 2012 genome research
#out-degree
ind=unique(net.con$V1)
out_degree=as.data.frame(matrix(NA,length(ind),2))
out_degree[,1]=as.character(unique(net.con$V1))
for (i in 1:length(ind)) {
  out_degree[i,2]=sum(net.con$V1==ind[i])
}
colnames(out_degree)=c("FBgn","out_de")
#hist(as.numeric(paste(out_degree[,2])))

#in-degree
ind=unique(net.con$V2)
in_degree=as.data.frame(matrix(NA,length(ind),2))
in_degree[,1]=as.character(unique(net.con$V2))
for (i in 1:length(ind)) {
  in_degree[i,2]=sum(net.con$V2==ind[i])
}
colnames(in_degree)=c("FBgn","in_de")
#hist(as.numeric(paste(in_degree[,2])))

conn_use=merge(in_degree,out_degree,by = "FBgn",all = T)
conn_use[which(is.na(conn_use$out_de)),3]=0
conn_use$conn=conn_use$in_de+conn_use$out_de

pleio_measure=merge(conn_use[,c(1,4)],tau_bg[,c(1,3)],by="FBgn")
pleio_measure$tau_male_4=1-pleio_measure$tau_male_4

####RNASeq data input and normalization####
dat_B_H=read.csv(file = "~/Dropbox (PopGen)/backup/Wei-Yun/monster_fullset_readcounts.csv",sep = ";")
colnames(dat_B_H)[1]="name"
dat_B_H_filtered=dat_B_H[apply(cpm(dat_B_H[,-1]),1,function(x){!any(x<=1)}),] # filter for lowly expressed genes
row.names(dat_B_H_filtered)=dat_B_H_filtered[,1]
dat_B_H_filtered=dat_B_H_filtered[,-1]

group_rep=paste(substr(colnames(dat_B_H_filtered),1,1),substr(colnames(dat_B_H_filtered),8,8),sep = "_")
group_rep[c(11,12,23,33,34)]="B"

y=DGEList(counts=dat_B_H_filtered,group = group_rep)
y=calcNormFactors(y)

dat_use=cpm(y)
b_idx=group_rep%in%"B"
dat_fc=log2(dat_use[,!b_idx]/apply(dat_use[,b_idx],1,mean)) # replicate-specific fold change
group_repH=paste(substr(colnames(dat_fc),1,1),substr(colnames(dat_fc),8,8),sep = "_")
test.rep=matrix(NA,nrow=10780,ncol = 4)
row.names(test.rep)=row.names(dat_fc)
# f value
for (i in 1:10780) {
  temp=data.frame(obs=as.numeric(dat_fc[i,]),trt=group_repH)
  fit=anova(lm(temp$obs~temp$trt))
  test.rep[i,1]=fit$`F value`[1]
  test.rep[i,2]=fit$`Pr(>F)`[1]
}
test.rep[,3]=p.adjust(test.rep[,2],method = "BH")
test.rep[,4]=row.names(test.rep)

colnames(test.rep)=c("F_value","p_value","adj_p","FBgn")
pleio_f=merge(test.rep[,c(1,4)],pleio_measure,by="FBgn")

####mean fc (magnitude of expression changes)####
dat_fc=as.data.frame(dat_fc)
dat_fc$FBgn=row.names(dat_fc)
fc_tau=merge(dat_fc,tau_bg,by="FBgn")
fc_tau$avg_fc=apply(fc_tau[,2:31],1,mean)
pleio_fc_f=merge(fc_tau[,c(1,34)],pleio_f,by="FBgn")


####ancestral variance####
# from Lai et al., 2021
dat.use=read.table(file = "~/Dropbox (PopGen)/backup/Wei-Yun/single_fly_B27_B28_H4_H9_filtered.csv",sep = ",",header = T,row.names = 1)
dat.use=dat.use[apply(cpm(dat.use),1,function(x){!(mean(x)<1)}),]

colnames(dat.use)

sample_ID=colnames(dat.use);sample_ID=as.data.frame(sample_ID)
sample_ID[,2]=strsplit2(sample_ID[,1],"_")[,2];sample_ID[,3]=strsplit2(sample_ID[,1],"_")[,4]
colnames(sample_ID)=c("sample","evo","pop")
sample_ID[,4]=paste0(sample_ID$evo,"_",sample_ID$pop)
colnames(sample_ID)=c("sample","evo","pop","evo_pop")
table(sample_ID$evo_pop)

y=DGEList(counts=dat.use,group = sample_ID$evo_pop)
y=calcNormFactors(y)

# SK: why did you do joint normalization?

#from 118-133, log-transformed the data to estimate the variance
 # count.use_old=dat.use[,1:43]
 # lib.size_old=colSums(count.use_old)
 # evo=strsplit2(colnames(count.use_old),"_")[,2]
 # y1=DGEList(counts=count.use_old,group = evo)
 # y1=calcNormFactors(y1)
 # count.use_new=dat.use[,44:81]
 # lib.size_new=colSums(count.use_new)
 # evo=strsplit2(colnames(count.use_new),"_")[,2]
 # y2=DGEList(counts=count.use_new,group = evo)
 # y2=calcNormFactors(y2)
 # 
 # dat_B_H_all = log(cpm(y))
 # dat_B_H_all_old=log(cpm(y1))
 # dat_B_H_all_new=log(cpm(y2))
 # dat_B_WY_old = dat_B_H_all_old[,which(substr(colnames(dat_B_H_all_old),3,3)=="B")]
 # dat_B_WY_new = dat_B_H_all_new[,which(substr(colnames(dat_B_H_all_new),3,3)=="B")]

y.b1=y
y.b1$counts=y.b1$counts[,y.b1$samples$group%in%"B_27"]
y.b1$samples=y.b1$samples[y.b1$samples$group%in%"B_27",]
DGE_b1=estimateDisp(y.b1)
cv_b1=apply(cpm(y.b1),1,function(x) sd(x)/mean(x))
logv_b1=apply(log(cpm(y.b1)),1,var)

y.b2=y
y.b2$counts=y.b2$counts[,y.b2$samples$group%in%"B_28"]
y.b2$samples=y.b2$samples[y.b2$samples$group%in%"B_28",]
DGE_b2=estimateDisp(y.b2)
cv_b2=apply(cpm(y.b2),1,function(x) sd(x)/mean(x))
logv_b2=apply(log(cpm(y.b2)),1,var)

anc_var_1=data.frame(FBgn=row.names(dat.use),anc_v=DGE_b1$tagwise.dispersion) # only repl. 1 was used
pleio_fc_f_anc=merge(anc_var_1,pleio_fc_f,by="FBgn")


####mean ancestral expression of each gene####
# from Lai et al., 2021 population mean in the ancestral population
dat_mean=data.frame(FBgn=row.names(dat_use),mean=apply(log(dat_use), 1, mean))
pleio_fc_f_anc_mean=merge(dat_mean,pleio_fc_f_anc,by="FBgn")
pleio_fc_f_anc_mean$tau_male_4=1-pleio_fc_f_anc_mean$tau_male_4


####all genes####
colnames(pleio_fc_f_anc_mean)=c("FBgn","mean","anc_v","fc","f-value","connectivity","ts_1-tau")
#write.table(pleio_fc_f_anc_mean,file = "pleio_evo_v2.txt",sep = "\t",col.names = T,row.names = F, quote = F)


####putatively adaptive genes####
pleio_use=read.table("./pleio_evo_v2.txt",sep = "\t",header = T,stringsAsFactors = F)

# reload data bc masked by previous script
dat_B_H=read.csv(file = "~/Dropbox (PopGen)/backup/Wei-Yun/monster_fullset_readcounts.csv",sep = ";")
colnames(dat_B_H)[1]="name"
dat_B_H_filtered=dat_B_H[apply(cpm(dat_B_H[,-1]),1,function(x){!any(x<=1)}),]

row.names(dat_B_H_filtered)=dat_B_H_filtered[,1]
dat_B_H_filtered=dat_B_H_filtered[,-1]
dat_B_H_filtered_log=log(dat_B_H_filtered)
group_rep=paste(substr(colnames(dat_B_H_filtered),1,1),substr(colnames(dat_B_H_filtered),8,8),sep = "_")
group_rep[c(11,12,23,33,34)]="B"

# DE for each rep and output the sig. gene ID
y=DGEList(counts = dat_B_H_filtered,group = group_rep)
ModelDesign_rep=model.matrix(~0+group_rep)
DGE_rep=estimateDisp(y,design = ModelDesign_rep,robust = T)
GLM_rep=glmFit(DGE_rep,design = ModelDesign_rep)

sig_ID_all=list()
for (i in 1:10){
  my.contrasts_m=rep(0,length(unique(group_rep)))
  my.contrasts_m[c(1,i+1)]=c(-1,1)
  LRT_m=glmLRT(GLM_rep,contrast = my.contrasts_m)
  res_m=LRT_m$table
  row.names(res_m)=row.names(dat_B_H_filtered)
  res_m$padj=p.adjust(res_m$PValue,method = "BH")
  sig_ID_all[[i]]=row.names(res_m[res_m$padj<0.05,])
}

# repl. frequency spectrum
idx=list()
for(i in 1:10){
  idx[[i]]=names(table(unlist(sig_ID_all)))[table(unlist(sig_ID_all))==i]
}

rep_spe=rep(0,dim(y)[1])
names(rep_spe)=rownames(y)
for(i in 1:10){
  rep_spe[idx[[i]]]=i
}

pleio_use$rep_spe=rep_spe[pleio_use$FBgn]
pleio_use_clean=na.omit(pleio_use)
pleio_use_no01=pleio_use_clean[!pleio_use_clean$rep_spe%in%c(0,1),]
#write.table(pleio_use_no01,file = "~/Dropbox (PopGen)/backup/Wei-Yun/project2/pleiotropy/pleio_putatively_adaptive_gene.txt",quote = F,sep = "\t")


