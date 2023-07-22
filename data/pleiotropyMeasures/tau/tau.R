library(data.table)
library(tidyverse)
library(gplots)
setwd("~/Dropbox (PopGen)/backup/Wei-Yun/pleiotropy/")
genefpkm <- fread("genefpkm.csv")
tissue <- fread("tissue.csv")

genefpkm %>% left_join(tissue) -> flyatlas2

flyatlas2 %>% filter(TissueName != "Whole body") -> flyatlas2

flyatlas2 %>%
  group_by(FBgn, Stage, Sex) %>%
  summarise(n_tissues = n()) -> n_tissues

flyatlas2 %>% left_join(n_tissues) -> flyatlas2

flyatlas2_male=subset(flyatlas2,subset = flyatlas2$Stage=="Adult" & flyatlas2$Sex=="Male")

flyatlas2_male %>%
  group_by(FBgn, Stage, Sex) %>%
  summarise(max_expr = max(FPKM)) -> max_expr

flyatlas2_male %>% left_join(max_expr) -> flyatlas2_male

flyatlas2_male %>%
  group_by(FBgn, Stage, Sex) %>%
  filter(FPKM == max_expr) %>%
  select(FBgn, Stage, Sex, TissueName, FPKM, max_expr) -> max_tissue

rename(max_tissue, max_tissue = TissueName) -> max_tissue

flyatlas2_male %>% left_join(select(max_tissue, FBgn, Stage, Sex, max_tissue)) -> flyatlas2_male


####pleiotropy log(FPKM)####
tau_cal2=function(x){
  x.max=max(x)
  tau=sum(1-x/x.max)/(length(x)-1)
}

tau_male_4=tapply(flyatlas2_male$FPKM,flyatlas2_male$FBgn,tau_cal2)
tau_male_4=as.data.frame(tau_male_4)
tau_male_4=cbind(row.names(tau_male_4),tau_male_4)
colnames(tau_male_4)[1]="FBgn"

# png(filename = "Distribution of tau estimate.png",height = 16,width = 20,units = "cm",res = 300,pointsize = 10)
# plot(density(tau_male_4,na.rm = T),col="red",main = "Distribution of tau estimate")
# legend("topleft",legend = c("all","adult male only"),col = c("black","red"),bty = "n",lty = c(1,1))
# dev.off()

flyatlas2_male_filter=subset(flyatlas2_male,subset = flyatlas2_male$FPKM!=0)
flyatlas2_male_filter %>%
  group_by(FBgn, Stage, Sex) %>%
  summarise(max_expr = max(FPKM)) -> max_expr
flyatlas2_male_filter %>% left_join(max_expr) -> flyatlas2_male_filter
flyatlas2_male_filter %>%
  group_by(FBgn, Stage, Sex) %>%
  filter(FPKM == max_expr) %>%
  select(FBgn, Stage, Sex, TissueName, FPKM, max_expr) -> max_tissue

rename(max_tissue, max_tissue = TissueName) -> max_tissue
flyatlas2_male_filter %>% left_join(select(max_tissue, FBgn, Stage, Sex, max_tissue)) -> flyatlas2_male_filter
test = unique(flyatlas2_male_filter[,c(1,18)])[!duplicated(unique(flyatlas2_male_filter[,c(1,18)])$FBgn),]
background=merge(test,tau_male_4,by = "FBgn")

#write.table(background,"~/Dropbox (PopGen)/backup/Wei-Yun/project2/pleiotropy/tau_flyatlas2_bg.txt",sep = "\t",quote = F,row.names = F)
