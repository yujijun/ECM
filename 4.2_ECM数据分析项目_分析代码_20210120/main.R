# Updata Time:2021.01.20 4:44pm
# Update Time:20210524
# Yujijun
 --------------------   
#   文件说明：
# 4.1_ECM数据分析项目_原始数据_20201231  # 此为原始数据文件【主要是质谱数据结果和RNA-seq结果】
# 4.2_ECM数据分析项目_分析代码_20210120  # 此为分析代码
# 4.3_ECM数据分析项目_输出结果_20210120  # 输出的中间过程文件存在在此【不包含图片】
# 4.4_ECM数据分析项目_图片输出_20210120  # output figure[include original figures and AI changed version]
# 4.5_ECM数据分析项目_conclusion_20201231  # Conclusion and disscussion powerpoint
 ---------------------

#### 01 load all package ####
library(ggplot2)
library(readxl)
library(plyr)
library(tidyverse)
library(export)
library(tidyfst)
library(VennDiagram)
library(gplots)
library(ggstatsplot)
library(RColorBrewer)
library(reshape2)
library(Rmisc)
library(pheatmap)
library(ggrepel)
library(readxl)
library(riverplot)
library(ggpubr)  
library(dplyr)
library(tidyr)
 
#### 02 Plot hyper-parameter ####
#plot title 16 bold;
 
if(T){
  mytheme <- theme_classic() + 
    theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black",face = "bold"), 
                   axis.text = element_text(size= 10,color = "black"),
                   axis.text.x = element_text(angle = 45, hjust = 1),
                   panel.grid=element_blank(),
                   legend.position= "top",
                   legend.box = "horizontal",
                   legend.text = element_text(size= 4),
                   legend.title= element_text(size= 4,face = "bold"),
                   text = element_text(size = 12,family = "Arial")) 
}
width = 3.35
width2 <- 7.01
height = 3.35
col_convert <- c("#B2182B","#2166AC")
color_4 <- c("#C58680","#B2182B","#87A3D2","#2166AC")
color_category <- c(brewer.pal(n = 8, name = "Pastel2")[1:6])
names(color_category) <- c("secreted_factors","proteoglycans","ECM-affiliated",
                           "regulators","glycoproteins","collagens")
figure_path <- "./4.4_ECM数据分析项目_图片输出_20210120/version8_last_V2_20210531/"
output_path <-"./4.3_ECM数据分析项目_输出结果_20210120/output_last/"

#### 03 Data input and cleaning ####
data <- read_excel("./4.1_ECM数据分析项目_原始数据_20201231/Table1 proteomics of DL-ECM vs Matrigel.xlsx")
#data <- read_excel("./4.1_ECM数据分析项目_原始数据_20201231/Annotation_combine.xlsx")
load("./4.1_ECM数据分析项目_原始数据_20201231/merge_orig_ref_filter.RData")
selectresult=subset(data,!(is.na(DL1)&is.na(DL2)&is.na(DL3)&is.na(Matrigel1)&is.na(Matrigel2)&is.na(Matrigel3)))
datafilter <- selectresult[!(selectresult$`Gene name` == "--" | is.na(selectresult$`Gene name`)), ]
test1<- count_dt(datafilter,`Gene name`)%>% filter_dt(n>1)
datadup= datafilter[which(datafilter$`Gene name` %in% test1$`Gene name`),]
p<- c('P02770','P09605','D3ZQ25','P04937','A0A0G2KAM4','A0A096MJI4','A0A0G2KAY3','D3ZN05','A0A0G2K484','A0A0G2K8C1','M0RBV9','A0A0G2JVG3','A0A0G2K781','Q5RKJ9','P16086','Q9EQT5','A0A0G2JSQ4','A0A140TAF0','Q9JJP9')
anno_filtered = datafilter[which(!datafilter$`Protein accession` %in% p),]
output_name <- paste0(output_path,"annotation_filtered.csv")
write_excel_csv(anno_filtered,output_name)
output_name <- paste0(output_path,"annotation_filtered.RData")
save(anno_filtered,file = output_name)

merge_orig_ref1<- count_dt(merge_orig_ref,`Gene name`)%>% filter_dt(n>1)
merge_orig_ref2 <- merge_orig_ref[which(merge_orig_ref$`Gene name` %in% merge_orig_ref1$`Gene name`),]
merge_orig_ref_filter <- merge_orig_ref[-c(67,76,96,101,134,169),]
choosen_column <- c(1,9:16,31)
merge_orig_ref_choosen <- merge_orig_ref_filter[,choosen_column]
merge_na <- merge_orig_ref_choosen[!(rowSums(is.na(merge_orig_ref_choosen)) > 7),]
rownames(merge_orig_ref_filter) <- merge_orig_ref_filter$`Gene name`
merge_orig_ref_filter_na <- merge_orig_ref_filter[merge_na$`Gene name`,]
output_name <- paste0(output_path,"merge_orig_ref_filter.csv")
write_excel_csv(merge_orig_ref_filter_na,output_name)
output_name <- paste0(output_path,"merge_orig_ref_filter.RData")
save(merge_orig_ref_filter_na,file = output_name)

#### 04 Visualization ####
## 1-Barplot for all genes ####
anno_filtered$EX_IN <- mapvalues(anno_filtered$`Subcellular localization`,
                                from = c("cytoplasm", "cytoplasm,mitochondria","cytoplasm,nucleus", 
                                  "cytoplasm,peroxisome","cytoskeleton","endoplasmic reticulum",
                                  "Golgi apparatus","mitochondria","nucleus","peroxisome",
                                  "extracellular","plasma membrane","endoplasmic reticulum,mitochondria",
                                  "extracellular,plasma membrane","mitochondria,nucleus"),
                                to=c("IN","IN","IN","IN","IN","IN","IN","IN","IN","IN","EX","EX","IN","EX","IN"))
anno_filtered[which(anno_filtered$`Gene name`=='Nup214'),]$EX_IN <- 'IN'
anno_filtered[which(anno_filtered$`Gene name`=='Lama3'),]$EX_IN <- 'EX'
###DL gene number###
DL_mat <- anno_filtered[,c("DL1","DL2","DL3","EX_IN")]
na.function <- function(DL_mat){
  DL_mat_na <- is.na(DL_mat)
  ratio_na <- rowSums(DL_mat_na)
  DL_mat_final <- DL_mat[ratio_na < 3,]
  return(DL_mat_final)
}
DL_mat_na <- na.function(DL_mat)
DL_table <- as.data.frame(table(DL_mat_na$EX_IN))
###M gene number###
ML_mat <- anno_filtered[,c("Matrigel1","Matrigel2","Matrigel3","EX_IN")]
na.function <- function(ML_mat){
  ML_mat_na <- is.na(ML_mat)
  ratio_na <- rowSums(ML_mat_na)
  ML_mat_final <- ML_mat[ratio_na < 3,]
  return(ML_mat_final)
}
ML_mat_na <- na.function(ML_mat)
Matrigel_table <- as.data.frame(table(ML_mat_na$EX_IN))
colnames(DL_table) <- c("Position","Freq")
colnames(Matrigel_table) <- c("Position","Freq")
DL_table$expr <- rep("DL-ECM",2)
Matrigel_table$expr <- rep("Matrigel",2)
All_table <- rbind(DL_table,Matrigel_table)
All_table$color <- paste0(All_table$expr,"_",All_table$Position)
## visulizaiton 
plot01_barplot <- ggplot(data = All_table,mapping = aes(x=expr,y=Freq,fill=color)) + 
  geom_bar(stat = "identity",width = 0.5) + 
  ylab("Average protein number") + 
  xlab("Experiment") + 
  scale_fill_manual(values = color_4)+  
  geom_text(aes(label=Freq), position = position_stack(0.5), color = "white")  + 
  guides(fill=guide_legend("Position")) + 
  mytheme
figure_name <- paste0(figure_path,"01_barplot_forallgene.pdf")
graph2pdf(x = plot01_barplot, file=figure_name,
          font = "Arial",
          width = width, 
          height = height/2, 
          bg = "transparent")

## 2-Vennplot about IN and EX####
DL_mat1 <- anno_filtered[,c("DL1","DL2","DL3","Gene name","EX_IN")]
na.function <- function(DL_mat1){
  DL_mat_na <- is.na(DL_mat1)
  ratio_na <- rowSums(DL_mat_na)
  DL_mat_final <- DL_mat1[ratio_na < 3,]
  return(DL_mat_final)
}
DL_mat_na1 <- na.function(DL_mat1)

ML_mat1 <- anno_filtered[,c("Matrigel1","Matrigel2","Matrigel3","Gene name","EX_IN")]
na.function <- function(ML_mat1){
  ML_mat_na <- is.na(ML_mat1)
  ratio_na <- rowSums(ML_mat_na)
  ML_mat_final <- ML_mat1[ratio_na < 3,]
  return(ML_mat_final)
}
ML_mat_na1 <- na.function(ML_mat1)

DL_IN <- DL_mat_na1$`Gene name`[which(DL_mat_na1$EX_IN == "IN")]
DL_EX <- DL_mat_na1$`Gene name`[which(DL_mat_na1$EX_IN == "EX")]
Matrigel_IN <- ML_mat_na1$`Gene name`[which(ML_mat_na1$EX_IN == "IN")]
Matrigel_EX <- ML_mat_na1$`Gene name`[which(ML_mat_na1$EX_IN == "EX")]

plot02_vennplo_IN <- venn.diagram(list(DL_IN=DL_IN,Matrigel_IN=Matrigel_IN),filename = NULL, 
                          height = 8.5, width = 8.5,resolution =300, 
                          units = "cm",main.fontfamily = "Arial",
                          imagetype="pdf", col="transparent",fill=color_4[c(2,4)],
                          alpha = 0.80, cex=0.9, cat.cex=0.9,force.unique = F)
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)
plot02_vennplot_IN <- grid.draw(plot02_vennplo_IN)
figure_name <- paste0(figure_path,"02_IN_Venn.pdf")
graph2pdf(x = plot02_vennplot_IN,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

plot02_vennplot_EX <- venn.diagram(list(DL_EX=DL_EX,Matrigel_EX=Matrigel_EX), filename = NULL, 
                          height = 8.5, width = 8.5,resolution =300, 
                          units = "cm", main.fontfamily = "Arial",
                          imagetype="pdf", col="transparent",fill=color_4[c(1,3)],
                          alpha = 0.80, cex=0.9, cat.cex=0.9,force.unique = F)
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)
plot02_vennplot_EX <- grid.draw(plot02_vennplot_EX)
figure_name <- paste0(figure_path,"02_EX_Venn.pdf")
graph2pdf(x = plot02_vennplot_EX, 
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")


## 3-Enrichment figure for all genes ####
# statistic the DL emrichment
data1 <- anno_filtered[,c("DL1","DL2","DL3","Matrigel1","Matrigel2","Matrigel3","Cellular Component")]
names(data1)[names(data1) == "Cellular Component"] <- "Cellular.Component"
m<-c("GO:0005576","GO:0005887","GO:0005856","GO:0005739","GO:0005737","GO:0043229","GO:0005777","GO:0042579","GO:0005768","GO:0005840","GO:0009986","GO:0005694","GO:0005783","GO:0032991","GO:0005634")
a<-data1[0,]
f<-function(x) sum(!is.na(x))
for (gene in m){
  test<- apply(data1[grep(gene,data1$`Cellular.Component`),],2,f)
  test$`Cellular.Component` <- gene
  a<-rbind(a,test)
}
a$DL <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, mean,na.rm=T)
a$DL_sd <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, sd,na.rm=T)
a$Matrigel <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, mean,na.rm=T)
a$Matrigel_sd <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, sd,na.rm=T)

a$DL <- a$DL/1482
a$DL_sd <- a$DL_sd/1482
a$Matrigel <- a$Matrigel/2004
a$Matrigel_sd <- a$Matrigel_sd/2004
a$CC <- mapvalues(a$`Cellular.Component`,
                  from = c("GO:0005576", "GO:0005887","GO:0005856","GO:0005739","GO:0005737","GO:0043229","GO:0005777","GO:0042579","GO:0005768","GO:0005840","GO:0009986","GO:0005694","GO:0005783","GO:0032991","GO:0005634"),
                  to=c("Extracellular","Plasma membrane","Cytoskeleton","Mitochondria","Cytoplasm","Other intracellular organelle","Peroxisome","Microbody","Endosome","Ribosome","Cell Surface","Chromosome","Endoplasmic Reticulum","Macromolecular Complex","Nucleus"))
data2 <- melt(a,
              id.vars = c("CC","DL_sd","Matrigel_sd"),
              measure.vars = c("DL","Matrigel"),
              variable.name = "condition",na.rm = T)
data2$se <- ifelse(data2$condition=="DL",data2$DL_sd,data2$Matrigel_sd)

test<-factor(a[sort(a$DL-a$Matrigel,index.return=TRUE,decreasing=TRUE)$ix,]$CC)
data2<-data2[order(factor(data2$CC, levels = test)),]
data2$CC <- factor(data2$CC, levels=unique(as.character(data2$CC)) )

plot03_CC_avg <- ggplot(data2, aes(x=CC, y=value, fill=condition)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  coord_flip() +
  ylab("Precentage") +
  xlab("") + 
  ggtitle("Cellular Component") +
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 12)) + 
  scale_fill_manual(values = color_4[c(2,4)]) + 
  mytheme + 
  theme(legend.position = "none")
figure_name <- paste0(figure_path,"03_Cellular.Component.avg.pdf")
graph2pdf(x=plot03_CC_avg,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height, bg = "transparent")

### cellular.component.total
a$DL <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, mean,na.rm=T)
a$DL_sd <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, sd,na.rm=T)
a$Matrigel <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, mean,na.rm=T)
a$Matrigel_sd <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, sd,na.rm=T)
a$CC <- mapvalues(a$`Cellular.Component`,
                  from = c("GO:0005576", "GO:0005887","GO:0005856","GO:0005739","GO:0005737","GO:0043229","GO:0005777","GO:0042579","GO:0005768","GO:0005840","GO:0009986","GO:0005694","GO:0005783","GO:0032991","GO:0005634"),
                  to=c("Extracellular","Plasma membrane","Cytoskeleton","Mitochondria","Cytoplasm","Other intracellular organelle","Peroxisome","Microbody","Endosome","Ribosome","Cell Surface","Chromosome","Endoplasmic Reticulum","Macromolecular Complex","Nucleus"))

data2 <- melt(a,
              id.vars = c("CC","DL_sd","Matrigel_sd"),
              measure.vars = c("DL","Matrigel"),
              variable.name = "condition",na.rm = T)
data2$se <- ifelse(data2$condition=="DL",data2$DL_sd,data2$Matrigel_sd)

test<-factor(a[sort(a$DL-a$Matrigel,index.return=TRUE,decreasing=TRUE)$ix,]$CC)
data2<-data2[order(factor(data2$CC, levels = test)),]
data2$CC <- factor(data2$CC, levels=unique(as.character(data2$CC)))

plot03.CC.total <- ggplot(data2, aes(x=`CC`, y=value, fill=condition)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  coord_flip() +
  ylab("Total protein number") +
  xlab("") + 
  ggtitle("Cellular Component") +
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  scale_fill_manual(values = color_4[c(2,4)]) + 
  mytheme + 
  theme(legend.position = "none")
figure_name <- paste0(figure_path,"03_Cellular.Component.total.pdf")
graph2pdf(x= plot03.CC.total,file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

## Biological Process of average
data1 <- anno_filtered[,c("DL1","DL2","DL3","Matrigel1","Matrigel2","Matrigel3","Biological Process")]
names(data1)[names(data1) == "Biological Process"] <- "Biological.Process"
m<-c("GO:0065007","GO:0009987","GO:0032502","GO:0008152","GO:0002376","GO:0000003","GO:0051179")
a<-data1[0,]
f<-function(x) sum(!is.na(x))
for (gene in m){
  test<- apply(data1[grep(gene,data1$`Biological.Process`),],2,f)
  test$`Biological.Process` <- gene
  a<-rbind(a,test)
}

a$DL <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, mean,na.rm=T)
a$DL_sd <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, sd,na.rm=T)
a$Matrigel <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, mean,na.rm=T)
a$Matrigel_sd <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, sd,na.rm=T)
a$DL <- a$DL/1482
a$DL_sd <- a$DL_sd/1482
a$Matrigel <- a$Matrigel/2004
a$Matrigel_sd <- a$Matrigel_sd/2004
a$BP <- mapvalues(a$`Biological.Process`,
                  from = c("GO:0065007","GO:0051179","GO:0032502","GO:0009987","GO:0008152","GO:0002376","GO:0000003"),
                  to=c("Biological regulation","Localization","Developmental process","Cellular process","Metabolic process","Immune system process","Reproduction"))

data2 <- melt(a,
              id.vars = c("BP","DL_sd","Matrigel_sd"),
              measure.vars = c("DL","Matrigel"),
              variable.name = "condition",na.rm = T)
data2$se <- ifelse(data2$condition=="DL",data2$DL_sd,data2$Matrigel_sd)
test<-factor(a[sort(a$DL-a$Matrigel,index.return=TRUE,decreasing=TRUE)$ix,]$BP)
data2<-data2[order(factor(data2$BP, levels = test)),]
data2$BP <- factor(data2$BP, levels=unique(as.character(data2$BP)) )

plot03.BP.precentage <- ggplot(data2, aes(x=`BP`, y=value, fill=condition)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  coord_flip() +
  ylab("Precentage") +
  xlab("") + 
  ggtitle("Biological Process") +
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(text = element_text(face = "bold")) + 
  scale_fill_manual(values = color_4[c(2,4)]) + mytheme + 
  theme(legend.position = "none")
figure_name <- paste0(figure_path,"03_Biological.Process.avg.pdf")
graph2pdf(x=plot03.BP.precentage,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height/2, 
          bg = "transparent")

#Biological Process of total 
a$DL <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, mean,na.rm=T)
a$DL_sd <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, sd,na.rm=T)
a$Matrigel <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, mean,na.rm=T)
a$Matrigel_sd <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, sd,na.rm=T)
a$BP <- mapvalues(a$`Biological.Process`,
                  from = c("GO:0065007","GO:0051179","GO:0032502","GO:0009987","GO:0008152","GO:0002376","GO:0000003"),
                  to=c("Biological regulation","Localization","Developmental process","Cellular process","Metabolic process","Immune system process","Reproduction"))
data2 <- melt(a,
              id.vars = c("BP","DL_sd","Matrigel_sd"),
              measure.vars = c("DL","Matrigel"),
              variable.name = "condition",na.rm = T)
data2$se <- ifelse(data2$condition=="DL",data2$DL_sd,data2$Matrigel_sd)

test<-factor(a[sort(a$DL-a$Matrigel,index.return=TRUE,decreasing=TRUE)$ix,]$BP)
data2<-data2[order(factor(data2$BP, levels = test)),]
data2$BP <- factor(data2$BP, levels=unique(as.character(data2$BP)) )

plot03_BP_total <- ggplot(data2, aes(x=`BP`, y=value, fill=condition)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  coord_flip() +
  ylab("Total protein number") +
  xlab("") + 
  ggtitle("Biological Process") +
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  scale_fill_manual(values = color_4[c(2,4)]) +
  mytheme
figure_name <- paste0(figure_path,"03_Biological.Process_total.pdf")
graph2pdf(x = plot03_BP_total,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

## MF average
data1 <- anno_filtered[,c("DL1","DL2","DL3","Matrigel1","Matrigel2","Matrigel3","Molecular Function")]
names(data1)[names(data1) == "Molecular Function"] <- "Molecular.Function"
####View(data1)
m<-c("GO:0003824","GO:0005198","GO:0060089","GO:0016209","GO:0005488")
a<-data1[0,]
f<-function(x) sum(!is.na(x))
for (gene in m){
  test<- apply(data1[grep(gene,data1$`Molecular.Function`),],2,f)
  test$`Molecular.Function` <- gene
  a<-rbind(a,test)
}

a$DL <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, mean,na.rm=T)
a$DL_sd <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, sd,na.rm=T)
a$Matrigel <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, mean,na.rm=T)
a$Matrigel_sd <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, sd,na.rm=T)
a$DL <- a$DL/1482
a$DL_sd <- a$DL_sd/1482
a$Matrigel <- a$Matrigel/2004
a$Matrigel_sd <- a$Matrigel_sd/2004
a$MF <- mapvalues(a$`Molecular.Function`,
                  from = c("GO:0060089","GO:0016209","GO:0005488","GO:0005198","GO:0003824"),
                  to=c("Molecular transducer activity","Antioxidant activity","Binding","Structural molecule activity","Catalytic activity"))
data2 <- melt(a,
              id.vars = c("MF","DL_sd","Matrigel_sd"),
              measure.vars = c("DL","Matrigel"),
              variable.name = "condition",na.rm = T)
####View(data2)
data2$se <- ifelse(data2$condition=="DL",data2$DL_sd,data2$Matrigel_sd)

test<-factor(a[sort(a$DL-a$Matrigel,index.return=TRUE,decreasing=TRUE)$ix,]$MF)
data2<-data2[order(factor(data2$MF, levels = test)),]
data2$MF <- factor(data2$MF, levels=unique(as.character(data2$MF)) )

plot03_MF_avg <- ggplot(data2, aes(x=`MF`, y=value, fill=condition)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  coord_flip() +
  ylab("Precentage") +
  xlab("") + 
  ggtitle("Molecular Function") +
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  scale_fill_manual(values = color_4[c(2,4)]) + 
  mytheme + 
  theme(legend.position = "none")
figure_name <- paste0(figure_path, "03_Molecular.Function_avg.pdf")
graph2pdf(x = plot03_MF_avg,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height/2,
          bg = "transparent")

a$DL <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, mean,na.rm=T)
a$DL_sd <- apply(data.frame(a$DL1,a$DL2,a$DL3), 1, sd,na.rm=T)
a$Matrigel <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, mean,na.rm=T)
a$Matrigel_sd <- apply(data.frame(a$Matrigel1,a$Matrigel2,a$Matrigel3), 1, sd,na.rm=T)
a$MF <- mapvalues(a$`Molecular.Function`,
                  from = c("GO:0060089","GO:0016209","GO:0005488","GO:0005198","GO:0003824"),
                  to=c("Molecular transducer activity","Antioxidant activity","Binding","Structural molecule activity","Catalytic activity"))

data2 <- melt(a,
              id.vars = c("MF","DL_sd","Matrigel_sd"),
              measure.vars = c("DL","Matrigel"),
              variable.name = "condition",na.rm = T)
####View(data2)
data2$se <- ifelse(data2$condition=="DL",data2$DL_sd,data2$Matrigel_sd)

test<-factor(a[sort(a$DL-a$Matrigel,index.return=TRUE,decreasing=TRUE)$ix,]$MF)
data2<-data2[order(factor(data2$MF, levels = test)),]
data2$MF <- factor(data2$MF, levels=unique(as.character(data2$MF)) )

plot03_MF_total <- ggplot(data2, aes(x=`MF`, y=value, fill=condition)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  coord_flip() +
  ylab("Total protein number") +
  xlab("") + 
  ggtitle("Molecular Function") +
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  scale_fill_manual(values = color_4[c(2,4)]) + 
  mytheme
figure_name <- paste0(figure_path,"03_Molecular.Function_total.pdf")
graph2pdf(x = plot03_MF_total,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")
## 4-Pie plot for categories####
DL_mat <- merge_orig_ref_filter[,c("DL1","DL2","DL3","Category")]
na.function <- function(DL_mat){
  DL_mat_na <- is.na(DL_mat)
  ratio_na <- rowSums(DL_mat_na)
  DL_mat_final <- DL_mat[ratio_na < 3,]
  return(DL_mat_final)
}
DL_mat_naomit <- na.function(DL_mat)
plot04_pie_DL <- ggstatsplot::ggpiestats(DL_mat_naomit, x = 'Category', 
                        results.subtitle = F, #标题中不显示统计结果
                        factor.levels = c("ECM Regulators", "ECM Glycoproteins", "Secreted Factors","ECM-affiliated Proteins", "Proteoglycans","Collagens"),#设置标签的名称
                        slice.label = 'percentage', #设置饼图中标签的类型（默认percentage：百分比, 还有counts：频数, both : 频数和百分比都显示）
                        perc.k = 2, #百分比的小数位数为2
                        direction = 1, #1为顺时针方向，-1为逆时针方向
                        palette = 'Pastel2', #设置调色板
                        title = NULL,
                        legend.title = NULL) + 
  mytheme
graph2ppt(file=paste0(figure_path,"05_D_pie_v2.pptx"),
          width = width, 
          height = height, 
          bg = "transparent")
figure_name <- paste0(figure_path,"05_D_pie_v2.pdf")
graph2pdf(x=plot04_pie_DL,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

ML_mat <- merge_orig_ref_filter[,c("Matrigel1","Matrigel2","Matrigel3","Category")]
na.function <- function(ML_mat){
  ML_mat_na <- is.na(ML_mat)
  ratio_na <- rowSums(ML_mat_na)
  ML_mat_final <- ML_mat[ratio_na < 3,]
  return(ML_mat_final)
}
ML_mat_naomit <- na.function(ML_mat)
plot04_pie_M <- ggstatsplot::ggpiestats(ML_mat_naomit, 
                        x = 'Category', 
                        results.subtitle = F, #标题中不显示统计结果
                        factor.levels = c("ECM Regulators", "ECM Glycoproteins", "Secreted Factors","ECM-affiliated Proteins", "Proteoglycans","Collagens"),#设置标签的名称
                        slice.label = 'percentage', #设置饼图中标签的类型（默认percentage：百分比, 还有counts：频数, both : 频数和百分比都显示）
                        perc.k = 2, #百分比的小数位数为2
                        direction = 1, #1为顺时针方向，-1为逆时针方向
                        palette = 'Pastel2', #设置调色板
                        title = NULL) + 
  theme(legend.position = "none")
graph2ppt(file=paste0(figure_path,"05_M_pie.pptx"))
figure_name <- paste0(figure_path,"05_M_pie.pdf")
graph2pdf(x=plot04_pie_M,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = width, 
          bg = "transparent")

## 4.1-Pie plot for categories_v2####
DL_mat <- merge_orig_ref_filter_na[,c("DL1","DL2","DL3","Category")]
na.function <- function(DL_mat){
  DL_mat_na <- is.na(DL_mat)
  ratio_na <- rowSums(DL_mat_na)
  DL_mat_final <- DL_mat[ratio_na < 3,]
  return(DL_mat_final)
}

DL_mat_naomit <- na.function(DL_mat)
ggstatsplot::ggpiestats(DL_mat_naomit, 'Category', 
                        results.subtitle = F, #标题中不显示统计结果
                        factor.levels = c("ECM Regulators", "ECM Glycoproteins", "Secreted Factors","ECM-affiliated Proteins", "Proteoglycans","Collagens"),#设置标签的名称
                        slice.label = 'percentage', #设置饼图中标签的类型（默认percentage：百分比, 还有counts：频数, both : 频数和百分比都显示）
                        perc.k = 2, #百分比的小数位数为2
                        direction = 1, #1为顺时针方向，-1为逆时针方向
                        palette = 'Pastel2', #设置调色板
                        title = 'Percent of Category for D')#设置标题 
#graph2ppt(file="./04_Output/Version3/04_D_pie.pptx")
figure_name <- paste0(figure_path,"04_D_pie.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

ML_mat <- merge_orig_ref_filter_na[,c("Matrigel1","Matrigel2","Matrigel3","Category")]
na.function <- function(ML_mat){
  ML_mat_na <- is.na(ML_mat)
  ratio_na <- rowSums(ML_mat_na)
  ML_mat_final <- ML_mat[ratio_na < 3,]
  return(ML_mat_final)
}
ML_mat_naomit <- na.function(ML_mat)
ggstatsplot::ggpiestats(ML_mat_naomit, 'Category', 
                        results.subtitle = F, #标题中不显示统计结果
                        factor.levels = c("ECM Regulators", "ECM Glycoproteins", "Secreted Factors","ECM-affiliated Proteins", "Proteoglycans","Collagens"),#设置标签的名称
                        slice.label = 'percentage', #设置饼图中标签的类型（默认percentage：百分比, 还有counts：频数, both : 频数和百分比都显示）
                        perc.k = 2, #百分比的小数位数为2
                        direction = 1, #1为顺时针方向，-1为逆时针方向
                        palette = 'Pastel2', #设置调色板
                        title = 'Percent of Category for M')#设置标题
figure_name <- paste0(figure_path,"04_M_pie_v2.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

library(dplyr)
data_group<- group_by(DL_mat_naomit, Category)
DL_count<- summarise(group_by(DL_mat_naomit, Category),count = n())
DL_count$name ="DL"
data_group<- group_by(ML_mat_naomit, Category)
ML_count<- summarise(group_by(ML_mat_naomit, Category),count = n())
ML_count$name ="ML"

ML_DL_count <- rbind(ML_count,DL_count)
#ML_DL_count$Category <- factor(rev(ML_DL_count$Category),ordered = T)

ML_DL_count<-ML_DL_count[order(factor(ML_DL_count$Category, levels = c("Collagens","ECM Glycoproteins","Proteoglycans","ECM-affiliated Proteins","ECM Regulators","Secreted Factors"))),]
ML_DL_count$Category <- factor(ML_DL_count$Category, levels=unique(as.character(ML_DL_count$Category)) )

ggplot(ML_DL_count, aes(x=name, y=count, fill=Category))+
  geom_bar(width = 0.5, stat = "identity")+
  coord_flip()+
  scale_fill_manual(name="Category",values = rev(c("#afdac7","#edc7dd","#cad3e5","#f8cbaa","#e2eec6","#f8ecac")))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.x = element_text(face = "bold",size = 10))+ 
  mytheme + theme(legend.position = "none")
figure_name <- paste0(figure_path,"05_M_DL_barplot_v2.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width2, 
          height = height, 
          bg = "transparent")


## 4.2-another pie plot for categories_v3 ####
DL_mat <- anno_filtered[,c("DL1","DL2","DL3","Subcellular localization")]
na.function <- function(DL_mat){
  DL_mat_na <- is.na(DL_mat)
  ratio_na <- rowSums(DL_mat_na)
  DL_mat_final <- DL_mat[ratio_na < 3,]
  return(DL_mat_final)
}
DL_mat_naomit <- na.function(DL_mat)
data<-DL_mat_naomit[,"Subcellular localization"]
data[which(data$`Subcellular localization`=='extracellular'),]$`Subcellular localization` <- 'Extracellular'
data[which(data$`Subcellular localization`%in% 
             c('cytoplasm',
               'cytoskeleton',
               'mitochondria',
               'endoplasmic reticulum',
               'peroxisome',
               'Golgi apparatus',
               'cytoplasm,peroxisome',
               'cytoplasm,mitochondria',
               'endoplasmic reticulum,mitochondria')),]$`Subcellular localization` <- 'Cytoplasm'
data[which(data$`Subcellular localization`=='nucleus'),]$`Subcellular localization` <- 'Nucleus'
data[which(data$`Subcellular localization`=='plasma membrane'),]$`Subcellular localization` <- 'Piasmamembrane'
data[which(data$`Subcellular localization`%in% 
             c("cytoplasm,nucleus",
               "mitochondria,nucleus",
               "extracellular,plasma membrane")),]$`Subcellular localization` <- 'Others'
data_group<- group_by(data, `Subcellular localization`)
data_GroupByID<- summarise(data_group,count = n())
data_GroupByID$per<- round(data_GroupByID$count/sum(data_GroupByID$count),3)*100
data_GroupByID$label <- paste(as.character(data_GroupByID$per),"%", data_GroupByID$`Subcellular localization`, sep = "")
data_GroupByID$`Subcellular localization` <- factor(data_GroupByID$`Subcellular localization`)
#data_GroupByID$label <- factor(data_GroupByID$label,levels =c("15% Extracellular","4% Piasmamembrane","30.1% Nucleus","45.5% Cytoplasm","5.4% Others"), ordered = T)
DL_count<- data.frame(count=data_GroupByID$count,name ='DL',label = data_GroupByID$`Subcellular localization`)
ggplot(data_GroupByID, aes(x = 2, y = count, fill = label)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  theme(axis.text = element_blank()) + 
  theme(axis.ticks = element_blank()) + 
  coord_polar(theta = "y", start = 0)+
  labs(x = '', y = '', title = '') + 
  theme_void()+
  theme(legend.position = "none") +
  scale_fill_manual(name="Percentage",values = c('#B25751','#5A80B8','#A1BA66','#7C659E','#65AAC3')) +
  xlim(0.5, 2.5) 

figure_name <- paste0(figure_path,"04_D_pie_subcell1.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width/2, 
          height = height/2, 
          bg = "transparent")
#graph2ppt(file="./04_Output/Version3/10_D_pie_subcell.pptx")
# graph2pdf(file="./04_Output/Version3/10_D_pie_subcell1.pdf",aspectr=2, font = "Arial",
#           width = 7, height = 5, bg = "transparent")

ML_mat <- anno_filtered[,c("Matrigel1","Matrigel1","Matrigel1","Subcellular localization")]
na.function <- function(ML_mat){
  ML_mat_na <- is.na(ML_mat)
  ratio_na <- rowSums(ML_mat_na)
  ML_mat_final <- ML_mat[ratio_na < 3,]
  return(ML_mat_final)
}
ML_mat_naomit <- na.function(ML_mat)
data<-ML_mat_naomit[,"Subcellular localization"]
data[which(data$`Subcellular localization`=='extracellular'),]$`Subcellular localization` <- 'Extracellular'
data[which(data$`Subcellular localization` %in% 
             c('cytoplasm',
               'cytoskeleton',
               'mitochondria',
               'endoplasmic reticulum',
               'peroxisome',
               'Golgi apparatus',
               'cytoplasm,peroxisome',
               'cytoplasm,mitochondria',
               'endoplasmic reticulum,mitochondria')),]$`Subcellular localization` <- 'Cytoplasm'
data[which(data$`Subcellular localization`=='nucleus'),]$`Subcellular localization` <- 'Nucleus'
data[which(data$`Subcellular localization`=='plasma membrane'),]$`Subcellular localization` <- 'Piasmamembrane'
data[which(data$`Subcellular localization` %in% 
             c("cytoplasm,nucleus",
               "mitochondria,nucleus",
               "extracellular,plasma membrane")),]$`Subcellular localization` <- 'Others'

data_group<- group_by(data, `Subcellular localization`)
data_GroupByID<- summarise(data_group,count = n())

data_GroupByID$per<- round(data_GroupByID$count/sum(data_GroupByID$count),3)*100
data_GroupByID$label <- paste(as.character(data_GroupByID$per),"% " , data_GroupByID$`Subcellular localization`, sep = "")
data_GroupByID$`Subcellular localization` <- factor(data_GroupByID$`Subcellular localization`)
#data_GroupByID$label <- factor(data_GroupByID$label,levels =c("12.4% Extracellular","4% Piasmamembrane","33.8% Nucleus","44.4% Cytoplasm","5.5% Others"), ordered = T)
ML_count<- data.frame(count=data_GroupByID$count,name ='ML',label = data_GroupByID$`Subcellular localization`)
ggplot(data_GroupByID, aes(x = 2, y = count, fill = label)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  theme(axis.text = element_blank()) + 
  theme(axis.ticks = element_blank()) + 
  coord_polar(theta = "y", start = 0)+
  labs(x = '', y = '', title = '') + 
  theme_void()+
  theme(legend.position = "none") +
  scale_fill_manual(name="Percentage",values = c('#B25751','#5A80B8','#A1BA66','#7C659E','#65AAC3')) +
  xlim(0.5, 2.5)
figure_name <- paste0(figure_path,"04_M_pie_subcell1.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width/2, 
          height = height/2, 
          bg = "transparent")
# graph2ppt(file="./04_Output/Version3/10_M_pie_subcell.pptx")
# graph2pdf(file="./04_Output/Version3/10_M_pie_subcell1.pdf",aspectr=2, font = "Arial",
#           width = 7, height = 5, bg = "transparent")

#Barplot
ML_DL_count <- rbind(ML_count,DL_count)
ggplot(ML_DL_count, aes(x=name, y=count, fill=label))+
  geom_bar(width = 0.5, stat = "identity")+
  scale_fill_manual(values = c("#786699",'#A55D53','#6280B2','#75A9BF','#A6BA6F')) +
  coord_flip()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.x = element_text(face = "bold",size = 10))  + mytheme 
  #theme(legend.position = "none")
figure_name <- paste0(figure_path,"04_M_L_pie_subcell1_withlegend.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height/2, 
          bg = "transparent")
# graph2pdf(file="./04_Output/Version3/10_M_L_pie_subcell.pdf",aspectr=2, font = "Arial",
#           width = 7, height = 3, bg = "transparent")
## 5-Protein abundance for categories####
data1 <- merge_orig_ref_filter[,c("DL1","DL2","DL3","Matrigel1","Matrigel2","Matrigel3","Category","Subcellular localization")]
data1<- data1[which(data1$`Subcellular localization`=='extracellular'),]
data1$DL <- apply(data.frame(data1$DL1,data1$DL2,data1$DL3), 1, mean,na.rm=T)
data1$Matrigel <- apply(data.frame(data1$Matrigel1,data1$Matrigel2,data1$Matrigel3), 1, mean,na.rm=T)
data2 <- melt(data1,
              id.vars = "Category",
              measure.vars = c("DL","Matrigel"),
              variable.name = "condition",na.rm = T)

data3 <- summarySE(data2, measurevar="value", groupvars=c("condition","Category"))
data3$Category <-factor(data3$Category,levels =rev(c('Collagens','ECM Glycoproteins','Proteoglycans','ECM Regulators','ECM-affiliated Proteins','Secreted Factors')), ordered = F)
####View(data3)
plot05_abundance <- ggplot(data3, aes(x=Category, y=value, fill=condition)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  coord_flip() +
  scale_fill_manual(values = alpha(color_4[c(2,4)],alpha = 0.3)) + 
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(legend.position="right")+ 
  guides(fill=guide_legend("Condition")) + 
  mytheme
figure_name <- paste0(figure_path,"06_Extracellular_protein_abundance.pdf")
graph2pdf(x = plot05_abundance,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")




## 6-Top10 genes in each group ####
# See it in v2_top10_yjj.R / v2_top10_mm.R
## 7-Heatmap figure for DEgenes ####
load("4.3_ECM数据分析项目_输出结果_20210120/output_v2/merge_orig_ref_filter.RData")
choosen_column <- c(1,9:16,31)
merge_orig_ref_choosen <- merge_orig_ref_filter_na[,choosen_column]
merge_na_heatmap <- merge_orig_ref_choosen[,c(2:7,10)]
merge_na_heatmap <- merge_na_heatmap[order(merge_na_heatmap$Category),]
merge_na_heatmap <- merge_na_heatmap[c(1:87,148:156,88:147,157:167),]
merge_na_heatmap$Category <- factor(merge_na_heatmap$Category, levels = c("Collagens","ECM Glycoproteins","Proteoglycans",
                                                                          "ECM Regulators","ECM-affiliated Proteins","Secreted Factors"))

annotation_row = data.frame(
  category = merge_na_heatmap$Category
)
rownames(annotation_row) <- rownames(merge_na_heatmap)
annotation_col = data.frame(
  Condition = factor(c(rep("DL",3),rep("Matrigel",3))))
rownames(annotation_col) <- colnames(merge_na_heatmap)[1:6]
#source("~/Documents/Code_library/Script/scRNAseq/scRNA_Color.R")
make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}
heatmap_col <- c("#05030d","#f9f998")
plot7.1_wholeheatmap <- pheatmap(merge_na_heatmap[,1:6],cluster_cols = F,
         cluster_rows = F,
         annotation_col = annotation_col,
         show_rownames = T,fontsize_row = 5,
         annotation_row = annotation_row,
         color = colorRampPalette(c("white","#67001F"))(15),
         gaps_row = c(23,87,96,134,156),
         cellwidth = 20,cellheight = 5,
         #main = "Heatmap of gene expression",
         na_col = "grey",
         border_color = "black",
         labels_row = make_bold_names(merge_na_heatmap[,1:6],rownames,rownames(merge_na_heatmap))
         )
figure_name <- paste0(figure_path,"7.1_wholeheatmap.pdf")
graph2pdf(x = plot7.1_wholeheatmap,
          file=figure_name,
          font = "Arial",
          width = width*2, 
          height = height*3, 
          bg = "transparent")

plot7.2_collagens <- pheatmap(merge_na_heatmap[1:23,1:6],cluster_cols = F,
         cluster_rows = F,
         show_rownames = T,fontsize_row = 6,
         color = colorRampPalette(c("white","#67001F"))(15),
         cellwidth = 20,cellheight = 6,
         border_color = "black",
         na_col = "grey",
         labels_row = make_bold_names(merge_na_heatmap[1:23,1:6],rownames,rownames(merge_na_heatmap[1:23,])),
         labels_col = make_bold_names(merge_na_heatmap[1:23,1:6],colnames,colnames(merge_na_heatmap[,1:6]))) 
figure_name <- paste0(figure_path,"7.2_pheatmap_collagens.pdf")
graph2pdf(x = plot7.2_collagens,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height,
          bg = "transparent")

#heatmap_anno <- read_excel("4.1_ECM数据分析项目_原始数据_20201231/Fig5核心分类_ECM糖蛋白与蛋白聚糖_生物学过程分类_v2.0_20210521.xlsx")
plot7.3_Glyco <- pheatmap(merge_na_heatmap[24:87,1:6],cluster_cols = F,
         cluster_rows = F,
         #annotation_col = annotation_col,
         show_rownames = T,fontsize_row = 6,
         #fontsize = 10,
         #annotation_row = annotation_row,
         #color = col_heatmap2[4:20],
         color = colorRampPalette(c("white","#67001F"))(15),
         # gaps_row = c(23,86,124,145,154),
         cellwidth = 20,cellheight = 6,
         border_color = "black",
         #main = "Heatmap of Collagens",
         na_col = "grey",
         labels_row = make_bold_names(merge_na_heatmap[24:87,1:6],rownames,rownames(merge_na_heatmap[1:23,])),
         labels_col = make_bold_names(merge_na_heatmap[24:87,1:6],colnames,colnames(merge_na_heatmap[,1:6]))) 
figure_name <- paste0(figure_path,"7.3_pheatmap_Glycoproteins.pdf")
graph2pdf(x = plot7.3_Glyco,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height*2, 
          bg = "transparent")

plot7.4_Proteoglycans <- pheatmap(merge_na_heatmap[88:96,1:6],cluster_cols = F,
         cluster_rows = F,
         #annotation_col = annotation_col,
         show_rownames = T,fontsize_row = 6,
         #fontsize = 10,
         #annotation_row = annotation_row,
         #color = col_heatmap2[4:20],
         color = colorRampPalette(c("white","#67001F"))(15),
         # gaps_row = c(23,86,124,145,154),
         cellwidth = 20,cellheight = 6,
         border_color = "black",
         #main = "Heatmap of Collagens",
         na_col = "grey",
         labels_row = make_bold_names(merge_na_heatmap[88:96,1:6],rownames,rownames(merge_na_heatmap[1:23,])),
         labels_col = make_bold_names(merge_na_heatmap[88:96,1:6],colnames,colnames(merge_na_heatmap[,1:6]))) 
figure_name <- paste0(figure_path,"7.4_pheatmap_Proteoglycans.pdf")
graph2pdf(x = plot7.4_Proteoglycans,
          file=figure_name,
          font = "Arial",
          width = width, 
          height = height/2,
          bg = "transparent")


## 8-Volcano plot for DEgenes ####
source("~/Documents/Code_library/Script/Visualization/volcanoPlot.R")
load("4.3_ECM数据分析项目_输出结果_20210120/output_v2/annotation_filtered.RData")
Orig_dataset_omit.na <- anno_filtered[!is.na(anno_filtered$`DL/Matrigel P value`),]
Orig_dataset_omit.na <- anno_filtered[!is.na(anno_filtered$`DL/Matrigel Ratio`),]
Orig_dataset_omit.na$logFc <- log2(Orig_dataset_omit.na$`DL/Matrigel Ratio`)
data <- Orig_dataset_omit.na
data$logFC <- data$logFc
data$pval <- data$`DL/Matrigel P value`
data$P.Value <- data$pval
data$gene <- data$`Gene name`
data_subset <- data[data$`Gene name` %in% unique(merge_orig_ref_filter$`Gene name`),]
data$adj.P.Val <- data$P.Value
data_subset$adj.P.Val <- data_subset$P.Value
data_subset <- data[data$`Gene name` %in% unique(merge_orig_ref_filter$`Gene name`),]
data_subset <- data_subset %>% subset(-log10(adj.P.Val) > -log10(0.001) & abs(logFC) > 5)
rownames(data_subset) <- data_subset$`Gene name`
data_subset_list <- list()
for (i in 1:6){
  gene_i <- rownames(merge_na_heatmap[which(merge_na_heatmap$Category == levels(merge_na_heatmap$Category)[i]),])
  data_subset_list[[i]] <- data_subset[which(data_subset$`Gene name` %in% gene_i),]
  print(dim(data_subset_list[[i]]))
}
subset_color <- c("#f8ecac","#edc7dd","#f8cbaa","#cad3e5","#e2eec6","#afdac7")
plot08_volcano <- ggplot(data=data, aes(x=logFC, y =-log10(adj.P.Val),color=significant)) +
  geom_point(alpha=0.8, size=0.2,col="black")+
  geom_point(data=subset(data, logFC > 1 & adj.P.Val < 0.05),alpha=0.8, size=0.2,col="red")+
  geom_point(data=subset(data, logFC < -1&adj.P.Val < 0.05),alpha=0.8, size=0.2,col="blue")+
  labs(x="log2 (fold change)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.5)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.5)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylab("-log10(p.val)") + 
  theme(axis.title = element_text(face = "bold")) + 
  geom_point(data=data_subset_list[[1]],alpha=0.85, size=1.5,col=subset_color[1])+
  geom_point(data=data_subset_list[[2]],alpha=0.85, size=1.5,col=subset_color[2])+
  geom_point(data=data_subset_list[[3]],alpha=0.85, size=1.5,col=subset_color[3])+
  geom_point(data=data_subset_list[[4]],alpha=0.85, size=1.5,col=subset_color[4])+
  geom_point(data=data_subset_list[[5]],alpha=0.85, size=1.5,col=subset_color[5])+
  geom_point(data=data_subset_list[[6]],alpha=0.85, size=1.5,col=subset_color[6])+
  geom_text_repel(data=data_subset,aes(label=gene),col="black",alpha = 0.8,size=3) + 
  mytheme
figure_name <- paste0(figure_path,"8_volcano.pdf")
graph2pdf(x = plot08_volcano,
          file=figure_name,
          font = "Arial",
          width = width*2,
          height = height*2, 
          bg = "transparent")
  
## 9-Ribbon plot for DLsignal but notMsignal####
sheet2 <- readxl::read_xlsx("./4.1_ECM数据分析项目_原始数据_20201231/Table2 Biological functions of proteins out of Matrisome in DL-ECM.xlsx",sheet = 2)
sheet2 <- sheet2[!is.na(sheet2$`Gene name`),]
sheet2 <- sheet2[-grep("--",sheet2$`Gene name`),]
sheet2 <- sheet2[,c(1,3)]

regu_cell <- readxl::read_xlsx("4.1_ECM数据分析项目_原始数据_20201231/Table2 Biological functions of proteins out of Matrisome in DL-ECM.xlsx",sheet = 3)
regu_expr <- readxl::read_xlsx("4.1_ECM数据分析项目_原始数据_20201231/Table2 Biological functions of proteins out of Matrisome in DL-ECM.xlsx",sheet = 4)
transfer <- readxl::read_xlsx("4.1_ECM数据分析项目_原始数据_20201231/Table2 Biological functions of proteins out of Matrisome in DL-ECM.xlsx",sheet = 5)
sign_transduct <- readxl::read_xlsx("4.1_ECM数据分析项目_原始数据_20201231/Table2 Biological functions of proteins out of Matrisome in DL-ECM.xlsx",sheet = 6)
metabolism <- readxl::read_xlsx("4.1_ECM数据分析项目_原始数据_20201231/Table2 Biological functions of proteins out of Matrisome in DL-ECM.xlsx",sheet = 7)
biosynthesis <- readxl::read_xlsx("4.1_ECM数据分析项目_原始数据_20201231/Table2 Biological functions of proteins out of Matrisome in DL-ECM.xlsx",sheet = 8)

regu_cell$class <- rep("Cell_regulation",nrow(regu_cell))
regu_expr$class <- rep("Expression_regulation",nrow(regu_expr))
transfer$class <- rep("Transfer",nrow(transfer))
sign_transduct$class <- rep("Signal_transduction",nrow(sign_transduct))
metabolism$class <- rep("Metabolism",nrow(metabolism))
biosynthesis$class <- rep("Biosynthesis",nrow(biosynthesis))

colnames(regu_cell)[1] <- "function"
colnames(regu_expr)[1] <- "function"
colnames(transfer)[1] <- "function"
colnames(sign_transduct)[1] <- "function"
colnames(metabolism)[1] <- "function"
colnames(biosynthesis)[1] <- "function"
All_df <- do.call("rbind", list(regu_cell, regu_expr, transfer,
                                sign_transduct,metabolism,biosynthesis))
All_df <- All_df[,c(5,1,2)]
colnames(All_df) <- c("class","func","Protein_name")
All_df_count <- All_df %>% 
  group_by(class,func) %>% 
  count()
All_df_class <- All_df %>%
  group_by((class)) %>% 
  count()
# nrow(regu_cell) + nrow(regu_expr) +nrow(transfer) + 
# nrow(sign_transduct) + nrow(metabolism) + nrow(biosynthesis)

#create node and edge:
nodes <- data.frame(
  ID = c(All_df_class$`(class)`,All_df_count$func),
  x = c(rep(1,nrow(All_df_class)),rep(2,nrow(All_df_count))),
  col=c(cluster_col[1:6],rep(NA,35)),
  labels = c(All_df_class$`(class)`,All_df_count$func),
  stringsAsFactors = F
)

edges <- All_df_count
colnames(edges) <- c("N1","N2","Value")
nodes <- as.data.frame(nodes)
edges <- as.data.frame(edges)
r <- makeRiver( nodes, edges )
op <- par(cex=0.8)
plot09_riverplot_v1 <- riverplot( r,srt=0,lty=0,textcel=1,
           nodewidth = 1,node_margin = 0.1,
           cex=0.2,plot_area = 0.9)
figure_name <- paste0(figure_path,"9_riverplot_v1.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width2, 
          height = height*2, 
          bg = "transparent")

edges_first <- All_df_class %>% add_column(Gene = rep("Genes_onlyinDL",6),.before=1)
edges_first <- edges_first %>% add_column(ID = 1:6,.after = 3)
colnames(edges_first) <- colnames(edges)
edges_all <- rbind(edges_first,edges)
edges_all$ID <- 1:nrow(edges_all)

nodes_all <- nodes %>% add_row(ID="Genes_onlyinDL",x=1,col="yellow",labels="Genes_onlyinDL",.before=1)
nodes_all$x <- c(1,rep(2,6),rep(3,35))

edges_all <- as.data.frame(edges_all)
nodes_all <- as.data.frame(nodes_all)


r_all <- makeRiver(nodes_all, edges_all)
op <- par(cex=0.8)
plot09.1_riverplot_v2 <- riverplot( r_all,srt=0,lty=0,textcel=1,
           nodewidth = 1,node_margin = 0.1,
           cex=0.2,plot_area = 0.9) + 
  scale_colour_manual(values = my36colors)
figure_name <- paste0(figure_path,"9.1_riverplot_v2.pdf")
graph2pdf(
          file=figure_name,
          font = "Arial",
          width = width2, 
          height = height*2, 
          bg = "transparent")

## different set
library(RColorBrewer)
display.brewer.all()
set3= brewer.pal(n=6,name = "Set3")
nodes_all$col[2:7] <- set3
r_all <- makeRiver(nodes_all, edges_all)
op <- par(cex=0.8)
plot09.1_riverplot_v3 <- riverplot( r_all,srt=0,lty=0,textcel=1,
           nodewidth = 1,node_margin = 0.1,
           cex=0.2,plot_area = 0.9)
figure_name <- paste0(figure_path,"9.1_riverplot_v3.pdf")
graph2pdf(
  file=figure_name,
  font = "Arial",
  width = width2, 
  height = height*2, 
  bg = "transparent")

## 10-GeomSplitViolin in Glycoproteins####

ECM_Glyco_function <- read_excel("./4.1_ECM数据分析项目_原始数据_20201231/Table3 Biological Function of ECM Glycoproteins.xlsx")
ECM_Glyco_function <- ECM_Glyco_function[,c(1,2)]
colnames(ECM_Glyco_function) <- c("feature","genename")
ECM_Glyco_function$feature[grep(pattern = "TGF",ECM_Glyco_function$feature)] <- "TGFβ"
#load("/Users/yujijun/Documents/work/4_ECM数据分析项目_20201120_/4.3_ECM数据分析项目_输出结果_20210120/output_v1/category_genename_deletena.RData")
ECM_Glyco_expre <- merge_na %>% filter(Category == "ECM Glycoproteins")

ECM_Glyco_all <- right_join(x = ECM_Glyco_function,
                            y = ECM_Glyco_expre,
                            by = c("genename" = "Gene name"))
ECM_Glyco_all <- ECM_Glyco_all[!is.na(ECM_Glyco_all$feature),1:8]

#### wigth to long
ECM_Glyco_all_long <- ECM_Glyco_all %>%  
  pivot_longer(
    cols = DL1:Matrigel3,
    names_to = "sample",
    values_to = "expression"
  )
ECM_Glyco_all_long$group <- "Matrigel"
ECM_Glyco_all_long$group[grep(pattern = "DL",x=ECM_Glyco_all_long$sample)] <- "DL"
ECM_Glyco_all_long <- ECM_Glyco_all_long[!is.na(ECM_Glyco_all_long$expression),]
ECM_Glyco_all_long$feature <- str_to_title(ECM_Glyco_all_long$feature)
ECM_Glyco_all_long_summary <- summarySE(ECM_Glyco_all_long, 
                                        measurevar="expression", 
                                        groupvars=c("feature","genename","group"))

source("/Users/yujijun/Documents/Code_library/Script/Visualization/geom_split_violin.R")
ggboxplot(ECM_Glyco_all_long, x = "feature", y = "expression", fill = "group",
          palette = c("#00AFBB", "#E7B800"))
# 自行调整下面的参数

plot10_split_violin <- ggplot(ECM_Glyco_all_long,aes(x= feature,y= expression,fill= group))+
  geom_split_violin(trim= F,color="white",scale = "area") + #绘制分半的小提琴图
  geom_point(data = ECM_Glyco_all_long_summary,aes(x= feature, y= expression),pch=19,
             position=position_dodge(0.5),size= 1)+ #绘制均值为点图
  # geom_errorbar(data = ECM_Glyco_all_long_summary,aes(ymin = expression-ci, ymax= expression+ci), 
  #               width= 0.05, 
  #               position= position_dodge(0.5), 
  #               color="black",
  #               alpha = 0.8,
  #               size= 0.5) +
  scale_fill_manual(values = color_4[c(2,4)])+ 
  #scale_fill_manual(values = col_convert) + 
  labs(y=("gene expression"),x=NULL,title = "ECM Glygoproteins") + 
  theme_bw() + 
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "anova",
                     label.y = max(ECM_Glyco_all_long$expression),
                     hide.ns = T) + mytheme + 
  theme(axis.text.x = element_text(angle =270,face = "bold",hjust = 0)) +
  theme(legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 12,face = "bold")) + 
  xlab("Molecular Function") + 
  ylab("Gene Expression") + 
  theme(plot.title = element_blank())
  
  #mytheme
figure_name <- paste0(figure_path,"10_split_violin.pdf")
graph2pdf(x = plot10_split_violin,
          file=figure_name,
          font = "Arial",
          width = width2*2/3, 
          height = height*2, 
          bg = "transparent")
output_name <- paste0(output_path,"01_ECM_Glyco_all_long_20210518_v1.RData")
save(ECM_Glyco_all_long,file = output_name)
output_name <- paste0(output_path,"01_ECM_Glyco_all_long_summary_20210518_v1.RData")
save(ECM_Glyco_all_long_summary,file = output_name)

## 11-Bubble Plot####
load("/Users/yujijun/Documents/work/4_ECM数据分析项目_20201120_/4.3_ECM数据分析项目_输出结果_20210120/output_v2/merge_orig_ref_filter.RData")
choosen_column <- c(1,9:17,31)
merge_orig_ref_choosen <- merge_orig_ref_filter_na[,choosen_column]
merge_na_heatmap <- merge_orig_ref_choosen[,c(1:7,10,11)]
merge_na <- merge_orig_ref_choosen[!(rowSums(is.na(merge_orig_ref_choosen)) > 7),]

rownames(merge_na_heatmap) <- merge_na$`Gene name`
merge_na_heatmap <- merge_na_heatmap[order(merge_na_heatmap$Category),]
merge_na_heatmap <- merge_na_heatmap[c(1:87,148:156,88:147,157:167),]
#merge_na_heatmap$Category <- factor(merge_na_heatmap$Category, levels = c("Collagens","ECM Glycoproteins","Proteoglycans",

merge_na_hehe <- merge_na_heatmap[1:23,]
merge_na_hehe <-  separate(merge_na_hehe,"Gene name", c("Gene name", "num"), "a")
merge_na_hehe <- merge_na_hehe[,c(1,2:8)]
merge_na_hehe_long <- merge_na_hehe %>%  pivot_longer(cols = DL1:Matrigel3,names_to = "sample",values_to = "expression")
merge_na_hehe_long$num <- NULL
merge_na_hehe_long <- merge_na_hehe_long[!is.na(merge_na_hehe_long$expression),]
colnames(merge_na_hehe_long) <- c("matrigeltype","sample","expression")
merge_na_hehe_long_summary <- merge_na_hehe_long %>%  summarySE(measurevar = "expression",
                                                                groupvars = c("matrigeltype","sample"))

my36colors <-c('#53A85F', '#F1BB72', '#F3B1A0', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175','#E5D2DD', '#D6E7A3')
# ggplot(merge_na_hehe_long_summary,aes(x=sample,y=expression)) + 
#   geom_point(aes(size=N,colour=matrigeltype)) + 
#   scale_colour_manual(values=my36colors)+
#   scale_size(rang=c(3,10))
#graph2ppt(file="/03 discussion/Bubble_Plotformatrigel.pptx")

plot11_bubble <- ggplot(merge_na_hehe_long_summary,aes(x=sample,y=expression)) + 
  geom_point(aes(size=N,colour=matrigeltype)) + 
  scale_colour_manual(values=my36colors)+
  scale_size(rang=c(3,10))+
  xlab("") +
  ylab("protein expression")+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.x = element_text(colour="grey20",size=24,angle=270,face="bold"),
        axis.text.y = element_text(colour="grey20",size=24,face="bold",hjust = 2),
        axis.title = element_text(face = "bold",size = 24)) 
#graph2ppt(file="/03 discussion/Bubble_Plotformatrigel1.pptx")
figure_name <- paste0(figure_path,"11_bubble_plot.pdf")
graph2pdf(x = plot11_bubble,
          file=figure_name,
          font = "Arial",
          width = width2, 
          height = width2, 
          bg = "transparent")


### 12-String Enrichment result visualization ####
library(ggpubr)
ECM.collagen.KEGG <- read.delim(file = "./4.3_ECM数据分析项目_输出结果_20210120/output_last/ECM_enrichment.KEGG_20210527.tsv")
ECM.collagen.KEGG <- ECM.collagen.KEGG[c(5,6,11),] 
ECM.collagen.KEGG$group <- "KEGG"
ECM.collagen.GO.BP <- read.delim(file = "./4.3_ECM数据分析项目_输出结果_20210120/output_last/ECM_Collagens_enrichment.Process_20210527.tsv")
ECM.collagen.GO.BP <- ECM.collagen.GO.BP[c(11,12,24,32),]
ECM.collagen.GO.BP$group <- "GO.BP"
ECM.collagen.GO.MF <- read.delim(file = "./4.3_ECM数据分析项目_输出结果_20210120/output_last/ECM_enrichment.Function_20210527.tsv")
ECM.collagen.GO.MF <- ECM.collagen.GO.MF[2,]
ECM.collagen.GO.MF$group <- "GO.MF"
ECM.collagen.Uniprot <- read.delim(file = "./4.3_ECM数据分析项目_输出结果_20210120/output_last/ECM_enrichment.Keyword_Uniprot_20210527.tsv")
ECM.collagen.Uniprot <- ECM.collagen.Uniprot[10,]
ECM.collagen.Uniprot$group <- "Uniprot"
ECM.collagen.enrichall <- bind_rows(ECM.collagen.KEGG,
                                    ECM.collagen.GO.BP,
                                    ECM.collagen.GO.MF,
                                    ECM.collagen.Uniprot)
colnames(ECM.collagen.enrichall)[6] <- "FDR"
plot12_String_enrichresult <- 
  ggbarplot(ECM.collagen.enrichall,
          x = "term.description",
          y = "observed.gene.count",
          fill = "FDR",
          color = "group",
          rotate = TRUE,
          #palette = "jco",#杂志jco的配色 
          sort.val = "asc",#上升排序,区别于desc，具体看图演示 
          x.text.angle=90,
          size = 1) + 
  scale_color_manual(values = my36colors[c(1,2,6,4)]) + 
  mytheme + 
  theme(axis.text = element_text(size = 6)) + 
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'))  + 
  theme(axis.text = element_text(size=10)) + 
  xlab("Function Description") + 
  ylab("Obversed Protein Number") +
  theme(legend.text = element_text(size = 8,face = "bold"),
        legend.title = element_text(size = 10,face = "bold")) 
  
figure_name <- paste0(figure_path,"12_String_enrichresult.pdf") 
graph2ppt(x = plot12_String_enrichresult,
          file=paste0(figure_path,"12_String_enrichresult.pptx"))
graph2pdf(x = plot12_String_enrichresult,
          file=figure_name,
          font = "Arial",
          width = width2, 
          height = height, 
          bg = "transparent")

### 13-RNAseq Enrichment result visualization #### 
library(readr)
library(ggplot2)
library(forcats)
library(enrichplot)
library(clusterProfiler)
library(clusterProfiler.dplyr)
load("./4.3_ECM数据分析项目_输出结果_20210120/output_last/RNAseq_allEnrich_20210528.RData")
base_path <- "/Users/yujijun/Documents/work/4_ECM数据分析项目_20201120_/4.1_ECM数据分析项目_原始数据_20201231/测序数据_version3"
setwd(base_path)
setwd('./1.deglist/DL_d5vsM_d5/')
DLvsM_d5 <- read.delim("DL_d5vsM_d5_deg.xls")
DLvsM_d5 <- DLvsM_d5 %>% filter(padj < 0.05,log2FoldChange > 1)
gene.df <- bitr(DLvsM_d5$gene_id, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
GO.MF <- enrichGO(gene.df$ENTREZID,ont = "MF",OrgDb = org.Hs.eg.db)
GO.CC <- enrichGO(gene.df$ENTREZID,ont = "CC",OrgDb = org.Hs.eg.db)
GO.BP <- enrichGO(gene.df$ENTREZID,ont = "BP",OrgDb = org.Hs.eg.db)
KEGG <- enrichKEGG(gene.df$ENTREZID,organism = "hsa")
save(GO.MF,GO.CC,GO.BP,KEGG, file = "4.3_ECM数据分析项目_输出结果_20210120/output_last/RNAseq_allEnrich_20210528.RData")
# GO.MF.test <- filter(GO.MF, p.adjust < .05, qvalue < 0.2)
# GO.MF.result <- GO.MF.test@result
GO.BP.test <- filter(GO.BP, p.adjust < .05, qvalue < 0.2)
GO.BP.result <- GO.BP.test@result
GO.BP.result <- GO.BP.result[grep(pattern = "vascu|tyrosine kinase activity|angiogenesis|collagen|MAP kinase activity|epidermal growth factor receptor signaling pathway|fibroblast growth factor receptor signaling pathway|regulation of ERK1 and ERK2 cascade|epithelial cell proliferation involved in liver morphogenesis|hepat|liver",GO.BP.result$Description),]
#GO.BP.result <- GO.BP.result[grep("negetive")]
rownames(GO.BP.result) <- 1:nrow(GO.BP.result)
GO.BP.result <- GO.BP.result[-c(1,2,3,5,14,30,32),]
GO.BP.test@result <- GO.BP.result
GO.BP.tyrosine <- GO.BP.result[grep("tyrosine",GO.BP.result$Description),] 

#write_tsv(KEGG.result,file="./4.3_ECM数据分析项目_输出结果_20210120/output_last/KEGG.test.20210528.tsv")
#write_tsv(KEGG.result,file = "4.3_ECM数据分析项目_输出结果_20210120/output_last/RNAseq_KEGG.result.tsv")
y <- mutate(GO.BP.test, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
Plot13_GO.BP <- ggplot(y, showCategory = 40, 
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("Rich Factor") +
  ylab(NULL) + 
  #ggtitle("Enriched GO Biological Process") + 
  mytheme + 
  theme(axis.text = element_text(size = 8)) + 
  theme(legend.position = c(0.8,0.4)) + 
  theme(axis.text = element_text(face = "bold",size = 8)) +
  ggtitle("Biological Process in GO") + 
  theme(plot.title = element_text(face = "bold",size = 12)) + 
  theme(legend.text = element_text(face = "bold",size = 8),
        legend.title = element_text(face = "bold",size = 10)) 

figure_name <- paste0(figure_path,"13_GO.BP.pdf")
graph2pdf(x = Plot13_GO.BP,
          file=figure_name,
          font = "Arial",
          width = width2*1.2, 
          height = height*2, 
          bg = "transparent")

KEGG.test <- filter(KEGG, p.adjust < .05, qvalue < 0.2)
KEGG.result <- KEGG.test@result
deleted_row <- c(5,9,11,12,13,16,17,19,20,24,26:36,40,43,44) -1 
KEGG.result <- KEGG.result[-deleted_row,]
KEGG.test@result <- KEGG.result
y <- mutate(KEGG.test, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
Plot13_KEGG <- ggplot(y, showCategory = 40, 
                       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("Rich Factor") +
  ylab(NULL) + 
  #ggtitle("Enriched GO Biological Process") + 
  mytheme + 
  theme(axis.text = element_text(size = 8)) + 
  theme(legend.position = c(0.8,0.4)) + 
  theme(axis.text = element_text(face = "bold",size = 8)) +
  ggtitle("Enrichment of KEGG") + 
  theme(plot.title = element_text(face = "bold",size = 12)) + 
  theme(legend.text = element_text(face = "bold",size = 8),
        legend.title = element_text(face = "bold",size = 10)) 
figure_name <- paste0(figure_path,"13_KEGG.pdf")
graph2pdf(x = Plot13_KEGG,
          file=figure_name,
          font = "Arial",
          width = width2, 
          height = height*1.5, 
          bg = "transparent")
