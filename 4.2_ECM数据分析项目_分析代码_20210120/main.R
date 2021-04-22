# Updata Time:2021.01.20 4:44pm
# Update Time:20210521
# Yujijun

#### load all package ####
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

col_convert <- c("#B2182B","#2166AC")
#### Data cleaning ####
data <- read_excel("./01_Original_data/Annotation_combine.xlsx")
load("./01_Original_data/merge_orig_ref_filter.RData")
selectresult=subset(data,!(is.na(DL1)&is.na(DL2)&is.na(DL3)&is.na(Matrigel1)&is.na(Matrigel2)&is.na(Matrigel3)))
datafilter <- selectresult[!(selectresult$`Gene name` == "--" | is.na(selectresult$`Gene name`)), ]
test1<- count_dt(datafilter,`Gene name`)%>% filter_dt(n>1)
datadup= datafilter[which(datafilter$`Gene name` %in% test1$`Gene name`),]
p<- c('P02770','P09605','D3ZQ25','P04937','A0A0G2KAM4','A0A096MJI4','A0A0G2KAY3','D3ZN05','A0A0G2K484','A0A0G2K8C1','M0RBV9','A0A0G2JVG3','A0A0G2K781','Q5RKJ9','P16086','Q9EQT5','A0A0G2JSQ4','A0A140TAF0','Q9JJP9')
anno_filtered = datafilter[which(!datafilter$`Protein accession` %in% p),]
write_excel_csv(anno_filtered,"4.3_ECM数据分析项目_输出结果_20210120/output_v2/annotation_filtered.csv")
save(anno_filtered,file = "4.3_ECM数据分析项目_输出结果_20210120/output_v2/annotation_filtered.RData")

merge_orig_ref1<- count_dt(merge_orig_ref,`Gene name`)%>% filter_dt(n>1)
merge_orig_ref2 <- merge_orig_ref[which(merge_orig_ref$`Gene name` %in% merge_orig_ref1$`Gene name`),]
merge_orig_ref_filter <- merge_orig_ref[-c(67,76,96,101,134,169),]
choosen_column <- c(1,9:16,31)
merge_orig_ref_choosen <- merge_orig_ref_filter[,choosen_column]
merge_na <- merge_orig_ref_choosen[!(rowSums(is.na(merge_orig_ref_choosen)) > 7),]
rownames(merge_orig_ref_filter) <- merge_orig_ref_filter$`Gene name`
merge_orig_ref_filter_na <- merge_orig_ref_filter[merge_na$`Gene name`,]
write_excel_csv(merge_orig_ref_filter_na,"4.3_ECM数据分析项目_输出结果_20210120/output_v2/merge_orig_ref_filter.csv")
save(merge_orig_ref_filter_na,file = "4.3_ECM数据分析项目_输出结果_20210120/output_v2/merge_orig_ref_filter.RData")

#### 01-Barplot for all genes####
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
DL_table$expr <- rep("DL",2)
Matrigel_table$expr <- rep("Matrigel",2)
All_table <- rbind(DL_table,Matrigel_table)
color_1 <- c("#B0C4DE","#4682B4","#CDBE70","#8B814C")
color_2 <- c("#E9967A","#FA8072","#CDBE70","#8B814C")
color_3 <- c("#C1CDC1","#838B83","#CDBE70","#8B814C")
color_4 <- c("#CDC1C5","#8B8386","#CDBE70","#8B814C")
color_5 <- c("#CDC1C5","#8B8386","#E9967A","#FA8072")
All_table$color <- paste0(All_table$expr,"_",All_table$Position)
## visulizaiton 
ggplot(data = All_table,mapping = aes(x=expr,y=Freq,fill=color)) + 
  geom_bar( stat = "identity",width = 0.5) + 
  ylab("Average protein number") + 
  scale_fill_manual(values = color_1)+  
  geom_text(aes(label=Freq), position = position_stack(0.5), color = "white") + 
  theme(text = element_text(face = "bold")) + 
  guides(fill=guide_legend("Position")) + 
  theme(text = element_text(face = "bold")) +  
  theme(axis.text = element_text(face = "bold",size = 10)) + 
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank())+
  theme_classic()
graph2ppt(file="01_barplot_for_all_gene.pptx")
graph2pdf(file="./04_Output/Version2/01_barplot_for_all_gene.pdf",aspectr=2, font = "Arial",
          width = 7, height = 5, bg = "transparent")

#### 02-Vennplot about IN and EX####
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

venn.plot <- venn.diagram(list(DL_IN=DL_IN,Matrigel_IN=Matrigel_IN),filename = NULL, height = 450, width = 450,resolution =300, 
                          imagetype="png", col="transparent",fill=c("#4682B4","#8B814C"),alpha = 0.7, cex=0.9, cat.cex=0.9,force.unique = F)

grid.draw(venn.plot)
graph2ppt(file="02_DL_Venn_1.pptx")
graph2pdf(file="./04_Output/Version2/02_IN_Venn.pdf",aspectr=2, font = "Arial",
          width = 5, height = 5, bg = "transparent")


venn.plot <- venn.diagram(list(DL_EX=DL_EX,Matrigel_EX=Matrigel_EX), filename = NULL, height = 450, width = 450,resolution =300, 
                          imagetype="png", col="transparent",fill=c("#B0C4DE","#CDBE70"),
                          alpha = 0.80, cex=0.9, cat.cex=0.9,force.unique = F)
grid.draw(venn.plot)
graph2ppt(file="02_M_Venn_1.pptx")
graph2pdf(file="04_Output/Version2/02_EX_Venn.pdf",aspectr=2, font = "Arial",
          width = 5, height = 5, bg = "transparent")

#### 03-Enrichment figure for all genes#### 
# statistic the DL emrichment
data1 <- anno_filtered[,c("DL1","DL2","DL3","Matrigel1","Matrigel2","Matrigel3","Cellular Component")]
names(data1)[names(data1) == "Cellular Component"] <- "Cellular.Component"
View(data1)
m<-c("GO:0005576","GO:0005887","GO:0005856","GO:0005739","GO:0005737","GO:0043229","GO:0005777","GO:0042579","GO:0005768","GO:0005840","GO:0009986","GO:0005694","GO:0005783","GO:0032991","GO:0005634")
a<-data1[0,]
f<-function(x) sum(!is.na(x))
for (gene in m){
  test<- apply(data1[grep(gene,data1$`Cellular.Component`),],2,f)
  test$`Cellular.Component` <- gene
  a<-rbind(a,test)
}
View(a)
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

library(reshape2)
data2 <- melt(a,
              id.vars = c("CC","DL_sd","Matrigel_sd"),
              measure.vars = c("DL","Matrigel"),
              variable.name = "condition",na.rm = T)
data2$se <- ifelse(data2$condition=="DL",data2$DL_sd,data2$Matrigel_sd)

library(ggplot2)
ggplot(data2, aes(x=`CC`, y=value, fill=condition)) + 
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
  theme(legend.position="right") + 
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 15)) + 
  theme(text = element_text(face = "bold")) + 
  scale_fill_manual(values = color_1[c(2,4)])
graph2ppt(file="05_Cellular.Component1.pptx")
graph2pdf(file="04_Output/Version2/05_Cellular.Component.pdf",aspectr=2, font = "Arial",
          width = 7, height = 10, bg = "transparent")

##
data1 <- anno_filtered[,c("DL1","DL2","DL3","Matrigel1","Matrigel2","Matrigel3","Biological Process")]
names(data1)[names(data1) == "Biological Process"] <- "Biological.Process"
View(data1)
m<-c("GO:0065007","GO:0009987","GO:0032502","GO:0008152","GO:0002376","GO:0000003","GO:0051179")
a<-data1[0,]
f<-function(x) sum(!is.na(x))
for (gene in m){
  test<- apply(data1[grep(gene,data1$`Biological.Process`),],2,f)
  test$`Biological.Process` <- gene
  a<-rbind(a,test)
}
View(a)

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


library(reshape2)
data2 <- melt(a,
              id.vars = c("BP","DL_sd","Matrigel_sd"),
              measure.vars = c("DL","Matrigel"),
              variable.name = "condition",na.rm = T)
View(data2)
data2$se <- ifelse(data2$condition=="DL",data2$DL_sd,data2$Matrigel_sd)

library(ggplot2)
ggplot(data2, aes(x=`BP`, y=value, fill=condition)) + 
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
  theme(legend.position="right") + 
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 15)) + 
  theme(text = element_text(face = "bold")) + 
  scale_fill_manual(values = color_1[c(2,4)])
#graph2ppt(file="05_Biological.Process1.pptx")
graph2pdf(file="04_Output/Version2/05_Biological.Process.pdf",aspectr=2, font = "Arial",
          width = 7, height = 5, bg = "transparent")


##
data1 <- anno_filtered[,c("DL1","DL2","DL3","Matrigel1","Matrigel2","Matrigel3","Molecular Function")]
names(data1)[names(data1) == "Molecular Function"] <- "Molecular.Function"
View(data1)
m<-c("GO:0003824","GO:0005198","GO:0060089","GO:0016209","GO:0005488")
a<-data1[0,]
f<-function(x) sum(!is.na(x))
for (gene in m){
  test<- apply(data1[grep(gene,data1$`Molecular.Function`),],2,f)
  test$`Molecular.Function` <- gene
  a<-rbind(a,test)
}
View(a)

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



library(reshape2)
data2 <- melt(a,
              id.vars = c("MF","DL_sd","Matrigel_sd"),
              measure.vars = c("DL","Matrigel"),
              variable.name = "condition",na.rm = T)
View(data2)
data2$se <- ifelse(data2$condition=="DL",data2$DL_sd,data2$Matrigel_sd)

library(ggplot2)
ggplot(data2, aes(x=`MF`, y=value, fill=condition)) + 
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
  theme(legend.position="right") + 
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 15)) + 
  theme(text = element_text(face = "bold")) + 
  scale_fill_manual(values = color_1[c(2,4)]) 
#graph2ppt(file="05_Molecular.Function1.pptx")
graph2pdf(file="04_Output/Version2/05_Molecular.Function.pdf",aspectr=2, font = "Arial",
          width = 7, height = 4, bg = "transparent")


#### 04-Pie plot for categories####
DL_mat <- merge_orig_ref_filter[,c("DL1","DL2","DL3","Category")]
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
graph2ppt(file="03_D_pie.pptx")
graph2pdf(file="03_D_pie.pdf",aspectr=2, font = "Arial",
          width = 7, height = 5, bg = "transparent")


ML_mat <- merge_orig_ref_filter[,c("Matrigel1","Matrigel2","Matrigel3","Category")]
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
graph2ppt(file="03_M_pie.pptx.pptx")
graph2pdf(file="03_M_pie.pdf",aspectr=2, font = "Arial",
          width = 7, height = 5, bg = "transparent")

#### 05-Protein abundance for categories####
data1 <- merge_orig_ref_filter[,c("DL1","DL2","DL3","Matrigel1","Matrigel2","Matrigel3","Category","Subcellular localization")]
data1<- data1[which(data1$`Subcellular localization`=='extracellular'),]

#View(data1)
data1$DL <- apply(data.frame(data1$DL1,data1$DL2,data1$DL3), 1, mean,na.rm=T)
data1$Matrigel <- apply(data.frame(data1$Matrigel1,data1$Matrigel2,data1$Matrigel3), 1, mean,na.rm=T)

library(reshape2)
data2 <- melt(data1,
              id.vars = "Category",
              measure.vars = c("DL","Matrigel"),
              variable.name = "condition",na.rm = T)
#View(data2)
#install.packages("Rmisc")
library(Rmisc)

data3 <- summarySE(data2, measurevar="value", groupvars=c("condition","Category"))
View(data3)
library(ggplot2)
ggplot(data3, aes(x=Category, y=value, fill=condition)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  coord_flip() +
  #xlab("") +
  #ylab("") + 
  #ggtitle("extracellular protein abundance") +
  scale_fill_manual(values = color_1[c(2,4)]) + 
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(legend.position="right")+ 
  guides(fill=guide_legend("Condition")) + 

#graph2ppt(file="04_extracellular protein abundance.pptx")
graph2pdf(file="04_Output/Version2/04_Extracellular protein abundance.pdf",aspectr=2, font = "Arial",
          width = 7, height = 5, bg = "transparent")


#### 06-Top10 expression of six category ####
library(ggplot2)
library(export)
data1 <- merge_orig_ref_filter[,c("DL1","DL2","DL3","Matrigel1","Matrigel2","Matrigel3","Category","Gene name")]
data1$DL <- apply(data.frame(data1$DL1,data1$DL2,data1$DL3), 1, mean,na.rm=T)
data1$Matrigel <- apply(data.frame(data1$Matrigel1,data1$Matrigel2,data1$Matrigel3), 1, mean,na.rm=T)

data2 = data1[which(data1$Category =='ECM Regulators'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:10),]
data2$DL <- round(data2$DL,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(colour="#668B8B", fill = '#96CDCD', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = DL,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$DL)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2ppt(file="06_DL_ECM Regulators2.pptx")
graph2pdf(file="06_DL_ECM Regulators.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

datann <- data1[which(data1$`Gene name`%in% c('Tgm2','Loxl1','Adamtsl4','F13b','Ambp','Serpinb6a','Hrg','Plg','Htra1','Lox')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#668B8B",fill = '#d7f5f5', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

graph2pdf(file="06_M_ECM Regulators_samegene.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

data2 = data1[which(data1$Category =='ECM Regulators'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#668B8B",fill = '#BBFFFF', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

graph2ppt(file="06_M_ECM Regulators1.pptx")
graph2pdf(file="06_M_ECM Regulators.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")


data2 = data1[which(data1$Category =='Secreted Factors'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:8),]
data2$DL <- round(data2$DL,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(colour="#698B22", fill = '#9ACD32', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = DL,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$DL)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2ppt(file="06_DL_Secreted Factors1.pptx")
graph2pdf(file="06_DL_Secreted Factors.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

datann <- data1[which(data1$`Gene name`%in% c('Megf6','Egfl7','S100a10','Pf4','Angptl6','Inhbc','Inhbe','Ccl9')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2 <- data2[c(1:3),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#698B22",fill = '#ddfc9d', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2pdf(file="06_M_Secreted Factors_samegene.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

data2 = data1[which(data1$Category =='Secreted Factors'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:6),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#698B22",fill = '#C0FF3E', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2ppt(file="06_M_Secreted Factors1.pptx")
graph2pdf(file="06_M_Secreted Factors.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")


data2 = data1[which(data1$Category =='ECM-affiliated Proteins'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:10),]
data2$DL <- round(data2$DL,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(colour="#1c1c30", fill = "#292985", stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = DL,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$DL)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2ppt(file="06_DL_ECM-affiliated Proteins1.pptx")
graph2pdf(file="06_DL_ECM-affiliated Proteins.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")



datann <- data1[which(data1$`Gene name`%in% c('Frem2','Anxa11','Anxa2','Anxa6','Frem1','Grem1','C1qtnf6','Anxa5','Lgals9','Anxa1')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2 <- data2[c(1:9),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#1818a3",fill = '#81b0f7', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2pdf(file="06_M_ECM-affiliated Proteins_samegene.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")



data2 = data1[which(data1$Category =='ECM-affiliated Proteins'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#1818a3",fill = '#0066FF', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2ppt(file="06_M_ECM-affiliated Proteins1.pptx")
graph2pdf(file="06_M_ECM-affiliated Proteins.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

data2 = data1[which(data1$Category =='Proteoglycans'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:9),]
data2$DL <- round(data2$DL,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(colour="#153e4a", fill = '#007799', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = DL,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$DL)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2ppt(file="06_DL_Proteoglycans1.pptx")
graph2pdf(file="06_DL_Proteoglycans.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

datann <- data1[which(data1$`Gene name`%in% c('Lum','Prg2','Prelp','Vcan','Prg3','Bgn','Prg4','Dcn','Aspn')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2 <- data2[c(1:8),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#0f8cba",fill = '#00BBFF', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.colour = 'black', size=0.5, linetype='solid')
graph2pdf(file="06_M_Proteoglycans_samegene.pdfy = element_line(",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")


data2 = data1[which(data1$Category =='Proteoglycans'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:8),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#0f8cba",fill = '#00BBFF', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2ppt(file="06_M_Proteoglycans1.pptx")
graph2pdf(file="06_M_Proteoglycans.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")


data2 = data1[which(data1$Category =='ECM Glycoproteins'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:10),]
data2$DL <- round(data2$DL,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(colour="#8B8B00", fill = '#CDCD00', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = DL,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$DL)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2ppt(file="06_DL_ECM Glycoproteins1.pptx")
graph2pdf(file="06_DL_ECM Glycoproteins.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

datann <- data1[which(data1$`Gene name`%in% c('Papln','Emilin1','Fbn1','Dpt','Postn','Ltbp4','Tgfbi','Lama3','Ecm1','Tnc')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#8B8B00",fill = '#FFFF00', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2pdf(file="06_M_ECM Glycoproteins_samegene.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

data2 = data1[which(data1$Category =='ECM Glycoproteins'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#8B8B00",fill = '#FFFF00', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2ppt(file="06_M_ECM Glycoproteins1.pptx")
graph2pdf(file="06_M_ECM Glycoproteins.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

data2 = data1[which(data1$Category =='Collagens'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:10),]
data2$DL <- round(data2$DL,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(colour="#8B636C", fill = '#CD919E', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = DL,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$DL)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2ppt(file="06_DL_Collagens1.pptx")
graph2pdf(file="06_DL_Collagens.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

datann <- data1[which(data1$`Gene name`%in% c('Col2a1','Col1a1','Col6a1','Col3a1','Col5a2','Col14a1','Col6a2','Col6a5','Col5a3','Col4a2')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#8B636C",fill = '#FFB5C5', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2pdf(file="06_M_Collagens_samegene.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")


data2 = data1[which(data1$Category =='Collagens'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(colour="#8B636C",fill = '#FFB5C5', stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2ppt(file="06_M_Collagens1.pptx")
graph2pdf(file="06_M_Collagens.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")


#### 07.1-Heatmap figure for DEgenes ####
load("4.3_ECM数据分析项目_输出结果_20210120/output_v2/merge_orig_ref_filter.RData")
choosen_column <- c(1,9:16,31)
merge_orig_ref_choosen <- merge_orig_ref_filter_na[,choosen_column]
merge_na_heatmap <- merge_orig_ref_choosen[,c(2:7,10)]
rownames(merge_na_heatmap) <- merge_na$`Gene name`
merge_na_heatmap <- merge_na_heatmap[order(merge_na_heatmap$Category),]
merge_na_heatmap <- merge_na_heatmap[c(1:87,148:156,88:147,157:167),]
merge_na_heatmap$Category <- factor(merge_na_heatmap$Category, levels = c("Collagens","ECM Glycoproteins","Proteoglycans",
                                                                          "ECM Regulators","ECM-affiliated Proteins","Secreted Factors"))
library(pheatmap)
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
pheatmap(merge_na_heatmap[,1:6],cluster_cols = F,
         cluster_rows = F,
         annotation_col = annotation_col,
         show_rownames = F,fontsize_row = 5,
         annotation_row = annotation_row,
         color = colorRampPalette(c("white","#67001F"))(15),
         gaps_row = c(23,87,96,134,156),
         cellwidth = 20,cellheight = 3,
         #main = "Heatmap of gene expression",
         na_col = "grey",
         border_color = "black",
         labels_row = make_bold_names(merge_na_heatmap[,1:6],rownames,rownames(merge_na_heatmap))
         )
#graph2ppt(file="04_Output_figure/version3/07_Heatmap2.pptx")
graph2pdf(file="04_Output_figure/version3/08_pheatmap_all_v1.pdf",aspectr=2, font = "Arial",
          width = 5, height = 12, bg = "transparent")

#### 07.2-Part of heatmap for DEgenes####
pheatmap(merge_na_heatmap[1:23,1:6],cluster_cols = F,
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
         labels_row = make_bold_names(merge_na_heatmap[1:23,1:6],rownames,rownames(merge_na_heatmap[1:23,])),
         labels_col = make_bold_names(merge_na_heatmap[1:23,1:6],colnames,colnames(merge_na_heatmap[,1:6]))) 
graph2pdf(file="04_Output_figure/Version3/08_pheatmap_collagens_v1.pdf",aspectr=2, font = "Arial",
          width = 3, height = 4, bg = "transparent")

heatmap_anno <- read_excel("4.1_ECM数据分析项目_原始数据_20201231/Fig5核心分类_ECM糖蛋白与蛋白聚糖_生物学过程分类_v2.0_20210521.xlsx")
pheatmap(merge_na_heatmap[24:87,1:6],cluster_cols = F,
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
graph2pdf(file="04_Output_figure/Version3/08_pheatmap_ECM_Glycoproteins_v1.pdf",aspectr=2, font = "Arial",
          width = 3, height = 7, bg = "transparent")

pheatmap(merge_na_heatmap[88:96,1:6],cluster_cols = F,
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
graph2pdf(file="04_Output_figure/Version3/08_pheatmap_Proteoglycans.pdf",aspectr=2, font = "Arial",
          width = 3, height = 2, bg = "transparent")


#### 08-Volcano plot for DEgenes ####
library(ggplot2)
library(ggrepel)
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

# ggplot(data=data, aes(x=logFC, y =-log10(pval))) +
#   ## 三个部分分别画点
#   geom_point(data=subset(data,abs(data$logFC) <= 1),aes(size=abs(logFC)),color="black",alpha=0.1) +
#   geom_point(data=subset(data,data$P.Value<0.05 & data$logFC > 1),aes(size=abs(logFC)),color="red",alpha=0.2) +
#   geom_point(data=subset(data,data$P.Value<0.05 & data$logFC < -1),aes(size=abs(logFC)),color="green",alpha=0.2) +
#   ## 画线
#   geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
#   geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
#   ## 主题
#   theme_bw()+
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"))+
#   labs(x="log2 (fold change)",y="-log10 (q-value)")+
#   theme(plot.title = element_text(hjust = 0.5))+
#   theme(legend.position='none')+
#   ## 标签
#   #geom_text_repel(data=subset(data, abs(logFC) > 6), aes(label=gene),col="black",alpha = 0.8)
#   geom_text_repel(data=data_subset, aes(label=gene),col="black",alpha = 0.8)

#换一种风格
library(ggplot2)
library(ggrepel)
# data <- allDiff
# data$gene <- data$gene_name
# data$significant <- as.factor(data$P.Value<0.05 & abs(data$logFC) > 0.5)
# data$adj.P.Val <- data$pval
# data$gene <- rownames(data)
data$adj.P.Val <- data$P.Value
data_subset$adj.P.Val <- data_subset$P.Value
data_subset <- data[data$`Gene name` %in% unique(merge_orig_ref_filter$`Gene name`),]
data_subset <- data_subset %>% subset(-log10(adj.P.Val) > -log10(0.01) & abs(logFC) > 1)
rownames(data_subset) <- data_subset$`Gene name`
data_subset_list <- list()
for (i in 1:6){
  gene_i <- rownames(merge_na_heatmap[which(merge_na_heatmap$Category == levels(merge_na_heatmap$Category)[i]),])
  data_subset_list[[i]] <- data_subset[which(data_subset$`Gene name` %in% gene_i),]
  print(dim(data_subset_list[[i]]))
}
subset_color <- c("#f8ecac","#edc7dd","#f8cbaa",
                  "#cad3e5","#e2eec6","#afdac7")
ggplot(data=data, aes(x=logFC, y =-log10(adj.P.Val),color=significant)) +
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
  #labs(title = "Volcano PLot of genes from matrisomeDB") + 
  #theme(plot.title = element_text(hjust=0.5,face = "bold",vjust = 1,size = 20)) + 
  theme(axis.title = element_text(face = "bold")) + 
  theme(text = element_text(face = "bold")) + 
  # geom_point(data=subset(data, abs(logFC) >= 6 & adj.P.Val < 0.05),alpha=0.8, size=3,col="green")+
  # geom_text_repel(data=subset(data, abs(logFC) > 6 & adj.P.Val < 0.05),
  #                 aes(label=gene),col="black",alpha = 0.8)
  geom_point(data=data_subset_list[[1]],alpha=0.85, size=1.5,col=subset_color[1])+
  geom_point(data=data_subset_list[[2]],alpha=0.85, size=1.5,col=subset_color[2])+
  geom_point(data=data_subset_list[[3]],alpha=0.85, size=1.5,col=subset_color[3])+
  geom_point(data=data_subset_list[[4]],alpha=0.85, size=1.5,col=subset_color[4])+
  geom_point(data=data_subset_list[[5]],alpha=0.85, size=1.5,col=subset_color[5])+
  geom_point(data=data_subset_list[[6]],alpha=0.85, size=1.5,col=subset_color[6])+
  geom_text_repel(data=data_subset,aes(label=gene),col="black",alpha = 0.8,size=3)

#coord_fixed(ratio = 2)
#library(export)
## 导成PPT可编辑的格式
graph2pdf(file="04_Output_figure/Version3/09_VolcanoPlot.pdf",aspectr=2, font = "Arial",
          width = 10, height = 8, bg = "transparent")
  
#### 09-Ribbon plot for DLsignal but notMsignal####
library(readxl)
sheet2 <- readxl::read_xlsx("./01_Original_data/DL_noMl.xlsx",sheet = 2)
sheet2 <- sheet2[!is.na(sheet2$`Gene name`),]
sheet2 <- sheet2[-grep("--",sheet2$`Gene name`),]
sheet2 <- sheet2[,c(1,3)]

regu_cell <- readxl::read_xlsx("01_Original_data/DL_noMl.xlsx",sheet = 3)
regu_expr <- readxl::read_xlsx("01_Original_data/DL_noMl.xlsx",sheet = 4)
transfer <- readxl::read_xlsx("01_Original_data/DL_noMl.xlsx",sheet = 5)
sign_transduct <- readxl::read_xlsx("01_Original_data/DL_noMl.xlsx",sheet = 6)
metabolism <- readxl::read_xlsx("01_Original_data/DL_noMl.xlsx",sheet = 7)
biosynthesis <- readxl::read_xlsx("01_Original_data/DL_noMl.xlsx",sheet = 8)

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
riverplot( r,srt=0,lty=0,textcel=1,
           nodewidth = 1,node_margin = 0.1,
           cex=0.2,plot_area = 0.9)
graph2pdf(file="04_Output_figure/Version3/10_riverPlot.pdf",aspectr=2, font = "Arial",
          width = 10, height = 8, bg = "transparent")

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
riverplot( r_all,srt=0,lty=0,textcel=1,
           nodewidth = 1,node_margin = 0.1,
           cex=0.2,plot_area = 0.9)

graph2pdf(file="04_Output_figure/Version3/10_riverPlot_v2.pdf",aspectr=2, font = "Arial",
          width = 10, height = 8, bg = "transparent")
library(RColorBrewer)
display.brewer.all()
set3= brewer.pal(n=6,name = "Set3")
nodes_all$col[2:7] <- set3
r_all <- makeRiver(nodes_all, edges_all)
op <- par(cex=0.8)
riverplot( r_all,srt=0,lty=0,textcel=1,
           nodewidth = 1,node_margin = 0.1,
           cex=0.2,plot_area = 0.9)
graph2pdf(file="04_Output_figure/Version3/10_riverPlot_v3.pdf",aspectr=2, font = "Arial",
          width = 10, height = 8, bg = "transparent")

#col <- c("#2F4F4F","#8B2323","#8B3E2F","#8B3626")


