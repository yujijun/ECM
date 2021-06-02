## library ####
library(ggplot2)
library(export)

## hyper-parameter ####
if(T){
  mytheme <- theme_classic() + 
    theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
          axis.title = element_text(size = 12,color ="black",face = "bold"), 
          axis.text = element_text(size= 8,color = "black"),
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
figure_path <- "./4.4_ECM数据分析项目_图片输出_20210120/version8_last_v2_20210531//"
output_path <-"./4.3_ECM数据分析项目_输出结果_20210120/output_last/"
color_category <- c(brewer.pal(n = 8, name = "Pastel2")[1:6])
names(color_category) <- c("secreted_factors","proteoglycans","ECM-affiliated",
                           "regulators","glycoproteins","collagens")

## visualization ####
### ECM regulator in DL ####
data1 <- merge_orig_ref_filter[,c("DL1","DL2","DL3","Matrigel1","Matrigel2","Matrigel3","Category","Gene name")]
data1$DL <- apply(data.frame(data1$DL1,data1$DL2,data1$DL3), 1, mean,na.rm=T)
data1$Matrigel <- apply(data.frame(data1$Matrigel1,data1$Matrigel2,data1$Matrigel3), 1, mean,na.rm=T)

data2 = data1[which(data1$Category =='ECM Regulators'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:10),]
data2$DL <- round(data2$DL,1)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
fill_color <- color_category["regulators"]
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(fill = fill_color, stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = DL,x=`Gene name`,size=10),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$DL)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none") 
figure_name <- paste0(figure_path,"06_DL_ECM Regulators.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

### ECM regulator of same genes in Matrigel ####
datann <- data1[which(data1$`Gene name`%in% c('Tgm2','Loxl1','Adamtsl4','F13b','Ambp','Serpinb6a','Hrg','Plg','Htra1','Lox')),]
datann$Matrigel <- round(datann$Matrigel,2)
datann$'Gene name' <-factor(datann$'Gene name', levels = rev(c('Tgm2','Loxl1','F13b','Adamtsl4','Ambp','Serpinb6a','Hrg','Plg','Htra1','Lox')),ordered = F)
#data2 <- datann[rank=(datann$Matrigel,decreasing = T),]
#datann1$Matrigel <- round(datann1$Matrigel,2)
#datann1$`Gene name` <- reorder(datann1$`Gene name`,datann$Matrigel)
fill_color <- alpha(color_category["regulators"],alpha = 0.5)
ggplot(data=datann, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = fill_color, stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(datann$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none") 
figure_name <- paste0(figure_path,"06_M_ECM Regulators_samegene.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")


### ECM regulator in Matrigel ####
# data2 = data1[which(data1$Category =='ECM Regulators'),]
# data2 <- data2[order(data2$Matrigel,decreasing = T),]
# data2 <- data2[c(1:10),]
# data2$Matrigel <- round(data2$Matrigel,2)
# data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
# fill_color <- alpha(color_category["regulators"],alpha = 0.5)
# ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
#   geom_bar(fill = '#c7659e', stat="identity",position=position_dodge(0.7),width=0.5) +
#   geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
#   scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
#   coord_flip()+
#   xlab("") +
#   ylab("") +
#   theme(panel.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         plot.background = element_rect(fill = "transparent",colour = NA))+
#   theme_classic() +
#   theme(axis.line.x = element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
#   theme(axis.text.y = element_text(face = "bold",size = 8)) + 
#   theme(legend.position = "none") 
# figure_name <- paste0(figure_path,"06_M_ECM Regulators.pdf")
# graph2pdf(file=figure_name,
#           font = "Arial",
#           width = width, 
#           height = height, 
#           bg = "transparent")

### Secreted Factors in DL ####
data2 = data1[which(data1$Category =='Secreted Factors'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:8),]
data2$DL <- round(data2$DL,1)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
fill_color <- color_category["secreted_factors"]
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(fill = color_category["secreted_factors"], stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none") 
figure_name <- paste0(figure_path,"06_DL_Secreted Factors.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

### Secreted factors of same genes in Matrigel ####
datann <- data1[which(data1$`Gene name`%in% c('Megf6','Egfl7','S100a10','Pf4','Inhbe','Inhbc','Angptl6','Ccl9')),]
#datann[,c('Megf6','Egfl7','S100a10','Pf4','Angptl6','Inhbc','Inhbe','Ccl9')]
#loc = match(c('Megf6','Egfl7','S100a10','Pf4','Angptl6','Inhbc','Inhbe','Ccl9'),datann$`Gene name`)
#datann1<-datann[loc,]
datann$Matrigel <- round(datann$Matrigel,2)
#datann1$`Gene name` <- reorder(datann1$`Gene name`,datann1$Matrigel)
datann$`Gene name`<-factor(datann$`Gene name`,levels =rev(c('Megf6','Egfl7','S100a10','Pf4','Inhbe','Inhbc','Angptl6','Ccl9')), ordered = F)
#?factor
fill_color <- alpha(color_category["secreted_factors"],alpha=0.5)
ggplot(data=datann, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = fill_color, stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(datann$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none") 
figure_name <- paste0(figure_path,"06_M_Secreted Factors_samegene.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

# data2 = data1[which(data1$Category =='Secreted Factors'),]
# data2 <- data2[order(data2$Matrigel,decreasing = T),]
# data2 <- data2[c(1:6),]
# data2$Matrigel <- round(data2$Matrigel,2)
# data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
# ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
#   geom_bar(fill = '#74a690', stat="identity",position=position_dodge(0.7),width=0.5) +
#   geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
#   scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
#   coord_flip()+
#   xlab("") +
#   ylab("") +
#   theme(panel.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         plot.background = element_rect(fill = "transparent",colour = NA))+
#   theme_classic() +
#   theme(axis.line.x = element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
#   theme(axis.text.y = element_text(face = "bold",size = 8)) + 
#   theme(legend.position = "none") 
# figure_name <- paste0(figure_path,"06_M_Secreted Factors.pdf")
# graph2pdf(file=figure_name,
#           font = "Arial",
#           width = width, 
#           height = height, 
#           bg = "transparent")

### ECM-affiliated Proteins ####
data2 = data1[which(data1$Category =='ECM-affiliated Proteins'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:10),]
data2$DL <- round(data2$DL,1)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
fill_color <- color_category["ECM-affiliated"]
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(fill = fill_color, stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none") 
figure_name <- paste0(figure_path,"06_DL_ECM-affiliated Proteins.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

### ECM-affiliated protein  in Matrigel####
datann <- data1[which(data1$`Gene name`%in% c('Frem2','Anxa11','Anxa2','Anxa6','Frem1','Grem1','C1qtnf6','Anxa5','Lgals9','Anxa1')),]
datann$Matrigel <- round(datann$Matrigel,2)
datann$`Gene name`<-factor(datann$`Gene name`,levels =rev(c('Frem2','Anxa11','Anxa2','Anxa6','Frem1','Grem1','C1qtnf6','Anxa5','Lgals9','Anxa1')), ordered = F)

#datann$Matrigel <- round(datann$Matrigel,2)
#datann1$`Gene name` <- reorder(datann1$`Gene name`,datann1$Matrigel)
fill_color <- alpha(color_category["ECM-affiliated"],alpha=0.5)
ggplot(data=datann, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = fill_color, stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(datann$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none") 
figure_name <- paste0(figure_path,"06_M_ECM-affiliated Proteins_samegene.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")


# data2 = data1[which(data1$Category =='ECM-affiliated Proteins'),]
# data2 <- data2[order(data2$Matrigel,decreasing = T),]
# data2 <- data2[c(1:10),]
# data2$Matrigel <- round(data2$Matrigel,2)
# data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
# ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
#   geom_bar(fill = '#5a719e', stat="identity",position=position_dodge(0.7),width=0.5) +
#   geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
#   scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
#   coord_flip()+
#   xlab("") +
#   ylab("") +
#   theme(panel.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         plot.background = element_rect(fill = "transparent",colour = NA))+
#   theme_classic() +
#   theme(axis.line.x = element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
#   theme(axis.text.y = element_text(face = "bold",size = 8)) + 
#   theme(legend.position = "none") 
# figure_name <- paste0(figure_path,"06_M_ECM-affiliated Proteins.pdf")
# graph2pdf(file=figure_name,
#           font = "Arial",
#           width = width, 
#           height = height, 
#           bg = "transparent")

### proteoglycans in DL ####
data2 = data1[which(data1$Category =='Proteoglycans'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:9),]
data2$DL <- round(data2$DL,1)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
fill_color <- color_category["proteoglycans"]
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(fill = fill_color, stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none") 
figure_name <- paste0(figure_path,"06_DL_Proteoglycans.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")


### proteoglycans in Matrigel ####
datann <- data1[which(data1$`Gene name`%in% c('Lum','Prg2','Prelp','Vcan','Prg3','Bgn','Prg4','Dcn','Aspn')),]
datann$Matrigel <- round(datann$Matrigel,2)
datann$`Gene name`<-factor(datann$`Gene name`,levels =rev(c('Lum','Prg2','Prelp','Vcan','Prg3','Bgn','Prg4','Dcn','Aspn')), ordered = F)
fill_color <- alpha(color_category["proteoglycans"],alpha=0.5)
ggplot(data=datann, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = fill_color, stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(datann$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none") 
figure_name <- paste0(figure_path,"06_M_Proteoglycans_samegene.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

# data2 = data1[which(data1$Category =='Proteoglycans'),]
# data2 <- data2[order(data2$Matrigel,decreasing = T),]
# data2 <- data2[c(1:8),]
# data2$Matrigel <- round(data2$Matrigel,2)
# data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
# ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
#   geom_bar(fill = '#d18956', stat="identity",position=position_dodge(0.7),width=0.5) +
#   geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
#   scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
#   coord_flip()+
#   xlab("") +
#   ylab("") +
#   theme(panel.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         plot.background = element_rect(fill = "transparent",colour = NA))+
#   theme_classic() +
#   theme(axis.line.x = element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
#   theme(axis.text.y = element_text(face = "bold",size = 8)) + 
#   theme(legend.position = "none") 
# figure_name <- paste0(figure_path,"06_M_Proteoglycans.pdf")
# graph2pdf(file=figure_name,
#           font = "Arial",
#           width = width, 
#           height = height, 
#           bg = "transparent")


###  glycoproteins ####
data2 = data1[which(data1$Category =='ECM Glycoproteins'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:10),]
data2$DL <- round(data2$DL,1)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
fill_color <- color_category["glycoproteins"]
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(fill = fill_color, stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none") 
figure_name <- paste0(figure_path,"06_DL_ECM Glycoproteins.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

### glycoproteins in Matrigel ####
datann <- data1[which(data1$`Gene name`%in% c('Papln','Tgfbi','Postn','Ltbp4','Lama3','Fbn1','Emilin1','Dpt','Ecm1','Tnc')),]
datann$Matrigel <- round(datann$Matrigel,2)
datann$`Gene name`<-factor(datann$`Gene name`,levels =rev(c('Papln','Tgfbi','Postn','Ltbp4','Lama3','Fbn1','Emilin1','Dpt','Ecm1','Tnc')), ordered = F)
fill_color <- alpha(color_category["glycoproteins"],alpha = 0.5)
ggplot(data=datann, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = fill_color, stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(datann$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none")
figure_name <- paste0(figure_path,"06_M_ECM Glycoproteins_samegene.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")


# data2 = data1[which(data1$Category =='ECM Glycoproteins'),]
# data2 <- data2[order(data2$Matrigel,decreasing = T),]
# data2 <- data2[c(1:10),]
# data2$Matrigel <- round(data2$Matrigel,2)
# data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
# ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
#   geom_bar(fill = '#a4b381', stat="identity",position=position_dodge(0.7),width=0.5) +
#   geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
#   scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
#   coord_flip()+
#   xlab("") +
#   ylab("") +
#   theme(panel.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         plot.background = element_rect(fill = "transparent",colour = NA))+
#   theme_classic() +
#   theme(axis.line.x = element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
#   theme(axis.text.y = element_text(face = "bold",size = 8)) + 
#   theme(legend.position = "none") 
# figure_name <- paste0(figure_path,"06_M_ECM Glycoproteins.pdf")
# graph2pdf(file=figure_name,
#           font = "Arial",
#           width = width, 
#           height = height, 
#           bg = "transparent")


### collagens in DL ####
data2 = data1[which(data1$Category =='Collagens'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:10),]
data2$DL <- round(data2$DL,)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
fill_color <- color_category["collagens"]
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(fill = fill_color, stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none") 
figure_name <- paste0(figure_path,"06_DL_Collagens.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")

### collagen in Matrigel ####
datann <- data1[which(data1$`Gene name`%in% c('Col2a1','Col1a1','Col6a1','Col5a2','Col3a1','Col14a1','Col6a2','Col6a5','Col5a3','Col4a2')),]
datann$Matrigel <- round(datann$Matrigel,2)
datann$`Gene name`<-factor(datann$`Gene name`,levels =rev(c('Col2a1','Col1a1','Col6a1','Col5a2','Col3a1','Col14a1','Col6a2','Col6a5','Col5a3','Col4a2')), ordered = F)
fill_color  <- alpha(color_category["collagens"],alpha = 0.5)
ggplot(data=datann, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = fill_color, stat="identity",position=position_dodge(0.7),width=0.5) +
  geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(datann$Matrigel)+0.25)) +
  coord_flip()+
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.y = element_text(face = "bold",size = 8)) + 
  theme(legend.position = "none") 
figure_name <- paste0(figure_path,"06_M_Collagens_samegene.pdf")
graph2pdf(file=figure_name,
          font = "Arial",
          width = width, 
          height = height, 
          bg = "transparent")


# data2 = data1[which(data1$Category =='Collagens'),]
# data2 <- data2[order(data2$Matrigel,decreasing = T),]
# data2 <- data2[c(1:10),]
# data2$Matrigel <- round(data2$Matrigel,2)
# data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
# ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
#   geom_bar(fill = '#cfc17a', stat="identity",position=position_dodge(0.7),width=0.5) +
#   geom_text(aes(label = Matrigel,x=`Gene name`),vjust=0.2,hjust=-0.1)+
#   scale_y_continuous(expand = c(0,0),limits = c(0,max(data2$Matrigel)+0.25)) +
#   coord_flip()+
#   xlab("") +
#   ylab("") +
#   theme(panel.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         plot.background = element_rect(fill = "transparent",colour = NA))+
#   theme_classic() +
#   theme(axis.line.x = element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
#   theme(axis.text.y = element_text(face = "bold",size = 8)) + 
#   theme(legend.position = "none") 
# figure_name <- paste0(figure_path,"06_M_Collagens.pdf")
# graph2pdf(file=figure_name,
#           font = "Arial",
#           width = width, 
#           height = height, 
#           bg = "transparent")

