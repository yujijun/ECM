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
  geom_bar(fill = '#96CDCD', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_DL_ECM Regulators2.pptx")
graph2pdf(file="06_DL_ECM Regulators.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")

datann <- data1[which(data1$`Gene name`%in% c('Tgm2','Loxl1','Adamtsl4','F13b','Ambp','Serpinb6a','Hrg','Plg','Htra1','Lox')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#d7f5f5', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

graph2pdf(file="06_M_ECM Regulators_samegene.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

data2 = data1[which(data1$Category =='ECM Regulators'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#BBFFFF', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_DL_ECM Regulators2.pptx")
graph2pdf(file="06_M_ECM Regulators.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")



data2 = data1[which(data1$Category =='Secreted Factors'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:8),]
data2$DL <- round(data2$DL,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(fill = '#9ACD32', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_DL_Secreted Factors1.pptx")
graph2pdf(file="06_DL_Secreted Factors.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")

datann <- data1[which(data1$`Gene name`%in% c('Megf6','Egfl7','S100a10','Pf4','Angptl6','Inhbc','Inhbe','Ccl9')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2 <- data2[c(1:3),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#ddfc9d', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2pdf(file="06_M_Secreted Factors_samegene.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

data2 = data1[which(data1$Category =='Secreted Factors'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:6),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#C0FF3E', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_M_Secreted Factors1.pptx")
graph2pdf(file="06_M_Secreted Factors.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")


data2 = data1[which(data1$Category =='ECM-affiliated Proteins'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:10),]
data2$DL <- round(data2$DL,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(fill = "#292985", stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_DL_ECM-affiliated Proteins1.pptx")
graph2pdf(file="06_DL_ECM-affiliated Proteins.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")



datann <- data1[which(data1$`Gene name`%in% c('Frem2','Anxa11','Anxa2','Anxa6','Frem1','Grem1','C1qtnf6','Anxa5','Lgals9','Anxa1')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2 <- data2[c(1:9),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#81b0f7', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2pdf(file="06_M_ECM-affiliated Proteins_samegene.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")



data2 = data1[which(data1$Category =='ECM-affiliated Proteins'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#0066FF', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_M_ECM-affiliated Proteins1.pptx")
graph2pdf(file="06_M_ECM-affiliated Proteins.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")

data2 = data1[which(data1$Category =='Proteoglycans'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:9),]
data2$DL <- round(data2$DL,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(fill = '#007799', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_DL_Proteoglycans1.pptx")
graph2pdf(file="06_DL_Proteoglycans.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")

datann <- data1[which(data1$`Gene name`%in% c('Lum','Prg2','Prelp','Vcan','Prg3','Bgn','Prg4','Dcn','Aspn')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2 <- data2[c(1:8),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#00BBFF', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2pdf(file="06_M_Proteoglycans_samegene.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")


data2 = data1[which(data1$Category =='Proteoglycans'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:8),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#00BBFF', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_M_Proteoglycans1.pptx")
graph2pdf(file="06_M_Proteoglycans.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")


data2 = data1[which(data1$Category =='ECM Glycoproteins'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:10),]
data2$DL <- round(data2$DL,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(fill = '#CDCD00', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_DL_ECM Glycoproteins1.pptx")
graph2pdf(file="06_DL_ECM Glycoproteins.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")

datann <- data1[which(data1$`Gene name`%in% c('Papln','Emilin1','Fbn1','Dpt','Postn','Ltbp4','Tgfbi','Lama3','Ecm1','Tnc')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#FFFF00', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2pdf(file="06_M_ECM Glycoproteins_samegene.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")

data2 = data1[which(data1$Category =='ECM Glycoproteins'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#FFFF00', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_M_ECM Glycoproteins1.pptx")
graph2pdf(file="06_M_ECM Glycoproteins.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")

data2 = data1[which(data1$Category =='Collagens'),]
data2 <- data2[order(data2$DL,decreasing = T),]
data2 <- data2[c(1:10),]
data2$DL <- round(data2$DL,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$DL)
ggplot(data=data2, aes(x=`Gene name`, y=DL)) +
  geom_bar(fill = '#CD919E', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_DL_Collagens1.pptx")
graph2pdf(file="06_DL_Collagens.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")

datann <- data1[which(data1$`Gene name`%in% c('Col2a1','Col1a1','Col6a1','Col3a1','Col5a2','Col14a1','Col6a2','Col6a5','Col5a3','Col4a2')),]
data2 <- datann[order(datann$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#FFB5C5', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
graph2pdf(file="06_M_Collagens_samegene.pdf",aspectr=2, font = "Arial",
          width = 10, height = 2.5, bg = "transparent")


data2 = data1[which(data1$Category =='Collagens'),]
data2 <- data2[order(data2$Matrigel,decreasing = T),]
data2 <- data2[c(1:10),]
data2$Matrigel <- round(data2$Matrigel,2)
data2$`Gene name` <- reorder(data2$`Gene name`,data2$Matrigel)
ggplot(data=data2, aes(x=`Gene name`, y=Matrigel)) +
  geom_bar(fill = '#FFB5C5', stat="identity",position=position_dodge(0.7),width=0.5) +
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
  theme(axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  theme(axis.text.y = element_text(face = "bold",size = 10)) + 
  theme(legend.position = "none") 
#graph2ppt(file="06_M_Collagens1.pptx")
graph2pdf(file="06_M_Collagens.pdf",aspectr=2, font = "Arial",
          width = 5, height = 2.0, bg = "transparent")
