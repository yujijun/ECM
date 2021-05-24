library(VennDiagram)
library(gplots)
library(readxl)
library(export)

setwd('/1.deglist/DL_d5vsRBE/')
a <- read_excel("DL_d5vsRBE_deg.xlsx",sheet=1)
a_1 <- a[which(a$padj<0.05 & abs(a$log2FoldChange)>1),]
a_2 <- a_1[1]

setwd('/1.deglist/DL_d10vsRBE')
b <- read_excel("DL_d10vsRBE_deg.xlsx",sheet=1)
b_1 <- b[which(b$padj<0.05 & abs(b$log2FoldChange)>1),]
b_2 <- b_1[1]

setwd('/1.deglist/DL_d15vsRBE/')
c <- read_excel("DL_d15vsRBE_deg.xlsx",sheet=1)
c_1 <- c[which(c$padj<0.05 & abs(c$log2FoldChange)>1),]
c_2 <- c_1[1]

setwd('/1.deglist/M_d5vsRBE')
d <- read_excel("M_d5vsRBE_deg.xlsx",sheet=1)
d_1 <- d[which(d$padj<0.05 & abs(d$log2FoldChange)>1),]
d_2 <- d_1[1]

setwd('/1.deglist/M_d10vsRBE/')
e <- read_excel("M_d10vsRBE_deg.xlsx",sheet=1)
e_1 <- e[which(e$padj<0.05 & abs(e$log2FoldChange)>1),]
e_2 <- e_1[1]

setwd('/1.deglist/M_d15vsRBE/')
f <- read_excel("M_d15vsRBE_deg.xlsx",sheet=1)
f_1 <- f[which(f$padj<0.05 & abs(f$log2FoldChange)>1),]
f_2 <- f_1[1]

setwd('/1.deglist/DL_d5vsM_d5/')
g <- read_excel("DL_d5vsM_d5_deg.xlsx",sheet=1)
g_1 <- g[which(g$padj<0.05 & abs(g$log2FoldChange)>1),]
g_2 <- g_1[1]

setwd('/1.deglist/DL_d10vsM_d10/')
h <- read_excel("DL_d10vsM_d10_deg.xlsx",sheet=1)
h_1 <- h[which(h$padj<0.05 & abs(h$log2FoldChange)>1),]
h_2 <- h_1[1]

setwd('/1.deglist/DL_d15vsM_d15/')
i <- read_excel("DL_d15vsM_d15_deg.xlsx",sheet=1)
i_1 <- i[which(i$padj<0.05 & abs(i$log2FoldChange)>1),]
i_2 <- i_1[1]

####venn plot####

a_1 <- a_1$gene_id
d_1 <- d_1$gene_id
Venn_list_IN <- list(DL_d5vsRBE = a_1,M_d5vsRBE = d_1)
venn.plot <- venn.diagram(x= Venn_list_IN, filename = NULL, height = 450, width = 450,resolution =300, 
                          imagetype="png", col="transparent",fill=c("#FCEBA0","#D7A3F0"),
                          alpha = 0.50, cex=0.8, cat.cex=0.8,force.unique = F)
grid.draw(venn.plot)
graph2ppt(file="/1.deglist/DL_d5vsRBE_&_M_d5vsRBE.pptx")

b_1 <- b_1$gene_id
e_1 <- e_1$gene_id
Venn_list_IN <- list(DL_d10vsRBE = b_1,M_d10vsRBE = e_1)
venn.plot <- venn.diagram(x= Venn_list_IN, filename = NULL, height = 450, width = 450,resolution =300, 
                          imagetype="png", col="transparent",fill=c("#FCEBA0","#D7A3F0"),
                          alpha = 0.50, cex=0.8, cat.cex=0.8,force.unique = F)
grid.draw(venn.plot)
graph2ppt(file="/1.deglist/DL_d10vsRBE_&_M_d10vsRBE.pptx")

c_1 <- c_1$gene_id
f_1 <- f_1$gene_id
Venn_list_IN <- list(DL_d15vsRBE = c_1,M_d15vsRBE = f_1)
venn.plot <- venn.diagram(x= Venn_list_IN, filename = NULL, height = 450, width = 450,resolution =300, 
                          imagetype="png", col="transparent",fill=c("#FCEBA0","#D7A3F0"),
                          alpha = 0.50, cex=0.8, cat.cex=0.8,force.unique = F)
grid.draw(venn.plot)
graph2ppt(file="/1.deglist/DL_d15vsRBE_&_M_d15vsRBE.pptx")

g_1 <- g_1$gene_id
h_1 <- h_1$gene_id
Venn_list_IN <- list(DL_d5vsM_d5 = g_1,DL_d10vsM_d10 = h_1)
venn.plot <- venn.diagram(x= Venn_list_IN, filename = NULL, height = 450, width = 450,resolution =300, 
                          imagetype="png", col="transparent",fill=c("#FCEBA0","#D7A3F0"),
                          alpha = 0.50, cex=0.8, cat.cex=0.8,force.unique = F)
grid.draw(venn.plot)
graph2ppt(file="/1.deglist/DL_d5vsM_d5_&_DL_d10vsM_d10.pptx")

g_1 <- g_1$gene_id
i_1 <- i_1$gene_id
Venn_list_IN <- list(DL_d5vsM_d5 = g_1,DL_d15vsM_d15 = i_1)
venn.plot <- venn.diagram(x= Venn_list_IN, filename = NULL, height = 450, width = 450,resolution =300, 
                          imagetype="png", col="transparent",fill=c("#FCEBA0","#D7A3F0"),
                          alpha = 0.50, cex=0.8, cat.cex=0.8,force.unique = F)
grid.draw(venn.plot)
graph2ppt(file="/1.deglist/DL_d5vsM_d5_&_DL_d15vsM_d15.pptx")

h_1 <- h_1$gene_id
i_1 <- i_1$gene_id
Venn_list_IN <- list(DL_d10vsM_d10 = g_1,DL_d15vsM_d15 = i_1)
venn.plot <- venn.diagram(x= Venn_list_IN, filename = NULL, height = 450, width = 450,resolution =300, 
                          imagetype="png", col="transparent",fill=c("#FCEBA0","#D7A3F0"),
                          alpha = 0.50, cex=0.8, cat.cex=0.8,force.unique = F)
grid.draw(venn.plot)
graph2ppt(file="/1.deglist/DL_d10vsM_d10_&_DL_d15vsM_d15.pptx")

####valcano plot####
library(ggplot2)
a$change <- as.factor(
  ifelse(
    a$padj<0.05 & abs(a$log2FoldChange)>1,
    ifelse(a$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)
valcano <- ggplot(data=a, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=20) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) + 
  ggtitle("DL_d5vsRBE_deg") + 
  scale_color_manual(name="", values=c("red", "green", "grey"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
valcano
graph2ppt(file="/1.deglist/DL_d5vsRBE_deg_valcano.pptx")

b$change <- as.factor(
  ifelse(
    b$padj<0.05 & abs(b$log2FoldChange)>1,
    ifelse(b$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)
valcano <- ggplot(data=b, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=20) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) + 
  ggtitle("DL_d10vsRBE_deg") + 
  scale_color_manual(name="", values=c("red", "green", "grey"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
valcano
graph2ppt(file="/1.deglist/DL_d10vsRBE_deg_valcano.pptx")

c$change <- as.factor(
  ifelse(
    c$padj<0.05 & abs(c$log2FoldChange)>1,
    ifelse(c$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)
valcano <- ggplot(data=c, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=20) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) + 
  ggtitle("DL_d15vsRBE_deg") + 
  scale_color_manual(name="", values=c("red", "green", "grey"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
valcano
graph2ppt(file="/1.deglist/DL_d15vsRBE_deg_valcano.pptx")


d$change <- as.factor(
  ifelse(
    d$padj<0.05 & abs(d$log2FoldChange)>1,
    ifelse(d$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)
valcano <- ggplot(data=d, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=20) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) + 
  ggtitle("M_d5vsRBE_deg") + 
  scale_color_manual(name="", values=c("red", "green", "grey"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
valcano
graph2ppt(file="/1.deglist/M_d5vsRBE_deg_valcano.pptx")

e$change <- as.factor(
  ifelse(
    e$padj<0.05 & abs(e$log2FoldChange)>1,
    ifelse(e$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)
valcano <- ggplot(data=e, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=20) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) + 
  ggtitle("M_d10vsRBE_deg") + 
  scale_color_manual(name="", values=c("red", "green", "grey"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
valcano
graph2ppt(file="/1.deglist/M_d10vsRBE_deg_valcano.pptx")

f$change <- as.factor(
  ifelse(
    f$padj<0.05 & abs(f$log2FoldChange)>1,
    ifelse(f$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)
valcano <- ggplot(data=f, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=20) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) + 
  ggtitle("M_d15vsRBE_deg") + 
  scale_color_manual(name="", values=c("red", "green", "grey"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
valcano
graph2ppt(file="/1.deglist/M_d15vsRBE_deg_valcano.pptx")


g$change <- as.factor(
  ifelse(
    g$padj<0.05 & abs(g$log2FoldChange)>1,
    ifelse(g$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)
valcano <- ggplot(data=g, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=20) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) + 
  ggtitle("DL_d5vsM_d5_deg") + 
  scale_color_manual(name="", values=c("red", "green", "grey"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
valcano
graph2ppt(file="/1.deglist/DL_d5vsM_d5_deg_valcano.pptx")

h$change <- as.factor(
  ifelse(
    h$padj<0.05 & abs(h$log2FoldChange)>1,
    ifelse(h$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)
valcano <- ggplot(data=h, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=20) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) + 
  ggtitle("DL_d10vsM_d10_deg") + 
  scale_color_manual(name="", values=c("red", "green", "grey"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
valcano
graph2ppt(file="/1.deglist/DL_d10vsM_d10_deg_valcano.pptx")

i$change <- as.factor(
  ifelse(
    i$padj<0.05 & abs(i$log2FoldChange)>1,
    ifelse(i$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)
valcano <- ggplot(data=i, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=20) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) + 
  ggtitle("DL_d15vsM_d15_deg") + 
  scale_color_manual(name="", values=c("red", "green", "grey"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
valcano
graph2ppt(file="/1.deglist/DL_d15vsM_d15_deg_valcano.pptx")


####bar plot####

library(dplyr)
a$change <- as.factor(
  ifelse(
    a$padj<0.05 & abs(a$log2FoldChange)>1,
    ifelse(a$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)

b$change <- as.factor(
  ifelse(
    b$padj<0.05 & abs(b$log2FoldChange)>1,
    ifelse(b$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)

c$change <- as.factor(
  ifelse(
    c$padj<0.05 & abs(c$log2FoldChange)>1,
    ifelse(c$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)

d$change <- as.factor(
  ifelse(
    d$padj<0.05 & abs(d$log2FoldChange)>1,
    ifelse(d$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)

e$change <- as.factor(
  ifelse(
    e$padj<0.05 & abs(e$log2FoldChange)>1,
    ifelse(e$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)

f$change <- as.factor(
  ifelse(
    f$padj<0.05 & abs(f$log2FoldChange)>1,
    ifelse(f$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)

g$change <- as.factor(
  ifelse(
    g$padj<0.05 & abs(g$log2FoldChange)>1,
    ifelse(g$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)

h$change <- as.factor(
  ifelse(
    h$padj<0.05 & abs(h$log2FoldChange)>1,
    ifelse(h$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)

i$change <- as.factor(
  ifelse(
    i$padj<0.05 & abs(i$log2FoldChange)>1,
    ifelse(i$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)

data_group<- group_by(a,change)
data_GroupByID1<- summarise(data_group,count = n())
data_GroupByID1[data_GroupByID1$change == 'NoDiff',]$count = sum(data_GroupByID1$count)-data_GroupByID1[data_GroupByID1$change == 'NoDiff',]$count
data_GroupByID1$type = 'A'


data_group<- group_by(b,change)
data_GroupByID2<- summarise(data_group,count = n())
data_GroupByID2[data_GroupByID2$change == 'NoDiff',]$count = sum(data_GroupByID2$count)-data_GroupByID2[data_GroupByID2$change == 'NoDiff',]$count
data_GroupByID2$type = 'B'

data_group<- group_by(c,change)
data_GroupByID3<- summarise(data_group,count = n())
data_GroupByID3[data_GroupByID2$change == 'NoDiff',]$count = sum(data_GroupByID3$count)-data_GroupByID3[data_GroupByID2$change == 'NoDiff',]$count
data_GroupByID3$type = 'C'

data_group<- group_by(d,change)
data_GroupByID4<- summarise(data_group,count = n())
data_GroupByID4[data_GroupByID2$change == 'NoDiff',]$count = sum(data_GroupByID4$count)-data_GroupByID4[data_GroupByID2$change == 'NoDiff',]$count
data_GroupByID4$type = 'D'

data_group<- group_by(e,change)
data_GroupByID5<- summarise(data_group,count = n())
data_GroupByID5[data_GroupByID2$change == 'NoDiff',]$count = sum(data_GroupByID5$count)-data_GroupByID5[data_GroupByID2$change == 'NoDiff',]$count
data_GroupByID5$type = 'E'

data_group<- group_by(f,change)
data_GroupByID6<- summarise(data_group,count = n())
data_GroupByID6[data_GroupByID2$change == 'NoDiff',]$count = sum(data_GroupByID6$count)-data_GroupByID6[data_GroupByID2$change == 'NoDiff',]$count
data_GroupByID6$type = 'F'

data_group<- group_by(g,change)
data_GroupByID7<- summarise(data_group,count = n())
data_GroupByID7[data_GroupByID2$change == 'NoDiff',]$count = sum(data_GroupByID7$count)-data_GroupByID7[data_GroupByID2$change == 'NoDiff',]$count
data_GroupByID7$type = 'G'

data_group<- group_by(h,change)
data_GroupByID8<- summarise(data_group,count = n())
data_GroupByID8[data_GroupByID2$change == 'NoDiff',]$count = sum(data_GroupByID8$count)-data_GroupByID8[data_GroupByID2$change == 'NoDiff',]$count
data_GroupByID8$type = 'H'

data_group<- group_by(i,change)
data_GroupByID9<- summarise(data_group,count = n())
data_GroupByID9[data_GroupByID2$change == 'NoDiff',]$count = sum(data_GroupByID9$count)-data_GroupByID9[data_GroupByID2$change == 'NoDiff',]$count
data_GroupByID9$type = 'I'

data_groupALL  = rbind(data_GroupByID1,data_GroupByID2,data_GroupByID3,data_GroupByID4,data_GroupByID5,data_GroupByID6,data_GroupByID7,data_GroupByID8,data_GroupByID9)

library(ggplot2)
ggplot(data = data_groupALL, mapping = aes(x = factor(type), y = count, fill = change)) + 
  geom_bar(stat = 'identity', position = 'dodge') + scale_fill_brewer(palette = 'Accent')+
  scale_fill_discrete(breaks=c("Up", "Down", "NoDiff"),labels=c("Up", "Down", "all"))
graph2ppt(file="/1.deglist/barplot_for_diffgenecount.pptx")


####KEGG/GO_term_analysis####
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("openxlsx")
BiocManager::install("enrichplot")
BiocManager::install("xml2")
BiocManager::install("clusterProfiler")
#BiocManager::install("org.Rn.eg.db") 
BiocManager::install("org.Hs.eg.db") 
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
install.packages("tweenr")

library(BiocManager)
BiocManager::valid()

library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)
library(enrichplot)

#GO
ego_BP <- enrichGO(a_2[[1]], org.Hs.eg.db,ont='BP',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_BP),file="/1.deglist/DL_d5vsRBE_GO_BP.csv")
ego_CC <- enrichGO(a_2[[1]], org.Hs.eg.db,ont='CC',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_CC),file="/1.deglist/DL_d5vsRBE_GO_CC.csv")
ego_MF <- enrichGO(a_2[[1]], org.Hs.eg.db,ont='MF',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_MF),file="/1.deglist/DL_d5vsRBE_GO_MF.csv")
#KEGG
a_3<-select(org.Hs.eg.db, keys=as.vector(a_2[[1]]), columns=c("ENSEMBL","UNIGENE","ENTREZID","CHR","GO","GENENAME"), keytype="ENSEMBL")
kegg <- enrichKEGG(a_3$ENTREZID, organism = "hsa",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod="fdr",minGSSize = 3, maxGSSize = 2000)
write.csv(summary(kegg),file="DL_d5vsRBE_KEGG.csv")


#GO
ego_BP <- enrichGO(b_2[[1]], org.Hs.eg.db,ont='BP',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_BP),file="/1.deglist/DL_d10vsRBE_GO_BP.csv")
ego_CC <- enrichGO(b_2[[1]], org.Hs.eg.db,ont='CC',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_CC),file="/1.deglist/DL_d10vsRBE_GO_CC.csv")
ego_MF <- enrichGO(b_2[[1]], org.Hs.eg.db,ont='MF',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_MF),file="/1.deglist/DL_d10vsRBE_GO_MF.csv")
#KEGG
b_3<-select(org.Hs.eg.db, keys=as.vector(b_2[[1]]), columns=c("ENSEMBL","UNIGENE","ENTREZID","CHR","GO","GENENAME"), keytype="ENSEMBL")
kegg <- enrichKEGG(b_3$ENTREZID, organism = "hsa",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod="fdr",minGSSize = 3, maxGSSize = 2000)
write.csv(summary(kegg),file="DL_d10vsRBE_KEGG.csv")


#GO
ego_BP <- enrichGO(c_2[[1]], org.Hs.eg.db,ont='BP',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_BP),file="/1.deglist/DL_d15vsRBE_GO_BP.csv")
ego_CC <- enrichGO(c_2[[1]], org.Hs.eg.db,ont='CC',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_CC),file="/1.deglist/DL_d15vsRBE_GO_CC.csv")
ego_MF <- enrichGO(c_2[[1]], org.Hs.eg.db,ont='MF',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_MF),file="/1.deglist/DL_d15vsRBE_GO_MF.csv")
#KEGG
c_3<-select(org.Hs.eg.db, keys=as.vector(c_2[[1]]), columns=c("ENSEMBL","UNIGENE","ENTREZID","CHR","GO","GENENAME"), keytype="ENSEMBL")
kegg <- enrichKEGG(c_3$ENTREZID, organism = "hsa",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod="fdr",minGSSize = 3, maxGSSize = 2000)
write.csv(summary(kegg),file="DL_d15vsRBE_KEGG.csv")

#GO
ego_BP <- enrichGO(d_2[[1]], org.Hs.eg.db,ont='BP',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_BP),file="/1.deglist/M_d5vsRBE_GO_BP.csv")
ego_CC <- enrichGO(d_2[[1]], org.Hs.eg.db,ont='CC',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_CC),file="/1.deglist/M_d5vsRBE_GO_CC.csv")
ego_MF <- enrichGO(d_2[[1]], org.Hs.eg.db,ont='MF',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_MF),file="/1.deglist/M_d5vsRBE_GO_MF.csv")
#KEGG
d_3<-select(org.Hs.eg.db, keys=as.vector(d_2[[1]]), columns=c("ENSEMBL","UNIGENE","ENTREZID","CHR","GO","GENENAME"), keytype="ENSEMBL")
kegg <- enrichKEGG(d_3$ENTREZID, organism = "hsa",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod="fdr",minGSSize = 3, maxGSSize = 2000)
write.csv(summary(kegg),file="M_d5vsRBE_KEGG.csv")


#GO
ego_BP <- enrichGO(e_2[[1]], org.Hs.eg.db,ont='BP',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_BP),file="/1.deglist/M_d10vsRBE_GO_BP.csv")
ego_CC <- enrichGO(e_2[[1]], org.Hs.eg.db,ont='CC',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_CC),file="/1.deglist/M_d10vsRBE_GO_CC.csv")
ego_MF <- enrichGO(e_2[[1]], org.Hs.eg.db,ont='MF',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_MF),file="/1.deglist/M_d10vsRBE_GO_MF.csv")
#KEGG
e_3<-select(org.Hs.eg.db, keys=as.vector(e_2[[1]]), columns=c("ENSEMBL","UNIGENE","ENTREZID","CHR","GO","GENENAME"), keytype="ENSEMBL")
kegg <- enrichKEGG(e_3$ENTREZID, organism = "hsa",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod="fdr",minGSSize = 3, maxGSSize = 2000)
write.csv(summary(kegg),file="M_d10vsRBE_KEGG.csv")

#GO
ego_BP <- enrichGO(f_2[[1]], org.Hs.eg.db,ont='BP',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_BP),file="/1.deglist/M_d15vsRBE_GO_BP.csv")
ego_CC <- enrichGO(f_2[[1]], org.Hs.eg.db,ont='CC',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_CC),file="/1.deglist/M_d15vsRBE_GO_CC.csv")
ego_MF <- enrichGO(f_2[[1]], org.Hs.eg.db,ont='MF',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_MF),file="/1.deglist/M_d15vsRBE_GO_MF.csv")
#KEGG
f_3<-select(org.Hs.eg.db, keys=as.vector(f_2[[1]]), columns=c("ENSEMBL","UNIGENE","ENTREZID","CHR","GO","GENENAME"), keytype="ENSEMBL")
kegg <- enrichKEGG(f_3$ENTREZID, organism = "hsa",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod="fdr",minGSSize = 3, maxGSSize = 2000)
write.csv(summary(kegg),file="M_d15vsRBE_KEGG.csv")

#GO
ego_BP <- enrichGO(g_2[[1]], org.Hs.eg.db,ont='BP',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_BP),file="/1.deglist/DL_d5vsM_d5_GO_BP.csv")
ego_CC <- enrichGO(g_2[[1]], org.Hs.eg.db,ont='CC',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_CC),file="/1.deglist/DL_d5vsM_d5_GO_CC.csv")
ego_MF <- enrichGO(g_2[[1]], org.Hs.eg.db,ont='MF',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_MF),file="/1.deglist/DL_d5vsM_d5_GO_MF.csv")
#KEGG
g_3<-select(org.Hs.eg.db, keys=as.vector(g_2[[1]]), columns=c("ENSEMBL","UNIGENE","ENTREZID","CHR","GO","GENENAME"), keytype="ENSEMBL")
kegg <- enrichKEGG(g_3$ENTREZID, organism = "hsa",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod="fdr",minGSSize = 3, maxGSSize = 2000)
write.csv(summary(kegg),file="DL_d5vsM_d5_KEGG.csv")

#GO
ego_BP <- enrichGO(h_2[[1]], org.Hs.eg.db,ont='BP',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_BP),file="/1.deglist/DL_d10vsM_d10_GO_BP.csv")
ego_CC <- enrichGO(h_2[[1]], org.Hs.eg.db,ont='CC',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_CC),file="/1.deglist/DL_d10vsM_d10_GO_CC.csv")
ego_MF <- enrichGO(h_2[[1]], org.Hs.eg.db,ont='MF',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_MF),file="/1.deglist/DL_d10vsM_d10_GO_MF.csv")
#KEGG
h_3<-select(org.Hs.eg.db, keys=as.vector(h_2[[1]]), columns=c("ENSEMBL","UNIGENE","ENTREZID","CHR","GO","GENENAME"), keytype="ENSEMBL")
kegg <- enrichKEGG(h_3$ENTREZID, organism = "hsa",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod="fdr",minGSSize = 3, maxGSSize = 2000)
write.csv(summary(kegg),file="DL_d10vsM_d10_KEGG.csv")

#GO
ego_BP <- enrichGO(i_2[[1]], org.Hs.eg.db,ont='BP',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_BP),file="/1.deglist/DL_d15vsM_d15_GO_BP.csv")
ego_CC <- enrichGO(i_2[[1]], org.Hs.eg.db,ont='CC',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_CC),file="/1.deglist/DL_d15vsM_d15_GO_CC.csv")
ego_MF <- enrichGO(i_2[[1]], org.Hs.eg.db,ont='MF',keyType="ENSEMBL",pvalueCutoff=1,qvalueCutoff=1)
write.csv(summary(ego_MF),file="/1.deglist/DL_d15vsM_d15_GO_MF.csv")
#KEGG
i_3<-select(org.Hs.eg.db, keys=as.vector(i_2[[1]]), columns=c("ENSEMBL","UNIGENE","ENTREZID","CHR","GO","GENENAME"), keytype="ENSEMBL")
kegg <- enrichKEGG(i_3$ENTREZID, organism = "hsa",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod="fdr",minGSSize = 3, maxGSSize = 2000)
write.csv(summary(kegg),file="DL_d15vsM_d15_KEGG.csv")

#####barplot of BP####
setwd("/1.deglist/for rna-seq of ecm")
library(DOSE)
data<- read.csv("DL_d5vsM_d5_GO_BP.csv")
data1 <- data[order(data$pvalue),]
data2 <- data1[1:10,]
count <- as.numeric(unlist(strsplit(data2$GeneRatio,"/284",fixed=T)))
data3 <- data.frame(data2[,3],count,data2[,7])
colnames(data3) <- c("Description","count","pvalue")
p <- ggplot(data=data3,aes(x=Description,y=count,fill=pvalue))
p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(fill='transparent',color='gray'),
                 axis.text.y=element_text(color="black",size=12))
p3 <- p2 + ylim(0,30) + scale_fill_gradient(low="red",high="blue")
p4 <- p3 + scale_x_discrete(limits=rev(data3[,1])) +labs(x="",y="",title="Biological_Process")
png("D5_Biological_Process_enrich.png",width=800,height=600, units = "px")
print(p4)
dev.off()

####Bubble Plot####

load("/ECM-main_test_yjj/4.3_ECM数据分析项目_输出结果_20210120/output_v2/merge_orig_ref_filter.RData")
choosen_column <- c(1,9:17,31)
merge_orig_ref_choosen <- merge_orig_ref_filter_na[,choosen_column]
merge_na_heatmap <- merge_orig_ref_choosen[,c(1:7,10,11)]
merge_na <- merge_orig_ref_choosen[!(rowSums(is.na(merge_orig_ref_choosen)) > 7),]

rownames(merge_na_heatmap) <- merge_na$`Gene name`
merge_na_heatmap <- merge_na_heatmap[order(merge_na_heatmap$Category),]
merge_na_heatmap <- merge_na_heatmap[c(1:87,148:156,88:147,157:167),]
#merge_na_heatmap$Category <- factor(merge_na_heatmap$Category, levels = c("Collagens","ECM Glycoproteins","Proteoglycans",
#                                                                          "ECM Regulators","ECM-affiliated Proteins","Secreted Factors"))

library(dplyr)
library(tidyr)
merge_na_hehe <- merge_na_heatmap[1:23,]
merge_na_hehe <-  separate(merge_na_hehe,"Gene name", c("Gene name", "num"), "a")
groupdata<-dplyr::group_by(merge_na_hehe,`Gene name`) %>% dplyr::summarize(DL1=sum(DL1>0,na.rm = T),DL2=sum(DL2>0,na.rm = T),DL3=sum(DL3>0,na.rm = T),
                                                                           Matrigel1=sum(Matrigel1>0,na.rm = T),Matrigel2=sum(Matrigel2>0,na.rm = T),Matrigel3=sum(Matrigel3>0,na.rm = T))
library(reshape2)
data_melt = melt(data = groupdata,id.vars=c("Gene name"),variable.name="cat",value.name="num")
data_melt$y = 1

####venn plot####
setwd('//03-ECM/03 discussion/')
aa <- read_excel("Colltype.xlsx",sheet=1)

dl <- aa$dl
m <- aa$M
venn.plot <- venn.diagram(x= Venn_list_IN, filename = NULL, height = 450, width = 600,resolution =300, 
                          imagetype="png", col="transparent",fill=c("#67001F","#756d70"),
                          alpha = 0.50, cex=1.0, cat.cex=1.0,force.unique = F,na = "none" )
grid.draw(venn.plot)
graph2ppt(file="/03 discussion/Colltype.pptx")

# Testing existence & format of the file
file.exists("DL_d15vsM_d15_deg.xlsx")
format_from_ext("DL_d15vsM_d15_deg.xlsx")
format_from_signature("DL_d15vsM_d15_deg.xlsx")

####Bubble Plot####
#### convert width into long ####
merge_na_hehe <- merge_na_hehe[,c(1,2:8)]
merge_na_hehe_long <- merge_na_hehe %>%  pivot_longer(cols = DL1:Matrigel3,names_to = "sample",values_to = "expression")
merge_na_hehe_long$num <- NULL
merge_na_hehe_long <- merge_na_hehe_long[!is.na(merge_na_hehe_long$expression),]
library(Rmisc)
colnames(merge_na_hehe_long) <- c("matrigeltype","sample","expression")
merge_na_hehe_long_summary <- merge_na_hehe_long %>%  summarySE(measurevar = "expression",
                                                                groupvars = c("matrigeltype","sample"))

my36colors <-c('#53A85F', '#F1BB72', '#F3B1A0', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175','#E5D2DD', '#D6E7A3')

library(ggplot2)
ggplot(merge_na_hehe_long_summary,aes(x=sample,y=expression)) + 
  geom_point(aes(size=N,colour=matrigeltype)) + 
  scale_colour_manual(values=my36colors)+
  scale_size(rang=c(3,10))
graph2ppt(file="/03 discussion/Bubble_Plotformatrigel.pptx")

ggplot(merge_na_hehe_long_summary,aes(x=sample,y=expression)) + 
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
        axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"))
graph2ppt(file="/03 discussion/Bubble_Plotformatrigel1.pptx")
