
############ fig 1b ####################

library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
library(dittoSeq)
library(readxl)

df <- as.data.frame(read_excel('data.xlsx'))
head(df)

mycol3 <- c4a('pu_bu',12)

ggplot(df, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = node,
               label = node)) +
  geom_sankey(flow.alpha = 0.5,
              flow.fill = 'grey', 
              flow.color = 'grey', 
              node.fill = mycol3, 
              smooth = 8,
              width = 0.3) +
  geom_sankey_text(size = 5,
                   color = "black")+
  theme_void() +
  theme(legend.position = 'none')

#########   fig 1c ########################

library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(org.Mm.eg.db)
library(GseaVis)

########## Duplicate code has been omitted and is shown here only once

setwd("DPP4/GSEA")
genelist_input <- as.data.frame(read_excel("DESeq2_result.xlsx"))

head(genelist_input)

genename <- genelist_input$gene_id

gene_map <- select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID")) #将SYMBOL格式的ID换成ENTREZ格式的ID。

non_duplicates_idx <- which(duplicated(gene_map$SYMBOL) == FALSE)

gene_map <- gene_map[non_duplicates_idx, ] 

colnames(gene_map)[1]<-"gene_id" 

temp<-inner_join(gene_map,genelist_input,by = "gene_id")

temp<-temp[,-1]

temp<-na.omit(temp)

temp$log2FoldChange<-sort(temp$log2FoldChange,decreasing = T) 

geneList = temp[,2] 

names(geneList) = as.character(temp[,1])

library(future)
plan("multisession", workers =4) 
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="ALL")   #使用GSEA进行GO富集分析

geneSetID = c('GO:0008239')

gseaNb(object = Go_gseresult,
       geneSetID = geneSetID,
       subPlot = 2,
       termWidth = 45,
       #addPval = T,
       #curveCol = c('#FF162D'), #ES折线图颜色更改
       #curveCol = c('#CD5555'), #ES折线图颜色更改
       #htCol = c("blue", "red"), #热图条颜色更改
       pDigit = 2,
       pvalSize = 5,
       legend.position = c(0.8,0.8),
       #addGene = gene,
       #pvalX = 0.3,pvalY = 0.2
)



####### fig 1e ##################

ToothGrowth <- as.data.frame(read_excel("data.xlsx"))

ToothGrowth$dose<-factor(ToothGrowth$dose)
head(ToothGrowth)

library(Rmisc)
ToothGrowth2 <- summarySE(ToothGrowth, measurevar = "log2(value)",groupvars = c("group","dose"))

ggplot(ToothGrowth,aes(x=group,y=`log2(value)`,fill=group))+
  geom_jitter(size=2 ,aes(fill=group),color="black",width = 0.08,shape=21)+
  facet_wrap(~dose,nrow = 1,ncol = 9)+  
  scale_y_continuous(limits = c(0,15),breaks = seq(0, 15, len = 5)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se,colour = "red"),
                width = 0.1,position = position_dodge(width = 0.6),cex=0.6)+ 
  stat_summary(fun = mean, geom = "point", size=1, color = "red")+
  # geom_boxplot(position=position_dodge(1))+ 
  # change color 
  scale_fill_manual(values = c("black","black"))+

  # geom_signif(comparisons = list(c("Normal","OA")),
  #              map_signif_level=T, 
  #            tip_length=0, 
  #           y_position = 14, 
  #          size=1, textsize = 3, 
  #         test = "t.test")+
  
theme(legend.position = "none",
      axis.title.x = element_blank(), 
      #axis.title.y = element_blank()  , 
      axis.text = element_text(color="black", size=10),
      #axis.line = element_line(linewidth = 1),
      axis.text.y = element_text(color="black",size = 10),
      axis.text.x = element_text(color="black",size = 10,angle = 0,face="plain"),
      axis.title.y = element_text(color="black",size = 10),
      axis.ticks.y = element_line(size=1),
      axis.ticks.x = element_blank(),
      plot.background = element_blank(),
      plot.title = element_text(color="black",size = 20,hjust = 0.5))


########### fig 2b #####################

BRCA_Match_DEG <- as.data.frame(read_excel("deseq2.xlsx"))

BRCA_Match_DEG$log10FDR <- -log10(BRCA_Match_DEG$padj)

colnames(BRCA_Match_DEG)[1] <- "gene_name"

BRCA_Match_DEG <- BRCA_Match_DEG %>% 
  mutate(DEG = case_when(
    logFC > 0.2 & padj < 0.05 ~ "Upregulated",   
    logFC < -0.2 & padj < 0.05 ~ "Downregulated", 
    TRUE ~ "Unchanged"                  
  ))

up <- BRCA_Match_DEG$gene_name[BRCA_Match_DEG$padj < 0.05 & BRCA_Match_DEG$logFC > 0.2]

down <- BRCA_Match_DEG$gene_name[BRCA_Match_DEG$padj < 0.05 & BRCA_Match_DEG$logFC < -0.2]

BRCA_Match_DEG$label <- ifelse(BRCA_Match_DEG$gene_name %in% c("S1PR1","LRP1","THAP2","FBXO36","TP53","CLDN1","IFI27","DMBT1","RTP4","HOXD4","BST2") ,
                               BRCA_Match_DEG$gene_name,"")


library(ggplot2)
library(ggrepel)
library(ggprism)

ggplot(BRCA_Match_DEG, aes(x = logFC, y=log10FDR, colour=DEG)) +
  geom_point(alpha=0.85, size=2) +  
  #geom_point(colour = "black", shape = 21, stroke = 0.5)+
  scale_color_manual(values=c('steelblue','gray','#DE4A4F')) + 
  xlim(c(-4, 4)) +
  geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.5) +
  ### 需要修改
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.5) + 
  labs(x="log2FoldChange", y="-log10pval") +  
  ggtitle("OENC vs OEDPP4") +
  
  theme_prism(border = T,  
              base_size = 12)+ 
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"))+
  geom_label_repel(data = BRCA_Match_DEG, aes(label = label),
                   size = 4,box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.color = "black",
                   show.legend = FALSE, max.overlaps = 10000) 

############  fig 2c ###########################
dat <- read_excel("data_dif1.xlsx")

diff_genes_venn <- dat$gene_id[dat$pval < 0.05 & abs(dat$logFC) > 0.1]

common_genes_MS <- diff_genes_venn

gene.MS <- bitr(common_genes_MS,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

gene_list_MS <- as.list(gene.MS$ENTREZID)

ego_MS <- enrichGO(gene = gene_list_MS,
                   OrgDb = org.Hs.eg.db,
                   ont = "all",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff =1, 
                   qvalueCutoff =1,
                   readable = TRUE)
GO <- ego_MS@result

y <- c(
  "regulation of cellular senescence",
  "aging",
  "cellular senescence",
  "regulation of response to DNA damage stimulus",
  "signal transduction in response to DNA damage",
  "DNA damage checkpoint signaling",
  "cellular response to oxidative stress",
  "response to oxidative stress",
  "inflammatory cell apoptotic process",
  "response to interferon-beta",
  "response to hypoxia"
)

MF <- GO[GO$Description %in% y, ]

res <- as.data.frame(read_excel("go_up_DEGs.xlsx")) 

res <- res %>%
  mutate(geneID = sapply(strsplit(as.character(geneID), "/"), function(x) paste(head(x, 4), collapse = "/")))

enrich <- res %>% 
  group_by(Cluster) %>% 
  top_n(n = nrow(res), wt = -pvalue) %>%   # 选择 Cluster下 n= 几 就有几个词条
  filter(Cluster %in% c("pathway"))

dt <- as.data.frame(enrich)

dt$Description <- str_replace(dt$Description, "^(.)", ~ toupper(.x))

x <- y
x <- str_replace(x, "^(.)", ~ toupper(.x))

dt$Description <- factor(dt$Description, levels = x)

dt <- dt[order(dt$Description), ]

options(scipen = 999)

dt$pvalue <- as.numeric(dt$pvalue)

mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), # 在自定义主题中去掉 y 轴通路标签:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  axis.ticks.y = element_line(color = "black"),  ## 给y轴添加刻度线
  #axis.title.y = element_text(color = "#7DA0C5"), ### 给y轴标题添加颜色 
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)

ggplot(data = dt, aes(x = -log10(pvalue), y = rev(Description),fill = Cluster)) +
  scale_fill_manual(values =c('#D33240')) +  
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0), limits = c(NA, 8)) + 
  labs(x = "-log10(pvalue)", y = "", title = " ") +
  geom_text(size=4, aes(x = 0.05, label = Description), hjust = 0) + 
  geom_text(size=4, aes(x = 0.05, label = geneID), hjust = 0, vjust = 2.5, 
            #color=rep(c('#6bb9d2', '#d55640'),
            color=rep(c('black'),
                      length.out = nrow(dt))) +
  theme_classic() + 
  mytheme +
  NoLegend()


############# fig 2d ######################
RNA =as.data.frame(read_excel("DPP4_RMA_seq.xlsx"))
protein=as.data.frame(read_excel("DPP4_protein.xlsx"))

dim(RNA)
head(RNA)
dim(protein)
head(protein)

combine= merge(RNA,protein,
               by.x="Symbol",
               by.y="Symbol",
               suffixes = c("_RNA","_Protein") ,
               all.x=FALSE,
               all.y=FALSE)

dim(combine)

combine <- as.data.frame(read_excel("RNA_protein.xlsx"))
data <- data.frame(combine[c(1,3,4,5,6,7,8,9)])

data$part <- case_when(abs(data$log2FC_RNA) >= 1 & abs(data$log2FC_Protein) >= 1 ~ "part1379",
                       abs(data$log2FC_RNA) < 1 & abs(data$log2FC_Protein) > 1 ~ "part28",
                       abs(data$log2FC_RNA) > 1 & abs(data$log2FC_Protein) < 1 ~ "part46",
                       abs(data$log2FC_RNA) < 1 & abs(data$log2FC_Protein) < 1 ~ "part5")
head(data)

data$sig <- ifelse(data$Symbol %in% c("DPP4","MMP3","MMP13","COL2A1","COL6A1","ATR","ATF6","FBXO31","STAT3","COL1A1","FN1","COL6A2","GPX4",
                                      "IL1RAP","TGFB2","TP53","CDKN1A","MMP3"), 
                   data$Symbol, "no")

top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2

mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="gray30",size = 12),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.ticks = element_line(size = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))

ggplot(data, aes(x=log2FC_RNA, y=log2FC_Protein)) +
  geom_point(aes(log2FC_RNA,log2FC_Protein),shape = 21, color="grey", size=5,fill = "#FFF0F5")+
  geom_point(data = filter(data, sig != "no"),aes(log2FC_RNA,log2FC_Protein),shape = 21, color="grey", size=5,fill = "red")+
  scale_color_manual(values = c("TRUE" = "#B40000", "FALSE" = "grey")) +
  geom_text_repel(data = filter(data, sig != "no"), aes(label=sig), size=3, 
                  box.padding = 0.35, point.padding = unit(0.3, "lines"),
                  segment.color = 'grey50', max.overlaps = 10000) +  #  文本能自动避让  不加线
  geom_hline(yintercept = c(-0,0),
             size = 0.8,
             color = "grey40",
             lty = "dashed")+
  geom_vline(xintercept = c(-0,0),
             size = 0.8,
             color = "grey40",
             lty = "dashed")+
  labs(title = " ",
       x = "Log2 Fold Change (RNA-seq)",
       y = "Log2 Fold Change (proteomics)")+mytheme


############# fig 2e   ###################
dat_draw <- as.data.frame(read_excel("dat_draw_GSVA_final_mod.xlsx"))

dat_draw$text_color <- ifelse(grepl("Cellular senescense|Response to oxidative stress|Signal transduction in response to DNA damage|
                                    |Cell cycle DNA replication|Extracellular matrix constituent secretion|Regulation of collagen metabolic process", 
                                    dat_draw$pathway, ignore.case = TRUE), "#DE4A4F", "black")

dat_draw$text_color2 <- ifelse(grepl("Cell cycle DNA replication|Extracellular matrix constituent secretion|Regulation of collagen metabolic process", 
                                     dat_draw$pathway, ignore.case = TRUE), "steelblue", "black")

ggplot(dat_draw)+
  geom_col(aes(reorder(pathway, log), log, fill = color))+
  scale_fill_manual(values = c("steelblue","#DE4A4F"))+
  annotate("segment", x = 0, xend = 20, y = 0, yend = 0)+
  theme_classic()+
  ylim(-4,4)+
  coord_flip()+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_blank(),
  )+
  #ylab("-log10(P Value) of GSVA score")+
  geom_text(data = dat_draw[1:13, ], aes(x = pathway, y = -0.1, label = pathway, color = text_color), hjust = 1,family = "Arial" , size = 4) +
  geom_text(data = dat_draw[14:19, ], aes(x = pathway, y = 0.1, label = pathway, color = text_color2), hjust = 0,family = "Arial" , size = 4) +
  ggtitle(" ") +  # title
  scale_x_discrete(expand = expansion(add = c(0, 1.5))) +
  scale_color_identity()



######### fig 7k  #################
library(readxl)
library(ggvenn)


setwd("DPP4/venn_heatmap")
OENC <- read_excel("deseq_oenc_oedpp4.xlsx")
OENC <- OENC[!is.na(OENC$padj), ]

Drug <- read_excel("deseq_oedpp4_drug.xlsx")
Drug <- Drug[!is.na(Drug$padj), ]

OEDPP4vsdrug_up <- Drug$gene_id[Drug$padj < 0.05 & Drug$log2FoldChange > 0.5]

OEDPP4vsdrug_down <- Drug$gene_id[Drug$padj < 0.05 & Drug$log2FoldChange < -0.5]

common_genes <- Reduce(intersect, list(OENCvsOEDPP4_up, OEDPP4vsdrug_down))
common_genes

common_genes1 <- Reduce(intersect, list(OENCvsOEDPP4_down, OEDPP4vsdrug_up))
common_genes1

x_common <- list(`OEDPP4 up` = OENCvsOEDPP4_up, `Drug down` = OEDPP4vsdrug_down )

x_common <- list(`OEDPP4 down` = OENCvsOEDPP4_down, `Drug up` = OEDPP4vsdrug_up )

ggvenn(x_common,
       columns = NULL,
       show_elements = F,
       show_percentage = F,
       fill_color = c("#EEF8FA","#F2E7E7"),
       fill_alpha = 0.5,
       stroke_color = "black" ,
       stroke_alpha = 0,
       stroke_size = 1,  # 不添加 黑框 
       stroke_linetype = "solid",
       set_name_color = "black",
       set_name_size = 10,
       text_color = "black",
       text_size = 20,
       label_sep = ",")


###########    fig 7l      #################
library(ComplexHeatmap)

result <- read_excel("3_oenc_oedpp4_drug.xlsx")
result <- as.data.frame(result)

result <- subset(result, result[,1] != "-")

result <- result[!duplicated(result[, 1]), ]
rownames(result)=result[,1]


result <- result[,-1]

data=log2(result+1)
dat=t(scale(t(data)))
dat[dat>1]=1  
dat[dat<(-1)]= -1
dim(dat)
head(dat)

result <- dat

mito <- c(common_genes, common_genes1)

common_genes <- mito[mito %in% rownames(result)]
selected_expression_data <- result[common_genes, ]

selected_expression_data <- na.omit(selected_expression_data)

my_breaks <- c( -1, 0, 1)

library(RColorBrewer)
my_palette <- colorRampPalette(c("#1E90FF", "white", "#FF0000"))(100)

col_anno <- rep(c("GroupA", "GroupB","GroupC"), each = 3)

col_anno_colors = c("GroupA" = "#8D00FF", "GroupB" = "#8BFD00", "GroupC" = "#00FFFE")

h1<-Heatmap(selected_expression_data, 
                  name = "Expression", 
                  col = structure(my_palette, breaks = my_breaks),
                  column_title = NULL, 
                  width = unit(6, "cm"),
                  height = unit(25,"cm"),
                  show_row_names = F,
                  column_order = colnames(selected_expression_data),  
                  column_gap = unit(0, "mm"),    
                  cluster_rows = T, 
                  rect_gp = gpar(col = "white", lwd = 0),  
                  show_row_dend = F,
                  show_column_dend = F,
                  row_names_side = c("right"), 
                  row_names_rot = 0,   
                  column_title_side = c("top"),
                  show_column_names = F,     
                  column_names_centered = F,    
                  show_heatmap_legend = F ,
                  heatmap_legend_param = list(   
                    at = c(-1, 0, 1),
                    grid_width = unit(3, "mm"),
                    # labels = c("low", "zero", "high"),
                    title = " ",
                    legend_height = unit(2, "cm"),   
                    legend_direction = "vertical"  ,
                    title_position = "lefttop-rot" 
                  ))   

h1

genes <- c("IL32",
           "PTGS2",
           "CXCL2",
           "ICAM1",
           "IL4I1",
           "CCL2",
           "COL1A1",
           "IL1A",
           "CDKN1A",
           "IL1B",
           "ADAM19",
           "TP53",
           "MMP3",
           "MMP14",
           "ADAMTS4",
           "CXCL8",
           "MMP13","AMD1",
           "FBXO5", "RRM1", "GPX4",
           "TK1"
)
genes <- as.data.frame(genes)


a<- h1 + rowAnnotation(link = anno_mark(at = which(rownames(selected_expression_data) %in% genes$genes), 
                                              labels = genes$genes, labels_gp = gpar(fontsize = 10)))

draw(a, align_heatmap_legend = "heatmap_top") 



######## fig 7m and 7n ####################
library(enrichplot)

gene.MS <- bitr(common_genes_MS,
                fromType = "SYMBOL", 
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

gene_list_MS <- as.list(gene.MS$ENTREZID)

ego_MS <- enrichKEGG(gene = diff_entrez$ENTREZID,
                        organism = "hsa", 
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        minGSSize = 10,
                        maxGSSize = 500)

#ENTREZ----symbol：
ego_MS <- setReadable(kk_MS,  
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID")

KEGG <- ego_MS@result

ora_pt <- pairwise_termsim(ego_MS)

set.seed(456)

emapplot(ora_pt
         ,showCategory = y
         ,color = "p.adjust" 
         ,shadowtext = TRUE 
         ,repel = T 
         ,node_label = "category" 
         ,layout.params = list(layout = NULL 
                               ,coords = NULL 
         )

         ,edge.params = list(show = TRUE 
                             , min = 0.2) 
         ,cex.params = list(category_node = 1 
                            , category_label = 1 
                            , line = 1 
         )
         
         ,hilight.params = list(category = NULL
                                ,alpha_hilight = 1
                                ,alpha_no_hilight = 0.3)
         ,cluster.params = list(cluster = T 
                                ,method = stats::kmeans 
                                ,n = 4 
                                , legend = F 
                                ,label_style = "shadowtext"
                                ,label_words_n = 4 
                                ,label_format = 30
         )
)+ 
  scale_fill_continuous(low = "#CD3333", high = "#FF4040", name = "p.adjust")+
  theme(legend.position = "right", legend.direction = "vertical")













