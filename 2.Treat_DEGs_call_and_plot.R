rm(list = ls())

library(scales)
library(ggsci)
library(dplyr)
library(ggpubr)
library(DESeq2) 
library(RColorBrewer)
library(amap)
library(scatterplot3d)
library(pcaMethods)
library(dendextend)
library(ggpubr)
library(ggsci)
library(ggcorrplot)
library(corrplot)
library(pheatmap)
library(gridExtra)  
library(UpSetR)
library(R.matlab)
library(VennDiagram)
library(ggrepel)
#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)#
pal2<-pal_jama("default",alpha = 1)(7)#
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)
pal5 <- pal_npg("nrc", alpha=0.5)(10)
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
show_col(ppCor_all)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]

expression_matrix0 <-read.table(file="/share/home/chenwei/LZJ_SCARE/File1_Wound_Time_final_count.txt",sep="\t", row.names=1,header =T)
head(expression_matrix0);dim(expression_matrix0)

###STEP2 reading meta data
colData <-read.csv(file="/share/home/chenwei/LZJ_SCARE/final_meta_data.csv",sep=",", row.names=1,header =T)

colData$Treat<-factor(colData$Treat,levels = c("CTRL","PBS"))
colData$Day<-factor(colData$Day,levels = c("D3","D7","D20"))
colData$Day_Treat<-factor(colData$Day_Treat,levels = c("D3_CTRL","D3_PBS","D7_CTRL","D7_PBS","D20_CTRL","D20_PBS" ))
head(colData);tail(colData)

###########################################
###################### for day  in special treat 
Treat_tag <-"final";average_count<-3
Treat_group<-"PBS_vs_CTRL"

table(colData$Treat,colData$Day)
##       D3 D7 D20
#   CTRL  3  3   3
#   PBS   3  3   3

expression_matrix1<-expression_matrix0[,colData$Mouse_ID]
dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix1, colData = colData, design = ~ Day_Treat)
#data pre_cleaning
dds1 <- dds0[ rowSums(counts(dds0)) > 1, ]
dds2 <- estimateSizeFactors(dds1)#data transform:sample number < 30 -> rlog
#step seven: count matrix generation
dds3 <- DESeq(dds2)
#log2
normalized_counts_log2 <- log2(counts(dds3, normalized=T)+1)
normalized_counts_mad2 <- apply(normalized_counts_log2, 1, mad)
normalized_counts_log2 <- normalized_counts_log2[order(normalized_counts_mad2, decreasing=T), ]


#################two group intergroup comparison#####################
##expression data save
Day3_target_genes <- c("Mmp3", "Mmp9", "Mmp10", "Mmp12", "Mmp13", "Rasgrf2", "Fgfbp1", "Traf1", "Igf2bp2", "Igfbp2", "Igfbp3", "Igfbp5", "Kras", "Fgfr2", "Fgfr4", "Mras", "Map3k12", "Fgf7", "Map4k2", "Map4k", "Dusp2", "Rasgrp1",  "Mapk13", "Map4k4", "Map3k3", "Map3k20", "Rasgrp3", "Mapkapk3", "Mapk9", "Fgf16", "Mapk11", "Fgfr3")
Day7_target_genes <- c("Mmp3", "Mmp9", "Mmp10", "Mmp12", "Mmp13", "Igf2", "Rasgrp3", "Fgf1", "Fgf7", "Fgfr4", "Igfbp2", "Igfbp3")
Day20_target_genes <- c("Mmp3", "Mmp9", "Mmp10", "Mmp12", "Mmp13", "Igf2bp3", "Rasgrf2", "Mapk7", "Fgfbp1", "Traf1", "Igf2", "Map3k2", "Kras", "Fgfr2", "Fgfr4", "Mras", "Map3k12", "Map4k1",  "Fgf7", "Map4k2", "Mapk1", "Map4k", "Dusp2", "Rasgrp1",  "Mapk13", "Map4k4", "Map3k3", "Map3k20", "Rasgrp3", "Rps6ka5", "Mapkapk3", "Mapk9", "Fgf1", "Mapk11", "Fgfr3" )
target_genes_list<-c(list(Day3_target_genes),list(Day7_target_genes),list(Day20_target_genes))
names(target_genes_list)<-c("D3","D7","D20")

Treat_select<-c("PBS","CTRL")
ann_colors_list = list(Treat=c(CTRL=ppCor[2],PBS=ppCor[1]),Class=c(Up_DEGs=ppCor[3],Down_DEGs=ppCor[4]))
for(Day_tag in c("D3","D7","D20")){
  #Day_tag <-"D20"
  colData_day<-colData[which(colData$Day == Day_tag),c("Mouse_ID","Treat")]
  head(colData_day)
  expression_matrix1<-expression_matrix0[,colData_day$Mouse_ID]
  evaluation.table<-expression_matrix1
  evaluation.table[evaluation.table !=0]<- 1
  gene_detected_number<-colSums(evaluation.table)
  gene_detected_number<-data.frame(gene_detected_number)
  colnames(gene_detected_number) <- "gene_number"
  
  gene_detected_number2<-merge(gene_detected_number,colData_day,by=0)
  gene_detected_number2<-gene_detected_number2[order(gene_detected_number2$Treat,gene_detected_number2$gene_number,decreasing = F),]
  gene_detected_number2$Row.names<-factor(gene_detected_number2$Row.names,levels=as.character(gene_detected_number2$Row.names))
  head(gene_detected_number2)
  plot_gene_detected_number<-ggplot(data=gene_detected_number2, mapping=aes(x=Row.names,y=gene_number))+scale_fill_manual(values=ppCor_all)+
    geom_text(stat="identity",aes(label=gene_number), color="black", size=2,position=position_stack(1.05))+
    theme_classic()+labs(x="",y="Number of genes",title=paste0("Number of gene detected with average count:",average_count))+guides(fill=guide_legend(ncol=1)) +
    theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
          axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
  count_plot<-plot_gene_detected_number+geom_bar(aes(fill=Treat),stat="identity",width=0.8)
  count_plot
  ggsave(file=paste0("/share/home/chenwei/LZJ_SCARE/intergroup/",Treat_group,"_",Day_tag,"_",Treat_tag,"_in_RNA_library.pdf"),count_plot,width = 5, height =5)
  
  #evaluation global effect
  head(expression_matrix1)
  head(colData_day)
  dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix1, colData = colData_day, design = ~ Treat)
  #data pre_cleaning
  dds1 <- dds0[ rowSums(counts(dds0)) > 1, ]
  dds2 <- estimateSizeFactors(dds1)
  dds3 <- DESeq(dds2)
  ###################
  rld <- rlog(dds2, blind = FALSE);vsd <- vst(dds2, blind = FALSE)
  #Visual sample: principal component analysis.
  plotPCA(rld,intgroup=c("Treat"))
  rld_PCA<-plotPCA(rld,intgroup=c("Treat"),returnData = TRUE)
  percentVar<-round(100*attr(rld_PCA,"percentVar"),1)
  plot_rld_PCA0<-ggplot(data=rld_PCA, mapping=aes(x=PC1,y=PC2,colour = Treat))+
    geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
    xlab(paste0("PC1: ",percentVar[1],"% variance"))+ ylab(paste0("PC2: ",percentVar[2],"% variance"))+
    labs(title ="Distribution of samples \n rld Deseq2 inner method")+scale_color_manual(values=c(ppCor))+
    theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                     axis.title = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),legend.position = "right")
  plot_rld_PCA0
  plot_rld_PCA1<-plot_rld_PCA0 + geom_text(aes(label=name), color="black", size=4,angle = 0)
  plot_rld_PCA1
  
  plot_PCA13<-grid.arrange(plot_rld_PCA0, plot_rld_PCA1, ncol = 2)  
  ggsave(file=paste0("/share/home/chenwei/LZJ_SCARE/intergroup/PCA_rld_Treat_text_",Treat_group,"_",Day_tag,"_",Treat_tag,".pdf"),plot_PCA13,width = 10, height =5)
  
  
  #样本距离 层次聚类,相关性分析 和主成分分析： RNA-Seq分析第一步通常是评估样本间的总体相似度
  data_normal<-log2(counts(dds2,normalized=TRUE)+1)
  anno1 <- as.data.frame(colData(dds2)[, c("Treat","Treat")])
  tag<-"normal"
  #cluster 
  dists <- dist(t(data_normal),method = "euclidean") 
  hc <- hcluster(dists, method="pearson")
  row_dend <- as.dendrogram(hc)
  colData_day
  lable_cor<-rep(ppCor[2],nrow(colData_day))
  lable_cor[which(hc$order %in% c(1:as.numeric(table(colData_day$Treat))[1]))] <- ppCor[1]
  labels_colors(row_dend) <- lable_cor
  pdf(file=paste0("/share/home/chenwei/LZJ_SCARE/intergroup/pearson_cluster_Treat_",Treat_group,"_",Day_tag,"_",Treat_tag,".pdf"),width = 10, height =6)
  plot(row_dend, main = paste0("samples similarity \n ",tag,"(pearson)"))
  dev.off()
  
  ##distance
  sampleDists <- dist(t(data_normal))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
  heat_plot<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main =paste0(tag,"_trans"))
  pdf(file=paste0("/share/home/chenwei/LZJ_SCARE/intergroup/distance_heatmap_",Treat_group,"_",Day_tag,"_",Treat_tag,".pdf"),width = 10, height =9)
  print(heat_plot)
  dev.off()
  
  #perform PCA analysis 
  PCA_data1 <- t(data_normal)
  set.seed(19921010)
  md <- prep(PCA_data1, scale="none", center=TRUE)
  resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 5)
  summary(resPPCA) 
  variation<-data.frame(summary(resPPCA))
  df <-as.data.frame(scores(resPPCA))
  PCA_data2<-merge(colData_day,df,by=0)
  
  plot_PCA0<-ggplot(data=PCA_data2, mapping=aes(x=PC1,y=PC2))+
    xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
    ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
    labs(title =paste0("Distribution of samples \n ",tag))+scale_color_manual(values=c(ppCor))+
    theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),panel.grid=element_blank(),plot.title=element_text(hjust=0.5),axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                     axis.text.y = element_text(size=15),legend.position = "right")
  
  plot_PCA1<-plot_PCA0+geom_point(aes(colour =Treat),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
  plot_PCA2<-plot_PCA1+geom_text(stat="identity",aes(label=Row.names), color="black", size=4)
  plot_PCA3<-grid.arrange(plot_PCA1, plot_PCA2, ncol = 2)  
  ggsave(paste0("/share/home/chenwei/LZJ_SCARE/intergroup/PCA_1_2",Treat_group,"_",Day_tag,"_",Treat_tag,".pdf"),plot_PCA3,width=10, height=5)
  
  plot_PCA0<-ggplot(data=PCA_data2, mapping=aes(x=PC2,y=PC3))+
    xlab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
    ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
    labs(title =paste0("Distribution of samples \n ",tag))+scale_color_manual(values=c(ppCor))+
    theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),panel.grid=element_blank(),plot.title=element_text(hjust=0.5),axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                     axis.text.y = element_text(size=15),legend.position = "right")
  plot_PCA1<-plot_PCA0+geom_point(aes(colour =Treat),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
  plot_PCA2<-plot_PCA1+geom_text(stat="identity",aes(label=Row.names), color="black", size=4)
  plot_PCA23<-grid.arrange(plot_PCA1, plot_PCA2, ncol = 2)  
  ggsave(paste0("/share/home/chenwei/LZJ_SCARE/intergroup/PCA_2_3",Treat_group,"_",Day_tag,"_",Treat_tag,".pdf"),plot_PCA23,width=10, height=5)
  
  plot_PCA0<-ggplot(data=PCA_data2, mapping=aes(x=PC1,y=PC3))+
    xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
    ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
    labs(title =paste0("Distribution of samples \n ",tag))+scale_color_manual(values=c(ppCor))+
    theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),panel.grid=element_blank(),plot.title=element_text(hjust=0.5),axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                     axis.text.y = element_text(size=15),legend.position = "right")
  plot_PCA1<-plot_PCA0+geom_point(aes(colour =Treat),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
  plot_PCA2<-plot_PCA1+geom_text(stat="identity",aes(label=Row.names), color="black", size=4)
  plot_PCA13<-grid.arrange(plot_PCA1, plot_PCA2,ncol = 2)  
  ggsave(paste0("/share/home/chenwei/LZJ_SCARE/intergroup/PCA_1_3",Treat_group,"_",Day_tag,"_",Treat_tag,".pdf"),plot_PCA13,width=10, height=5)
  
  ##############
  ###call DEGs
  sampleA <-Treat_select[1]; sampleB<-Treat_select[2]
  contrastV <- c("Treat",sampleA, sampleB)
  res1 <- results(dds3, contrast=contrastV) 
  summary(res1)
  
  ##给DESeq2的原始输出结果增加样品平均表达信息，使得结果更容易理解和解析。
  # 获得第一组数据均值
  baseA <- counts(dds3, normalized=TRUE)[, colData(dds3)$Treat == sampleA]
  if (is.vector(baseA)){
    baseMeanA <- as.data.frame(baseA)
  } else {
    baseMeanA <- as.data.frame(rowMeans(baseA,na.rm = T))
  }
  colnames(baseMeanA) <-paste0(sampleA,"_mean")
  head(baseMeanA) 
  
  # 获得第二组数据均值
  baseB <- counts(dds3, normalized=TRUE)[, colData(dds3)$Treat == sampleB]
  if (is.vector(baseB)){
    baseMeanB <- as.data.frame(baseB)
  } else {
    baseMeanB <- as.data.frame(rowMeans(baseB,na.rm = T))
  }
  colnames(baseMeanB) <-paste0(sampleB,"_mean")
  head(baseMeanB)
  # 结果组合
  res2 <- cbind(baseMeanA, baseMeanB, as.data.frame(res1));head(res2) 
  res3 <- cbind(ID=rownames(res2), as.data.frame(res2));head(res3)
  
  # 校正后p-value为NA的复制为1
  res3$padj[is.na(res3$padj)] <- 1;res3$pvalue[is.na(res3$pvalue)] <- 1
  # 按pvalue排序, 把差异大的基因放前面
  res3 <- res3[order(res3$pvalue,decreasing = F),]
  head(res3);tail(res3) ;dim(res3)
  #whole RNA_data combined
  head(as.data.frame(res3));dim(as.data.frame(res3))#22231     9
  head(normalized_counts_log2);dim(normalized_counts_log2)
  
  merge_data <- merge(normalized_counts_log2,as.data.frame(res3),by="row.names",sort=FALSE)
  merge_data <- merge_data[order(merge_data$pvalue,decreasing = F),]
  rownames(merge_data)<- merge_data$Row.names
  merge_data<-merge_data[,-1]
  file <- paste0("/share/home/chenwei/LZJ_SCARE/intergroup/file2_",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
  write.table(merge_data, file=file, sep="\t", quote=F, row.names=F,col.names=T)
  
  table(merge_data$pvalue<0.05)  
  table(merge_data$pvalue<0.01) 
  table(merge_data$pvalue<0.05 & abs(merge_data$log2FoldChange)>=0.585)   
  table(merge_data$padj<0.1)  
  table(merge_data$padj<0.05)   
  table(merge_data$padj<0.05 & abs(merge_data$log2FoldChange)>=0.585)
  
  #差异基因的进一步筛选
  # p_value<0.05 & Foldchange >=1.5
  res_de_up <- subset(merge_data, pvalue<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
  res_de_dw <- subset(merge_data, pvalue<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
  nrow(res_de_up);nrow(res_de_dw)
  DEG_information<-rbind(res_de_up,res_de_dw)
  file <- paste0("/share/home/chenwei/LZJ_SCARE/intergroup/file3_",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
  write.table(as.data.frame(DEG_information), file=file, sep="\t", quote=F, row.names=T, col.names=T)
  
  #plot 
  head(merge_data) 
  merge_data$log10_pvalue<-(-log10(merge_data$pvalue))
  merge_data$threshold = as.factor(ifelse(merge_data$pvalue < 0.05 & abs(merge_data$log2FoldChange)>= 0.585, ifelse(merge_data$log2FoldChange>= 0.585 ,'Up','Down'),'NoSignifi'))
  table(merge_data$threshold)
  merge_data$threshold<-factor(merge_data$threshold,level=c("Up","Down","NoSignifi"))
  range(merge_data$log2FoldChange)#-5.242734  6.042089
  range( merge_data$log10_pvalue)# 0.000000 9.492076
  # merge_data[which(merge_data$log10_pvalue > 30),]$log10_pvalue <-30
  plot_vocano<-ggplot(data = merge_data, aes(x = log2FoldChange, y = log10_pvalue, colour=threshold)) +
    geom_point(alpha=0.4, size=3.5) +
    scale_color_manual(values=c("red","blue","grey")) +
    scale_x_continuous(limits=c(-ceiling(max(abs(merge_data$log2FoldChange))),ceiling(max(abs(merge_data$log2FoldChange)))),
                       breaks=seq(-ceiling(max(abs(merge_data$log2FoldChange))),ceiling(max(abs(merge_data$log2FoldChange))), 2))+
    #xlim(-ceiling(max(abs(merge_data$log2FoldChange)))-1,ceiling(max(abs(merge_data$log2FoldChange)))+1)+
    ylim(0,ceiling(max(merge_data$log10_pvalue)))+
    geom_vline(xintercept=c(-0.585,0.585),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x=paste0("log2(Foldchange (",sampleA,"/",sampleB,") in gene expression level"),y="-log10 (p value)",
         title=paste(Treat_group,"_",Day_tag,"_",sampleA,"_vs._",sampleB,":Up_DEGs: ",length(which(merge_data$threshold == "Up"))," & Down_DEG: ",length(which(merge_data$threshold == "Down"))," [p value<0.05 & FoldChange > 1.5]",sep="")) +
    theme_bw()+
    theme(panel.border = element_rect(colour="black",fill=NA),panel.grid.major =element_blank())+
    theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5), legend.position="right", legend.title = element_blank(),axis.line = element_line(colour="black"))+
    theme(axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
  plot_vocano
  ggsave(plot_vocano,file=paste0("/share/home/chenwei/LZJ_SCARE/intergroup/",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,"_vocano_DEGs.pdf"), width = 8, height = 8)
  ggsave(plot_vocano,file=paste0("/share/home/chenwei/LZJ_SCARE/intergroup/",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,"_vocano_DEGs.png"), width = 8, height = 8)
  
  plot_vocano2<-plot_vocano+theme(legend.position='none')+labs(x="",y="",title="")+
    theme(plot.title=element_blank(),axis.text.y = element_blank(), axis.text.x = element_blank(),axis.text = element_blank(),legend.text = element_blank(),legend.title = element_blank(),axis.title=element_blank())
  plot_vocano2
  ggsave(plot_vocano2,file=paste0("/share/home/chenwei/LZJ_SCARE/intergroup/",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,"_vocano_DEGs_notext.pdf"), width = 8, height = 8)
  ggsave(plot_vocano2,file=paste0("/share/home/chenwei/LZJ_SCARE/intergroup/",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,"_vocano_DEGs_notext.png"), width = 8, height = 8)
  
  # 添加目标基因标签
  target_genes<-as.character(unlist(target_genes_list[Day_tag]))
  merge_data$GENE_SYMBOL<-rownames(merge_data)
  merge_data$is_target <- merge_data$GENE_SYMBOL %in% target_genes
  
  Targetsig<-merge_data[merge_data$is_target,]
  dim(Targetsig)
  Target_down<-Targetsig[which(Targetsig$threshold== "Down"),]
  Target_Up<-Targetsig[which(Targetsig$threshold== "Up"),]
  
  plot_vocano3<-plot_vocano+geom_label_repel(data=Target_Up,aes(log2FoldChange,log10_pvalue,label=GENE_SYMBOL),max.overlaps =100,
                               nudge_x = 1.5,nudge_y = -0.8,color = "white",alpha = 0.9,point.padding = 0.5, size = 2.5,
                               fill = "#96C93D", segment.size = 0.5,segment.color = "grey50",direction = "y", hjust = 0.5) +
    geom_label_repel(data = Target_down,aes(log2FoldChange,log10_pvalue,label=GENE_SYMBOL),
                     nudge_x = -3, nudge_y = 0.2,color = "white",alpha = 0.9,point.padding = 0.5, size = 2.5,
                     fill = "#9881F5",segment.size = 0.5, segment.color = "grey50",direction = "y",hjust = 0.5)#+  scale_x_continuous(limits = NULL)
  plot_vocano3
  
  ggsave(plot_vocano3,file=paste0("/share/home/chenwei/LZJ_SCARE/intergroup/",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,"_vocano_DEGs_lable.pdf"), width = 8, height = 8)
  ggsave(plot_vocano3,file=paste0("/share/home/chenwei/LZJ_SCARE/intergroup/",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,"_vocano_DEGs_lable.png"), width = 8, height = 8)
  
  #绘制各组特异高表达基因热图
 
  res_de_up <- subset(merge_data, pvalue<0.05&log2FoldChange>=0.585 )
  res_de_dw <- subset(merge_data, pvalue<0.05&log2FoldChange<=(-1)*0.585 )
  de_id_whole = rbind(res_de_up, res_de_dw)
  nrow(de_id_whole);head(de_id_whole)
  
  ##select target Day sample
  A_list<-colData_day[which(colData_day$Treat == sampleA),]$Mouse_ID
  B_list<-colData_day[which(colData_day$Treat == sampleB),]$Mouse_ID 
  
  mat<-de_id_whole[,c(A_list,B_list)]
  #anno1 <- data.frame(sample=colnames(mat))
  #anno1$Day<-ifelse(grepl("MCDAN",anno1$sample),"MCDAN","sIUGR")
  #rownames(anno1)<-anno1$sample
  annotation_col <-data.frame(Treat=colData_day[c(A_list,B_list),]$Treat)
  rownames(annotation_col) = colnames(mat)
  
  annotation_row = data.frame(Class = factor(rep(c("Up_DEGs","Down_DEGs"),c(nrow(res_de_up),nrow(res_de_dw)))))
  rownames(annotation_row) = rownames(mat)
  
  ann_colors =ann_colors_list
  range(mat)# 0.00000 18.32706
  bk <- c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01)) 
  # 做热图： 
  heat_plot<-pheatmap(mat, scale = "row", 
                      # color = colorRampPalette(colors = c("blue","white","red"))(100),
                      color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), 
                      cluster_row =FALSE,cluster_col =FALSE,                    
                      legend_breaks=seq(-3,3,1),  breaks=bk,
                      annotation_col = annotation_col,annotation_row=annotation_row,
                      annotation_colors = ann_colors, 
                      gaps_row = nrow(res_de_up),gaps_col =length(A_list),cutree_col = 2,
                      show_rownames = F,
                      main=paste0(Day_tag,"_",sampleA,"_vs_",sampleB,"_DEGs (pvalue <0.05& FoldChange>=1.5"))
  pdf(file=paste0("/share/home/chenwei/LZJ_SCARE/intergroup/heatmap_",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,"_DEGs_p005_FC1.5.pdf"),width = 7, height =7)
  print(heat_plot)
  dev.off()
}

##相互关系探讨
#selected DEGs
# p_value<0.05 & Foldchange >=1.5
DEG_list_up<-list();DEG_list_dw<-list()
for(Day_tag in c("D3","D7","D20")){
  #Day_tag <-"D14"
  file_DEG <- paste0("/share/home/chenwei/LZJ_SCARE/intergroup/file3_",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
  DEG_file<-read.table(file=file_DEG,sep="\t", row.names=1,header =T)
  
  DEG_up <-  rownames(DEG_file[which(DEG_file$log2FoldChange> 0),])
  DEG_dw <-  rownames(DEG_file[which(DEG_file$log2FoldChange< 0),])
  DEG_list_up<-c(DEG_list_up,list(DEG_up))
  DEG_list_dw<-c(DEG_list_dw,list(DEG_dw))
}

#for THREE Day 
tag<- "Three list of Up DEGs:pvalue0.05 & FC1.5"
mainname<-tag
venn1 <-venn.diagram(list(D7_up=DEG_list_up[[1]],D14_up= DEG_list_up[[2]],D30_up=DEG_list_up[[3]]),
                     alpha=c(0.5,0.5,0.5),lwd=1,lty=1,col="white",fill=ppCor[c(5:6,10)], 
                     cex = 1.5,cat.col=ppCor[c(5:6,10)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                     cat.fontface=4, cat.cex = 1.5, main=mainname, 
                     main.cex = 2,  filename = NULL)

tag<- "Three list of Down DEGs:pvalue0.05 & FC1.5"
mainname<-tag
venn2 <-venn.diagram(list(D7_dw=DEG_list_dw[[1]],D14_dw= DEG_list_dw[[2]],D30_dw=DEG_list_dw[[3]]),
                     alpha=c(0.5,0.5,0.5),
                     lwd=1,lty=1,col="white",fill=ppCor[c(5:6,10)], 
                     cex = 1.5,cat.col=ppCor[c(5:6,10)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                     cat.fontface=4, cat.cex = 1.5, main=mainname, 
                     main.cex = 2, 
                     # main.fontface = 2, main.fontfamily = 3, 
                     filename = NULL)

pdf(paste0("/share/home/chenwei/LZJ_SCARE/intergroup/heatmap_",Treat_group,"_",Treat_tag,"_venn_three_lists_DEGs.pdf"),width = 8,height=8)
grid.newpage(); 
grid.draw(venn1)
grid.newpage();
grid.draw(venn2)
dev.off()

#upset plot for all eight lists of DEGs
listinput<-c(DEG_list_up,DEG_list_dw)
names(listinput)<-c("D3_Up","D7_Up","D20_Up","D3_Down","D7_Down","D20_Down")

upset_dataframe<-as.data.frame(fromList(listinput))
dim(upset_dataframe)# 3939    6
upset_dataframe[1:10,1:5]
upset_dataframe_rowSums<-rowSums(upset_dataframe)
range(upset_dataframe_rowSums)#1:3
upset_dataframe_colSums<-colSums(upset_dataframe)
range(upset_dataframe_colSums)#359 1965
upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]
upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400)

pdf(paste0("/share/home/chenwei/LZJ_SCARE/intergroup/",Treat_group,"_",Treat_tag,"_two_group_upsetplot_six_lists_interDay_DEGs.pdf"),width = 10,height=8)

upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400,
      sets=c("D3_Up","D7_Up","D20_Up","D3_Down","D7_Down","D20_Down"),
      keep.order = TRUE,
      query.legend = "top",
      sets.bar.color = "brown",
      #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
      shade.color="pink",
      #      matrix.color="purple",
      order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
      point.size = 3,line.size = 1.3,
      mainbar.y.label = "DMRs number Intersections", sets.x.label = "DMRs number per subset",
      mb.ratio = c(0.60, 0.40),
      text.scale = c(1.5, 1.5,1.2,1.5,1.5,1),
      show.numbers = 'yes'#,=queries = list(list2_q1,list2_q2,list2_q3,list2_q4,list2_q5,list2_q6)
)
dev.off()

