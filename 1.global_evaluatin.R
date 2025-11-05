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

##############
##############
Expression_matrix<-read.table(file="/share/home/chenwei/LZJ_SCARE/expression_matrix.txt",sep="\t", row.names=1,header =T)
head(Expression_matrix);dim(Expression_matrix)

count.table2<-Expression_matrix[,c("gene_symbol",paste0("read_count_CON_",c("d3_1","d3_2","d3_3","d7_1","d7_2","d7_3","d20_1","d20_2","d20_3")), 
                                   paste0("read_count_PBS_",c("d3_1","d3_2","d3_3","d7_1","d7_2","d7_3","d20_1","d20_2","d20_3")))]
rownames(count.table2)<-count.table2$gene_symbol;count.table2<-count.table2[,-1]
colnames(count.table2)<-c(paste0("CON_",c("d3_1","d3_2","d3_3","d7_1","d7_2","d7_3","d20_1","d20_2","d20_3")),paste0("PBS_",c("d3_1","d3_2","d3_3","d7_1","d7_2","d7_3","d20_1","d20_2","d20_3")))
head(count.table2);dim(count.table2)


###STEP2 reading meta data
colData<-data.frame(Mouse_ID=c(paste0("CON_",c("d3_1","d3_2","d3_3","d7_1","d7_2","d7_3","d20_1","d20_2","d20_3")), paste0("PBS_",c("d3_1","d3_2","d3_3","d7_1","d7_2","d7_3","d20_1","d20_2","d20_3"))),
                    Treat=c(rep("CTRL",9),rep("PBS",9)),
                    Day=c(rep(c("D3","D7","D20"),each=3),rep(c("D3","D7","D20"),each=3)))
rownames(colData)<-colData$Mouse_ID

colData$Treat<-factor(colData$Treat,levels = c("CTRL","PBS"))
colData$Day<-factor(colData$Day,levels = c("D3","D7","D20"))
colData<-colData[order(colData$Day,colData$Treat,decreasing = F),]
colData$Day_Treat<-paste0(colData$Day,"_",colData$Treat)
colData$Day_Treat<-factor(colData$Day_Treat,levels = c("D3_CTRL","D3_PBS","D7_CTRL","D7_PBS","D20_CTRL","D20_PBS" ))
head(colData);tail(colData)
write.table(colData,file="/share/home/chenwei/LZJ_SCARE/final_meta_data.csv",sep=",", quote=F, row.names=T,col.names=T)

tag_name<-"All_sample" #test line
expression_matrix0<-count.table2[,rownames(colData)]
# 转换单元格值为整数
expression_matrix0[,] <- lapply(expression_matrix0[,], function(x) as.integer(x))
all(is.integer(expression_matrix0))
head(expression_matrix0);dim(expression_matrix0)#  18713    18


##考虑消除低表达量基因的影响
average_count <-3
print(paste0("gene with expressioncount: ",average_count))
expression_matrix<-expression_matrix0[which(rowSums(expression_matrix0 > 0)>= average_count),];dim(expression_matrix)#17679    18
#expression_matrix<-expression_matrix0[which(rowSums(expression_matrix0)>= average_count * ncol(expression_matrix0)),]
write.table(expression_matrix,file="/share/home/chenwei/LZJ_SCARE/File1_Wound_Time_final_count.txt",sep="\t", quote=F, row.names=T,col.names=T)

#step two : evaluation gene detective number
evaluation.table<-expression_matrix
evaluation.table[evaluation.table !=0]<- 1
gene_detected_number<-colSums(evaluation.table)
gene_detected_number<-data.frame(gene_detected_number)
colnames(gene_detected_number) <- "gene_number"

gene_detected_number2<-merge(gene_detected_number,colData,by=0)
gene_detected_number2<-gene_detected_number2[order(gene_detected_number2$Day_Treat,gene_detected_number2$gene_number,decreasing = F),]
gene_detected_number2$Row.names<-factor(gene_detected_number2$Row.names,levels=as.character(gene_detected_number2$Row.names))
head(gene_detected_number2)
plot_gene_detected_number<-ggplot(data=gene_detected_number2, mapping=aes(x=Row.names,y=gene_number))+scale_fill_manual(values=ppCor_all)+
  geom_text(stat="identity",aes(label=gene_number), color="black", size=2,position=position_stack(1.05))+
  theme_classic()+labs(x="",y="Number of genes",title=paste0("Number of gene detected in samples:",average_count))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))

plot_gene_detected_number1<-plot_gene_detected_number+geom_bar(aes(fill=Day_Treat),stat="identity",width=0.8)
plot_gene_detected_number2<-plot_gene_detected_number+geom_bar(aes(fill=Day),stat="identity",width=0.8)
plot_gene_detected_number3<-plot_gene_detected_number+geom_bar(aes(fill=Treat),stat="identity",width=0.8)
count_plot<-ggarrange(plot_gene_detected_number1, plot_gene_detected_number2, plot_gene_detected_number3,labels = c("A", "B", "C"),ncol = 1, nrow = 3)
count_plot
ggsave(file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/figure1_",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,"_in_RNA_library.pdf"),count_plot,width = 7, height =9)

#evaluation global effect
head(expression_matrix)
head(colData)

dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix, colData = colData, design = ~ Day+Treat)
#data pre_cleaning
dds1 <- dds0[ rowSums(counts(dds0)) > 1, ]
dds2 <- estimateSizeFactors(dds1)#data transform:sample number < 30 -> rlog
#step three : count matrix generation
dds3 <- DESeq(dds2)

#show transform result
rld <- rlog(dds2)
vsd <- vst(dds2)
df <- bind_rows(
  as.data.frame(log2(counts(dds2, normalized=TRUE)[, 2:3]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as.data.frame(assay(rld)[, 2:3]) %>% mutate(transformation = "rlog"),
  as.data.frame(assay(vsd)[, 2:3]) %>% mutate(transformation = "vst"))
colnames(df)[1:2] <- c("sample2","sample3")  
trans_plot<-ggplot(df, aes(x = sample2, y = sample3)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)
trans_plot
ggsave(file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,"_trans_three_type_in_RNA_library.pdf"),trans_plot,width = 15, height =5)

##the normalized counts (normalized for library size, which is the total number of gene counts per sample, while accounting for library composition)
normalized_counts <- counts(dds3, normalized=TRUE)
#根据基因在不同的样本中表达变化的差异程度mad值对数据排序，差异越大的基因排位越前。
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.table(normalized_counts,file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,"_deseq2_normalized_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)

#log2
normalized_counts_log2 <- log2(counts(dds3, normalized=T)+1)
normalized_counts_mad2 <- apply(normalized_counts_log2, 1, mad)
normalized_counts_log2 <- normalized_counts_log2[order(normalized_counts_mad2, decreasing=T), ]
write.table(normalized_counts,file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,"_deseq2_log2_normalized_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)

# other log transform
rlog_dds3_Mat <- assay(rld)
rlog_dds_mad3 <- apply(rlog_dds3_Mat, 1, mad)
rlog_dds3_Mat <- rlog_dds3_Mat[order(rlog_dds_mad3, decreasing=T), ]
write.table(normalized_counts,file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,"_deseq2_rlog_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)
# 只在Linux下能运行，目的是填补表格左上角内容 #system(paste("sed -i '1 s/^/ID\t/'", "ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.xls"))

#Visual sample: principal component analysis.
plotPCA(rld,intgroup=c("Day_Treat"))
rld_PCA<-plotPCA(rld,intgroup=c("Day_Treat"),returnData = TRUE)
percentVar<-round(100*attr(rld_PCA,"percentVar"),1)
#rld_PCA$Day_Treat<-factor(rld_PCA$Day_Treat,levels=c("D7_sham","D7_saline","D7_hgf","D7_luc","D14_sham","D14_saline","D14_hgf","D14_luc","D30_sham","D30_saline","D30_hgf","D30_luc"))

plot_rld_PCA0<-ggplot(data=rld_PCA, mapping=aes(x=PC1,y=PC2,colour = Day_Treat))+
  geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+ ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  labs(title ="Distribution of samples \n rld Deseq2 inner method")+scale_color_manual(values=c(ppCor))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),legend.position = "right")
plot_rld_PCA0
ggsave(file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/PCA_rld_Day_Treat_",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),plot_rld_PCA0,width = 10, height =10)
plot_rld_PCA1<-plot_rld_PCA0 + geom_text(aes(label=name), color="black", size=4,angle = 0)
plot_rld_PCA1
ggsave(file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/PCA_rld_Day_Treat_text_",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),plot_rld_PCA1,width = 10, height =10)

plotPCA(rld,intgroup=c("Day"))
rld_PCA<-plotPCA(rld,intgroup=c("Day"),returnData = TRUE)
percentVar<-round(100*attr(rld_PCA,"percentVar"),1)
rld_PCA$Day<-factor(rld_PCA$Day,levels=c("D3","D7","D20"))
plot_rld_PCA2<-ggplot(data=rld_PCA, mapping=aes(x=PC1,y=PC2,colour = Day))+
  geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+ ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  labs(title ="Distribution of samples \n rld Deseq2 inner method")+scale_color_manual(values=c(ppCor))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),legend.position = "right")
plot_rld_PCA2
ggsave(file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/PCA_rld_Day_",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),plot_rld_PCA2,width = 10, height =10)

plotPCA(rld,intgroup=c("Treat"))
rld_PCA<-plotPCA(rld,intgroup=c("Treat"),returnData = TRUE)
percentVar<-round(100*attr(rld_PCA,"percentVar"),1)
rld_PCA$Treat<-factor(rld_PCA$Treat,levels=c("CTRL","PBS"))
plot_rld_PCA3<-ggplot(data=rld_PCA, mapping=aes(x=PC1,y=PC2,colour = Treat))+
  geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+ ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  labs(title ="Distribution of samples \n rld Deseq2 inner method")+scale_color_manual(values=c(ppCor))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),legend.position = "right")
plot_rld_PCA3
ggsave(file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/PCA_rld_Day_",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),plot_rld_PCA3,width = 10, height =10)

#样本距离 层次聚类,相关性分析 和主成分分析： RNA-Seq分析第一步通常是评估样本间的总体相似度
data_not_normal<-log2(counts(dds2,normalized=F)+1)
data_normal<-log2(counts(dds2,normalized=TRUE)+1)
data_rld<-assay(rld)
data_vsd<-assay(vsd)
data_list<-c(list(data_not_normal),list(data_normal),list(data_rld),list(data_vsd))

data_name<-c("not_normal","normal","rld","vsd")
anno1 <- as.data.frame(colData(dds2)[, c("Day","Treat","Day_Treat")])

for(i in 2:3){
#i=3
tag<-data_name[i]
print(tag)
data_raw_bin <-data_list[[i]]

#cluster 
dists <- dist(t(data_raw_bin),method = "euclidean") 
hc <- hcluster(dists, method="pearson")
row_dend <- as.dendrogram(hc)

lable_cor<-rep(ppCor[12],18)
lable_cor[which(hc$order %in% c(1:3))] <- ppCor[1]
lable_cor[which(hc$order %in% c(4:6))] <- ppCor[2]
lable_cor[which(hc$order %in% c(7:9))] <- ppCor[3]
lable_cor[which(hc$order %in% c(10:12))] <- ppCor[4]
lable_cor[which(hc$order %in% c(13:15))] <- ppCor[5]

labels_colors(row_dend) <- lable_cor

pdf(file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/pearson_cluster_Treat_Day_",tag,"_",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),width = 10, height =6)
plot(row_dend, main = paste0("samples similarity \n ",tag,"(pearson)"))
dev.off()

lable_cor<-rep(ppCor[2],18)
lable_cor[which(hc$order %in% c(1:3,7:9,13:15))] <- ppCor[1]
labels_colors(row_dend) <- lable_cor

pdf(file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/pearson_cluster_Day_",tag,"_",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),width = 10, height =6)
plot(row_dend, main = paste0("samples similarity \n ",tag,"(pearson)"))
dev.off()

lable_cor<-rep(ppCor[3],18)
lable_cor[which(hc$order %in% c(1:6))] <- ppCor[1]
lable_cor[which(hc$order %in% c(7:12))] <- ppCor[2]
labels_colors(row_dend) <- lable_cor
pdf(file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/pearson_cluster_Treat_",tag,"_",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),width = 10, height =6)
plot(row_dend, main = paste0("samples similarity \n ",tag,"(pearson)"))
dev.off()

##distance
sampleDists <- dist(t(data_raw_bin))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
heat_plot<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main =paste0(tag,"_trans"))
pdf(file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/distance_heatmap_",tag,"_",tag_name,"_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),width = 10, height =9)
print(heat_plot)
dev.off()

## correlationship analysis
gene_exp_corr<-cor(data_raw_bin,method="pearson",use = "p") #默认是pearson法，method="spearman"
gene_exp_p.mat <- cor_pmat(data_raw_bin,method="pearson",use = "p") # 还可以计算下p值
sample_names<-colnames(gene_exp_corr)
gene_exp_corr2<-gene_exp_corr[sample_names,sample_names]
gene_exp_p.mat2<-gene_exp_p.mat[sample_names,sample_names]
write.table(gene_exp_corr2,file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/R_pearson_correlationship_data_",tag,"_",tag_name,"_coeffector_among_gene_expression_gene_detected_with_sample_number_no_less_than_",average_count,".txt"),sep="\t", quote=F, row.names=T,col.names=T)
write.table(gene_exp_p.mat2,file=paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/R_pearson_correlationship_data_",tag,"_",tag_name,"_pvalue_among_gene_expression_gene_detected_with_sample_number_no_less_than_",average_count,".txt"),sep="\t", quote=F, row.names=T,col.names=T)

##plot heatmap
range(gene_exp_corr2)#  0.9642436 1.0000000
col=colorRampPalette(c("navy", "white", "firebrick3")) #设置颜色

range(gene_exp_corr2)
corheatmap3<-ggcorrplot(round(gene_exp_corr2,2),outline.color = "grey",ggtheme = ggplot2::theme_bw,
                        colors = c("#6D9EC1", "white", "#E46726"),title = "Corralation:pearson",lab = F)+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0.9,limits=c(0.8,1),breaks=c(0.8,0.85,0.9,0.95,1))
#     scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red', midpoint = 0.95,limits=c(0.9,1),breaks=c(0.9,0.925,0.95,0.975,1))
corheatmap3
ggsave(paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/heatmap_pearson_correlationship_data_",tag,"_",tag_name,"_pvalue_among_gene_expression_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),corheatmap3,width=10, height=10)

#perform PCA analysis 
PCA_data1 <- t(data_raw_bin)
set.seed(19921010)
md <- prep(PCA_data1, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 5)
summary(resPPCA) 
variation<-data.frame(summary(resPPCA))
df <-as.data.frame(scores(resPPCA))
PCA_data2<-merge(colData,df,by=0)
#PCA_data2$Group<-factor(PCA_data2$Group,levels=c("CTRL","HA","HAI"))

plot_PCA0<-ggplot(data=PCA_data2, mapping=aes(x=PC1,y=PC2))+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title =paste0("Distribution of samples \n ",tag))+scale_color_manual(values=c(ppCor))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),panel.grid=element_blank(),plot.title=element_text(hjust=0.5),axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
plot_PCA1<-plot_PCA0+geom_point(aes(colour =Day_Treat),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
plot_PCA2<-plot_PCA0+geom_point(aes(colour =Day),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
plot_PCA3<-plot_PCA0+geom_point(aes(colour =Treat),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
plot_PCA4<-plot_PCA1+geom_text(stat="identity",aes(label=Row.names), color="black", size=4)

plot_PCA5<-grid.arrange(plot_PCA1, plot_PCA2, plot_PCA3, plot_PCA4, ncol = 2)  
ggsave(paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/PCA_1_2",tag,"_",tag_name,"_pvalue_among_gene_expression_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),plot_PCA5,width=15, height=15)

plot_PCA0<-ggplot(data=PCA_data2, mapping=aes(x=PC2,y=PC3))+
  xlab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title =paste0("Distribution of samples \n ",tag))+scale_color_manual(values=c(ppCor))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),panel.grid=element_blank(),plot.title=element_text(hjust=0.5),axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
plot_PCA1<-plot_PCA0+geom_point(aes(colour =Day_Treat),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
plot_PCA2<-plot_PCA0+geom_point(aes(colour =Day),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
plot_PCA3<-plot_PCA0+geom_point(aes(colour =Treat),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
plot_PCA4<-plot_PCA1+geom_text(stat="identity",aes(label=Row.names), color="black", size=4)

plot_PCA523<-grid.arrange(plot_PCA1, plot_PCA2, plot_PCA3, plot_PCA4, ncol = 2)  
ggsave(paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/PCA_2_3",tag,"_",tag_name,"_pvalue_among_gene_expression_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),plot_PCA523,width=10, height=10)

plot_PCA0<-ggplot(data=PCA_data2, mapping=aes(x=PC1,y=PC3))+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title =paste0("Distribution of samples \n ",tag))+scale_color_manual(values=c(ppCor))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),panel.grid=element_blank(),plot.title=element_text(hjust=0.5),axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
plot_PCA1<-plot_PCA0+geom_point(aes(colour =Day_Treat),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
plot_PCA2<-plot_PCA0+geom_point(aes(colour =Day),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
plot_PCA3<-plot_PCA0+geom_point(aes(colour =Treat),stat= "identity",size=4,alpha=0.5,show.legend = TRUE)
plot_PCA4<-plot_PCA1+geom_text(stat="identity",aes(label=Row.names), color="black", size=4)

plot_PCA513<-grid.arrange(plot_PCA1, plot_PCA2, plot_PCA3, plot_PCA4, ncol = 2)  
ggsave(paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/PCA_1_3",tag,"_",tag_name,"_pvalue_among_gene_expression_gene_detected_with_sample_number_no_less_than_",average_count,".pdf"),plot_PCA513,width=10, height=10)


color_list<-rep(ppCor[2],18)
color_list[which(hc$order %in% c(1:3,7:9,13:15))] <- ppCor[1]
pdf(paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/D3_PCA_",tag,"_",tag_name,"_pvalue_among_gene_expression_gene_detected_with_sample_number_no_less_than_",average_count,"_1.pdf"),width = 10, height =10)
plot3d <- with(PCA_data2, scatterplot3d(PC1, PC2, PC3, color =color_list,cex.symbols = 1.2,pch = 16, font.lab = 2, font.axis = 2))
legend(plot3d$xyz.convert(60,20,20),pch=16,legend=levels(PCA_data2$Treat), col =ppCor,cex = rel(1.1), bty = 'n',yjust=0.5, xjust = 0, horiz = F)
dev.off()

color_list<-rep(ppCor[3],18)
color_list[which(hc$order %in% c(1:6))] <- ppCor[1]
color_list[which(hc$order %in% c(7:12))] <- ppCor[2]
pdf(paste0("/share/home/chenwei/LZJ_SCARE/1.global_evaluation/D3_PCA_",tag,"_",tag_name,"_pvalue_among_gene_expression_gene_detected_with_sample_number_no_less_than_",average_count,"_2.pdf"),width = 10, height =10)
plot3d <- with(PCA_data2, scatterplot3d(PC1, PC2, PC3, color =color_list,cex.symbols = 1.2,pch = 16, font.lab = 2, font.axis = 2))
legend(plot3d$xyz.convert(50,20,20),pch=16,legend=levels(PCA_data2$Day), col =ppCor,cex = rel(1.1), bty = 'n',yjust=0.5, xjust = 0, horiz = F)
dev.off()
}
