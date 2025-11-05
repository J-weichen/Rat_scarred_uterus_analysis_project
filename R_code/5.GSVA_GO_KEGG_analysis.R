rm(list = ls())  

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(msigdbr)
library(GSEABase)
library(GSVA)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

##Step1： Load data
##加载sample*Gene symbol矩阵

expression_matrix <-read.table(file="/share/home/chenwei/Scarred_uterus/3.analysis_result/1.global_evaluation/Final_remian_sample_gene_detected_with_sample_number_no_less_than_3_deseq2_normalized_count.txt",sep="\t", row.names=1,header =T)
expression_matrix0 <- log2(expression_matrix+1)
head(expression_matrix0);dim(expression_matrix0)#20279    40

###STEP2 reading meta data
colData <-read.csv(file="/share/home/chenwei/Scarred_uterus/3.analysis_result/final_meta_data.csv",sep=",", row.names=1,header =T)
colData$Treat<-factor(colData$Treat,levels = c("sham","saline","luc","hgf"))
colData$Day<-factor(colData$Day,levels = c("D7","D14","D30"))
colData$Day_Treat<-factor(colData$Day_Treat,levels = c("D7_sham","D7_saline","D7_luc","D7_hgf","D14_sham","D14_saline","D14_luc","D14_hgf","D30_sham","D30_saline","D30_luc","D30_hgf"))
head(colData);dim(colData)

###3. GSVA analysis
##3.1 Gene set data download (usually GO and [or] KEGG are selected for GSVA analysis)
#GVSA analysis with GO data sets
GO_df_all <- msigdbr(species="Mus musculus", category="C5")# Homo sapiens or Mus musculus
GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
head(GO_df)
go_list <- split(GO_df$gene_symbol, GO_df$gs_name) #按照gs_name给gene_symbol分

##3.3 GSVA分析
GSVA.res <- gsva(as.matrix(expression_matrix0), gset.idx.list=go_list,min.sz>1, max.sz=Inf,verbose=TRUE, method="gsva",kcdf="Gaussian",parallel.sz = parallel::detectCores())#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
head(GSVA.res)
dim(GSVA.res)#10395    40
saveRDS(GSVA.res,file="/share/home/chenwei/LZJ_SCARE/GSVA_result/GSVA_GSVA_res_sample_remained_gene_detected_with_sample_number_no_less_than_3.rds")


GSVA.res[rownames(GSVA.res)[grepl(toupper("healing"),rownames(GSVA.res))],]#"REACTOME_PYROPTOSIS"
#gsva.df <- data.frame(Genesets=rownames(GSVA.res), GSVA.res, check.names = F)
##randomly select 30 terms
gsva_d = GSVA.res[sample(nrow(GSVA.res),30),]
pheatmap::pheatmap(gsva_d, show_colnames = T,  scale = "row",angle_col = "45",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

##4. 差异分析（Mann-Whitney U检验,也被称为Wilcoxon秩和检验）
head(colData)
group_info<- colData[,c("Day_Treat","Mouse_ID")]
colnames(group_info)<-c("Group","sample")
group_info$group <- as.numeric(factor(group_info$Group))  
results <- data.frame()
# 两两比较  
for (i in 1:(length(levels(group_info$Group))-1)) { 
  #i=1
  for (j in (i + 1):length(levels(group_info$Group))) { 
    # j=2
    # 获取属于当前组
    Class1<-unique(as.character(group_info[group_info$group == i,]$Group))
    Class2<-unique(as.character(group_info[group_info$group == j,]$Group))
    print(Class1);print(Class2)
    
    # 提取对应样本的数据  
    class1_samples <- group_info$sample[group_info$group == i]  
    class2_samples <- group_info$sample[group_info$group == j]
    
    ## 创建一个空的数据框用于存储结果
    result_df <- data.frame(Geneset = character(), P_Value = numeric(), stringsAsFactors =FALSE)
    for(geneset in rownames(GSVA.res)){
      
      class1_data <- GSVA.res[geneset, class1_samples]  
      class2_data <- GSVA.res[geneset, class2_samples]  
      #执行Wilcoxon秩和检验
      wilcox_result <- wilcox.test(class1_data, class2_data, paired = FALSE)
      #提取p值
      p_value<-wilcox_result$p.value
      #将结果添加到结果数据框
      result_df<- bind_rows(result_df,data.frame(Geneset = geneset,P_Value= p_value))
    }
    result_df$Adjusted_P_Value <- p.adjust(result_df$P_Value, method ="BH")
    result_df$compare_group<-paste0(Class1,"-",Class2)
    results<-rbind(results,result_df)
  }  
}  
dim(results)
results_sig<-results[which(results$P_Value<0.05),]
table(results_sig$compare_group)
###D14_saline-D30_luc  D7_hgf-D14_saline     D7_hgf-D30_luc  D7_luc-D14_saline       D7_luc-D30_luc     D7_luc-D7_hgf   D7_sham-D7_luc  D7_luc-D30_hgf  D7_luc-D30_sham  
###        567               3535               2604               2211               3284                1979                 1               1                    1 
results_sig[grepl(toupper("healing"),results_sig$Geneset),]
#result_df <- results[which(results$compare_group =="D7_luc-D7_hgf"),] %>% top_n(-20, wt =P_Value)
write.table(results,file="/share/home/chenwei/LZJ_SCARE/GSVA_result/GSVA_GO_sample_remained_gene_detected_with_sample_number_no_less_than_3.txt",sep="\t", quote=F, row.names=T,col.names=T)
###########KEGG##########
KEGG_df_all <- msigdbr(species="Mus musculus", category="C2")# Homo sapiens or Mus musculus
KEGG_df <- dplyr::select(KEGG_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
#GO_df <- GO_df[KEGG_df$gs_subcat!="HPO",]
head(KEGG_df)
KEGG_list <- split(KEGG_df$gene_symbol, KEGG_df$gs_name) #按照gs_name给gene_symbol分

##3.3 GSVA分析
GSVA_KEGG.res <- gsva(as.matrix(expression_matrix0), gset.idx.list=KEGG_list,min.sz>1, max.sz=Inf,verbose=TRUE, method="gsva",kcdf="Gaussian",parallel.sz = parallel::detectCores())#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
head(GSVA_KEGG.res)
dim(GSVA_KEGG.res)#6365   40
saveRDS(GSVA_KEGG.res,file="/share/home/chenwei/LZJ_SCARE/GSVA_result/GSVA_GSVA_KEGG_res_sample_remained_gene_detected_with_sample_number_no_less_than_3.rds")

####
rownames(GSVA_KEGG.res)[grepl(toupper("healing"),rownames(GSVA_KEGG.res))]#"WP_BURN_WOUND_HEALING"
#gsva.df <- data.frame(Genesets=rownames(GSVA_KEGG.res), GSVA_KEGG.res, check.names = F)
##randomly select 30 terms
gsva_d = GSVA_KEGG.res[sample(nrow(GSVA_KEGG.res),30),]
pheatmap::pheatmap(gsva_d, show_colnames = T,  scale = "row",angle_col = "45",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

##4. 差异分析（Mann-Whitney U检验,也被称为Wilcoxon秩和检验）
head(colData)
group_info<- colData[,c("Day_Treat","Mouse_ID")]
colnames(group_info)<-c("Group","sample")
group_info$group <- as.numeric(factor(group_info$Group))  
results <- data.frame()
# 两两比较  
for (i in 1:(length(levels(group_info$Group))-1)) { 
  #i=1
  for (j in (i + 1):length(levels(group_info$Group))) { 
    # j=2
    # 获取属于当前组
    Class1<-unique(as.character(group_info[group_info$group == i,]$Group))
    Class2<-unique(as.character(group_info[group_info$group == j,]$Group))
    print(Class1);print(Class2)
    
    # 提取对应样本的数据  
    class1_samples <- group_info$sample[group_info$group == i]  
    class2_samples <- group_info$sample[group_info$group == j]
    
    ## 创建一个空的数据框用于存储结果
    result_df <- data.frame(Geneset = character(), P_Value = numeric(), stringsAsFactors =FALSE)
    for(geneset in rownames(GSVA_KEGG.res)){
      
      class1_data <- GSVA_KEGG.res[geneset, class1_samples]  
      class2_data <- GSVA_KEGG.res[geneset, class2_samples]  
      #执行Wilcoxon秩和检验
      wilcox_result <- wilcox.test(class1_data, class2_data, paired = FALSE)
      #提取p值
      p_value<-wilcox_result$p.value
      #将结果添加到结果数据框
      result_df<- bind_rows(result_df,data.frame(Geneset = geneset,P_Value= p_value))
    }
    result_df$Adjusted_P_Value <- p.adjust(result_df$P_Value, method ="BH")
    result_df$compare_group<-paste0(Class1,"-",Class2)
    results<-rbind(results,result_df)
  }  
}  
dim(results)
results_sig<-results[which(results$P_Value<0.05),]
table(results_sig$compare_group)
###D14_saline-D30_luc  D7_hgf-D14_saline     D7_hgf-D30_luc  D7_luc-D14_saline     D7_luc-D30_luc      D7_luc-D7_hgf 
###       299               2212               1625               1228               1891               1139
results_sig[grepl(toupper("healing"),results_sig$Geneset),]
write.table(results,file="/share/home/chenwei/LZJ_SCARE/GSVA_result/GSVA_KEGG_sample_remained_gene_detected_with_sample_number_no_less_than_3.txt",sep="\t", quote=F, row.names=T,col.names=T)


##############
GSVA.res<-readRDS(file="/share/home/chenwei/LZJ_SCARE/GSVA_result/GSVA_GSVA_res_sample_remained_gene_detected_with_sample_number_no_less_than_3.rds")
file_target<- "/share/home/chenwei/LZJ_SCARE/GSVA_result/GSVA_GO_sample_remained_gene_detected_with_sample_number_no_less_than_3.txt"
results<-read.table(file=file_target,sep="\t", row.names=1,header =T)
results_sig<-results[which(results$P_Value<0.05),]

##5.GSVA可视化
Group<-c("D7_luc","D7_hgf")
ann_colors <- list(group = c(D7_luc ="blue",D7_hgf = "firebrick"))

Group<-c("D14_luc","D14_hgf")
ann_colors <- list(group = c(D14_luc ="blue",D14_hgf = "firebrick"))

Group<-c("D30_luc","D30_hgf")
ann_colors <- list(group = c(D30_luc ="blue",D30_hgf = "firebrick"))


sample_list<-colData[which(colData$Day_Treat %in% Group),]$Mouse_ID
##提取目标基因集在每个样本中的表达矩阵
row_names <- unique(results_sig[grepl(toupper("angiogenesis"),results_sig$Geneset),]$Geneset)
#row_names <- unique(results_sig[grepl(toupper("healing"),results_sig$Geneset),]$Geneset)
gsva_selected <- GSVA.res[row_names,sample_list]

##5.1 设置注释信息
annotation_col <- data.frame(colData[which(colData$Day_Treat %in% Group),]$Day_Treat)
rownames(annotation_col) <- sample_list
colnames(annotation_col)<-"Group"


##5.3 pheatmap
plot1 <- pheatmap::pheatmap(gsva_selected,cluster_cols=F,cluster_rows=T,scale ="row", 
                            show_colnames=T,how_rownames=T,
                            annotation_col=annotation_col, annotation_colors=ann_colors,
                            #fontsize_row=8,fontsize_col=8,
                            #height=8,width=10,cellwidth=8,cellheight=8,treeheight_row=20,treeheight_col=20,
                            color=colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100),border_color="black",
                            gaps_col=length(colData[which(colData$Day_Treat == Group[1]),]$Mouse_ID),
                            main ="GSVA level of selected GO term :scale by row")
plot2 <- pheatmap::pheatmap(gsva_selected,cluster_cols=F,cluster_rows=T,#scale ="row", 
                            show_colnames=T,how_rownames=T,
                            annotation_col=annotation_col, annotation_colors=ann_colors,
                            color = colorRampPalette(c("navy","pink","orange","firebrick3"))(50),
                            #color=colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100),border_color="black",
                            gaps_col=length(colData[which(colData$Day_Treat == Group[1]),]$Mouse_ID),
                            main ="GSVA level of selected GO term :scale by row")

pdf(paste0("/share/home/chenwei/LZJ_SCARE/GSVA_result/heatmap_GO_woud_healing_",Group[1],"_vs_",Group[2],".pdf"),width = 10, height =10)
print(plot1);print(plot2)
dev.off()

