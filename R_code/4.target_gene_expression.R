rm(list = ls())  
options(stringsAsFactors = F)
library(org.Rn.eg.db)
library(GO.db)
library(KEGG.db)

library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggsci)
library(DESeq2)
library(pheatmap)
library(reshape2)
library(ggsci)
library(ggpubr)
library(scales)

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
#n=5;barplot(rep(1,n), col=colorRampPalette(colors = c('red', 'white'))( n ))
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]
show_col(ppCor_all2)

expression_matrix0 <-read.table(file="/share/home/chenwei/LZJ_SCARE/File1_Wound_Time_final_count.txt",sep="\t", row.names=1,header =T)
head(expression_matrix0);dim(expression_matrix0)

###STEP2 reading meta data
colData <-read.csv(file="/share/home/chenwei/LZJ_SCARE/final_meta_data.csv",sep=",", row.names=1,header =T)
colData$Treat<-factor(colData$Treat,levels = c("CTRL","PBS"))
colData$Day<-factor(colData$Day,levels = c("D3","D7","D20"))
colData$Day_Treat<-factor(colData$Day_Treat,levels = c("D3_CTRL","D3_PBS","D7_CTRL","D7_PBS","D20_CTRL","D20_PBS" ))
colData$pair_id<-as.character(unlist(lapply(strsplit(as.character(colData$Mouse_ID),"_"), function(x) paste0(x[2],"_",x[3]))))
colData$pair_id<-factor(colData$pair_id,levels =unique(colData$pair_id))
colData<-colData[order(colData$pair_id,colData$Day,colData$Treat,decreasing = F),]
head(colData);tail(colData)
str(colData)

expression_matrix1<-expression_matrix0[,colData$Mouse_ID]
dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix1, colData = colData, design = ~ pair_id+Treat)
dds1 <- dds0[ rowSums(counts(dds0)) > 1, ]
dds2 <- estimateSizeFactors(dds1)#data transform:sample number < 30 -> rlog
data_normal<-log2(counts(dds2,normalized=T)+1)

########################################
###gene_list for select_MAPK_Gene
genelist_name <-"select_MAPK_Gene"
Day3_target_genes <- c("Mmp3", "Mmp9", "Mmp10", "Mmp12", "Mmp13", "Rasgrf2", "Fgfbp1", "Traf1", "Igf2bp2", "Igfbp2", "Igfbp3", "Igfbp5", "Kras", "Fgfr2", "Fgfr4", "Mras", "Map3k12", "Fgf7", "Map4k2", "Map4k1", "Dusp2", "Rasgrp1", "Mapk13", "Map4k4", "Map3k3", "Map3k20", "Rasgrp3", "Mapkapk3", "Mapk9", "Fgf16", "Mapk11", "Fgfr3" )
Day7_target_genes <- c("Mmp3", "Mmp9", "Mmp10", "Mmp12", "Mmp13", "Igf2", "Rasgrp3", "Fgf1", "Fgf7", "Fgfr4", "Igfbp2", "Igfbp3")
Day20_target_genes <- c("Mmp3", "Mmp9", "Mmp10", "Mmp12", "Mmp13", "Igf2bp3", "Rasgrf2", "Mapk7", "Fgfbp1", "Traf1", "Igf2", "Map3k2", "Kras", "Fgfr2", "Fgfr4", "Mras", "Map3k12", "Map4k1",  "Fgf7", "Map4k2", "Mapk1", "Map4k1", "Dusp2", "Rasgrp1",  "Mapk13", "Map4k4", "Map3k3", "Map3k20", "Rasgrp3", "Rps6ka5", "Mapkapk3", "Mapk9", "Fgf1", "Mapk11", "Fgfr3" )

Commom_gene<-Reduce(intersect, list(Day3_target_genes, Day7_target_genes, Day20_target_genes))
final_gene<-c(Commom_gene,"Igfbp2", "Igfbp3", "Igfbp4","Igfbp5","Irs2","Fgfr2","Fgfr3")  #c(reference_genes0,"Igf1","Fgf2")

genelist_name <-"new_select_MAPK_Gene"
final_gene<-c("Mmp12", "Mmp10", "Mmp3", "Mmp13", "Mmp7", "Mmp11", "Mmp20", "Mmp19", "Mmp2", "Mmp9", "Mmp14", "Mmp23", "Fgf11", "Fgfr2", "Fgfr4", "Fgf7", "Fgf16", "Fgfr3", "Igf2", "Igfbp3", "Igfbp5", ",gfbp2", "Igf2bp2", "Igfbp4", "Mapk7", "Mapk13", "Mapk8ip3", "Mapkapk3", "Mapk9", "Mapk6", "Mapk11", "Traf1", "Kras", "Mras", "Dusp2", "Rasgrp1", "Rasgrf2", "Rasgrp3")
target_gene<-final_gene#c("Igf1","Fgf2")
Over_DEGs<-rownames(data_normal)[which(rownames(data_normal) %in% unique(target_gene))];length(Over_DEGs)#114
Over_DEGs2<-intersect(final_gene,Over_DEGs)

Gene_expression<-data_normal[which(rownames(data_normal)  %in% Over_DEGs2),]
Gene_expression <- Gene_expression[rowSums(Gene_expression) > 0, ]
Gene_expression[,1:4];dim(Gene_expression)# 114  18
range(Gene_expression)#0.00000 16.72701

#for RNA expression
head(colData)
annotation_col<-colData[order(colData$Day_Treat,decreasing = T),c("Treat","Day")]
ann_colors = list( Treat=c(CTRL=ppCor[2], PBS=ppCor[1]),Day=c(D3=ppCor[3], D7=ppCor[4],D20=ppCor[5]))#Day_Treat=c(D3_CTRL=ppCor[1],D3_PBS=ppCor[2],D7_CTRL=ppCor[3],D7_PBS=ppCor[4],D20_CTRL=ppCor[5],D20_PBS=ppCor[6]),
Gene_expression<-Gene_expression[Over_DEGs2,rownames(annotation_col)]

plot_heat<-pheatmap(Gene_expression, cluster_row =T,cluster_col =FALSE,scale ="row", na_col = "black", 
                    clustering_distance_rows ="euclidean",#correlation
                    show_rownames = T,show_colnames = T,
                    gaps_col =c(3,6,9,12,15) ,
                    annotation_col = annotation_col,#annotation_row=annotation_row,
                    annotation_colors = ann_colors, 
                    treeheight_col = 20, #treeheight_row = 30, 
                    border=FALSE,
                    color = colorRampPalette(c("navy","white","red"))(50),
                    main ="Gene expression level of significant differential genes(ttest_data):Scale value",angle_col ="45")

pdf(file=paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",genelist_name,"_heatmap_gene_cluster.pdf"),width = 8, height =8)
print(plot_heat)
dev.off()
plot_heat<-pheatmap(Gene_expression, cluster_row =F,cluster_col =FALSE,scale ="row", na_col = "black", 
                    clustering_distance_rows ="euclidean",#correlation
                    show_rownames = T,show_colnames = T,
                    gaps_col =c(3,6,9,12,15) ,
                    annotation_col = annotation_col,#annotation_row=annotation_row,
                    annotation_colors = ann_colors, 
                    treeheight_col = 20, #treeheight_row = 30, 
                    border=FALSE,
                    color = colorRampPalette(c("navy","white","red"))(50),
                    main ="Gene expression level of significant differential genes(ttest_data):Scale value",angle_col ="45")

pdf(file=paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",genelist_name,"_heatmap_gene_nocluster.pdf"),width = 8, height =8)
print(plot_heat)
dev.off()

###################
data_plot<-merge(t(Gene_expression),colData,by=0);data_plot<-data_plot[,-1]
data_plot2 <- melt(data_plot,variable.name="gene",value.name = "gene_expression",id.vars = c("Mouse_ID","Treat","Day","Day_Treat","pair_id"))
data_plot2$gene<-factor(data_plot2$gene,levels = Over_DEGs2)

data_plot2<-data_plot2[order(data_plot2$gene,data_plot2$pair_id,data_plot2$Day,data_plot2$Treat,decreasing = F),]

###Treat
head(data_plot2);range(data_plot2$gene_expression)
#my_comparisons <- list(c(Group[1],Group[2]),c(Group[1],Group[3]),c(Group[2],Group[3]))
###个体分开绘制geom_violin ggboxplot
head(data_plot2)
stat_data_ttest<-compare_means(gene_expression ~ Treat,group.by = c("Day","gene"), data = data_plot2, method = "t.test", paired = TRUE)
stat_data_wtest<-compare_means(gene_expression ~ Treat,group.by = c("Day","gene"),data = data_plot2, method = "wilcox.test", paired = TRUE)  
write.table(stat_data_ttest,file= paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",genelist_name,"_PBS_vs_CTRL_stat_data_paired_ttest_model.txt"),sep="\t", quote=F, row.names=T,col.names=T)
write.table(stat_data_wtest,file= paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",genelist_name,"_PBS_vs_CTRL_stat_data_paired_wtest_model.txt"),sep="\t", quote=F, row.names=T,col.names=T)

##plot for target gene
#reference_genes0<-c("Actb","Gapdh","Ppib")
###Treat###个体分开绘制geom_violin ggboxplot
head(data_plot2);range(data_plot2$gene_expression)
#my_comparisons <- list(c(Group[1],Group[2]),c(Group[1],Group[3]),c(Group[2],Group[3]))
#data_plot2$pair_id<-as.character(unlist(lapply(strsplit(as.character(data_plot2$Mouse_ID),"_"), function(x) paste0(x[2],"_",x[3]))))
compare_means(gene_expression ~ Treat,group.by = c("Day","gene"), data = data_plot2, method = "t.test", paired = TRUE)
compare_means(gene_expression ~ Treat,group.by = c("Day","gene"),data = data_plot2, method = "wilcox.test", paired = TRUE) 

P_mean_site<-ggpaired(data_plot2,x = "Day_Treat", y = "gene_expression",id = "pair_id",  color = "Treat",palette = ppCor[2:1],  line.color = "gray", line.size = 0.4) +
  facet_wrap( ~gene, scales = "free_y",nrow=4)+   stat_summary(aes(group = Treat),fun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9))+
  stat_compare_means(aes(group = Treat),paired = TRUE,method = "t.test", label = "p.signif",hide.ns = TRUE,color="red",
                     comparisons = list(
                       c(levels(data_plot2$Day_Treat)[1], levels(data_plot2$Day_Treat)[2]),  # 第一个Day的Treat比较
                       c(levels(data_plot2$Day_Treat)[3], levels(data_plot2$Day_Treat)[4]),  
                       c(levels(data_plot2$Day_Treat)[5], levels(data_plot2$Day_Treat)[6])),# 第二个Day的Treat比较
                     label.y = max(data_plot2$gene_expression)) +
  labs(title="mRNA expression level for PBS related genes(paired t.test)", x ="sample_group", y ="log2(normalized count+1)")+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),
        axis.title.y = element_text(size=10,colour = "black",face = "bold"),
        axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))+#+stat_summary(afun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9)) 
  ylim(0,max(data_plot2$gene_expression)+1)
P_mean_site
ggsave(paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",genelist_name,"_PBS_vs_CTRL_stat_data_paired_ttest_box_plot.pdf"),P_mean_site,width=17, height=14)

#ggsave(paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",genelist_name,"_PBS_vs_CTRL_stat_data_paired_ttest_box_plot.pdf"),P_mean_site,width=32, height=14)

########################################
###gene list for commom DEGs
# p_value<0.05 & Foldchange >=1.5
##DEGs reading
Treat_select<-c("PBS","CTRL")
Treat_tag <-"final";average_count<-3;Treat_group<-"PBS_vs_CTRL"
sampleA <-Treat_select[1]; sampleB<-Treat_select[2]
DEG_list_up<-list();DEG_list_dw<-list()
for(Day_tag in c("D3","D7","D20")){
  #Day_tag <-"D14"
  file_DEG <- paste0("/share/home/chenwei/LZJ_SCARE/intergroup/file3_paired_",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
  DEG_file<-read.table(file=file_DEG,sep="\t", row.names=1,header =T)
  
  DEG_up <-  rownames(DEG_file[which(DEG_file$log2FoldChange> 0),])
  DEG_dw <-  rownames(DEG_file[which(DEG_file$log2FoldChange< 0),])
  DEG_list_up<-c(DEG_list_up,list(DEG_up))
  DEG_list_dw<-c(DEG_list_dw,list(DEG_dw))
}


Commom_UP_gene<-Reduce(intersect, list(DEG_list_up[[1]], DEG_list_up[[2]], DEG_list_up[[3]]))
Commom_dw_gene<-Reduce(intersect, list(DEG_list_dw[[1]], DEG_list_dw[[2]], DEG_list_dw[[3]]))
length(Commom_UP_gene);length(Commom_dw_gene)#833 611

genelist_name<-"Commom_UP_down";target_gene<-c(Commom_UP_gene,Commom_dw_gene)
Over_DEGs<-rownames(data_normal)[which(rownames(data_normal) %in% unique(target_gene))];length(Over_DEGs)#114
Over_DEGs2<-Over_DEGs[match(target_gene,Over_DEGs)]

Gene_expression<-data_normal[which(rownames(data_normal)  %in% Over_DEGs2),]
Gene_expression <- Gene_expression[rowSums(Gene_expression) > 0, ]
Gene_expression[,1:4];dim(Gene_expression)# 114  18
range(Gene_expression)#0.00000 16.72701

#for RNA expression
head(colData)
annotation_col<-colData[order(colData$Day_Treat,decreasing = T),c("Treat","Day")]
ann_colors = list( Treat=c(CTRL=ppCor[2], PBS=ppCor[1]),Day=c(D3=ppCor[3], D7=ppCor[4],D20=ppCor[5]))#Day_Treat=c(D3_CTRL=ppCor[1],D3_PBS=ppCor[2],D7_CTRL=ppCor[3],D7_PBS=ppCor[4],D20_CTRL=ppCor[5],D20_PBS=ppCor[6]),
Gene_expression<-Gene_expression[Over_DEGs2,rownames(annotation_col)]


plot_heat<-pheatmap(Gene_expression, cluster_row =F,cluster_col =FALSE,scale ="row", na_col = "black", 
                    clustering_distance_rows ="euclidean",#correlation
                    show_rownames = F,show_colnames = T,
                    gaps_col =c(3,6,9,12,15) ,
                    annotation_col = annotation_col,#annotation_row=annotation_row,
                    annotation_colors = ann_colors, 
                    treeheight_col = 20, #treeheight_row = 30, 
                    border=FALSE,
                    color = colorRampPalette(c("navy","white","red"))(50),
                    main ="Gene expression level of common DEGs:Scale value",angle_col ="45")

pdf(file=paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",genelist_name,"_heatmap_gene_nocluster.pdf"),width = 8, height =8)
print(plot_heat)
dev.off()

##################待修改 plot for genes in target GO term or KEGG
##DEGs reading
Treat_select<-c("PBS","CTRL")
Treat_tag <-"final";average_count<-3;Treat_group<-"PBS_vs_CTRL"
sampleA <-Treat_select[1]; sampleB<-Treat_select[2]
# p_value<0.05 & Foldchange >=1.5
DEG_final<-data.frame()
for(Day_tag in c("D3","D7","D20")){
  #Day_tag <-"D14"
  file_DEG <- paste0("/share/home/chenwei/LZJ_SCARE/intergroup/file3_paired_",Treat_group,"_",Day_tag,"_",Treat_tag,"_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
  DEG_file<-read.table(file=file_DEG,sep="\t", row.names=1,header =T)
  
  DEG_up <-  rownames(DEG_file[which(DEG_file$log2FoldChange> 0),])
  DEG_dw <-  rownames(DEG_file[which(DEG_file$log2FoldChange< 0),])
  DEGs_temp<-data.frame(DEGs=c(DEG_up,DEG_dw),class =c(rep("DEG_up",length(DEG_up)),rep("DEG_dw",length(DEG_dw))),Day=rep(Day_tag,length(c(DEG_up,DEG_dw))))
  
  DEG_final<-rbind(DEG_final,DEGs_temp)
}
############# 获取所有 GO term 及其对应的基因
go_genes <- as.list(org.Rn.egGO2ALLEGS)
#go_terms <- as.list(org.Rn.egGO)
go_descriptions <- as.list(Term(GO.db::GOTERM))
go_id_list<-names(go_descriptions[grepl("MAPK", toupper(go_descriptions))])
go_info <- data.frame(GO_term = character(), Description = character(), Gene_symbol = character(), stringsAsFactors = FALSE)
# 遍历每个 GO term，提取名称、描述和 gene symbol
for (go_id in go_id_list) {
  # 获取 GO term 名称和描述
  #go_id<-"GO:0051015"
  go_description <- go_descriptions[[go_id]]
  print(go_description)
  gene_symbols <- mget(go_genes[[go_id]], org.Rn.egSYMBOL, ifnotfound = NA)
  gene_symbols <- as.character(unlist(gene_symbols))
  gene_symbols <- gene_symbols[!is.na(gene_symbols)]
  
    # 如果 gene_symbols 依然有值，添加到数据框
    if (length(gene_symbols) > 0) {
      go_info <- rbind(go_info, data.frame(GO_term = go_id, Description = go_description, Gene_symbol = gene_symbols, stringsAsFactors = FALSE))
    }
  }
head(go_info)
table(go_info$Description)
names(table(go_info$Description))[table(go_info$Description) !=  47817]

# 获取所有 KEGG pathway 及其对应的基因
kegg_genes <- as.list(org.Rn.egPATH2EG)
kegg_terms <- as.list(KEGGPATHID2NAME)
KEGG_id_list<-names(kegg_terms[grepl("MAPK", toupper(kegg_terms))])
kegg_info <- data.frame(Pathway_ID = character(),Pathway_Name = character(),Description = character(),Gene_symbol = character(), stringsAsFactors = FALSE)
for (pathway_id in KEGG_id_list) {
  # pathway_id<-"04010" 
  pathway_name <- kegg_terms[[pathway_id]]
  gene_symbols <- mget(kegg_genes[[pathway_id]], org.Rn.egSYMBOL, ifnotfound = NA)
  gene_symbols <- as.character(unlist(gene_symbols))
  gene_symbols <- gene_symbols[!is.na(gene_symbols)]
  
  if (length(gene_symbols) > 0) {
    kegg_info <- rbind(kegg_info, data.frame(Pathway_ID = pathway_id, Pathway_Name = pathway_name, Gene_symbol = gene_symbols, stringsAsFactors = FALSE))
}
  }
head(kegg_info)
names(table(kegg_info$Pathway_Name))[table(kegg_info$Pathway_Name) !=  47817]


########PLOT####################
pathway_name<-"MAPK signaling pathway"
target_PATHWAY_set <-kegg_info[which(kegg_info$Pathway_Name == pathway_name),]
Over_DEGs<-unique(DEG_final[which(DEG_final$DEGs %in% unique(target_PATHWAY_set$Gene_symbol)),]$DEGs);length(Over_DEGs)#114
Gene_expression<-data_normal[which(rownames(data_normal)  %in% Over_DEGs),]
Gene_expression <- Gene_expression[rowSums(Gene_expression) > 0, ]
Gene_expression[,1:4];dim(Gene_expression)# 114  18
range(Gene_expression)#0.00000 16.72701

#for RNA expression
head(colData)
annotation_col<-colData
ann_colors = list(Day_Treat=c(D3_CTRL=ppCor[1],D3_PBS=ppCor[2],D7_CTRL=ppCor[3],D7_PBS=ppCor[4],D20_CTRL=ppCor[5],D20_PBS=ppCor[6]),
                  Treat=c(CTRL=ppCor[2], PBS=ppCor[1]),Day=c(D3=ppCor[1], D7=ppCor[2],D20=ppCor[3]))
Gene_expression<-Gene_expression[Over_DEGs,rownames(annotation_col)]

plot_heat<-pheatmap(Gene_expression, cluster_row =T,cluster_col =FALSE,scale ="row", na_col = "black", 
                    clustering_distance_rows ="euclidean",#correlation
                    show_rownames = T,show_colnames = T,
                    gaps_col =c(3,6,9,12,15) ,
                    annotation_col = annotation_col[,2:4],#annotation_row=annotation_row,
                    annotation_colors = ann_colors, 
                    treeheight_col = 20, #treeheight_row = 30, 
                    border=FALSE,
                    color = colorRampPalette(c("lightblue","white","red"))(50),
                    main ="Gene expression level of significant differential genes(ttest_data):Scale value",angle_col ="45")

pdf(file=paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",Treat_group,"_",pathway_name,"_heatmap_gene.pdf"),width = 8, height =16)
print(plot_heat)
dev.off()

###################
data_plot<-merge(t(Gene_expression),colData,by=0);data_plot<-data_plot[,-1]
head(data_plot)
data_plot2 <- melt(data_plot,variable.name="gene",value.name = "gene_expression",id.vars = c("Mouse_ID","Treat","Day","Day_Treat"))

###Treat
head(data_plot2);range(data_plot2$gene_expression)
#my_comparisons <- list(c(Group[1],Group[2]),c(Group[1],Group[3]),c(Group[2],Group[3]))
###个体分开绘制geom_violin ggboxplot
head(data_plot2)
stat_data_ttest<-compare_means(gene_expression ~ Treat,group.by = c("Day","gene"), data = data_plot2, method = "t.test")
stat_data_wtest<-compare_means(gene_expression ~ Treat,group.by = c("Day","gene"),data = data_plot2, method = "wilcox.test")  
write.table(stat_data_ttest,file= paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",Treat_group,"_",pathway_name,"_stat_data_ttest_model.txt"),sep="\t", quote=F, row.names=T,col.names=T)
write.table(stat_data_wtest,file= paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",Treat_group,"_",pathway_name,"_stat_data_wtest_model.txt"),sep="\t", quote=F, row.names=T,col.names=T)

##plot for target gene
reference_genes0<-c("Actb","Gapdh","Ppib")

final_gene<-c(reference_genes0,"Mmp3","Mmp9")
Gene_expression<-data_normal[which(rownames(data_normal)  %in% final_gene),]
data_plot<-merge(t(Gene_expression),colData,by=0);data_plot<-data_plot[,-1]
data_plot2 <- melt(data_plot,variable.name="gene",value.name = "gene_expression",id.vars = c("Mouse_ID","Treat","Day","Day_Treat"))

###Treat
head(data_plot2);range(data_plot2$gene_expression)
#my_comparisons <- list(c(Group[1],Group[2]),c(Group[1],Group[3]),c(Group[2],Group[3]))
###个体分开绘制geom_violin ggboxplot
head(data_plot2)
compare_means(gene_expression ~ Treat,group.by = c("Day","gene"), data = data_plot2, method = "t.test")
compare_means(gene_expression ~ Treat,group.by = c("Day","gene"),data = data_plot2, method = "wilcox.test")  
data_plot2$gene<-factor(data_plot2$gene,levels = final_gene)

P_mean_site<-ggboxplot(data_plot2, x = "Day", y = "gene_expression", color = "Treat",palette = ppCor[2:1], add = "jitter")+facet_wrap( ~gene, scales = "free_y",nrow=2)+ # 
  #ggviolin(data_plot2, x = "group", y = "gene_expression", color = "group",palette = ppCor,outlier.shape = NA)+#ylim(c(70,85))+, add = "jitter"
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=0,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),
        axis.title.y = element_text(size=10,colour = "black",face = "bold"),
        axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))+#+stat_summary(afun=mean, geom="point", shape=18, size=2, col="black", position = position_dodge(0.9)) 
  ylim(0,max(data_plot2$gene_expression)+1)

#P_mean_site2<-P_mean_site+stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif",hide.ns = TRUE,color="red") +labs(title="mRNA expression level for AMA related genes(wilcox.test)", x ="sample_group", y ="log2(count+1)")
#P_mean_site3<-P_mean_site+stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.signif",hide.ns = TRUE,color="red") +labs(title="mRNA expression level for AMA related genes(wilcox.test)", x ="sample_group", y ="log2(count+1)")
P_mean_site2<-P_mean_site+stat_compare_means(aes(group = Treat),method = "wilcox.test", label = "p.signif",hide.ns = TRUE,color="red",label.y = max(data_plot2$gene_expression)) +labs(title="mRNA expression level for PBS related genes(wilcox.test)", x ="sample_group", y ="log2(normalized count+1)")
P_mean_site3<-P_mean_site+stat_compare_means(aes(group = Treat),method = "t.test", label = "p.signif",hide.ns = TRUE,color="red",label.y = max(data_plot2$gene_expression)) +labs(title="mRNA expression level for PBS related genes(t.test)", x ="sample_group", y ="log2(normalized count+1)")
P_mean_site2
P_mean_site3

ggsave(paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",Treat_group,"_box_target_gene_plot_stat_wtest.pdf"),P_mean_site2,width=15, height=10)
ggsave(paste0("/share/home/chenwei/LZJ_SCARE/gene_expression/",Treat_group,"_box_target_gene_plot_stat_ttest.pdf"),P_mean_site3,width=15, height=10)
