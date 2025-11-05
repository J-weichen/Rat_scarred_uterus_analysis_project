rm(list = ls())  
options(stringsAsFactors = F)
#library(org.Hs.eg.db)
#library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(GseaVis)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggsci)
library(scales)
library(GOstats)
library(KEGGREST)

length(keggList("pathway",'rno'))
#update.packages("KEGG.db","org.Rn.eg.db")

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

###expression file reading
expression_matrix0 <-read.table(file="/share/home/chenwei/LZJ_SCARE/File1_Wound_Time_final_count.txt",sep="\t", row.names=1,header =T)
head(expression_matrix0);dim(expression_matrix0)
###STEP2 reading meta data
colData <-read.csv(file="/share/home/chenwei/LZJ_SCARE/final_meta_data.csv",sep=",", row.names=1,header =T)
colData$Treat<-factor(colData$Treat,levels = c("CTRL","PBS"))
colData$Day<-factor(colData$Day,levels = c("D3","D7","D20"))
colData$Day_Treat<-factor(colData$Day_Treat,levels = c("D3_CTRL","D3_PBS","D7_CTRL","D7_PBS","D20_CTRL","D20_PBS" ))
head(colData);tail(colData)

###目标通路hallmark genelist分析

############group select###############
Group<-c("D3_CTRL","D3_PBS")
Group<-c("D7_CTRL","D7_PBS")
Group<-c("D20_CTRL","D20_PBS")

#########################################
# 提取对应样本的数据  
class1_samples <- colData$Mouse_ID[colData$Day_Treat == Group[1]]  
class2_samples <- colData$Mouse_ID[colData$Day_Treat == Group[2]]

expression_matrix1<-expression_matrix0[,c(class1_samples,class2_samples)]
colData_select<-colData[which(colData$Day_Treat %in% Group),]
colData_select$Day_Treat<-factor(colData_select$Day_Treat,levels = Group)

colData_select$pair_id<-as.character(unlist(lapply(strsplit(as.character(colData_select$Mouse_ID),"_"), function(x) paste0(x[2],"_",x[3]))))
colData_select$pair_id<-factor(colData_select$pair_id,levels =unique(colData_select$pair_id))
colData_select<-colData_select[order(colData_select$pair_id,colData_select$Treat,decreasing = F),]
head(colData_select);tail(colData_select)
str(colData_select)

expression_matrix1<-expression_matrix1[,colData_select$Mouse_ID]
dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix1, colData = colData_select, design = ~ pair_id+Day_Treat)
dds1 <- dds0[ rowSums(counts(dds0)) > 1, ]
dds2 <- estimateSizeFactors(dds1)
dds3 <- DESeq(dds2)

sampleA <-Group[2]; sampleB<-Group[1]
contrastV <- c("Day_Treat",sampleA, sampleB)

res1 <- as.data.frame(results(dds3, contrast=contrastV))
summary(res1)

need_DEG <- res1[,c(2,5)] #选择log2FoldChange和pvalue（凑成数据框）
colnames(need_DEG) <- c('log2FoldChange','pvalue')
need_DEG$SYMBOL <- rownames(need_DEG)
head(need_DEG)

##for human#######################################
##### 创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）
geneList <- need_DEG$log2FoldChange
names(geneList) <- toupper(need_DEG$SYMBOL)
geneList <- sort(geneList, decreasing = T)   #从大到小排序

KEGG_MAPK_PATHWAY_set <- read.gmt('/share/home/chenwei/LZJ_SCARE/KEGG_MAPK_SIGNALING_PATHWAY.v2024.1.Hs.gmt')
ST_ERK1_ERK2_MAPK_PATHWAY_set <- read.gmt('/share/home/chenwei/LZJ_SCARE/ST_ERK1_ERK2_MAPK_PATHWAY.v2024.1.Hs.gmt')
BIOCARTA_ERK_PATHWAY_set <- read.gmt('/share/home/chenwei/LZJ_SCARE/BIOCARTA_ERK_PATHWAY.v2024.1.Hs.gmt')
head(KEGG_MAPK_PATHWAY_set)
head(ST_ERK1_ERK2_MAPK_PATHWAY_set)
head(BIOCARTA_ERK_PATHWAY_set)

backgene_set<-rbind(KEGG_MAPK_PATHWAY_set,ST_ERK1_ERK2_MAPK_PATHWAY_set,BIOCARTA_ERK_PATHWAY_set)
colnames(backgene_set)<-c("ont","gene")
head(backgene_set)
table(backgene_set$ont)
backgene_set$ont<-factor(backgene_set$ont,levels=c("KEGG_MAPK_SIGNALING_PATHWAY","ST_ERK1_ERK2_MAPK_PATHWAY", "BIOCARTA_ERK_PATHWAY"),ordered=TRUE)

target_GO_GSEA <- GSEA(geneList,TERM2GENE = backgene_set,  pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                       #minGSSize = 10,maxGSSize = 500,pvalueCutoff = 0.05,
                       pAdjustMethod = 'BH')  #指定 p 值校正方法
gseaplot2(target_GO_GSEA, color = colorspace::rainbow_hcl(4),geneSetID = rownames(target_GO_GSEA[1,]),title=(paste0(Group[2],"_vs_",Group[1],"_KEGG_MAPK_PATHWAY_set")), pvalue_table = TRUE)

P1<-gseaplot2(target_GO_GSEA, color = ppCor[1:3],geneSetID = c("KEGG_MAPK_SIGNALING_PATHWAY","ST_ERK1_ERK2_MAPK_PATHWAY", "BIOCARTA_ERK_PATHWAY"),title=(paste0(Group[2],"_vs_",Group[1],"_Three target pathway")),rel_heights = c(1.3, 0.3, 0.6), pvalue_table = T)
P1
ggsave(P1,file=paste0("/share/home/chenwei/LZJ_SCARE/GSEA_result/",Group[2],"_vs_",Group[1],"_human_gsea_MAPK_SIGNALING_PATHWAY.pdf"),width = 12,height =8,limitsize = F)


##for mouse#######################
geneList <- need_DEG$log2FoldChange
names(geneList) <- need_DEG$SYMBOL
geneList <- sort(geneList, decreasing = T)   #从大到小排序
head(geneList)

BIOCARTA_MAPK_PATHWAY_set <- read.gmt('/share/home/chenwei/LZJ_SCARE/BIOCARTA_MAPK_PATHWAY.v2024.1.Mm.gmt')
BIOCARTA_ERK_PATHWAY_set <- read.gmt('/share/home/chenwei/LZJ_SCARE/BIOCARTA_ERK_PATHWAY.v2024.1.Mm.gmt')
head(BIOCARTA_MAPK_PATHWAY_set)
head(BIOCARTA_ERK_PATHWAY_set)

backgene_set<-rbind(BIOCARTA_MAPK_PATHWAY_set,BIOCARTA_ERK_PATHWAY_set)
colnames(backgene_set)<-c("ont","gene")
head(backgene_set)
table(backgene_set$ont)
backgene_set$ont<-factor(backgene_set$ont,levels=c("BIOCARTA_MAPK_PATHWAY","BIOCARTA_ERK_PATHWAY"),ordered=TRUE)

target_GO_GSEA <- GSEA(geneList,TERM2GENE = backgene_set,  pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                       #minGSSize = 10,maxGSSize = 500,pvalueCutoff = 0.05,
                       pAdjustMethod = 'BH')  #指定 p 值校正方法
gseaplot2(target_GO_GSEA, color = colorspace::rainbow_hcl(4),geneSetID = rownames(target_GO_GSEA[1,]),title=(paste0(Group[2],"_vs_",Group[1],"_BIOCARTA_MAPK_PATHWAY_set")), pvalue_table = TRUE)


P1<-gseaplot2(target_GO_GSEA, color = ppCor[1:2],geneSetID = c("BIOCARTA_MAPK_PATHWAY","BIOCARTA_ERK_PATHWAY"),title=(paste0(Group[2],"_vs_",Group[1],"_two target pathway")),rel_heights = c(1.3, 0.3, 0.6), pvalue_table = T)
P1
ggsave(P1,file=paste0("/share/home/chenwei/LZJ_SCARE/GSEA_result/",Group[2],"_vs_",Group[1],"_mouse_gsea_MAPK_SIGNALING_PATHWAY.pdf"),width = 12,height =8,limitsize = F)


##############################################
###所有通路分析
## 物种设置
organism = 'rno'    #  人类'hsa' 小鼠'mmu'   
OrgDb = 'org.Rn.eg.db'#人类"org.Hs.eg.db" 小鼠"org.Mm.eg.db" 大鼠：org.Rn.eg.db

##4. GSEA分析（Mann-Whitney U检验,也被称为Wilcoxon秩和检验）
head(colData)
group_info<- colData[,c("Day","Treat","Mouse_ID")]
colnames(group_info)<-c("Day","Group","sample")

group_info$pair_id<-as.character(unlist(lapply(strsplit(as.character(group_info$sample),"_"), function(x) paste0(x[2],"_",x[3]))))
group_info$pair_id<-factor(group_info$pair_id,levels =unique(group_info$pair_id))
group_info<-group_info[order(group_info$pair_id,group_info$Group,decreasing = F),]

#group_info$group <- as.numeric(factor(group_info$Group))  
results <- data.frame()


# 两两比较  不考虑匹配
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
    
    expression_matrix1<-expression_matrix0[,c(class1_samples,class2_samples)]
    colData_select<-group_info[which(group_info$Group %in% c(Class1,Class2)),]
    colData_select$Group<-factor(colData_select$Group,levels = c(Class1,Class2))
    #head(colData_select)
    
    dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix1, colData = colData_select, design = ~ Group)
    dds1 <- dds0[ rowSums(counts(dds0)) > 1, ]
    dds2 <- estimateSizeFactors(dds1)
    dds3 <- DESeq(dds2)
    
    sampleA <-Class2; sampleB<-Class1
    contrastV <- c("Group",sampleA, sampleB)
    res1 <- as.data.frame(results(dds3, contrast=contrastV))
    summary(res1)
    need_DEG <- res1[,c(2,5)] #选择log2FoldChange和pvalue（凑成数据框）
    colnames(need_DEG) <- c('log2FoldChange','pvalue')
    need_DEG$SYMBOL <- rownames(need_DEG)
    ##### 创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）
    #####转化id 
    df <- bitr(rownames(need_DEG),fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = OrgDb) #人数据库org.Hs.eg.db 小鼠org.Mm.eg.db
    need_DEG <- merge(need_DEG, df, by='SYMBOL')  #按照SYMBOL合并注释信息
    geneList <- need_DEG$log2FoldChange
    names(geneList) <- need_DEG$ENTREZID
    geneList <- sort(geneList, decreasing = T)   #从大到小排序
    
    ##3. 利用clusterProfiler包进行GSEA富集（gseGO()和gseKEGG()函数可以很方便地对GO与KEGG通路进行GSEA， 再使用DOSE::setReadable转化id 。）
    ##### gsea富集 
    ####
    KEGG_kk_entrez <- gseKEGG(geneList= geneList,organism= organism,pvalueCutoff = 0.25)  #实际为padj阈值可调整 
    KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez,OrgDb=OrgDb,keyType='ENTREZID')#转化id               
    GO_kk_entrez <- gseGO(geneList= geneList,ont="BP", OrgDb = OrgDb,keyType = "ENTREZID",pvalueCutoff = 0.25)   #实际为padj阈值可调整 # "BP"、"MF"和"CC"或"ALL"   
    GO_kk <- DOSE::setReadable(GO_kk_entrez,OrgDb=OrgDb,keyType='ENTREZID')#转化id 
    
    save(KEGG_kk_entrez, GO_kk_entrez, file = paste0("/share/home/chenwei/LZJ_SCARE/GSEA_result/",sampleA,"_vs_",sampleB,"_GSEA_result.RData"))
    
  }}

# 两两比较  考虑匹配
for (i in c("D3","D7","D20")) { 
  #i="D3"
  colData_select<-group_info[which(group_info$Day == i),]
  colData_select$Group<-factor(colData_select$Group,levels = c( "CTRL","PBS"))
  colData_select$pair_id<-factor(colData_select$pair_id,levels =unique(colData_select$pair_id))
  
  expression_matrix1<-expression_matrix0[,colData_select$sample]
    #head(colData_select)
  dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix1, colData = colData_select, design = ~ pair_id + Group)
  dds1 <- dds0[ rowSums(counts(dds0)) > 1, ]
  dds2 <- estimateSizeFactors(dds1)
  dds3 <- DESeq(dds2)
  
  Treat_select<-c("PBS","CTRL")
  sampleA <-Treat_select[1]; sampleB<-Treat_select[2]
  contrastV <- c("Group",sampleA,sampleB)
  res1 <- as.data.frame(results(dds3, contrast=contrastV))
  summary(res1)
  need_DEG <- res1[,c(2,5)] #选择log2FoldChange和pvalue（凑成数据框）
  colnames(need_DEG) <- c('log2FoldChange','pvalue')
  need_DEG$SYMBOL <- rownames(need_DEG)
  ##### 创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）
  #####转化id 
  df <- bitr(rownames(need_DEG),fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = OrgDb) #人数据库org.Hs.eg.db 小鼠org.Mm.eg.db
  need_DEG <- merge(need_DEG, df, by='SYMBOL')  #按照SYMBOL合并注释信息
  geneList <- need_DEG$log2FoldChange
  names(geneList) <- need_DEG$ENTREZID
  geneList <- sort(geneList, decreasing = T)   #从大到小排序
  
  ##3. 利用clusterProfiler包进行GSEA富集（gseGO()和gseKEGG()函数可以很方便地对GO与KEGG通路进行GSEA， 再使用DOSE::setReadable转化id 。）
  ##### gsea富集 
  ####
  KEGG_kk_entrez <- gseKEGG(geneList= geneList,organism= organism,pvalueCutoff = 0.25)  #实际为padj阈值可调整 
  KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez,OrgDb=OrgDb,keyType='ENTREZID')#转化id               
  GO_kk_entrez <- gseGO(geneList= geneList,ont="BP", OrgDb = OrgDb,keyType = "ENTREZID",pvalueCutoff = 0.25)   #实际为padj阈值可调整 # "BP"、"MF"和"CC"或"ALL"   
  GO_kk <- DOSE::setReadable(GO_kk_entrez,OrgDb=OrgDb,keyType='ENTREZID')#转化id 
  
  save(KEGG_kk_entrez, GO_kk_entrez, file = paste0("/share/home/chenwei/LZJ_SCARE/GSEA_result/paired_",i,"_",sampleA,"_vs_",sampleB,"_GSEA_result.RData"))
}

###########################
##GSEA富集结果可视化
####GSEA的可视化主要是GSEA结果图、 gsearank plot和ridgeplot山脊图。
####同样也可以进行其他可视化如barplot、dotplot、cnetplot等等
Group<-c("CTRL","PBS")
Day<-"D20"
Day<-"D3"
Day<-"D7"

load(paste0("/share/home/chenwei/LZJ_SCARE/GSEA_result/paired_",Day,"_PBS_vs_CTRL_GSEA_result.RData")) 
KEGG_kk_entrez;GO_kk_entrez
KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez,OrgDb=OrgDb,keyType='ENTREZID')#转化id               
GO_kk <- DOSE::setReadable(GO_kk_entrez,OrgDb=OrgDb,keyType='ENTREZID')#转化id 

##4.1 gseaplot GSEA结果图下面
##一般认为|NES|>1，NOM pvalue<0.05，FDR（padj）<0.25的通路是显著富集的；
##还可以从结果中细分出上下调通路单独绘图，以下代码仅展示KEGG通路富集结果的上调通路。

##选取富集结果
kk_gse <- KEGG_kk;kk_gse_entrez <- KEGG_kk_entrez
row_names <- unique(kk_gse_entrez[grepl(toupper("mapk"),toupper(kk_gse_entrez$Description)),]$Description)
ID_select<-unique(kk_gse_entrez[grepl(toupper("mapk"),toupper(kk_gse_entrez$Description)),]$ID)
row_names

KEGG_MAPK_signaling_pathway <- kk_gse[which(kk_gse$Description %in%  row_names),]
write.table(KEGG_MAPK_signaling_pathway,file=paste0("/share/home/chenwei/LZJ_SCARE/GSEA_result/paired_",Day,"_PBS_vs_CTRL_rat_KEGG_MAPK_signaling_pathway.txt"),sep="\t", quote=F, row.names=T,col.names=T)

#### 合并 GSEA通路 
gseap2 <- gseaplot2(kk_gse, ID_select,title = paste0("paired_",Day,"_PBS_vs_CTRL_KEGG MAPK signaling pathway"),color =ppCor[1:10],#GSEA线条颜色 
                    base_size = 20,rel_heights = c(1.5, 0.5, 1),subplots = 1:3, ES_geom = "line",#enrichment score用线还是用点"dot" 
                    pvalue_table = T) #显示pvalue等信息
gseap2
ggsave(gseap2,file=paste0("/share/home/chenwei/LZJ_SCARE/GSEA_result/paired_",Day,"_PBS_vs_CTRL_gsea_rat_KEGG_MAPK_signaling_PATHWAY.pdf"),width = 12,height =8,limitsize = F)


###条件筛选
kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.25 & abs(kk_gse$NES)>1]
kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]
write.table(kk_gse_cut,file=paste0("/share/home/chenwei/LZJ_SCARE/GSEA_result/paired_",Day,"_PBS_vs_CTRL_rat_All_significant_KEGG_pathway.txt"),sep="\t", quote=F, row.names=T,col.names=T)


##选取富集结果
GO_gse <- GO_kk;GO_gse_entrez <- GO_kk_entrez
head(GO_gse)

unique(GO_kk_entrez[grepl(toupper("mapk"),toupper(GO_kk_entrez$Description)),]$Description)
row_names <- unique(GO_kk_entrez[grepl(toupper("mapk"),toupper(GO_kk_entrez$Description)),]$Description)
ID_select<-unique(GO_kk_entrez[grepl(toupper("mapk"),toupper(GO_kk_entrez$Description)),]$ID)
#row_names <- unique(GO_kk_entrez[grepl(toupper("healing"),toupper(GO_kk_entrez$Description)),]$Description)
#ID_select<-unique(GO_kk_entrez[grepl(toupper("healing"),toupper(GO_kk_entrez$Description)),]$ID)
row_names
GO_MAPK_signaling_pathway<-GO_gse[which(GO_gse$Description %in%  row_names),]
write.table(GO_MAPK_signaling_pathway,file=paste0("/share/home/chenwei/LZJ_SCARE/GSEA_result/paired_",Day,"_PBS_vs_CTRL_rat_GO_MAPK_signaling_pathway.txt"),sep="\t", quote=F, row.names=T,col.names=T)

###条件筛选
GO_gse_cut <- GO_gse[GO_gse$pvalue<0.05 & GO_gse$p.adjust<0.25 & abs(GO_gse$NES)>1]
GO_gse_cut_down <- GO_gse_cut[GO_gse_cut$NES < 0,]
GO_gse_cut_up <- GO_gse_cut[GO_gse_cut$NES > 0,]
write.table(GO_gse_cut,file=paste0("/share/home/chenwei/LZJ_SCARE/GSEA_result/paired_",Day,"_PBS_vs_CTRL_rat_All_significant_GO_term.txt"),sep="\t", quote=F, row.names=T,col.names=T)

#选择展现NES前几个通路 
down_gsea <- GO_gse_cut_down[tail(order(GO_gse_cut_down$NES,decreasing = T),10),]
up_gsea <- GO_gse_cut_up[head(order(GO_gse_cut_up$NES,decreasing = T),10),]
diff_gsea <- GO_gse_cut[head(order(abs(GO_gse_cut$NES),decreasing = T),10),]

#### 经典的GSEA图 
gseap0 <- gseaplot2(GO_gse, ID_select,title = paste0("paired_",Day,"_PBS_vs_CTRL_GO MAPK cascade"),color =ppCor[1:length(row_names)],#GSEA线条颜色 
                    base_size = 20,rel_heights = c(1.5, 0.5, 1),subplots = 1:3, ES_geom = "line",#enrichment score用线还是用点"dot" 
                    pvalue_table = T) #显示pvalue等信息
gseap0
ggsave(gseap0,file=paste0("/share/home/chenwei/LZJ_SCARE/GSEA_result/paired_",Day,"_PBS_vs_CTRL_gsea_rat_GO_MAPK_cascade_term.pdf"),width = 12,height =8,limitsize = F)









######本次未跑#####

############
##4.3 ridgeplot山脊图
## ridgeplot
ridgep1<- ridgeplot(KEGG_kk_entrez,showCategory = 15,fill = "p.adjust",core_enrichment = TRUE, label_format = 30,orderBy = "NES",decreasing = F) 
ridgep1
ridgep2<- ridgeplot(GO_kk_entrez,showCategory = 15,fill = "p.adjust",core_enrichment = TRUE, label_format = 30,orderBy = "NES",decreasing = F) 
ridgep2
#ggsave(ridgep,filename = 'ridgeplot.pdf',width =10,height =10)

#选择展现NES前几个通路 
down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES,decreasing = T),10),]
up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T),10),]
diff_gsea <- kk_gse_cut[head(order(abs(kk_gse_cut$NES),decreasing = T),10),]

#### 经典的GSEA图 
up_gsea$Description
i=2
gseap1 <- gseaplot2(kk_gse,up_gsea$ID[i],title = up_gsea$Description[i],#标题
                    color = "red", base_size = 20,rel_heights = c(1.5, 0.5, 1),#副图的相对高
                    subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line", pvalue_table = T) #显示pvalue等信息
##ggsave(gseap1, filename = 'GSEA_up_1.pdf', width =10, height =8) 
#### 合并 GSEA通路 
gseap2 <- gseaplot2(kk_gse, ID_select,title = "UP_GSEA_all",color =ppCor[1:10],#GSEA线条颜色 
                    base_size = 20,rel_heights = c(1.5, 0.5, 1),subplots = 1:3, ES_geom = "line",#enrichment score用线还是用点"dot" 
                    pvalue_table = T) #显示pvalue等信息

##ggsave(gseap2, filename = "GSEA_up_all.pdf",width =12,height =12)

#####################
gseap01 <- gseaplot2(GO_gse, "GO:0042060",title = "wound healing",color ="red",#GSEA线条颜色 
                     base_size = 20,rel_heights = c(1.5, 0.5, 1),subplots = 1:3, ES_geom = "line",#enrichment score用线还是用点"dot" 
                     pvalue_table = T) #显示pvalue等信息
gseap01

up_gsea$Description
i=2
gseap1 <- gseaplot2(GO_gse,up_gsea$ID[i],title = up_gsea$Description[i],#标题
                    color = "red", base_size = 20,rel_heights = c(1.5, 0.5, 1),#副图的相对高
                    subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line", pvalue_table = T) #显示pvalue等信息
##ggsave(gseap1, filename = 'GSEA_up_1.pdf', width =10, height =8) 
#### 合并 GSEA通路 
gseap2 <- gseaplot2(GO_gse, up_gsea$ID,title = "UP_GSEA_all",color =ppCor[1:10],#GSEA线条颜色 
                    base_size = 20,rel_heights = c(1.5, 0.5, 1),subplots = 1:3, ES_geom = "line",#enrichment score用线还是用点"dot" 
                    pvalue_table = T) #显示pvalue等信息

###############
##提取目标通路的基因
row_names
#[1]  "positive regulation of MAPK cascade" "negative regulation of MAPK cascade"
pathway_of_interest <-"regulation of mapk"    
genes_in_pathway <- unlist(geneIdsByCategory(GO_gse, pathway_of_interest))  ## 提取该通路中的所有基因  


# 将基因ID转换为人类基因符号  
gene_symbols <- mapIds(org.Hs.eg.db, keys=genes_in_pathway, keytype="ENTREZID", column="SYMBOL")  

# 输出结果  
print(gene_symbols)  
