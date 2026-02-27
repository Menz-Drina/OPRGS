#GSVA----
library(GSVA)
ex <- exp
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex+1)
print("log2 transform finished")}else{print("log2 transform not needed")}
gsva_data <- GSE36389 
load("KEGG-84-metabolism.rdata")
gene_set <- geneset
params <- gsvaParam(as.matrix(gsva_data), gene_set, minSize = 1, 
                    maxSize = Inf, kcdf = "Gaussian", tau = 1, maxDiff = TRUE, 
                    absRanking = FALSE)

gsva_result <- gsva(params, verbose = TRUE, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE))
#GSVA-DEG----
library(limma)
library(stringr)
group_list = ifelse(str_detect(colnames(gsva_result),"GTEX"),"Normal","Cancer")
group_list = ifelse(str_detect(substr(colnames(gsva_result),14,16),"11A"),"Normal","Cancer")
group_list <- pd$group
group_list = factor(group_list,levels = c("Normal","Cancer"))
table(group_list)
design = model.matrix(~group_list)
fit = lmFit(gsva_result, design)
fit = eBayes(fit)
DEG = topTable(fit, coef = 2, number = Inf)
head(DEG)
DEG <- DEG[DEG$P.Value<0.05,]
GSE26712_DEG <- DEG
save.image(file = "GSE26712-GSVA.rdata")
save(GSE26712_DEG,file = "GSVA-GSE26712.rdata")
#The processing of other datasets is similar
rm(list = ls())
#GSVA-DEG-RRA----
library(clusterProfiler)
library(RobustRankAggreg)
library(pheatmap)
deg1 <- TCGA_GETX 
deg2 <- GSE66957  
deg3 <- GSE26712  
deg4 <- GPL570

get_up_deg1 <- deg1[deg1$change=="UP",]
get_up_deg1<- get_up_deg1[order(get_up_deg1$log2FoldChange,decreasing = T),]
get_up_deg2 <- deg2[deg2$change=="UP",]
get_up_deg2<- get_up_deg2[order(get_up_deg2$log2FoldChange,decreasing = T),]
get_up_deg3 <- deg3[deg3$change=="UP",]
get_up_deg3<- get_up_deg3[order(get_up_deg3$log2FoldChange,decreasing = T),]
get_up_deg4 <-deg4[deg4$change=="UP",]
get_up_deg4<- get_up_deg4[order(get_up_deg4$log2FoldChange,decreasing = T),]

get_up_deg1 <- rownames(get_up_deg1)
get_up_deg2 <- rownames(get_up_deg2)
get_up_deg3 <- rownames(get_up_deg3)
get_up_deg4 <- rownames(get_up_deg4)

glist_up=list(get_up_deg1,get_up_deg2,get_up_deg3,get_up_deg4)
ups=aggregateRanks(glist = glist_up,N=length(unique(unlist(glist_up))))
tmp_up <- as.data.frame(table(unlist(glist_up)))
ups$Freq <- tmp_up[match(ups$Name,tmp_up[,1]),2]
gs <- ups[ups$Score<0.05,1]
updat=data.frame(deg1=deg1[gs,'log2FoldChange'],
                 deg2=deg2[gs,'log2FoldChange'],
                 deg3=deg3[gs,'log2FoldChange'],
                 deg4=deg4[gs,'log2FoldChange'])
rownames(updat)=gs
head(updat)
colnames(updat) <- c('TCGA-GTEX','GSE26712','GSE66957','GPL570')
pheatmap(updat,display_numbers=T)

#ConsensusClusterPlus----
library(ConsensusClusterPlus)
library(tinyarray)
library(tidyverse)
##step1
exp <- cancer_vst
coo <- intersect(rownames(exp),Necropotosis)
exp <- exp[coo,rownames(meta)]
identical(rownames(meta),colnames(exp))
exp = as.matrix(exp)
df <-  sweep(exp,1, apply(exp,1,median,na.rm=T))
df1 <- df
##step2.1
maxK <-  6 
results <-  ConsensusClusterPlus(df1,
                                 maxK = maxK,
                                 reps = 500,
                                 pItem = 0.8,
                                 pFeature = 1,
                                 clusterAlg = "pam",
                                 seed = 985,
                                 title="Necropotosis",
                                 innerLinkage="complete",
                                 plot="pdf")

icl = calcICL(results,
              title="Necropotosis",
              plot="pdf")
##step2.2
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") 

for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}

optK = Kvec[which.min(PAC)]
optK
table(results[[optK]]$consensusClass)
Cluster = results[[optK]]$consensusClass
identical(names(Cluster),rownames(meta))
meta$Cluster = Cluster
library(survival)
library(survminer)
sfit <- survfit(Surv(time, event) ~ Cluster,
                data = meta)
ggsurvplot(sfit, surv.median.line = "hv", 
           legend.title = "Strata",
           legend.labs = c("Cluster = 1", "Cluster = 2"),
           pval = TRUE, 
           conf.int = TRUE,
           palette = c("#E7B800", "#2E9FDF"),  
           font.legend=18,
           font.title = c(20, "black"),
           font.subtitle = c(20, "black"),
           font.caption = c(20, "black"),
           font.x = c(20,"black"),
           font.y = c(20, "black"),
           font.tickslab = c(18, "black"),
           ggtheme = theme_bw())
##step4
library(Rtsne)
set.seed(2023)
tsne_out = Rtsne(t(exp),perplexity = 30)
pdat = data.frame(tsne_out$Y,factor(Cluster))
colnames(pdat) = c("Dim1","Dim2","Cluster")
library(ggplot2)
library(paletteer)
ggplot(pdat,aes(Dim1,Dim2))+
  geom_point(aes(Dim1,Dim2,fill = Cluster),shape = 21,color = "black")+
  stat_ellipse(aes(color = Cluster,fill = Cluster),
               geom = "polygon",
               alpha = 0.3,
               linetype = 2,
               level = 0.9)+
  scale_color_paletteer_d("RColorBrewer::Set3")+
  scale_fill_paletteer_d("RColorBrewer::Set3")+
  theme_classic()+
  theme(legend.position = "top", axis.text = element_text(size = 18,color="black"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=18))
#DEA----
suppressMessages(library(genefilter))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(TCGAbiolinks))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(ggplotify))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(AnnoProbe))
suppressMessages(library(tinyarray))
suppressMessages(library(GEOquery)) 
suppressMessages(library(ggplot2)) 
suppressMessages(library(ggstatsplot)) 
suppressMessages(library(reshape2))
suppressMessages(library(FactoMineR))
suppressMessages(library(factoextra))
suppressMessages(library(pheatmap))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(UpSetR))
suppressMessages(library(WGCNA))

counts1 <- SummarizedExperiment(counts) 
selectGenes <- function(counts, min.count=10, N=0.90){
  
  lib.size <- colSums(counts)
  MedianLibSize <- median(lib.size)
  CPM.Cutoff <- min.count / MedianLibSize*1e6
  CPM <- edgeR::cpm(counts,lib.size=lib.size)
  
  min.samples <- round(N * ncol(counts))
  
  f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
  flist <- genefilter::filterfun(f1)
  keep <- genefilter::genefilter(CPM, flist)
  return(keep)
}
keep.exprs <- selectGenes(assay(counts1), min.count=10, N=0.90)
myFilt <- counts1[keep.exprs,]
dim(myFilt)
exp <- counts[myFilt@NAMES,]
dat <- exp
library(stringr)
group_list = ifelse(str_detect(colnames(exp),"GTEX"),"Normal","Cancer")
group_list = factor(group_list,levels = c("Normal","Cancer"))
table(group_list)
cancer_type <- "TCGA-GTEX"
dat <- exp
exp=t(exp)
exp=as.data.frame(exp)
dat.pca <- PCA(exp , graph = FALSE)
gse_number <- cancer_type
## PCA
this_title <- paste0(gse_number,'-PCA')
p1 <- fviz_pca_ind(dat.pca,
                   geom.ind = "point", 
                   col.ind = group_list, 
                   palette = "Dark2",
                   addEllipses = TRUE, 
                   legend.title = "Groups")+
  ggtitle(this_title)+
  theme_ggstatsplot()+
  theme(plot.title = element_text(size=20,hjust = 0.5),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text = element_text(color="black",size=16),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18))
p1
exp <- dat
#Deseq2
colData <- data.frame(row.names =colnames(exp), 
                      condition=group_list)
if(!file.exists(paste0(cancer_type,"dd.Rdata"))){
  dds <- DESeqDataSetFromMatrix(
    countData = exp,
    colData = colData,
    design = ~ condition)
  dds <- DESeq(dds)
  save(dds,file = paste0(cancer_type,"dd.Rdata"))
}

res <- results(dds, contrast = c("condition",rev(levels(group_list))))
resOrdered <- res[order(res$padj),]
DEG <- as.data.frame(resOrdered)
head(DEG)

logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)))
k1 = (DEG$padj < 0.05)&(DEG$log2FoldChange <= -logFC_cutoff)
k2 = (DEG$padj < 0.05)&(DEG$log2FoldChange >= logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
DEG$symbol <- rownames(DEG)
View(DEG)
table(DEG$change)
TCGA_DEG <- DEG
save(DEG,group_list,file = paste0(cancer_type,"-Deseq2.Rdata"))
#Venn----
library(ggvenn)
a <- list(rownames(DEG1),rownames(DEG2))
p1 <- ggvenn(a, columns = c("TCGA+GTEx", "Oxidative Phosphorylation"))           
#Cox----
library(stringr)
library(AnnoProbe)
library(tinyarray)
data <- surv_cox(exprSet_hub1,meta1,cut.point = T)
#Mime1----
suppressMessages(library(Mime1))
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(stringi))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(beepr))
suppressMessages(library(data.table))
suppressMessages(library(pheatmap))
suppressMessages(library(ggsignif))
suppressMessages(library(RColorBrewer))
suppressMessages(library(future.apply))
suppressMessages(library(gplots))
suppressMessages(library(DESeq2))
suppressMessages(library(ggrepel))
suppressMessages(library(Rcpp))
suppressMessages(library(survivalsvm))
suppressMessages(library(dplyr))
suppressMessages(library(rms))
suppressMessages(library(pec))
suppressMessages(library(ggDCA))
suppressMessages(library(glmnet))
suppressMessages(library(foreign))
suppressMessages(library(regplot))
suppressMessages(library(randomForestSRC))
suppressMessages(library(timeROC))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(caret))
suppressMessages(library(regplot))
suppressMessages(library(gbm))
suppressMessages(library(tidyverse))
train <- exp
trainlist=list(train=train,test1=exp1,test3=exp3,test4=exp4,test5=exp5)
genelist=colnames(exp)[-c(1:3)]
res <- ML.Dev.Prog.Sig(train_data = trainlist$train,
                       list_train_vali_Data = trainlist,
                       unicox.filter.for.candi = F,
                       candidate_genes = genelist,
                       mode = 'all',nodesize =5,seed = 2025 )
list_train_vali_Data=trainlist
cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],order =names(list_train_vali_Data),width = 0.35)
model <- res$ml.res$`StepCox[both] + SuperPC`
fit <- res$riskscore$plsRcox
train <- fit$train
test1 <- fit$test1
test2 <- fit$test2
test3 <- fit$test3
res.cut <- surv_cutpoint(test2, time = "OS.time", event = "OS",
                         variables = c("RS") )


res.cat <- surv_categorize(res.cut)

fit <- survfit(Surv(OS.time, OS) ~ RS, data = res.cat)
ggsurvplot(fit, data = res.cat, surv.median.line = "hv", 
           legend.title = "Risk Group",
           legend.labs = c("High Risk", "Low Risk"),
           pval = TRUE, 
           conf.int = TRUE,
           palette = c("#E7B800", "#2E9FDF"),
           font.legend=16,
           font.title = c(20, "black"),
           font.subtitle = c(20, "black"),
           font.caption = c(20, "black"),
           font.x = c(20,"black"),
           font.y = c(20, "black"),
           font.tickslab = c(16, "black"),
           ggtheme = theme_bw())
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["train"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.2y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["train"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 2,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["train"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.4y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["train"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 4,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["train"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")
# SAME AS OTHER DATASETS

#Boxplot----
library(tinyarray)
library(ggpubr)
library(tidyverse)
library(ggplot2)
selection <- genes
dat = t(vst[,selection])
p1 <- draw_boxplot(dat,
                   meta$group,
                   method="wilcox.test",
                   size=8)+ 
  facet_wrap(~rows,scales = "free",nrow = 1)+
  theme(legend.position = "top", axis.text = element_text(size = 18,color="black"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=18),
        strip.text = element_text(size=18))
pdf("gene-group.pdf",width = 10.00, height=14.00) 
p1
dev.off()
#Meta----
library(meta)
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,
                                   optimal.model = "RSF + SuperPC",
                                   type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
meta_result <- meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))
meta_trimfill <- trimfill(meta_result)
summary(meta_trimfill)
funnel(meta_trimfill,
       studlab = TRUE,
       pch = ifelse(meta_trimfill$trimfill, 17, 19),  
       col = ifelse(meta_trimfill$trimfill, "red", "blue"),  
       main = "Trim-and-Fill Adjusted Funnel Plot",
       xlab = "HR",
       ylab = "Standard Error",
       level = 0.95,       
       contour = c(0.9, 0.95, 0.99),  
       col.contour = c("darkgreen", "green", "lightgreen"))
legend(1,0.2,       
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"),bty="n",
       fill = c("darkgreen", "green", "lightgreen")) 


loo_result <- metainf(meta_result,
                      pooled='common') 
forest(loo_result,
       layout = "RevMan5", 
       col.bg = "#1E90FF",   
       col.common="#DC143C",
       common=TRUE,
       digits=2,
       digits.tau=2,
       digits.I2=2,
       digits.pval=2)
#Contrast----
OV <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = T,
                              collected_sig_table = contrast,
                              list_input_data = list_train_vali_Data)
cindex_comp(OV,
            res1,
            model_name="RSF + SuperPC",
            dataset=names(list_train_vali_Data))

#Bubble Plot and Corelation-----
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),
                    collapse = ""))
}

randomColors <- replicate(200,randomColor())



corResult <- pdat
corResult$immune <-rownames(corResult)
corResult <- corResult%>%
  dplyr::select("immune", everything())
corResult1 <- corResult
corResult1 <- separate(corResult1,celltype,into = c('a','b',"c","d","e","f"),"_")

write.csv(corResult1,file = "corResult1.csv")
corResult1 <- fread("corResult1.csv")
corResult1 <- as.data.frame(corResult1)
rownames(corResult1) <- corResult1$V1
corResult1 <- corResult1[,-1]
corResult <- corResult1
corResult$Software=factor(corResult$Software,level=as.character(unique(corResult$Software[rev(order(as.character(corResult$Software)))])))
b=corResult[order(corResult$Software),]
b$immune=factor(b$immune,levels=rev(as.character(b$immune)))
colslabels=rep(hue_pal()(length(levels(b$Software))),table(b$Software))     

ggplot(data=b, aes(x=cor, y=immune, color=Software))+
  labs(x="Correlation coefficient",y="Immune cell")+
  geom_point(size=4.1)+
  theme(panel.background=element_rect(fill="white",size=1,color="black"),
        panel.grid=element_line(color="grey75",size=0.5),
        axis.ticks = element_line(size=0.5),
        axis.text.y = ggtext::element_markdown(colour=rev(colslabels)))

identical(rownames(immu_expr),rownames(TCGA1))
immu_expr$Risk_score <- TCGA1$RS
immu_expr$ID <- rownames(immu_expr)
df_riskscore <- immu_expr%>%
  dplyr::select(Risk_score,ID, everything())%>%
  pivot_longer(-c(1,2),names_to = "cell_type",values_to = "score")

ggplot(df_riskscore, aes(RS,score))+
  geom_point()+
  geom_smooth(method = "lm",color="blue")+
  stat_cor(method = "spearman",color="red")+
  theme_classic()+
  theme(legend.position = "top", axis.text = element_text(size = 18,color="black"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=18))+
  facet_wrap(~cell_type,scales = "free_y",ncol = 5)
#Estimate-----
suppressMessages(library(IOBR))
suppressMessages(library(ggplot2))
suppressMessages(library(linkET))
suppressMessages(library(tidyverse))
suppressMessages(library(corrplot))
suppressMessages(library(ggcor))
suppressMessages(library(ggnewscale))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(forcats))
im_estimate <- deconvo_tme(eset = exp1,
                           method = "estimate")
#Immunotherapy----
rm(list = ls())
library(superpc)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ggstatsplot)
model <- res$ml.res$`RSF + SuperPC`
data <- list(x=t(trainlist$train[,-c(1:3)]),
             y=trainlist$train$OS.time,
             censoring.status=trainlist$train$OS,
             featurenames=colnames(trainlist$train)[-c(1:3)])
expr <- eset
pdata <- phe
expr <- expr[genelist,]
expr <- t(expr)
test <- list(x=t(expr),
             y=pdata$OS,
             censoring.status=pdata$status,
             featurenames=colnames(expr))
pred_superpc <- superpc.predict(model[[1]],data,
                                test,
                                threshold = model[[2]]$thresholds[which.max(model[[2]][["scor"]][1,])],
                                n.components = 1)

pred_superpc <- as.numeric(pred_superpc$v.pred)
#
expr <- as.data.frame(expr)
expr$RS <- pred_superpc
identical(rownames(expr),rownames(pdata))
expr$OS.time <- pdata$OS
expr$OS <- pdata$status
#
res.cut <- surv_cutpoint(expr, time = "OS.time", event = "OS",
                         variables = c("RS") )

res.cat <- surv_categorize(res.cut)
phe <- pdata[,1:2]
phe <- na.omit(phe)
coo <- intersect(rownames(phe),rownames(expr))
phe <- phe[coo,]
expr <- expr[coo,]
res.cat <- res.cat[coo,]
identical(rownames(phe),rownames(expr))
identical(rownames(phe),rownames(res.cat))
phe$RS <- expr$RS
phe$riskgroup <- res.cat$RS
phe$Response <-pdata$response
fit <- survfit(Surv(OS.time, OS) ~ RS, data = res.cat)

ggsurvplot(fit, data = res.cat, surv.median.line = "hv", 
           legend.title = "Risk Group",
           legend.labs = c("High Risk", "Low Risk"),
           pval = TRUE, 
           conf.int = TRUE,
           palette = c("#E7B800", "#2E9FDF"),
           font.legend=16,
           font.title = c(20, "black"),
           font.subtitle = c(20, "black"),
           font.caption = c(20, "black"),
           font.x = c(20,"black"),
           font.y = c(20, "black"),
           font.tickslab = c(16, "black"),
           ggtheme = theme_bw())
#
subtype <- phe$riskgroup
Responder <- phe$Response
dat <- data.frame(subtype,Responder)
table(dat$subtype)
dat$subtype <- factor(dat$subtype,levels = c('low','high'))
dat = dplyr::count(dat,subtype,Responder)
dat = dat %>% group_by(subtype) %>% 
  summarise(Responder = Responder,n = n/sum(n))
table(dat$Responder)
dat$Responder = factor(dat$Responder,levels = c("Responder","Non-Responder"))
dat <- as.data.frame(dat)
f = chisq.test(table(subtype,Responder),correct = F)  #chisq #fisher.test
f = fisher.test(table(subtype,Responder)) 
label = paste("chisq.test P value =",round(f$p.value,3))
label

#"#0571B0","#a93a50","#7570B3","#1B9E77"
p1 <- ggplot(data = dat)+
  geom_bar(aes(x = subtype,y = n,
               fill = Responder),
           stat = "identity")+
  scale_fill_manual(values =  c("#a93a50","#0571B0"))+
  geom_label(aes(x = subtype,y = n,
                 label = scales::percent(n),
                 fill = Responder),
             color = "white",
             size = 4,label.size = 0,
             show.legend = F,
             position = position_fill(vjust = 0.5))+
  annotate("text", x = 1.1, y = 1.1, label = label,size = 4) +
  theme(plot.title = element_text(size=16,hjust = 0.5),
        axis.text = element_text(size = 16,color="black"),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title = element_text(size=16), 
        legend.text = element_text(size=16))+
  theme_bw(base_size = 20)+
  guides(label = "none")
p1
#TIDE----
View(paletteer::palettes_d_names)
library(ggstatsplot)
p1 <- ggbetweenstats(data=res,  
                     x = subtype,  
                     y = TIDE, 
                     title = "TIDE score",
                     xlab = "Subtype",
                     ylab = "TIDE score",
                     package = "ggsci",
                     palette = "category20_d3") 
pdf("TIDE score.pdf",width = 5.5, height=5.5) 
p1
dev.off()
p2 <- ggbetweenstats(data=res,  
                     x = subtype,  
                     y = IFNG, 
                     title = "Interferon-gamma score",
                     xlab = "Subtype",
                     ylab = "Interferon-gamma score",
                     package = "ggsci",
                     palette = "category20b_d3") 
pdf("IFNG score.pdf",width = 5.5, height=5.5) 
p2
dev.off()
p3 <- ggbetweenstats(data=res,  
                     x = subtype,  
                     y = Dysfunction, 
                     title = "T cell dysfunction score",
                     xlab = "Subtype",
                     ylab = "T cell dysfunction score",
                     package = "ggsci",
                     palette = "category20c_d3") 
pdf("T cell dysfunction score.pdf",width = 5.5, height=5.5) 
p3
dev.off()
p4 <- ggbetweenstats(data=res,  
                     x = subtype,  
                     y = Exclusion, 
                     title = "T cell exclusion score",
                     xlab = "Subtype",
                     ylab = "T cell exclusion score",
                     package = "ggsci",
                     palette = "hallmarks_dark_cosmic") 
pdf("T cell exclusion score.pdf",width = 5.5, height=5.5) 
p4
dev.off()
p5 <- ggbetweenstats(data=res,  
                     x = subtype,  
                     y = MDSC,
                     title = "Myeloid-derived suppressor cells score",
                     xlab = "Subtype",
                     ylab = "Myeloid-derived suppressor cells score") 

pdf("MDSC score.pdf",width = 5.5, height=5.5) 
p5
dev.off()
p6 <- ggbetweenstats(data=res,  
                     x = subtype,  
                     y = CAF,
                     title = "Cancer-associated fibroblasts score",
                     xlab = "Subtype",
                     ylab = "Cancer-associated fibroblasts score",
                     package = "ggsci",
                     palette = "hallmarks_light_cosmic") 

pdf("Cancer-associated fibroblasts score.pdf",width = 5.5, height=5.5) 
p6
dev.off()
p7 <- ggbetweenstats(data=res,  
                     x = subtype,  
                     y = `TAM M2`,
                     title = "M2-tumor-associated mocrophage score",
                     xlab = "Subtype",
                     ylab = "M2-tumor-associated mocrophage score",
                     package = "ggprism",
                     palette = "floral") 

pdf("TAM M2.pdf",width = 5.5, height=5.5) 
p7
dev.off()
p8 <- ggbetweenstats(data=res,
                     x = subtype, 
                     y = `MSI Expr Sig`,
                     title = "MSI score",
                     xlab = "Subtype",
                     ylab = "MSI score",
                     package = "ggprism",
                     palette = "muted_rainbow") 

pdf("MSI score.pdf",width = 5.5, height=5.5) 
p8
dev.off()

#TCIA----
library(ggstatsplot)
exprset <- t(cancer_vst) 
exprset <- as.data.frame(exprset)
sample_group <- TCGA$riskgroup
table(sample_group)
sample_group <- factor(sample_group,levels = c('Low','High'))
table(sample_group)
exprset$sample_subtypes <- sample_group
table(exprset$sample_subtypes)
rownames(exprset) <- substr(rownames(exprset),1,12)
TCIA <- TCIA[rownames(exprset),]
identical(rownames(TCIA),rownames(exprset))
TCIA$subtype <- exprset$riskgroup
TCIA$subtype <- exprset$sample_subtypes
table(TCIA$subtype)
TCIA$subtype <- factor(TCIA$subtype,levels = c('Low','High'))
TCIA$subtype <- factor(TCIA$subtype,levels = c('low','high'))
p1 <- ggbetweenstats(data=TCIA,  
                     x = subtype,  
                     y = `CTLA-4(-) & PD-1(-)`, 
                     title = "CTLA-4(-) & PD-1(-)",
                     xlab = "Subtype",
                     ylab = "CTLA-4(-) & PD-1(-)",
                     package = "ggsci",
                     palette = "category20_d3") 
pdf("CTLA-4(-) & PD-1(-).pdf",width = 5.5, height=5.5) 
p1
dev.off()
p2 <- ggbetweenstats(data=TCIA,  
                     x = subtype,  
                     y = `CTLA-4(-) & PD-1(+)`, 
                     title = "CTLA-4(-) & PD-1(+)",
                     xlab = "Subtype",
                     ylab = "CTLA-4(-) & PD-1(+)",
                     package = "ggsci",
                     palette = "category20b_d3") 
pdf("CTLA-4(-) & PD-1(+).pdf",width = 5.5, height=5.5) 
p2
dev.off()
p3 <- ggbetweenstats(data=TCIA,  
                     x = subtype,  
                     y = `CTLA-4(+) & PD-1(-)`, 
                     title = "CTLA-4(+) & PD-1(-)",
                     xlab = "Subtype",
                     ylab = "CTLA-4(+) & PD-1(-)",
                     package = "ggsci",
                     palette = "category20c_d3") 
pdf("CTLA-4(+) & PD-1(-).pdf",width = 5.5, height=5.5) 
p3
dev.off()
p4 <- ggbetweenstats(data=TCIA,  
                     x = subtype,  
                     y = `CTLA-4(+) & PD-1(+)`, 
                     title = "CTLA-4(+) & PD-1(+)",
                     xlab = "Subtype",
                     ylab = "CTLA-4(+) & PD-1(+)",
                     package = "ggsci",
                     palette = "hallmarks_dark_cosmic") 
pdf("CTLA-4(+) & PD-1(+).pdf",width = 5.5, height=5.5) 
p4
dev.off()

#Drug----
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggsci))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(oncoPredict))

res <- read.csv("DrugPredictions.csv")
drug <- c('paclitaxel',
          "Platin",
          "carboplatin",
          "oxaliplatin",
          "fluorouracil",
          "decitabine",
          "cyclophosphamide",
          "doxorubicin",
          "etoposide",
          "tamoxifen",
          "ifosfamide",
          "docetaxel",
          "olaparib",
          "gemcitabine")
res <- res[,-1]
res <- res[,drug]
res %>% 
  bind_cols(sample_group = sample_group) %>% 
  pivot_longer(1:14,names_to = "drugs",values_to = "ic50") %>% 
  ggplot(., aes(sample_group,ic50))+
  geom_boxplot(aes(fill=sample_group))+
  scale_fill_jama()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")+
  facet_wrap(vars(drugs),scales = "free_y",nrow = 2)+
  stat_compare_means()
#

#Violin----
library(ggplot2)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
#data was created from GSVA
im_ssgsea <- data
im_ssgsea %>% 
  mutate(sample_subtypes = factor(sample_subtypes)) %>% 
  pivot_longer(-c(ID,sample_subtypes), names_to = "cell_type",values_to = "value") %>% 
  ggplot(aes(cell_type,value,fill=sample_subtypes))+
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  theme_bw()+
  stat_compare_means(aes(group = sample_subtypes,label = ..p.signif..),
                     method = "kruskal.test",label.y = 18)+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,hjust = 1)
  )

ssgsea_list <- ssgsea_long%>% 
  split(.$cell_type) %>% 
  purrr::map(~ ggplot(., aes(cell_type,Score,fill=sample_subtypes))+
               geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
               theme_bw()+
               stat_compare_means(aes(group = sample_subtypes,label = ..p.signif..),
                                  method = "kruskal.test",label.y = 0.57)+
               theme(legend.position = "top",
                     axis.text.x = element_text(angle = 45,hjust = 1))+
               theme_bw())
paths <- paste0(names(ssgsea_list),".pdf")

pwalk(list(paths,ssgsea_list),ggsave,width=6,height=6)
#scRNA-seq----
suppressMessages(library(SingleR))
suppressMessages(library(ggsci))
suppressMessages(library(dplyr))
suppressMessages(library(future))
suppressMessages(library(Seurat))
suppressMessages(library(clustree))
suppressMessages(library(cowplot))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(stringr))
suppressMessages(library(qs))
suppressMessages(library(Matrix))
suppressMessages(library(GEOquery))
suppressMessages(library(magrittr))
suppressMessages(library(limma))
suppressMessages(library(harmony))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(S4Vectors))
suppressMessages(library(tibble))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(DESeq2))
suppressMessages(library(ggrepel))
suppressMessages(library(COSG))
suppressMessages(library(SCP))
suppressMessages(library(Nebulosa))
suppressMessages(library(fs))
suppressMessages(library(celldex))
suppressMessages(library(gplots))
suppressMessages(library(ggthemes))
suppressMessages(library(tinyarray))
suppressMessages(library(gridExtra))
suppressMessages(library(VennDetail))
suppressMessages(library(VennDiagram))
suppressMessages(library(scatterplot3d))
suppressMessages(library(singleseqgset))
suppressMessages(library(devtools))
suppressMessages(library(grid))
suppressMessages(library(ggnewscale))
suppressMessages(library(ggstatsplot))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(FactoMineR))
suppressMessages(library(factoextra))
suppressMessages(library(pheatmap))
suppressMessages(library(genefilter))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(TCGAbiolinks))
suppressMessages(library(edgeR))
suppressMessages(library(ggplotify))
suppressMessages(library(tidyr))
suppressMessages(library(AnnoProbe))
suppressMessages(library(GEOquery)) 
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(UpSetR))
suppressMessages(library(WGCNA))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(GSVA))
suppressMessages(library(presto))
suppressMessages(library(DoubletFinder))
suppressMessages(library(decontX))
suppressMessages(library(CellChat))
suppressMessages(library(infercnv))
suppressMessages(library(infercnvNGCHM))
suppressMessages(library(rjags))
dir='EC'
samples=list.files( dir )
samples 
sceList = lapply(samples,function(pro){ 
  print(pro)  
  tmp = Read10X(file.path(dir,pro )) 
  if(length(tmp)==2){
    ct = tmp[[1]] 
  }else{ct = tmp}
  sce =CreateSeuratObject(counts =  ct ,
                          project =  pro  ,
                          min.cells = 3,
                          min.features = 500 )
  return(sce)
}) 

do.call(rbind,lapply(sceList, dim))
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples  )
names(sce.all@assays$RNA@layers)
dim(sce.all[["RNA"]]$counts)
LayerData(sce.all, assay = "RNA", layer = "counts")
sce.all
sce.all <- JoinLayers(sce.all)
sce.all
dim(sce.all[["RNA"]]$counts)
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 
length(sce.all$orig.ident)
phe = sce.all@meta.data
table(phe$orig.ident)
phe$patient = phe$orig.ident
table(phe$patient)
phe$sample = phe$patient
table(phe$sample)
phe$sample <- ifelse(str_detect(phe$sample,"NORM"),"Normal","Cancer")
phe$sample = gsub("Cancer1|Cancer2|Cancer3|Cancer4|Cancer5|Cancer6|Cancer7", "Cancer", phe$sample) 
phe$sample = gsub("Norm1|Norm2|Norm3|Norm4|Norm5", "Normal", phe$sample) 
phe$group <- phe$sample
table(phe$group)
sce.all@meta.data = phe
save(sce.all,file = "OV-ORIGINAL.rdata")

load("OV-ORIGINAL.rdata")
sp='human'
dir.create("./1-QC")
setwd("./1-QC")
source('../scRNA_scripts/qc.R')
#qc
basic_qc <- function(input_sce){

  mito_genes=rownames(input_sce)[grep("^MT-", rownames(input_sce),ignore.case = T)] 
  print(mito_genes) 
  #input_sce=PercentageFeatureSet(input_sce, "^MT-", col.name = "percent_mito")
  input_sce=PercentageFeatureSet(input_sce, features = mito_genes, col.name = "percent_mito")
  fivenum(input_sce@meta.data$percent_mito)
  
  ribo_genes=rownames(input_sce)[grep("^Rp[sl]", rownames(input_sce),ignore.case = T)]
  print(ribo_genes)
  input_sce=PercentageFeatureSet(input_sce,  features = ribo_genes, col.name = "percent_ribo")
  fivenum(input_sce@meta.data$percent_ribo)

  Hb_genes=rownames(input_sce)[grep("^Hb[^(p)]", rownames(input_sce),ignore.case = T)]
  print(Hb_genes)
  input_sce=PercentageFeatureSet(input_sce,  features = Hb_genes,col.name = "percent_hb")
  fivenum(input_sce@meta.data$percent_hb)

  feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito",
             "percent_ribo", "percent_hb")
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
  p1 
  w=length(unique(input_sce$orig.ident))/3+5;w
  ggsave(filename="Vlnplot1.pdf",plot=p1,width = w,height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
    scale_y_continuous(breaks=seq(0, 100, 5)) +
    NoLegend()
  p2	
  w=length(unique(input_sce$orig.ident))/2+5;w
  ggsave(filename="Vlnplot2.pdf",plot=p2,width = w,height = 5)
  
  p3=FeatureScatter(input_sce, 
                    "nCount_RNA", 
                    "nFeature_RNA", 
                    group.by = "orig.ident", 
                    pt.size = 0.5)
  ggsave(filename="Scatterplot.pdf",plot=p3,width = 10,height = 10)
  #selected_c <- WhichCells(input_sce, expression = nFeature_RNA > 500 )
  selected_f <- rownames(input_sce)[Matrix::rowSums(input_sce@assays$RNA$counts > 0 ) > 3]
  input_sce.filt <- subset(input_sce, features = selected_f)
  input_sce.filt <- subset(input_sce.filt, subset = nFeature_RNA > 500 & nFeature_RNA < 8000)
  
  dim(input_sce) 
  dim(input_sce.filt) 
  
  
  C=subset(input_sce.filt,downsample=100)@assays$RNA$counts
  dim(C)
  C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
  
  most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
  
  pdf("TOP50_most_expressed_gene.pdf",width=8)
  boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
          cex = 0.1, las = 1, 
          xlab = "% total count per cell", 
          col = (scales::hue_pal())(50)[50:1], 
          horizontal = TRUE)
  dev.off()
  rm(C)
  selected_mito <- WhichCells(input_sce.filt, expression = percent_mito < 20)
  selected_ribo <- WhichCells(input_sce.filt, expression = percent_ribo > 3)
  selected_hb <- WhichCells(input_sce.filt, expression = percent_hb < 1 )
  length(selected_hb)
  length(selected_ribo)
  length(selected_mito)
  
  input_sce.filt <- subset(input_sce.filt, cells = selected_mito)
  #input_sce.filt <- subset(input_sce.filt, cells = selected_ribo)
  input_sce.filt <- subset(input_sce.filt, cells = selected_hb)
  dim(input_sce.filt)
  
  table(input_sce.filt$orig.ident) 
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
  w=length(unique(input_sce.filt$orig.ident))/3+5;w 
  ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered,width = w,height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
    NoLegend()
  w=length(unique(input_sce.filt$orig.ident))/2+5;w 
  ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered,width = w,height = 5) 
  return(input_sce.filt) 
  
}

#
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
save(sce.all.filt,file = "QC.rdata")
setwd('../')

rm(list = ls())
set.seed(2028)
load("1-QC/QC.rdata") 
table(sce.all.filt$orig.ident)

if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('../scRNA_scripts/harmony.R')
  sce.all.int = run_harmony(sce.all.filt)
  setwd('../')
}
#harmony
run_harmony <- function(input_sce){
  
  print(dim(input_sce))
  input_sce <- NormalizeData(input_sce, 
                             normalization.method = "LogNormalize",
                             scale.factor = 1e4) 
  input_sce <- FindVariableFeatures(input_sce)
  input_sce <- ScaleData(input_sce)
  input_sce <- RunPCA(input_sce, features = VariableFeatures(object = input_sce))
  seuratObj <- RunHarmony(input_sce, "orig.ident")
  names(seuratObj@reductions)
  seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                       reduction = "harmony")
  
  # p = DimPlot(seuratObj,reduction = "umap",label=T ) 
  # ggsave(filename='umap-by-orig.ident-after-harmony',plot = p)
  
  input_sce=seuratObj
  input_sce <- FindNeighbors(input_sce, reduction = "harmony",
                             dims = 1:15) 
  input_sce.all=input_sce
  for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
    input_sce.all=FindClusters(input_sce.all, #graph.name = "CCA_snn", 
                               resolution = res, algorithm = 1)
  }
  colnames(input_sce.all@meta.data)
  apply(input_sce.all@meta.data[,grep("RNA_snn",colnames(input_sce.all@meta.data))],2,table)
  
  p1_dim=plot_grid(ncol = 3, 
                   CellDimPlot(input_sce.all, 
                               group.by = "RNA_snn_res.0.01", 
                               reduction = "UMAP", 
                               label = T,
                               label.size = 4, 
                               label_repel = T, 
                               label_insitu = T,
                               label_point_size = 1, 
                               label_point_color =NA ,
                               label_segment_color = NA,
                               theme_use = ggplot2::theme_classic, 
                               theme_args = list(base_size = 13))+ 
                     ggtitle("louvain_0.01"), 
                   CellDimPlot(input_sce.all, 
                               group.by = "RNA_snn_res.0.05", 
                               reduction = "UMAP", 
                               label = T,
                               label.size = 4, 
                               label_repel = T, 
                               label_insitu = T,
                               label_point_size = 1, 
                               label_point_color =NA ,
                               label_segment_color = NA,
                               theme_use = ggplot2::theme_classic, 
                               theme_args = list(base_size = 13))+ 
                     ggtitle("louvain_0.05"), 
                   CellDimPlot(input_sce.all, 
                               group.by = "RNA_snn_res.0.1", 
                               reduction = "UMAP", 
                               label = T,
                               label.size = 4, 
                               label_repel = T, 
                               label_insitu = T,
                               label_point_size = 1, 
                               label_point_color =NA ,
                               label_segment_color = NA,
                               theme_use = ggplot2::theme_classic, 
                               theme_args = list(base_size = 13))+ 
                     ggtitle("louvain_0.1"))
  
  ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low.pdf",width = 21,height = 9)
  
  p1_dim=plot_grid(ncol = 3, 
                   CellDimPlot(input_sce.all, 
                               group.by = "RNA_snn_res.0.2", 
                               reduction = "UMAP", 
                               label = T,
                               label.size = 4, 
                               label_repel = T, 
                               label_insitu = T,
                               label_point_size = 1, 
                               label_point_color =NA ,
                               label_segment_color = NA,
                               theme_use = ggplot2::theme_classic, 
                               theme_args = list(base_size = 13))+ 
                     ggtitle("louvain_0.2"), 
                   CellDimPlot(input_sce.all, 
                               group.by = "RNA_snn_res.0.3", 
                               reduction = "UMAP", 
                               label = T,
                               label.size = 4, 
                               label_repel = T, 
                               label_insitu = T,
                               label_point_size = 1, 
                               label_point_color =NA ,
                               label_segment_color = NA,
                               theme_use = ggplot2::theme_classic, 
                               theme_args = list(base_size = 13))+ 
                     ggtitle("louvain_0.3"), 
                   CellDimPlot(input_sce.all, 
                               group.by = "RNA_snn_res.0.5", 
                               reduction = "UMAP", 
                               label = T,
                               label.size = 4, 
                               label_repel = T, 
                               label_insitu = T,
                               label_point_size = 1, 
                               label_point_color =NA ,
                               label_segment_color = NA,
                               theme_use = ggplot2::theme_classic, 
                               theme_args = list(base_size = 13))+ 
                     ggtitle("louvain_0.5"))
  
  ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_middle.pdf",width = 21,height = 9)
  
  p1_dim=plot_grid(ncol = 2, 
                   CellDimPlot(input_sce.all, 
                               group.by = "RNA_snn_res.0.8", 
                               reduction = "UMAP", 
                               label = T,
                               label.size = 4, 
                               label_repel = T, 
                               label_insitu = T,
                               label_point_size = 1, 
                               label_point_color =NA ,
                               label_segment_color = NA,
                               theme_use = ggplot2::theme_classic, 
                               theme_args = list(base_size = 13))+ 
                     ggtitle("louvain_0.8"), 
                   CellDimPlot(input_sce.all, 
                               group.by = "RNA_snn_res.1", 
                               reduction = "UMAP", 
                               label = T,
                               label.size = 4, 
                               label_repel = T, 
                               label_insitu = T,
                               label_point_size = 1, 
                               label_point_color =NA ,
                               label_segment_color = NA,
                               theme_use = ggplot2::theme_classic, 
                               theme_args = list(base_size = 13))+ 
                     ggtitle("louvain_1"))
  ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_high.pdf",width = 12,height = 6)
  
  p2_tree=clustree(input_sce.all@meta.data, 
                   prefix = "RNA_snn_res.")+
    theme_ggstatsplot()+
    theme(plot.title = element_text(size=20,hjust = 0.5),
          axis.title.x=element_text(size=20),
          axis.title.y=element_text(size=20),
          axis.text = element_text(color="black",size=16),
          legend.title = element_text(size=20),
          legend.text = element_text(size=18))
  
  ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf",width = 15,height = 9)
  table(input_sce.all@active.ident) 
  saveRDS(input_sce.all, "sce.all_int.rds")
  return(input_sce.all)
  
}
#
getwd()
dir.create('check-by-0.5')
setwd('check-by-0.5')
sel.clust = "RNA_snn_res.0.5"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
scRNA <- sce.all.int 
sce.all.int = readRDS('2-harmony/sce.all_int.rds')

pdf('RNA_snn_res.0.5_umap.pdf',width = 9,height = 6)
umap <- CellDimPlot(scRNA, 
                    group.by = sel.clust, 
                    reduction = "UMAP", 
                    label = F,
                    label.size = 4, 
                    label_repel = T, 
                    label_insitu = T,
                    label_point_size = 1, 
                    label_point_color =NA ,
                    label_segment_color = NA,
                    theme_use = ggplot2::theme_classic, 
                    theme_args = list(base_size = 16))
umap
dev.off()
pdf('RNA_snn_res.0.5_sample_umap.pdf',width = 9,height = 6)
sample_umap <- CellDimPlot(scRNA, 
                           group.by = "group", 
                           reduction = "UMAP", 
                           label = F,
                           label.size = 4, 
                           label_repel = T, 
                           label_insitu = T,
                           label_point_size = 1, 
                           label_point_color =NA ,
                           label_segment_color = NA,
                           theme_use = ggplot2::theme_classic, 
                           theme_args = list(base_size = 16))
sample_umap
dev.off()
pdf('RNA_snn_res.0.5_umap_sample_umap.pdf',width = 12,height = 4)
umap+sample_umap
dev.off()

cell_markers <- list()
cell_markers$B_plasma <- c("CD79A","JCHAIN")
cell_markers$T_cell <- c("CD3D","CD3E",'CD8A')
cell_markers$Epithelial_cell <- c("KRT18","EPCAM","CD24","KRT19")
cell_markers$Endothelial_cell <- c("PECAM1","CLDN5")
cell_markers$Cycle_cell <- c("MKI67","TOP2A")
cell_markers$Fibroblast <- c("DCN","OGN")
cell_markers$Monocyte <- c("CD14","C1QA")
cell_markers$SMC <- c("ACTA2","MYH11","TAGLN")

print(cell_markers)
p <- DotPlot(scRNA, 
             features = cell_markers, 
             assay='RNA',
             group.by = sel.clust,
             cols = "RdYlBu") + 
  ggtitle(paste0(sel.clust)) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p[["theme"]][["strip.text"]]$angle <- 90
p[["theme"]][["strip.text"]]$hjust <- 0
pdf(paste0(sel.clust,"_","Markers_dotplot8.pdf"), width=18, height = 6,bg="white")
p+umap
dev.off()


celltype=data.frame(ClusterID=0:10,
                    celltype= 0:10)



celltype[celltype$ClusterID %in% c(7,9),2]='B_plasma'
celltype[celltype$ClusterID %in% c( 0),2]='T_cell'
celltype[celltype$ClusterID %in% c(2,10),2]='Epithelial_cell'
celltype[celltype$ClusterID %in% c( 4 ),2]='Endothelial_cell'
celltype[celltype$ClusterID %in% c( 6 ),2]='Cycle_cell'
celltype[celltype$ClusterID %in% c( 1,5),2]='Fibroblast'
celltype[celltype$ClusterID %in% c( 3 ),2]='Monocyte'
celltype[celltype$ClusterID %in% c( 8),2]='SMC'

table(scRNA@meta.data$RNA_snn_res.0.5)
table(celltype$celltype)

scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),
                  'celltype'] <- celltype$celltype[i]}
table(scRNA@meta.data$celltype)


th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
pdf('umap_by_celltype-END.pdf',width = 9,height = 6)

celltype_umap <- CellDimPlot(scRNA, 
                             group.by = "celltype", 
                             reduction = "UMAP", 
                             label = F,#label = T
                             label.size = 4, 
                             label_repel = T, 
                             label_insitu = T,
                             label_point_size = 1, 
                             label_point_color =NA ,
                             label_segment_color = NA,
                             theme_use = ggplot2::theme_classic, 
                             theme_args = list(base_size = 16))

celltype_umap
dev.off()
saveRDS(scRNA, "sce_celltype.rds")

setwd('../')
phe=sce.all.int@meta.data
save(phe,file = 'phe.Rdata')

sce.all$sample <- sce.all$group
tb=table(sce.all$sample, sce.all$celltype)
head(tb)
balloonplot(tb)
#
bar_data <- as.data.frame(tb)

bar_per <- bar_data %>% 
  group_by(Var1) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per) 
col =c("#3176B7","#F78000","#3FA116","#CE2820","#9265C1",
       "#885649","#DD76C5","#BBBE00")#,"#41BED1"
col <- my36colors
colnames(bar_per)

p1 = ggplot(bar_per, aes(x = percent, y = Var1)) +
  geom_bar(aes(fill = Var2) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = " ", fill = NULL)+labs(x = 'Relative proportion(%)')+
  scale_fill_manual(values=col)+
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5))
p1

p2 = ggplot(bar_per, aes(x = Freq, y = Var1)) +
  geom_bar(aes(fill = Var2) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = " ", fill = NULL)+labs(x = 'Total cell number')+
  scale_fill_manual(values=col)+
  theme_classic()+
  theme(plot.title = element_text(size=12,hjust=0.5))

pdf(filename="prop1.pdf",width = 6,height = 6)
p1
dev.off()
pdf(filename="prop2.pdf",width = 6,height = 6)
p2
dev.off()
pdf(filename="prop3.pdf",width = 12,height = 6)
p1+p2
dev.off()
setwd('../')
#
sce.all=readRDS( "../3-Celltype/sce_celltype.rds")
pdf('Nebulosa_marker.pdf',width = 15,height = 9)
p2 <- plot_density(sce.all, 
                   coo,
                   size = 0.3,
                   ncol = 3)
p2
dev.off()
#Enrichment analysis----
suppressMessages(library(genefilter))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(TCGAbiolinks))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(ggplotify))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(AnnoProbe))
suppressMessages(library(tinyarray))
suppressMessages(library(GEOquery)) 
suppressMessages(library(ggplot2)) 
suppressMessages(library(ggstatsplot)) 
suppressMessages(library(reshape2))
suppressMessages(library(FactoMineR))
suppressMessages(library(factoextra))
suppressMessages(library(pheatmap))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(UpSetR))
suppressMessages(library(WGCNA))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(GSEABase))
suppressMessages(library(DOSE))
suppressMessages(library(ggstatsplot))
suppressMessages(library(ggsci))
suppressMessages(library(linkET))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(paletteer))
suppressMessages(library(Rtsne))

expr <- as.data.frame(cancer_vst)
HOPX_expr<-expr["KIF1A",]
HOPX_expr[,1:4]

cormatrixes<-function(x,y){
  tem<-linkET::correlate(t(x),t(y),engine="WGCNA")
  tem1<-as.data.frame(tem[[1]])
  tem1<-cbind(rownames(tem1),tem1)
  tem1_long<-reshape2::melt(tem1,value.name="correlation")
  tem2<-as.data.frame(tem[[2]])
  tem2<-cbind(rownames(tem2),tem2)
  tem2_long<-reshape2::melt(tem2,value.name="pvalue")
  result<-cbind(tem1_long,tem2_long$pvalue)
  names(result)<-c("v1","v2","correlation","pvalue")
  return(result)
}

identical(colnames(expr),colnames(HOPX_expr))
cor_res<-cormatrixes(HOPX_expr,expr)
cor_res<-na.omit(cor_res)

rownames(cor_res) <- cor_res$v2
colnames(cor_res) <- c("KIF1A","symbol","correlation","pvalue")

hopx_related_mrna<-cor_res%>%
  filter(correlation > 0.3,pvalue<0.05)%>%
  distinct(symbol)%>%
  pull(symbol)
length(hopx_related_mrna)
DEG <- cor_res
logFC_cutoff <- 0.3
k1 = (DEG$pvalue < 0.05)&(DEG$correlation <= -logFC_cutoff)
k2 = (DEG$pvalue < 0.05)&(DEG$correlation >= logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
df <- bitr(unique(DEG$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)#org.Hs.eg.db
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
gene_up= DEG[DEG$change=="UP",'ENTREZID']
gene_down=DEG[DEG$change == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
geneList=DEG$correlation
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = T)
###GO----
ego_up <- enrichGO(gene       = gene_up,
                   universe      = names(geneList),
                   OrgDb         = org.Hs.eg.db,#org.Hs.eg.db
                   ont           = "ALL",#"BP", "MF", and "CC"
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)
ego_up <- as.data.frame(ego_up)
go1 <- ego_up 
go1$p.adjust <- -log10(go1$p.adjust)
go1$GeneRatio <- ifelse(go1$ONTOLOGY=="BP",go1$Count/568*100,
                        ifelse(go1$ONTOLOGY=="CC", go1$Count/597*100, 
                               go1$Count/591*100))

go1<-go1[order(-go1$GeneRatio)[1:30],]

go1<-go1[order(-go1$GeneRatio),]
go1<-go1%>%
  group_by(ONTOLOGY)%>%
  slice(1:10)

p3 <- ggplot(go1,
             aes(Description,GeneRatio,fill=p.adjust))+  
  geom_bar(stat = "identity")+  
  coord_flip()+  
  scale_fill_gradient(low="blue",high="red")+  
  labs(x="",fill="-log10(adj.P)")+  
  ggtitle("GO enrichment")+  
  theme_bw()+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=15),        
        legend.title = element_text(size = 18),  
        legend.text = element_text(size = 14), 
        plot.title = element_text(size = 20,face = "bold",hjust = 0.5),
        strip.text = element_text(size=18))
p3+facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y',)

###KEGG----
kk_up <- enrichKEGG(gene         = gene_up,
                    organism     = 'hsa',
                    pvalueCutoff = 0.99)#0.05
kk_up <- as.data.frame(kk_up)


go1<-kk_up
go1$p.adjust <- -log10(go1$p.adjust)
go1$GeneRatio <- go1$Count/250*100
go1<-go1[order(-(go1$GeneRatio)),]
go1<-go1[order(-go1$GeneRatio)[1:20],]

p3 <- ggplot(go1,
             aes(Description,Count,fill=p.adjust))+  
  geom_bar(stat = "identity")+  
  coord_flip()+  
  scale_fill_gradient(low="blue",high="red")+  
  labs(x="",fill="-log10(adj.P)")+  
  ggtitle("KEGG enrichment")+  
  theme_bw()+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=15),        
        legend.title = element_text(size = 18),  
        legend.text = element_text(size = 14), 
        plot.title = element_text(size = 20,face = "bold",hjust = 0.5))
p3
###HALLMARK-GSEA----

geneset<-read.gmt("h.all.v2025.1.Hs.entrez.gmt")
geneset$term=str_remove(geneset$term,"HALLMARK_")
geneList=DEG$log2FoldChange
geneList=DEG$correlation
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = T)
egmt<-GSEA(geneList,TERM2GENE=geneset,verbose=F)
egmt2<-setReadable(egmt,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
kk_gse=egmt2 <- -
kk_gse1 <- as.data.frame(kk_gse)
down<-kk_gse[kk_gse@result$p.adjust<0.05 & kk_gse@result$NES< 0,];down$group=-1
up<-kk_gse[kk_gse@result$p.adjust<0.05 & kk_gse@result$NES>0,];up$group=1
up <- down
pro='gsea'
lapply(1:nrow(up),function(i){
  plot <- gseaplot2(egmt2,up$ID[i],
                    title=up$Description[i],pvalue_table=T)
  
  ggsave(paste0(pro,'_up_',
                gsub('/','-',up$Description[i]),
                '.pdf'),plot = plot,width=8,height=6) 
  
  
})

group_list <- ifelse(exp$KIF1A > median(exp$KIF1A), "High","Low")

lapply(1:nrow(up),function(i){
  cg=strsplit(up$core_enrichment[i],'/')[[1]]
  dat <- exp 
  col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
  top_annotation = HeatmapAnnotation(
    cluster = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
                         labels = c("Low","High"),
                         labels_gp = gpar(col = "white", 
                                          fontsize = 20)))
  plot <- Heatmap(t(scale(t(dat[cg,]))),name = " ",
                  col = col_fun,
                  cluster_columns = F,
                  top_annotation = top_annotation,
                  column_split = group_list,
                  show_heatmap_legend = T,
                  border = F,
                  show_column_names = F,
                  show_row_names = T,
                  column_title = up$Description[i],
                  row_names_gp = gpar(fontsize = 16),
                  column_title_gp = gpar(fontsize = 24)
  )
  plot <- as.ggplot(plot)
  ggsave(paste0(pro,'_up_pheatmap_',
                gsub('/','-',up$Description[i]),
                '.pdf'), plot,width=16,height=16)
})
