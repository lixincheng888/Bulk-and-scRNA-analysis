####Step 1、gene annotatio####
#gene annotation
# exp=read.csv("data.csv",sep="\t",header=T,check.names=F)#row data
# sample=read.csv("sample.csv",header=T)
#sample<-sample$sample
# exp<-na.omit(exp)
# 
# GPL6244_anno <-data.table::fread("GSE117261_family.soft",skip ="ID")
# ann <- GPL6244_anno %>% 
#   dplyr::select(ID,gene_assignment) %>% 
#   dplyr::filter(gene_assignment != "---") %>% 
#   separate(gene_assignment,c("drop","symbol"),sep="//") %>% 
#   dplyr::select(-drop)
# 
# write.csv(ann,file="GPL6244.csv")
# 
# exp$ID = as.character(exp$ID_REF)
# 
# 
# exp <- exp %>% 
#   inner_join(ID,by="ID") %>% 
#   dplyr::select(-ID) %>% 
#   dplyr::select(Gene.Symbol, everything()) %>% 
#   mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>%
#   arrange(desc(rowMean)) %>% 
#   distinct(Gene.Symbol,.keep_all = T) %>% 
#   dplyr::select(-rowMean)  
# 
# rownames(exp)=exp[,1]#
# exp=exp[,sample]
# write.csv(exp,file="exp.csv")

####Step 2、WGCNA analysis####
#1、 R package
# rm(list = ls())
# library(WGCNA)
# library(tinyarray)
#2、read data
#setwd("pathway")
#expr<-read.csv("exp.csv",header=T,check.names=F)#gene expression data

#clindata<-read.csv("clinaldata.csv",row.names = 1)#clinical data
#expr<-expr[,rownames(clindata)]#IPAH samples
#datExpr0 = t(expr[order(apply(expr, 1, var), decreasing = T)[1:round(0.75*nrow(expr))],])

#Expression matrix after cleaning

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(datExpr0), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 1.5)
abline(h = 70, col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = 70, minSize = 10)
table(clust)   
keepSamples = (clust==1)  
datExpr = datExpr0[keepSamples, ]  
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datExpr = as.data.frame(datExpr)
##2、clincal data 
library(stringr)
clindata<-clindata[rownames(datExpr),]#Eliminate outliers
colnames(clindata)
traitData = data.frame(Gender=as.numeric(as.factor(clindata$Gender)),
                       Age=as.numeric(clindata$Age),
                       Inflammatory_score=as.numeric(clindata$Inflammatory_score),
                       Ad_thickness=as.numeric(clindata$Adventitia_Fractional_thickness),
                       Me_thickness=as.numeric(clindata$Media_Fractional_thickness),
                       In_thickness=as.numeric(clindata$Intima_Fractional_thickness),
                       In_Me_thickness=as.numeric(clindata$I.M),
                       To_thickness=as.numeric(clindata$Total_thickness),
                       mPAP=as.numeric(clindata$mPAP),
                       PVR=as.numeric(clindata$PVR))
rownames(traitData)<-rownames(clindata)


femaleSamples = rownames(datExpr)

sampleTree2 = hclust(dist(datExpr), method = "average")

traitColors = numbers2colors(traitData, signed = FALSE)

#soft value validation
powers = c(1:10, seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

cex1 = 0.9
#png(file = "Soft threshold.png", width = 2000, height = 1500,res = 300)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))+theme(axis.title.x = element_text(size = 30))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))+theme(axis.title.x = element_text(size = 30))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

power = sft$powerEstimate

#One-step construction

net = blockwiseModules(datExpr, 
                       power = power,
                       TOMType = "unsigned", 
                       minModuleSize = 30, 
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,
                       deepSplit = 2 ,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "testTOM",
                       verbose = 3)

save(net,datTraits,file="one_step_net.Rdata")

load("one_step_net.Rdata")

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

table(mergedColors)

moduleLabels = net$colors

moduleColors = labels2colors(net$colors)

MEs = net$MEs;
geneTree = net$dendrograms[[1]];

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Merge Module
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 
sizeGrWindow(8,16)
# P value
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(1, 8.5, 1, 1));

library(RColorBrewer)
mycol <- colorRampPalette(c("#5088D0",'white','#a13037'))(20)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = mycol,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Model identification

Inflammatory_score= as.data.frame(datTraits$Inflammatory_score);
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, Inflammatory_score, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Inflammatory_score), sep="");
names(GSPvalue) = paste("p.GS.", names(Inflammatory_score), sep="");

#
module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module,"module"),
                   ylab = "Gene significance for Inflammatory_score",
                   main = paste(""),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# module gene

moduleColors = mergedColors
colorOrder = c("grey", unique(c(standardColors(50), unique(mergedColors))))
moduleLabels = match(moduleColors, colorOrder)-1

# Module
module_dataframe <- data.frame(gene_id=colnames(datExpr), 
                               module_name=paste0('module_', moduleLabels), 
                               module_color=moduleColors)

write.csv(module_dataframe, file = "module_dataframe_gene.csv")


####Step 3、DEGs analysis ####
library(limma)
library("impute")
library(dplyr)
library(tidyr)
library(sva)
setwd("working pathway")
exp<-read.csv(file="exp.csv",header = T)

dimnames=list(rownames(exp),colnames(exp))
exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=impute.knn(exp)
rt=mat$data

#normalize
pdf(file="rawBox.pdf")
boxplot(rt,col = "blue",xaxt = "n",outline = F)
dev.off() #
rt2=normalizeBetweenArrays(as.matrix(rt))
pdf(file="normalBox.pdf")
boxplot(rt2,col = "red",xaxt = "n",outline = F)
dev.off()

#data filtering

cutoff = .75
exp <- data.frame(rt2[which(apply(rt2, 1, function(x){length(which
                                                           (x!= 0))/length(x)}) >= cutoff),])


# group design

colnames(exp)

Type=c(rep("IPAH",24),rep("Control",21))
design <- model.matrix(~0+factor(Type))
names<-colnames(exp)
rownames(design)<-names
colnames(design) <- c("Control","IPAH")
fit <- lmFit(exp,design)
cont.matrix<-makeContrasts(IPAH-Control,levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)

write.csv(allDiff,file="limmaDEG.csv")

DEGs<-na.omit(allDiff)

diff<- DEGs[with(DEGs, (abs(DEGs$logFC)> 0.58 & adj.P.Val< 0.05 )), ]

write.csv(diff,file="DEGs.csv")
#GO and KEGG analysis

#The GO and KEGG analysis of the functions of DEGs related to inflammation were 
#performed using the "clusterProfiler" package (version 4.2.0).
# The enrichment results were visualized using Xiantao Academic Online website (https://www.xiantaozi.com/).

####Step 3、scRNA analysis####
rm(list=ls())
# R packages

library(Seurat)
library(metap)
library(ggplot2)
library(cowplot)
library(rhdf5)
library(glmGamPoi)
library(tidyverse)
library(patchwork)
library(limma)
library(stringr)
library(openxlsx)
library(celldex)
library(dplyr)
library(ggsci)
library(sctransform)
#setwd("")

#Normal samples
Normal<-list()

samples <- basename(list.files("Normal/",recursive = F))#sample pathway

for (sample in samples) {
  # 
  scrna_data <-  Read10X_h5(filename =str_c("Normal/", sample))
  # 
  seob1 <- CreateSeuratObject(
    counts = scrna_data,
    project = sample,
    min.cells = 3,  
    min.features = 200)
  seob1[['sample']] <- sample
  seob1[['status']]<-"Normal"
 Normal[[sample]] = seob1}

Normal <- merge(x = Normal[[1]], # 第一个样本
                y = Normal[-1], # 其他的样本
                add.cell.ids = names(IPAH)) # cell id 添加前缀


#IPAH samples

IPAH<-list()

samples <- basename(list.files("IPAH/",recursive = F))#sample pathway

for (sample in samples) {
  # 
  scrna_data <-  Read10X_h5(filename =str_c("IPAH/", sample))
  # 
  seob1 <- CreateSeuratObject(
    counts = scrna_data,
    project = sample,
    min.cells = 3,  
    min.features = 200)
  seob1[['sample']] <- sample
  seob1[['status']]<-"IPAH"
  Normal[[sample]] = seob1}

IPAH <- merge(x = IPAH[[1]], # 
               y = IPAH[-1], # 
add.cell.ids = names(IPAH)) # 


seob<-merge(x=Normal, y=IPAH)

seob[["percent.mt"]] <- PercentageFeatureSet(
  seob, pattern = "^MT-")

#
VlnPlot(seob, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by  = "sample")

#
# p1 <- FeatureScatter(seob, 
#                      feature1 = "nCount_RNA", 
#                      feature2 = "nFeature_RNA",
#                      group.by = "sample")
# 
# p2 <- FeatureScatter(seob, 
#                      feature1 = "nCount_RNA", 
#                      feature2 = "percent.mt",
#                      group.by = "sample")
# p1+p2+plot_layout(guides = "collect")
# DefaultAssay(sce2)<-'RNA'
#
seob <- subset(seob,
               subset = nFeature_RNA > 200 &
                 nFeature_RNA < 3000 &
                 percent.mt < 10)
# SCTransform

seob_list <- SplitObject(seob, split.by = "sample")

# 1. SCTransform

for(i in 1:length(seob_list)){
  seob_list[[i]] <- SCTransform(
    seob_list[[i]],
    method = "glmGamPoi",
    #vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
    verbose = FALSE)
}

# 2. Integration
## 
features <- SelectIntegrationFeatures(object.list = seob_list,
                                      nfeatures =3000) 
## 
seob_list <- PrepSCTIntegration(object.list = seob_list, 
                                anchor.features = features)
seob_list <- lapply(X = seob_list, FUN = RunPCA, features = features)
## anchors
anchors <- FindIntegrationAnchors(object.list = seob_list, 
                                  reference = c(3,7), 
                                  normalization.method = "SCT",
                                  dims = 1:30, 
                                  reduction = "rpca",
                                  k.anchor = 20,
                                  anchor.features = features)

## integration
rm(features)
seob <- IntegrateData(anchorset = anchors, 
                      dims = 1:30,
                      normalization.method = "SCT")

DefaultAssay(seob) <- "integrated"

#Dimensionality reduction analysis

#PCA
seob <- RunPCA(seob)

ElbowPlot(seob, ndims = 50)

#clustering analysis

seob <- FindNeighbors(seob,
                      dims = 1:30)


seob <- FindClusters(seob,
                     resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                     random.seed = 1)

cluster_umap <- plot_grid(ncol = 2, rel_heights = 1,
                          DimPlot(seob, reduction = "umap", group.by = "integrated_snn_res.0.1", label = T) & NoAxes(), 
                          DimPlot(seob, reduction = "umap", group.by = "integrated_snn_res.0.3", label = T)& NoAxes()
                          
)


# p1 <- DimPlot(seob, 
#               reduction = "pca",
#               group.by  = "seurat_clusters",
#               label = T)
# 
# p2 <- DimPlot(seob, 
#               reduction = "umap", 
#               group.by  = "seurat_clusters",
#               label = T)
# 
# p1  + p2+plot_layout(guides = "collect")& theme(legend.position = "top")

#DEGs of clusters

logFCfilter=0.25
adjPvalFilter=0.05
seob.markers <- FindAllMarkers(object = seob,
                               only.pos =T,#高表达
                               min.pct = 0.20,
                               logfc.threshold = logFCfilter)
sig.markers=seob.markers[(abs(seob.markers$avg_log2FC))>logFCfilter &(seob.markers$p_val_adj<adjPvalFilter),]

                              
top10 <- seob.markers %>% group_by(cluster) %>% top_n(n = 10, wt= avg_log2FC)

#write.csv(top10,file="top10DEGs.csv")

#cell marker
markers <- c( "CD68",
              "LUM","ACTA2","COL1A1",
              "VWF","PECAM1",
              "EPCAM",
              "CD3D",
              "NKG7",
              "TPSAB1", #Mast cell
              #"S100A8",
              "FCN1",        #mono
              "CD79A",      #B cells
              "MZB1",
              "MS4A6A", 
              "LYVE1")

library(scCustomize)
#top 5 gene
all_markers <- FindAllMarkers(object = seob)

top5_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 5, named_vector = FALSE,
                                    make_unique = TRUE)

Clustered_DotPlot(seurat_object = pbmc, features = top5_markers)

# Cell annotation
celltype=data.frame(ClusterID=0:24,
                    celltype='NA')
celltype[celltype$ClusterID %in% c(0,1,2,4,9,19),2]='Macrophages'  
celltype[celltype$ClusterID %in% c(5,13,24),2]='Low quality Macs'
celltype[celltype$ClusterID %in% c(22),2]='Low quality T cells'
celltype[celltype$ClusterID %in% c(11),2]='Fibroblasts' 
celltype[celltype$ClusterID %in% c(8),2]='SMCs' 
celltype[celltype$ClusterID %in% c(3,12),2]='VECs'
celltype[celltype$ClusterID %in% c(20),2]='LECs'
celltype[celltype$ClusterID %in% c(10,14,16,21,23),2]='Epithelial cells'
celltype[celltype$ClusterID %in% c(7),2]='NK cells'
celltype[celltype$ClusterID %in% c(0),2]='T cells' 
celltype[celltype$ClusterID %in% c(11),2]='B/Plasma cells'  
celltype[celltype$ClusterID %in% c(18),2]='B cells'
celltype[celltype$ClusterID %in% c(11),2]='Mast cells' 
celltype[celltype$ClusterID %in% c(21),2]='Proliferating cells' 
celltype[celltype$ClusterID %in% c(15),2]='Dendritic cells' 
celltype[celltype$ClusterID %in% c(6),2]='Monocytes'  

for(i in 1:nrow(celltype)){
  seob@meta.data[which(seob@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

table(seob@meta.data$celltype)



my36colors <-c( '#BD956A', '#585658',"#EDAA05","#2E8B57","#3C719E","#0086B3","#E2D200","#E7921E","#B40C1F","#B91D2E","#42ABC7",'#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E59CC4', '#AB3282', '#23452F',
               '#9FA3A8', '#E0D4CA', '#5F3D69',  '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')

sce2<-subset(seob, cellltype==c("Low quality T cells","Low quality Macs"),invert=T)

DefaultAssay(sce2) <- "RNA"
Idents(sce2) <- "celltype"

DimPlot(sce2, 
        label = T,
        pt.size = 0.5,
        group.by = "celltype",
        reduction = "umap", 
        cols = my36colors)+
        tidydr::theme_dr()+
       theme(panel.grid = element_blank())



VlnPlot(sce2, 
        features = markers,
        stacked=T,
        pt.size=0,
        cols = my36colors,
        direction = "vertical",
        x.lab = '', y.lab = '')+
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
#function
plot_density <- function(obj,
                            marker,
                            dim=c("TSNE","UMAP"),
                            size,
                            ncol=NULL
){
  require(ggplot2)
  require(ggrastr)
  require(Seurat)
  
  cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494'))
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c'))
  mypalette <- c(rev(cold(11)), warm(10))
  
  if(dim=="TSNE"){
    
    xtitle = "tSNE1"
    ytitle = "tSNE2"
    
  }
  
  if(dim=="UMAP"){
    
    xtitle = "UMAP1"
    ytitle = "UMAP2"
  }
  
  
  if(length(marker)==1){
    
    plot <- FeaturePlot(obj, features = marker)
    data <- plot$data
    
    
    if(dim=="TSNE"){
      
      colnames(data)<- c("x","y","ident","gene")
      
    }
    
    if(dim=="UMAP"){
      
      colnames(data)<- c("x","y","ident","gene")
    }
    
    
    #ggplot
    p <- ggplot(data, aes(x, y)) +
      geom_point_rast(shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = size) +
      geom_density_2d(data=data[data$gene>0,], 
                      aes(x=x, y=y), 
                      bins = 5, colour="black") +
      scale_fill_gradientn(colours = mypalette)+
      scale_colour_gradientn(colours = mypalette)+
      theme_bw()+ggtitle(marker)+
      labs(x=xtitle, y=ytitle)+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    return(p)
    
  }else{
    
    gene_list <- list()
    
    
    
    for (i in 1:length(marker)) {
      plot <- FeaturePlot(obj, features = marker[i])
      data <- plot$data
      
      
      if(dim=="TSNE"){
        
        colnames(data)<- c("x","y","ident","gene")
      }
      
      if(dim=="UMAP"){
        
        colnames(data)<- c("x","y","ident","gene")
      }
      
      gene_list[[i]] <- data
      names(gene_list) <- marker[i]
    }
    
    plot_list <- list()
    
    
    for (i in 1:length(marker)) {
      
      p <- ggplot(gene_list[[i]], aes(x, y)) +
        geom_point_rast(shape = 21, stroke=0.25,
                        aes(colour=gene, 
                            fill=gene), size = size) +
        geom_density_2d(data=gene_list[[i]][gene_list[[i]]$gene>0,], 
                        aes(x=x, y=y), 
                        bins = 5, colour="black") +
        scale_fill_gradientn(colours = mypalette)+
        scale_colour_gradientn(colours = mypalette)+
        theme_bw()+ggtitle(marker[i])+
        labs(x=xtitle, y=ytitle)+
        theme(
          plot.title = element_text(size=12, face="bold.italic", hjust = 0),
          axis.text=element_text(size=8, colour = "black"),
          axis.title=element_text(size=12),
          legend.text = element_text(size =10),
          legend.title=element_blank(),
          aspect.ratio=1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      plot_list[[i]] <- p
    }
    
    
    Seurat::CombinePlots(plot_list, ncol = ncol)
    
    
  }
  
  
}


#UMAP
plot_density(obj=sce2, 
                marker=c("CCL5", "GZMA"), 
                dim = "UMAP", 
                size =1, 
                ncol =2)

plot_density(obj=sce2, 
             marker=c("CXCL9"), 
             dim = "UMAP", 
             size =1, 
             ncol =1)

#Gene difference analysis
library(ggsci)
library(ggsignif)
library(tidyverse)

seu_object<-subset(sce2,celltype==c("NK cells","T cells"))

seu_object <- ScaleData(seu_object)

seu_object<-subset(seu_object,CCL5>0&GZMA>0)



seu_object=seu_object[,seu_object$celltype=="NK cells"] # T cells

data.frame(Group=seu_object$group,
           gene=sseu_object@assays[["RNA"]]@data["CCL5",])%>%  #GZMA
           ggplot(aes(x = Group, y= gene,fill =Group))+
           geom_violin(width = 0.8,trim=T,color="black",position = position_nudge(x = 0, y = 0),lwd=0.8,alpha=0.2) + 
           geom_point(shape=21,size=1,position = position_jitterdodge(),alpha=0.9,col="black") +
          scale_fill_manual(values = c("#F2300F", "#57C3F3")) +
           ggsignif::geom_signif(comparisons = list(c("Normal","IPAH")),map_signif_level = T)+
           theme_classic(base_size=14)+
           xlab("NK cells")+   #T cells
           ylab("CCL5")+       # GZMA
           theme(axis.text.y = element_text(size=12, colour = "black"))+
           theme(axis.text.x = element_text(size=12, colour = "black"))

    
####Step 5、iTALK analysis  ####

rm(list = ls())

setwd("working path")
library(circlize)
library(iTALK)
source("working path/sc_function.R")
load("working path/sce2.rdata")
sce3<-subset(sce2,celltype==c("T cells","NK cells","Fibroblasts","SMCs","VECs"))

#expression data

exp <- as.data.frame(t(as.matrix(sce3@assays$RNA@counts)))

exp$cell_type <- sce3@meta.data$celltype

exp$compare_group <- sce3@meta.data$group


highly_exprs_genes<-rawParse(exp,top_genes=21095,stats='mean')

comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_col<-structure(c("#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3"),names=unique(exp$cell_type))

#Cell interaction analysis

par(mfrow=c(1,2))
res<-NULL

for(comm_type in comm_list){
  res_cat<-FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
  
  NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
 
  LRPlot(res_cat[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],link.arr.width=res_cat$cell_to_mean_exprs[1:20])
  title(comm_type)
 
  res<-rbind(res,res_cat)
}

#sort
res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),]


#write.csv(res,file="net_inter.csv")

# Identify the major receptors for CCL5 on VECs, SMCs and fibroblasts

A=subset(res,ligand==c("CCL5"))

gene_col<-structure(c(rep('#CC3333',length(A[1:40,]$ligand)),
                      rep("#006699",length(A[1:40,]$receptor))),
                    names=c(A[1:40,]$ligand,
                            A[1:40,]$receptor))

cell_col <- structure(c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                        "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                        "#33A02C", "#B2DF8A", "#55A1B1"),
                      names=unique(A$cell_from))



pdf(file="circlize.pdf", width=8, height=6, pointsize=8)

iTALK::LRPlot(A[c(),],#CCL5: From (T cells and NK cells) To (SMCs, VECs and fibroblasts) 
              datatype='mean count',
              link.arr.lwd=A$cell_from_mean_exprs[c()],#CCL5: From (T cells and NK cells) To (SMCs, VECs and fibroblasts)
              link.arr.width=0.1,
              link.arr.col = 'grey20',
              print.cell = T,
              track.height_1=uh(1, "mm"),
              track.height_2 = uh(15, "mm"),
              text.vjust = "0.5cm",
              cell_col = cell_col[1:5])+
  
         legend("right",
         pch=20,
         legend=unique(A$cell_from),
         bty="n",
         col =cell_col[1:5],
         cex=1,pt.cex=3,
         border="black",
         title = 'Celltype') 

###https://github.com/Coolgenome/iTALK/blob/master/example/example_code.r



