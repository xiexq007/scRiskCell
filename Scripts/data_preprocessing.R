
# ------------------------------------------------------------------------------

# scRNA-seq meta-analysis of following published studies
# 1. Human Pancreas Analysis Program (HPAP)
# 2. Motakis (GSE221156)
# 3. Zhang (GSE195986) 
# 4. Segerstolpe et al.2016 (E-MATB-5061) 
# 5. Xin et al.2018 (GSE114297)

# ------------------------------------------------------------------------------

library(Seurat) 
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(clustree)
library(harmony)
library(SeuratWrappers)
library(patchwork)
library(cols4all)
library(cowplot)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# ------------------------------------------------------------------------------
# 1. create Seurat objects
# ------------------------------------------------------------------------------

# Human Pancreas Analysis Program (HPAP)
dir = list.dirs("./rawdata/hpap")[-1] 
names(dir) = list.files("./rawdata/hpap",recursive = F)
hpap = CreateSeuratObject(counts = Read10X(dir), min.cells = 3, min.features = 200)
dim(hpap) # 34097 genes 195134 cells
table(hpap@meta.data$orig.ident)  
hpap@meta.data$donor = str_split(hpap@meta.data$orig.ident, "-", simplify = TRUE)[,1]
hpap@meta.data$disease_status = str_split(hpap@meta.data$orig.ident, "-", simplify = TRUE)[,2]
hpap@meta.data$gender = str_split(hpap@meta.data$orig.ident, "-", simplify = TRUE)[,3]
hpap@meta.data$data_source = "hpap"

# saveRDS(hpap, "./rawdata/Seurat_objects/hpap.rds")

# Motakis (GSE221156)
dir = list.dirs("./rawdata/GSE221156/data")[-1] 
names(dir) = list.files("./rawdata/GSE221156/data",recursive = F)
motakis = CreateSeuratObject(counts = Read10X(dir), min.cells = 3, min.features = 200)
dim(motakis) # 31683 genes 322314 cells
table(motakis@meta.data$orig.ident)  
motakis@meta.data$donor = str_split(motakis@meta.data$orig.ident, "-", simplify = TRUE)[,1]
motakis@meta.data$disease_status = str_split(motakis@meta.data$orig.ident, "-", simplify = TRUE)[,2]
motakis@meta.data$gender = str_split(motakis@meta.data$orig.ident, "-", simplify = TRUE)[,3]
motakis@meta.data$data_source = "motakis"

# saveRDS(motakis, "./rawdata/Seurat_objects/motakis.rds")

# Zhang (GSE195986) 
samples = list.files("./rawdata/GSE195986")

scRNAlist = list()
for(i in 1:length(samples)){
  print(i)
  counts = read.table(samples[i], sep = "\t",header = T,row.names = 1)
  name = rownames(counts)
  rownames(counts) = gsub('_', '-', name)
  donor = sub("^[^-]+-(\\w+).*", "\\1", samples[i])
  scRNAlist[[i]] = CreateSeuratObject(counts, min.cells = 3, min.features = 200, project = donor)
}

zhang = merge(scRNAlist[[1]], add.cell.ids = sub("^[^-]+-(\\w+).*", "\\1", samples),
              y = c(scRNAlist[[2]], scRNAlist[[3]], scRNAlist[[4]], 
                    scRNAlist[[5]], scRNAlist[[6]], scRNAlist[[7]], 
                    scRNAlist[[8]], scRNAlist[[9]], scRNAlist[[10]],scRNAlist[[11]]))

zhang@meta.data$disease_status = ifelse(zhang@meta.data$orig.ident %in% c("HT2","HT3","HT5","HT6","HT7"), "ND",
                               ifelse(zhang@meta.data$orig.ident %in% c("HT1", "HT4"), "preT2D",
                                      ifelse(zhang@meta.data$orig.ident %in% c("T2D1","T2D2","T2D3","T2D4"), "T2D", NA)))
zhang@meta.data$donor = zhang@meta.data$orig.ident
zhang@meta.data$gender = ifelse(zhang@meta.data$donor %in% c("HT3","HT7"),"Female","Male")
zhang@meta.data$orig.ident = paste(zhang@meta.data$donor,zhang@meta.data$disease_status,zhang@meta.data$gender,sep = "-")

dim(zhang) # 17570 genes 55708 cells 
table(zhang@meta.data$orig.ident)  
zhang@meta.data$data_source = "zhang"
zhang = JoinLayers(zhang)
# saveRDS(zhang, "./rawdata/Seurat_objects/zhang.rds")

# Segerstolpe et al.2016 (E-MATB-5061)
data = read.csv("./rawdata/E-MATB-5061/E-MATB-5061-count.csv",header = T,row.names = 1)
ann = read.csv("./rawdata/E-MATB-5061/E-MATB-5061-match.csv",header = T)
rownames(data) = gsub('_', '-', rownames(data))

segerstolpe = CreateSeuratObject(data, min.cells = 3, min.features = 200, project = "segerstolpe")
ann = ann[match(rownames(segerstolpe@meta.data),ann$ERR),]
segerstolpe@meta.data$donor = ann$Donor
segerstolpe@meta.data$disease_status = ann$Group
segerstolpe@meta.data$gender = ann$sex
segerstolpe@meta.data$data_source = "segerstolpe"
segerstolpe@meta.data$orig.ident = paste(segerstolpe@meta.data$donor,segerstolpe@meta.data$disease_status,segerstolpe@meta.data$gender,sep = "-")

dim(segerstolpe) # 39173 genes 3318 cells 
table(segerstolpe@meta.data$orig.ident)  

# saveRDS(segerstolpe, "./rawdata/Seurat_objects/segerstolpe.rds")

# Xin et al.2018 (GSE114297)
dir = list.dirs("./rawdata/GSE114297/data")[-1] 
names(dir) = list.files("./rawdata/GSE114297/data",recursive = F)
xin = CreateSeuratObject(counts = Read10X(dir), min.cells = 3, min.features = 200)
dim(xin) # 20497 genes 20778 cells
table(xin@meta.data$orig.ident)  
xin@meta.data$donor = str_split(xin@meta.data$orig.ident, "-", simplify = TRUE)[,1]
xin@meta.data$disease_status = str_split(xin@meta.data$orig.ident, "-", simplify = TRUE)[,2]
xin@meta.data$gender = str_split(xin@meta.data$orig.ident, "-", simplify = TRUE)[,3]
xin@meta.data$data_source = "xin"

# saveRDS(xin, "./rawdata/Seurat_objects/xin.rds")

# ------------------------------------------------------------------------------
# 2. load data
# ------------------------------------------------------------------------------

# Human Pancreas Analysis Program (HPAP)
hpap = readRDS("./rawdata/Seurat_objects/hpap.rds")

# Motakis (GSE221156)
motakis = readRDS("./rawdata/Seurat_objects/motakis.rds")

# Zhang (GSE195986) 
zhang = readRDS("./rawdata/Seurat_objects/zhang.rds")

# Segerstolpe et al.2016 (E-MATB-5061)
segerstolpe = readRDS("./rawdata/Seurat_objects/segerstolpe.rds")

# Xin et al.2018 (GSE114297)
xin = readRDS("./rawdata/Seurat_objects/xin.rds")

dim(hpap) # 34097 genes 195134 cells
dim(motakis) # 31683 genes 322314 cells
dim(zhang) # 17570 genes 55708 cells
dim(segerstolpe) # 39173 genes 3318 cells
dim(xin) # 20497 genes 20778 cells

# ------------------------------------------------------------------------------
# 3. QC
# ------------------------------------------------------------------------------

objects = list(hpap,motakis,zhang,segerstolpe,xin)

for (i in 1:5) {
  
  print(i)
  objects[[i]][["percent.mt"]] = PercentageFeatureSet(objects[[i]], pattern = "^MT-")
  objects[[i]][["percent.ribo"]] = PercentageFeatureSet(objects[[i]], pattern = "^RPS|^RPL")
  VlnPlot(objects[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, raster = F)
  
}

#  Filter cells with genes < 200 and > 8000, mt percent > 25%
filtered_hpap = subset(objects[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25) 
filtered_motakis = subset(objects[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25) 
filtered_zhang = subset(objects[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25) 
filtered_segerstolpe = subset(objects[[4]], subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25) 
filtered_xin = subset(objects[[5]], subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25) 

dim(filtered_hpap) # 34097 genes 181895 cells
dim(filtered_motakis) # 31683 genes 271243 cells
dim(filtered_zhang) # 17570 genes 55468 cells
dim(filtered_segerstolpe) # 39173 genes 2143 cells
dim(filtered_xin) # 20497 genes 20623 cells

# ------------------------------------------------------------------------------
# 4. Normalize
# ------------------------------------------------------------------------------

filtered_hpap = NormalizeData(filtered_hpap, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_motakis = NormalizeData(filtered_motakis, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_zhang = NormalizeData(filtered_zhang, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_segerstolpe = NormalizeData(filtered_segerstolpe, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_xin = NormalizeData(filtered_xin, normalization.method = "LogNormalize", scale.factor = 10000)

# ------------------------------------------------------------------------------
# 5. Remove double positive hormonal cells in each dataset
# ------------------------------------------------------------------------------

# hpap
hormones = c("INS","GCG","PPY","SST","GHRL")
n = length(hormones)
hpap_ridges = list()
for (i in 1:5) {
  hpap_ridges[[i]] = RidgePlot(filtered_hpap,features=hormones[i],group.by = "data_source")
}

hpap_scatters = list()
for (i in 1:(n-1)) {
  for (j in (i+1):n){
    hpap_scatter = FeatureScatter(filtered_hpap, feature1 = hormones[i],feature2 = hormones[j], group.by = "data_source")+
      labs(title = NULL, caption = NULL)+
      theme(legend.position = "none",axis.text = element_blank(),axis.ticks = element_blank())
    hpap_scatters[[length(hpap_scatters)+1]] = hpap_scatter
  }
}

filtered_hpap = subset(filtered_hpap, subset = ((INS>6 & GCG>4.5)|(INS>6 & PPY>4.5)|(INS>6 & SST>4.5)|(INS>6 & GHRL>4)|(GCG>4.5 & PPY>4.5)|(GCG>4.5 & SST>4.5)|(GCG>4.5 & GHRL>4)|(PPY>4.5 & SST>4.5)|(PPY>4.5 & GHRL>4)|(SST>4.5 & GHRL>4)), invert = T)
dim(filtered_hpap) # 34097 genes 171085 cells

# motakis
motakis_ridges = list()
for (i in 1:5) {
  motakis_ridges[[i]] = RidgePlot(filtered_motakis,features=hormones[i],group.by = "data_source")
}

motakis_scatters = list()
for (i in 1:(n-1)) {
  for (j in (i+1):n){
    motakis_scatter = FeatureScatter(filtered_motakis, feature1 = hormones[i],feature2 = hormones[j], group.by = "data_source")+
      labs(title = NULL, caption = NULL)+
      theme(legend.position = "none",axis.text = element_blank(),axis.ticks = element_blank())
    motakis_scatters[[length(motakis_scatters)+1]] = motakis_scatter
  }
}

filtered_motakis = subset(filtered_motakis, subset = ((INS>6.5 & GCG>5)|(INS>6.5 & PPY>4.5)|(INS>6.5 & SST>5)|(INS>6.5 & GHRL>4)|(GCG>5 & PPY>4.5)|(GCG>5 & SST>5)|(GCG>5 & GHRL>4)|(PPY>4.5 & SST>5)|(PPY>4.5 & GHRL>4)|(SST>5 & GHRL>4)), invert = T)
dim(filtered_motakis) # 31683 genes 252097 cells

# zhang
zhang_ridges = list()
for (i in 1:5) {
  zhang_ridges[[i]] = RidgePlot(filtered_zhang,features=hormones[i],group.by = "data_source")
}

zhang_scatters = list()
for (i in 1:(n-1)) {
  for (j in (i+1):n){
    zhang_scatter = FeatureScatter(filtered_zhang, feature1 = hormones[i],feature2 = hormones[j], group.by = "data_source",pt.size = 0.00000)+
      labs(title = NULL, caption = NULL)+
      theme(legend.position = "none",axis.text = element_blank(),axis.ticks = element_blank())
    zhang_scatters[[length(zhang_scatters)+1]] = zhang_scatter
  }
}

filtered_zhang = subset(filtered_zhang, subset = ((INS>5.5 & GCG>5)|(INS>5.5 & PPY>5)|(INS>5.5 & SST>5)|(INS>5.5 & GHRL>5)|(GCG>5 & PPY>5)|(GCG>5 & SST>5)|(GCG>5 & GHRL>5)|(PPY>5 & SST>5)|(PPY>5 & GHRL>5)|(SST>5 & GHRL>5)), invert = T)
dim(filtered_zhang) # 17570 genes 50765 cells

# segerstolpe
segerstolpe_ridges = list()
for (i in 1:5) {
  segerstolpe_ridges[[i]] = RidgePlot(filtered_segerstolpe,features=hormones[i],group.by = "data_source")
}

segerstolpe_scatters = list()
for (i in 1:(n-1)) {
  for (j in (i+1):n){
    segerstolpe_scatter = FeatureScatter(filtered_segerstolpe, feature1 = hormones[i],feature2 = hormones[j], group.by = "data_source",pt.size = 0.00000)+
      labs(title = NULL, caption = NULL)+
      theme(legend.position = "none",axis.text = element_blank(),axis.ticks = element_blank())
    segerstolpe_scatters[[length(segerstolpe_scatters)+1]] = segerstolpe_scatter
  }
}

filtered_segerstolpe = subset(filtered_segerstolpe, subset = ((INS>6 & GCG>5)|(INS>6 & PPY>5.5)|(INS>6 & SST>5.5)|(INS>6 & GHRL>4)|(GCG>5 & PPY>5.5)|(GCG>5 & SST>5.5)|(GCG>5 & GHRL>4)|(PPY>5.5 & SST>5.5)|(PPY>5.5 & GHRL>4)|(SST>5.5 & GHRL>4)), invert = T)
dim(filtered_segerstolpe) # 39173 genes 2050 cells

# xin
xin_ridges = list()
for (i in 1:5) {
  xin_ridges[[i]] = RidgePlot(filtered_xin,features=hormones[i],group.by = "data_source")
}

xin_scatters = list()
for (i in 1:(n-1)) {
  for (j in (i+1):n){
    xin_scatter = FeatureScatter(filtered_xin, feature1 = hormones[i],feature2 = hormones[j], group.by = "data_source",pt.size = 0.00000)+
      labs(title = NULL, caption = NULL)+
      theme(legend.position = "none",axis.text = element_blank(),axis.ticks = element_blank())
    xin_scatters[[length(xin_scatters)+1]] = xin_scatter
  }
}

filtered_xin = subset(filtered_xin, subset = ((INS>6.5 & GCG>5.5)|(INS>6.5 & PPY>4)|(INS>6.5 & SST>4)|(INS>6.5 & GHRL>4)|(GCG>5.5 & PPY>4)|(GCG>5.5 & SST>4)|(GCG>5.5 & GHRL>4)|(PPY>4 & SST>4)|(PPY>4 & GHRL>4)|(SST>4 & GHRL>4)), invert = T)
dim(filtered_xin) # 20497 genes 19948 cells

# ------------------------------------------------------------------------------
# 6.  Combining datasets (hpap, motakis, zhang, segerstolpe, xin)
# ------------------------------------------------------------------------------

DefaultAssay(object = filtered_hpap) = "RNA"
panc_combined = merge(filtered_hpap, y = c(filtered_motakis,filtered_zhang,filtered_segerstolpe,filtered_xin), add.cell.ids = c("hpap", "motakis", "zhang", "segerstolpe","xin"), project = "AllPanc")
panc_combined 
# remove(hpap,motakis,zhang,segerstolpe,filtered_hpap,filtered_motakis,filtered_zhang,filtered_segerstolpe)
# saveRDS(panc_combined, "./rawdata/Seurat_objects/unintegrated_data.rds")

# ------------------------------------------------------------------------------
# 7.  Integrating datasets using harmony (hpap, motakis, zhang, segerstolpe, xin)
# ------------------------------------------------------------------------------

panc_combined = readRDS("./rawdata/Seurat_objects/unintegrated_data.rds")

# panc_combined = JoinLayers(panc_combined)
# panc_combined = split(panc_combined[["RNA"]], f = panc_combined$donor)
# panc_combined = CreateSeuratObject(panc_combined)

panc_combined = NormalizeData(panc_combined)

# Peek at highly variable genes
panc_combined = FindVariableFeatures(panc_combined, selection.method = "vst", nfeatures = 8000)

# Scaling
panc_combined = ScaleData(panc_combined)

# Perform PCA
panc_combined = RunPCA(panc_combined)

# Determining PCs to use downstream
ElbowPlot(panc_combined) # 15 PCs

# Visualize by data source before integrating
panc_combined = RunUMAP(panc_combined, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")

mycol= c("#E15759", "#4E79A7", "#76B7B2", "#F28E2B", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F","#7F3C8D", "#11A579", "#3969AC" ,"#F2B701" ,"#E73F74" ,"#80BA5A", "#E68310", "#008695", "#A5AA99")
DimPlot(panc_combined, reduction = "umap.unintegrated", group.by = "data_source")+
  scale_color_manual(values = mycol)+
  labs(title = "")+
  theme_void()
  
# Harmony integration
panc_integrated = IntegrateLayers(object = panc_combined, method = HarmonyIntegration,
                                  orig.reduction = "pca", new.reduction = "harmony",verbose = F)
# saveRDS(panc_integrated, "./rawdata/Seurat_objects/integrated_harm_data.rds")

# ------------------------------------------------------------------------------
# 8.  Clustering cells
# ------------------------------------------------------------------------------

panc_integrated = readRDS("./rawdata/Seurat_objects/integrated_harm_data.rds")
panc_integrated = FindNeighbors(panc_integrated, dims = 1:15, reduction = "harmony")
panc_integrated = FindClusters (panc_integrated, resolution = 0.9)

# resolution decision
seq = seq(0.9, 4.0, by = 0.5)
for (res in seq) {
  panc_integrated = FindClusters(panc_integrated, resolution = res)
}

clustree(panc_integrated, prefix = 'RNA_snn_res.') 

panc_integrated = RunUMAP(panc_integrated, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")

# Clustering visuals 
p = DimPlot(panc_integrated, reduction = "umap.harmony", group.by = "seurat_clusters",label = T)

# Visualize by batch after data integration
p1 = DimPlot(panc_integrated, reduction = "umap.harmony", group.by = "data_source")
p2 = DimPlot(panc_integrated, reduction = "umap.harmony", group.by = "gender")
p3 = DimPlot(panc_integrated, reduction = "umap.harmony", group.by = "donor") 
p4 = DimPlot(panc_integrated, reduction = "umap.harmony", group.by = "disease_status")

p1 | p2 | p4
p3 | p

# ------------------------------------------------------------------------------
# 9.  CellType annotation
# ------------------------------------------------------------------------------

panc_markers = c("GCG","TTR","INS","IAPP","SST","HHEX","LEPR",
                 "PPY","GHRL","KRT19","CFTR","KRT7","PRSS1","CPA1",
                 "CPA2","REG1A","COL1A1","PDGFRB","COL1A2","TIMP1","FLT1",
                 "VWF","PECAM1","CD68","CD74","CD163","TPSAB1","CPA3")

DotPlot(panc_integrated, features = panc_markers, group.by = "seurat_clusters")+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3497d8","white","#e23232"))+
  coord_flip()

new.cluster.ids = c("beta","alpha","beta","ductal","acinar","alpha","beta","alpha","alpha","stellate",
                    "acinar","beta","alpha","alpha","delta","beta","alpha","alpha","endothelial",
                    "beta","ductal","stellate","beta","acinar","stellate","PP/gamma","ductal","mast",
                    "endothelial","alpha","macrophage","delta","acinar","ductal","acinar","endothelial","beta",
                    "beta","stellate","alpha","endothelial","ductal","stellate")

names(new.cluster.ids) = levels(panc_integrated)
panc_integrated = RenameIdents(panc_integrated, new.cluster.ids)
panc_integrated@meta.data$cell_type = Idents(object = panc_integrated)

table(panc_integrated@meta.data$cell_type)
table(panc_integrated@meta.data$seurat_clusters)

DimPlot(panc_integrated, reduction = "umap.harmony", group.by = "cell_type")+
  scale_color_manual(values = mycol)+
  labs(title = "")+
  theme_void()

# Based on bimodality of cell type markers, 
# set threshold on log exp to remove counts below the threshold
# This is only for the purpose of cell annotation on featurePlots. Counts were not removed in the original gene expression matrix
panc_integrated = JoinLayers(panc_integrated)

panc_data = GetAssayData(object = panc_integrated, layer  = "data")
panc_data["INS",] = ifelse(panc_data["INS",]<7,0,panc_data["INS",])
panc_data["GCG",] = ifelse(panc_data["GCG",]<5.5,0,panc_data["GCG",])
panc_data["SST",] = ifelse(panc_data["SST",]<5,0,panc_data["SST",])
panc_data["PPY",] = ifelse(panc_data["PPY",]<5,0,panc_data["PPY",])
panc_data["CD68",] = ifelse(panc_data["CD68",]<3,0,panc_data["CD68",])

panc_filtered = SetAssayData(panc_integrated, layer = "data", new.data = panc_data)

INS = FeaturePlot(panc_filtered, features = "INS",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1)) 
GCG = FeaturePlot(panc_filtered, features = "GCG",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1)) 
PPY = FeaturePlot(panc_filtered, features = "PPY",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1)) 
SST = FeaturePlot(panc_filtered, features = "SST",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1)) 
CPA1 = FeaturePlot(panc_filtered, features = "CPA1",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1))
KRT19 = FeaturePlot(panc_filtered, features = "KRT19",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1))
PDGFRB = FeaturePlot(panc_filtered, features = "PDGFRB",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1))
FLT1 = FeaturePlot(panc_filtered, features = "FLT1",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1))
CD68 = FeaturePlot(panc_filtered, features = "CD68",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1))
TPSAB1 = FeaturePlot(panc_filtered, features = "TPSAB1",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1))


INS + GCG + PPY + SST + CPA1 + KRT19 + PDGFRB + FLT1 + CD68 + TPSAB1 + plot_layout(ncol = 5)

saveRDS(panc_integrated, "./rawdata/Seurat_objects/integrated_data.rds")

# ------------------------------------------------------------------------------
# 10.  CellType visuals 
# ------------------------------------------------------------------------------

# celltype cluster and cell number count
panc_integrated = readRDS("./rawdata/integrated_data.rds")

cell_counts = as.data.frame(table(panc_integrated$cell_type))
p_umap = DimPlot(panc_integrated, reduction = "umap.harmony", group.by = "cell_type", label = T)+
  scale_color_manual(values = mycol)+
  labs(title = "")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1),
    ) 
  
p_bar = ggplot(cell_counts, aes(x = reorder(Var1,Freq), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  coord_flip()+
  scale_fill_manual(values = mycol) +
  labs(x = "", y = "Cell Number") +
  theme_bw()+
  theme( 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1)
    )  

cell_legend = get_legend(p_umap)
p_umap = p_umap + theme(legend.position = "none")
p_umap + p_bar + cell_legend + plot_layout(widths = c(2,0.7,0.7))

# cell marker visuals
panc_markers = c("GCG","TTR","INS","IAPP","SST","HHEX","LEPR","PPY", "KRT19","CFTR","KRT7",
                 "PRSS1","CPA1", "CPA2","REG1A","COL1A1","PDGFRB","COL1A2","TIMP1","FLT1",
                 "VWF","PECAM1","TPSAB1","CPA3","CD68","CD74","CD163"
                )

mycol2 = c("#adc8e8","#ffbc78","#99e08b","#fe9896","#c5b2d6","#c49d93","#f8b8d2","#dadada","#ffe8a3","#95cfcd")

VlnPlot(panc_integrated, features = panc_markers, stack = TRUE,fill.by = "ident",cols = mycol)+
  labs(y = "")

DotPlot(panc_integrated, features = panc_markers, group.by = "cell_type")+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3497d8","white","#e23232"))

# heatmap: the top highly expressed genes in each of the cell clusters.
table(panc_integrated@meta.data$cell_type)
markers = FindAllMarkers(panc_integrated, logfc.threshold = 0.5, only.pos = TRUE,min.pct = 0.25)
markers = markers[which(markers$p_val_adj < 0.05),]
sig_markers = markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
head(sig_markers)
genes = unique(sig_markers$gene)

aver_exp = AverageExpression(panc_integrated,features = genes,group.by = 'cell_type',layer = 'data')
aver_exp = as.data.frame(aver_exp$RNA)
aver_exp = t(scale(t(aver_exp)))
aver_exp = aver_exp[sig_markers$gene,]

# heatmap visuals

mycol = colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))
cell_anno = data.frame(cell_anno = colnames(aver_exp), row.names = colnames(aver_exp))

cols = c("#afc8e9","#ffbc79","#99e08b","#ff9997","#c6b2d6","#dadada","#f7b7d3","#ffe8a3","#9fdbe6","#c49d94")
names(cols) = cell_anno$cell_anno

# col annotations
cell = data.frame(colnames(aver_exp))
colnames(cell) = 'cell'

col_anno = HeatmapAnnotation(df = cell,show_annotation_name = F,
                             gp = gpar(col = 'white', lwd = 2),col = list(cell = cols))

# row annotations
row_cols = setNames(rep(cols, each = 5), rownames(aver_exp))
row_anno = rowAnnotation(foo = anno_text(rownames(aver_exp),location = 0,just = "left",
                                         gp = gpar(fill = row_cols,col = "black",fontface = 'italic'),
                                         width = max_text_width(rownames(aver_exp))*1.25))

Heatmap(aver_exp,
        name = 'expression',
        col = mycol,
        cluster_columns = F,
        cluster_rows = F,
        column_names_side = c('top'),
        column_names_rot = 60,
        row_names_gp = gpar(fontsize = 12, fontface = 'italic'),
        rect_gp = gpar(col = "white", lwd = 1.5),
        top_annotation = col_anno) + row_anno

write.csv(markers,"./celltype_markers.csv")

