library(Seurat)
library(Matrix)
library(ggplot2)
library(utils)
library(dplyr)
library(qusage)
library(cowplot)

# import the data as Seurat S4 datasets
load_matrix <- function(file) {
  file_name <- paste0('../rna-seq/',file, "/filtered_feature_bc_matrix/")
  print(file_name)
  barcode.path <- paste0(file_name, "barcodes.tsv.gz")
  features.path <- paste0(file_name,"features.tsv.gz")
  matrix.path <- paste0(file_name,"matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names <- read.delim(features.path,
                              header = FALSE,
                              stringsAsFactors = FALSE
  )
  barcode.names <- read.delim(barcode.path,
                              header = FALSE,
                              stringsAsFactors = FALSE
  )
  colnames(mat) <- barcode.names$V1
  rownames(mat) <- feature.names$V2
  
  # setup data
  mat <- CreateSeuratObject(
    counts = mat,
    project = "Adipocyte",
    min.cells = 3,
    min.features = 200
  )
  mat$phenotype <- file  # set the phenotype
  return(mat)
}

# preprocess
down_stream <- function(data) {
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data <- subset(
    data, 
    subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 40
  )
  data <- NormalizeData(
    data, 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
  )
  data <- NormalizeData(data)
  data <- FindVariableFeatures(
    data, 
    selection.method = "vst", 
    nfeatures = 2000
  )
  all.genes <- rownames(data)
  data <- ScaleData(data, features = all.genes)  # scale the data
  data <- RunPCA(data, features = VariableFeatures(object = data))
  data <- JackStraw(data, num.replicate = 100)
  data <- ScoreJackStraw(data, dims = 1:20)
  return(data)
}

# clustering
cluster_vis <- function(data) {
  data <- FindNeighbors(data, dims = 1:20)
  data <- FindClusters(data, resolution = 0.7)
  data <- RunUMAP(data, dims = 1:20)
  return(data)
}

# plot the experssion of 'Adipoq', 'Cdh5', and 'Pdgfra'
analysis <- function(data) {
  data <- down_stream(data)
  data <- cluster_vis(data)
  VlnPlot(data, c('Adipoq', 'Cdh5', 'Pdgfra'))
  return(data)
}

# reperform the subset data
reperform <- function (data, name) {
  data <- GetAssay(data)@counts
  data <- CreateSeuratObject(
    counts = data, 
    project = "Adipocyte", 
    min.cells = 0, 
    min.features = 0
  )
  data$phenotype <- name
  data <- down_stream(data)
  data <- cluster_vis(data)
  return(data)
}

# merge the subset 
merge_downstream <- function (data1, data2) {
  data <- merge(data1, data2)
  data <- down_stream(data)
  data <- cluster_vis(data)
  return(data)
}

# find the markers between GFPpos and GFPneg
marker_table <- function(data, top_gene) {
  Idents(data) <- data$phenotype
  markers <- FindAllMarkers(
    data,
    only.pos = TRUE,
    min.pct = 0.2,
    logfc.threshold = 0.25
  )
  return(markers)   
}

# for four datasets
if (getOption('run.main', default=TRUE)) {
  # read data from 4 datasets
  pb <- txtProgressBar(style=3)
  file_name <- c('iBATneg', 'iBATpos', 'iWATneg', 'iWATpos')
  for (i in 1:length(file_name)) {
    assign(file_name[i], load_matrix(file_name[i]))
    setTxtProgressBar(pb, i/length(file_name))
  }
  close(pb)
  
  # for the iBATneg data
  iBATneg <- analysis(iBATneg)
  VlnPlot(iBATneg, c('Adipoq', 'Cdh5', 'Pdgfra'))
  iBATneg_subset <- subset(iBATneg, idents = c(0, 4, 8)) # select only the 0, 4 and 8
  
  # for the iBATpos data
  iBATpos <- analysis(iBATpos)
  VlnPlot(iBATpos,c('Adipoq', 'Cdh5', 'Pdgfra'))
  FeaturePlot(iBATpos, features = c("Adipoq"), cols = c("grey", "red"), pt.size = 1)
  DimPlot(iBATpos, label = T, pt.size = 1, label.size = 6)
  iBATpos_subset <- subset(iBATpos, idents = c(0, 1, 2, 3, 4, 5, 7)) # select only the 0, 1, 2, 3, 4, 5 and 7
  
  # for the iWATneg data
  iWATneg <- analysis(iWATneg)
  VlnPlot(iWATneg,c('Adipoq', 'Cdh5', 'Pdgfra'))
  DimPlot(iWATneg, label = T, pt.size = 1, label.size = 6)
  FeaturePlot(iWATneg, features = c("Adipoq"), cols = c("grey", "red"), pt.size = 1)
  FeaturePlot(iWATneg, features = c("Pdgfra"), cols = c("grey", "red"), pt.size = 1)
  iWATneg_subset <- subset(iWATneg, idents = c(4, 6))  #  select only the  4 and 6
  
  # for the iWATpos data
  iWATpos <- analysis(iWATpos)
  VlnPlot(iWATpos,c('Adipoq', 'Cdh5', 'Pdgfra'))
  iWATpos_subset <- subset(iWATpos, idents = c(1, 2, 6, 8)) # select only the 1, 2, 6 and 8
  
  # reanalysis the subset data 
  iBATneg_subset <- reperform(iBATneg_subset, 'iBATneg')
  iBATpos_subset <- reperform(iBATpos_subset, 'iBATpos')
  iWATneg_subset <- reperform(iWATneg_subset, 'iWATneg')
  iWATpos_subset <- reperform(iWATpos_subset, 'iWATpos')
  
  # merge the data
  iBAT <- merge_downstream(iBATneg_subset, iBATpos_subset)  # for iBAT cells
  iWAT <- merge_downstream(iWATneg_subset, iWATpos_subset)  # for iWAT cells
  
  # save the UMAP figure for iBAT and iWAT
  pdf('/path/to/iBAT_umap.pdf')
  DimPlot(iBAT, group.by = 'phenotype', pt.size = 1) + NoLegend()
  dev.off()
  
  pdf('/path/to/iWAT_umap.pdf')
  DimPlot(iWAT, group.by = 'phenotype', pt.size = 1) + NoLegend()
  dev.off()
  
  # Heat map
  iBAT_markers <- marker_table(iBAT)
  top30_B <- iBAT_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
  DoHeatmap(iBAT, features = top30_B$gene, group.by = "phenotype")
}