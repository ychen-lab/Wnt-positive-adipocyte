library(Signac)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(JASPAR2020)
library(TFBSTools)
library(ggsignif)
library(motifmatchr)
library(chromVAR)
library(EnsDb.Mmusculus.v79)

### --------- functions ------- ####
# import data
import_data <- function(path) {
  
  # count
  count_path <- paste0(
    '../scATAC-rawdata/',
    path,
    '/filtered_peak_bc_matrix.h5'
  )
  counts <- Read10X_h5(count_path)
  
  # meta
  meta_path <- paste0(
    '../scATAC-rawdata/',
    path,
    '/singlecell.csv'
  )
  meta <- read.csv(
    meta_path,
    header = T,
    row.names = 1
  )
  
  # Chromatinassay
  fragment_path <- paste0(
    '../scATAC-rawdata/',
    path,
    '/fragments.tsv.gz'
  )
  
  # create chromatin assay
  assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(':', '-'),
    genome = 'mm10',
    fragments = fragment_path,
    min.cells = 1
  )
  
  # Create Seurat
  data <- CreateSeuratObject(
    counts = assay,
    assay = 'peaks',
    project = 'ATAC',
    meta.data = meta
  )
  return(data)
}

# data preprocess and filtering
prep_fil <- function(data, peak_limit){
  Annotation(data) <- annotations
  data <- NucleosomeSignal(object = data)
  data$nucleosome_group <- ifelse(data$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  data <- TSSEnrichment(data, fast = FALSE)
  data$high.tss <- ifelse(data$TSS.enrichment > 2, 'High', 'Low')
  data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100
  data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments
  data_sub <- subset(
    x = data,
    subset = peak_region_fragments > 2000 &
      peak_region_fragments < peak_limit &
      pct_reads_in_peaks > 20 &
      blacklist_ratio < 0.025 &
      nucleosome_signal < 4 &
      TSS.enrichment > 2
  )
  return(data_sub)
}

# find clusters
find_cluster <- function(data) {
  data <- RunTFIDF(data)
  data <- FindTopFeatures(data, min.cutoff = 'q0')
  data <- RunSVD(data)
  data <- RunUMAP(data, dims = 2:30, reduction = 'lsi')
  data <- FindNeighbors(data,reduction = 'lsi', dims = 2:30)
  data <- FindClusters(data, algorithm = 3, resolution = 1.2, verbose = F)
  data.gene.activities <- GeneActivity(data)
  data[['RNA']] <- CreateAssayObject(counts = data.gene.activities)
  data <- NormalizeData(
    object = data,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(data$nCount_RNA)
  )
  DefaultAssay(data) <- 'RNA'
  return(data)
}

### ----- analysis ------ ###
# import data
data_name <- c('BATpos', 'BATneg', 'iWATpos', 'iWATneg')
for (i in 1:4) {
  data <- import_data(data_name[i])
  assign(data_name[i], data)
}

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# preporcess and filtering the data
BATpos <- prep_fil(BATpos, 40000)
BATneg <- prep_fil(BATneg, 40000)
iWATpos <- prep_fil(iWATpos, 60000)
iWATneg <- prep_fil(iWATneg, 90000)

# find adipocyte clusters
BATpos <- find_cluster(BATpos)
VlnPlot(BATpos, c("Adipoq"))
BATpos <- subset(BATpos, idents = c(0, 1, 2, 3, 4, 5, 6, 8))

BATneg <- find_cluster(BATneg)
VlnPlot(BATneg, c("Adipoq"))
BATneg <- subset(BATneg, idents = c(1, 4, 6))

iWATpos <- find_cluster(iWATpos)
VlnPlot(iWATpos, c("Adipoq"))
iWATpos <- subset(iWATpos, idents = c(0, 1, 3, 5, 7, 8))

iWATneg <- find_cluster(iWATneg)
VlnPlot(iWATneg, c("Adipoq"))
iWATneg <- subset(iWATneg, idents = c(0, 1, 3, 4, 6, 8))

# add information to identify dataset of origin
iBATpos$dataset <- 'iBATpos'
iBATneg$dataset <- 'iBATneg'
iWATpos$dataset <- 'iWATpos'
iWATneg$dataset <- 'iWATneg'

# merge all datasets, adding a cell ID to make sure cell names are unique
iBAT <- merge(
  x = iBATpos,
  y = list(iBATneg),
  add.cell.ids = c("iBATpos", "iBATneg")
)
iWAT <- merge(
  x = iWATpos,
  y = list(iWATneg),
  add.cell.ids = c("iWATpos", "iWATneg")
)


# normal analysis for iBAT
iBAT <- RunTFIDF(iBAT)
iBAT <- FindTopFeatures(iBAT, min.cutoff = 20)
iBAT <- RunSVD(iBAT)
iBAT <- RunUMAP(iBAT, dims = 1:50, reduction = "lsi")
iBAT <- FindNeighbors(
  object = iBAT,
  reduction = "lsi",
  dims = 2:30
)
iBAT <- FindClusters(
  object = iBAT,
  algorithm = 3,
  resolution = 0.7,
  verbose = FALSE
)
DimPlot(iBAT, group.by = 'dataset', pt.size = 1) + NoLegend()
DimPlot(iBAT, label = T, pt.size = 1, label.size = 6) + NoLegend()

DefaultAssay(iBAT) <- 'peaks'
Annotation(iBAT) <- annotations
iBAT.gene.activities <- GeneActivity(iBAT)

# add the gene activity matrix to the Seurat object as a new assay
iBAT[["RNA"]] <- CreateAssayObject(counts = iBAT.gene.activities)
iBAT <- NormalizeData(
  object = iBAT,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = median(iBAT$nCount_RNA)
)

DefaultAssay(iBAT) <- "peaks"
Idents(iBAT) <- "dataset"
iBAT.markers <- FindAllMarkers(
  iBAT,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25
)
iBAT_topmarkers <- iBAT.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_logFC)
write.csv(iBAT_topmarkers, "/path/to/iBAT_markers.csv")
saveRDS(iBAT, "/path/to/iBAT.RDS")

# ----------------- iWAT datasets -----------------------------------
# normal analysis for iWAT
iWAT <- RunTFIDF(iWAT)
iWAT <- FindTopFeatures(iWAT, min.cutoff = 20)
iWAT <- RunSVD(iWAT)
iWAT <- RunUMAP(
  iWAT,
  dims = 1:50,
  reduction = "lsi"
)
iWAT <- FindNeighbors(
  object = iWAT,
  reduction = "lsi",
  dims = 2:30
)
iWAT <- FindClusters(
  object = iWAT,
  algorithm = 3,
  resolution = 0.7,
  verbose = FALSE
)
DimPlot(iWAT, group.by = 'dataset', pt.size = 1) + NoLegend()
DimPlot(iWAT, label = T, pt.size = 1, label.size = 6) + NoLegend()

DefaultAssay(iWAT) <- 'peaks'
Annotation(iWAT) <- annotations
iWAT.gene.activities <- GeneActivity(iWAT)

# add the gene activity matrix to the Seurat object as a new assay
iWAT[["RNA"]] <- CreateAssayObject(counts = iWAT.gene.activities)
iWAT <- NormalizeData(
  object = iWAT,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = median(iWAT$nCount_RNA)
)

DefaultAssay(iWAT) <- "peaks"

pdf("/path/to/iWAT_dimplot.pdf")
DimPlot(iWAT, group.by = "dataset")
dev.off()

Idents(iWAT) <- "dataset"
iWAT.markers <- FindAllMarkers(
  iWAT,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25
)
iWAT_topmarkers <- iWAT.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_logFC)
write.csv(iWAT_topmarkers, "/path/to/iWAT_markers.csv")
saveRDS(iWAT, "/path/to/iWAT.RDS")

# load data
iBAT <- readRDS('/path/to/iBAT.RDS')
iWAT <- readRDS('/path/to/iWAT.RDS')

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(iBAT),
  pwm = pfm,
  genome = 'mm10',
  use.counts = FALSE
)

motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

iBAT <- SetAssayData(
  object = iBAT,
  assay = 'peaks',
  slot = 'motifs',
  new.data = motif
)
iBAT[["peaks"]]

iBAT <- RegionStats(object = iBAT, genome = BSgenome.Mmusculus.UCSC.mm10)

iBAT <- RunChromVAR(
  object = iBAT,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(iBAT) <- 'chromvar'

iBATpos.differential.activity <- FindMarkers(
  object = iBAT,
  ident.1 = 'iBATpos',
  ident.2 = 'iBATneg',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
MotifPlot(
  object = iBAT,
  motifs = head(rownames(iBATpos.differential.activity)),
  assay = 'peaks'
)

iBATneg.differential.activity <- FindMarkers(
  object = iBAT,
  ident.1 = 'iBATneg',
  ident.2 = 'iBATpos',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
MotifPlot(
  object = iBAT,
  motifs = head(rownames(iBATneg.differential.activity)),
  assay = 'peaks'
)

# Find all markers
iBAT.motifs <- FindAllMarkers(
  iBAT, only.pos = TRUE, 
  min.pct = 0.2, 
  logfc.threshold = 0.25
)

key <- iBAT.motifs %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

pdf('/path/to/iBAT.motifs.pdf', width = 15, height = 8)
VlnPlot(iBAT, key$gene, n = 5, pt.size = F)
dev.off()

FeaturePlot(iBAT, key$gene, n = 5,
            min.cutoff = 'q10',
            max.cutoff = 'q90',
            pt.size = 1)

MotifPlot(
  object = iBAT,
  motifs = key$gene,
  assay = 'peaks'
)

#### ---------- for iWAT --------------------- #############
# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(iWAT),
  pwm = pfm,
  genome = 'mm10',
  use.counts = FALSE
)

motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

iWAT <- SetAssayData(
  object = iWAT,
  assay = 'peaks',
  slot = 'motifs',
  new.data = motif
)

iWAT <- RegionStats(object = iWAT, genome = BSgenome.Mmusculus.UCSC.mm10)

iWAT <- RunChromVAR(
  object = iWAT,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(iWAT) <- 'chromvar'

# gene markers of different pheno
iWATpos.differential.activity <- FindMarkers(
  object = iWAT,
  ident.1 = 'iWATpos',
  ident.2 = 'iWATneg',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
MotifPlot(
  object = iWAT,
  motifs = head(rownames(iWATpos.differential.activity)),
  assay = 'peaks'
)


iWATneg.differential.activity <- FindMarkers(
  object = iWAT,
  ident.1 = 'iWATneg',
  ident.2 = 'iWATpos',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
MotifPlot(
  object = iWAT,
  motifs = head(rownames(iWATneg.differential.activity)),
  assay = 'peaks'
)

iWAT.motifs <- FindAllMarkers(iWAT, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
key <- iWAT.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

MotifPlot(
  object = iWAT,
  motifs = key$gene,
  assay = 'peaks'
)

FeaturePlot(
  object = iWAT,
  features = key$gene,
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  ncol = 3
)

