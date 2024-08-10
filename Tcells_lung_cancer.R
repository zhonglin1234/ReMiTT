
############################################################
#Load scRNA data of T cells from lung cancer
############################################################

require(openxlsx)
require(dplyr)
require(Seurat)

setwd('xx')
#############Smart data
lung_smart=read.xlsx("GSE99254/lung_smart_TCRtable.xlsx",sheet=1)
lung_names=gsub('-','.',lung_smart$Cell.Name)
lung_smart$Cell.Name=lung_names

#############Raw scRNA data
lung_scRNA=read.csv("/xx/GSE99254/GSE99254_NSCLC.TCell.S12346.count.txt",sep='\t',header=T)
lung_names=lung_names[which(lung_names %in% colnames(lung_scRNA))]
lung_scRNA=lung_scRNA[,c('symbol',lung_names)]
dim(lung_scRNA) 
lung_scRNA=lung_scRNA[which(!is.na(lung_scRNA$symbol)),]
rownames(lung_scRNA)=lung_scRNA[,'symbol']
lung_scRNA=lung_scRNA[,2:ncol(lung_scRNA)]  
dim(lung_scRNA) #23370 10096


#Seurat object
lung_scRNA.sparse=as(as.matrix(lung_scRNA),"sparseMatrix")
rm(lung_scRNA)
lung_seurat=CreateSeuratObject(counts=lung_scRNA.sparse,project='scRNA',min.cells = 20,min.features = 500, assay = 'RNA') 
rm(lung_scRNA.sparse)

#Normalize the data
lung_seurat <- NormalizeData(lung_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
#Find the highly variable genes
lung_seurat=FindVariableFeatures(lung_seurat, selection.method = "vst", nfeatures = 2000)
#Scale the data
all.genes <- rownames(lung_seurat)
lung_seurat <- ScaleData(lung_seurat, features = all.genes)
#Run PCA
lung_seurat <- RunPCA(lung_seurat, features = VariableFeatures(object = lung_seurat))

lung_seurat <- FindNeighbors(lung_seurat, dims = 1:30)
lung_seurat <- FindClusters(lung_seurat, resolution = .5)

head(Idents(lung_seurat), 5)

lung_seurat<- RunTSNE(lung_seurat, dims = 1:30)
DimPlot(lung_seurat, reduction = "tsne")

##Save the data
save(lung_seurat, file="xx/all_lung_cells_seurat.Rdata")
