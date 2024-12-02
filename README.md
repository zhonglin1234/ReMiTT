# T Cell Migrating Trail Identification Workflow

This repository contains functions and data (downloaded from Visium 10x website) for identifying T cell migrating trails with spatial RNA-seq data.

1. To run this tutorial, you **_MUST DOWNLOAD_** the **HDF5 matrix file** from the following link:
   [Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.h5](https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.h5)

2. Put the `Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.h5` file under the folder: **data/visium_ovarian**

After doing this, you can start running the codes in this tutorial.

3. For applying this method to other datasets, please use samples with abundant T cell infiltration. Migration trails can hardly be identified with samples with low-level T cell infiltration.



## Installation and Setup

Ensure you have all required R packages installed and your working directory is correctly set. Before running the code, source the necessary functions:
```r
source("functions/get_paths_git.R")
source("functions/path_downstream_git.R")
source("functions/packaged_functions_git.R")
source("functions/figures_tables_git.R")
```

## Step 1: Load Data

#### 1.1 Load CD4/CD8 T Cell Stage Markers
```r
load('data/cd4_8_cluster_markers.Rdata')
all.clus.markers <- unique(clus.markers[, 1])
```

#### 1.2 Load Cancer Markers
```r
load('data/cancer_markers.Rdata')
```

#### 1.3 Load GSEA Data
```r
c5 <- msigdbr(species = 'Homo sapiens', category = 'C5')
c5.go <- c5[which(c5$gs_subcat != 'HPO'), ]
```

## Step 2: Prepare Spatial Data

```r
dat1 <- data_preparation(
  folder = "data/visium_ovarian",  
  filenm = "Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.h5", 
  res = 0.2,
  loc_file = 'data/visium_ovarian/spatial/tissue_positions_list.csv',
  ex.list = c('PDCD1', 'GZMB', 'LAG3', 'HAVCR2', 'TIGIT', 'CD37', 'CTLA4'), #markers of T cell exhaustion, check if you have these in your dataset
  cell = c('CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B'), # T cell surface markers, check if you have these in your data set
  cell_cut = 5,
  cut.t = 1e-10,
  n.vgg = 3000,
  clustering = TRUE  ### We recommend using clustering = FALSE for other datasets
)

# Extract prepared data
sp.counts.r <- dat1[[1]]
sp.counts.norm <- dat1[[2]]
vgg <- dat1[[3]]
loc.raw <- dat1[[4]]
meta_dat <- dat1[[5]]
keep.ids <- dat1[[6]]
meta.cts <- dat1[[7]]
loc.cts <- dat1[[8]]
pt.cts <- dat1[[9]]
r.unit <- dat1[[10]]
sub.clusters <- dat1[[11]]
sp.counts.tpm <- dat1[[12]]
```

## Step 3: Visualize Clusters (optional)

#### 3.1 Inspect Seurat Clusters and identify potential clusters of tumor

```r
tmp.mean <- plot_seurat_clusters(
  loc = loc.raw,
  meta_dat = meta_dat,
  sp.counts.norm = sp.counts.norm,
  sp.counts.r = sp.counts.r,
  cancer.markers = cancer.markers$ovarian
)
tumor.clus <- order(tmp.mean[[2]][, 1], decreasing = TRUE)[1:3] - 1
```

#### 3.2 Visualize Tumor Regions and T Cell Spots

```r
plot_tumor(loc = loc.raw, meta_dat = meta_dat, tumor_cluster = tumor.clus)
points(loc[keep.ids, 'x'], loc[keep.ids, 'y'], col = alpha('green', 0.4), pch = 19, cex = 0.5)
```

## Step 4: Analyze Exhaustion Scores (optional)

#### 4.1 Check distribution of exhaustion Score
```r
plot_pt(loc = loc.raw, meta_dat = meta_dat)
```

#### 4.2 Exhaustion Score by Seurat clusters
```r
boxplot_pt(loc = loc.raw, meta_dat = meta_dat)
```

## Step 5: Identify Potential Trails

```r
potential.paths <- get_all_paths(
  loc.cts = loc.cts,
  pt.cts = pt.cts,
  cut.xydist = r.unit * 1.3,
  cut.sdist.n = 6,
  pt.dist.cuts.n = 11,
  min.length = 5,
  r.unit = r.unit,
  nr.graph = 1.3
)
```

## Step 6: Filter Trails by Correlation

```r
gpt.norm <- path_selection(
  paths = potential.paths,
  meta_dat = meta_dat,
  meta.cts = meta.cts,
  loc.cts = loc.cts,
  pt.cts = pt.cts,
  r.unit = r.unit,
  tcell_spots = keep.ids,
  sub.clusters = sub.clusters,
  sp.counts.norm = sp.counts.norm,
  tumor.clus = tumor.clus,
  cor.cut = 0.2,
  cell = c('CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B'), #markers of T cell exhaustion, check if you have these in your dataset
  ex.markers = c('PDCD1', 'GZMB', 'LAG3', 'HAVCR2', 'TIGIT', 'CD37', 'CTLA4') # T cell surface markers, check if you have these in your data set
)  ##These are the final algorithm-defined T cell migrating trails
```

## Step 7: Visualize Trails
```r
plot_all_shortest_paths(all.paths = gpt.norm, loc = loc.cts, length.cut = 5)
```

## Step 8: Systematic analysis (optional)

#### 8.1 Assign cluster label for each migration trails
```r
path.clus=unlist(lapply(gpt.norm,function(x)check_path_cluster(x,meta.cts=meta.cts)))
table(path.clus)
```

#### 8.2 Get 10000 control trail sets
```r
perm.pathset.norm.all=get_control_paths(paths=gpt.norm,meta.cts=meta.cts,nr=1.3,s.range1=7*r.unit,s.range2=3*r.unit,n.paths=90000,loc.cts=loc.cts,r.unit=r.unit,n.set=10000)
```

#### 8.3 Calculate the mean expression of each gene in each control set and the real migration trail set
```r

#10000 control sets
all.means.vgg=c()
for (i in 1:length(perm.pathset.norm.all)){
  pathset=perm.pathset.norm.all[[i]]
  pathset.spots=lapply(pathset,function(x) unlist(sub.clusters[x]))
  pathset.genes=lapply(pathset.spots,function(x) sp.counts.tpm[x,vgg])
  pathset.mean=do.call(rbind,lapply(pathset.genes,function(x) apply(x,2,mean)))
  all.means.vgg=c(all.means.vgg,list(pathset.mean))
}
controlset.means=do.call(rbind,lapply(all.means.vgg,function(x) apply(x, 2, mean))) #10000 control sets: 10000*2957

# Real set
gpt.mat.norm=get_gpt_gg1(dd=sp.counts.tpm,pathlist=gpt.norm,gg.markers=vgg) 
pathset.mean.all=apply(gpt.mat.norm,2,mean) #real pathset: 1*2957
```


#### 8.4 GSEA
```r
gg.all=vocano_plot1(pathset.mean=pathset.mean.all,controlset.means=controlset.means,plot=TRUE)
gg.all.enrich=enricher(gene=gg.all[[1]],TERM2GENE = c5.go[,c('gs_name','gene_symbol')])
enrich_barplot(gg.all.enrich,1:50)
```


