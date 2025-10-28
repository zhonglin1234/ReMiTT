# T Cell Migrating Trail Identification Workflow

This repository contains functions and data (downloaded from Visium 10x website) for identifying T cell migrating trails with Visium HD spatial RNA-seq data.

1. To run this tutorial, you **_MUST DOWNLOAD_** the **HDF5 matrix file** from the following link:
https://cf.10xgenomics.com/samples/spatial-exp/3.1.1/Visium_HD_Human_Lung_Cancer_Fixed_Frozen/Visium_HD_Human_Lung_Cancer_Fixed_Frozen_binned_outputs.tar.gz

2. Use the data in the folder of 'square_016um' and put this folder under: **data_HD**

3. For applying this method to other datasets, please use samples with abundant T cell infiltration. Migration trails can hardly be identified with samples with low-level T cell infiltration.


## Installation and Setup

Ensure you have all required R packages installed and your working directory is correctly set. Before running the code, source the necessary functions:
```r
source("functions_HD/path_downstream_HD_git.R")
source("functions_HD/packaged_functions_HD_git.R")
```

## Step 1: Load Data

#### 1.1 Load CD4/CD8 T Cell Stage Markers
```r
load('data_HD/cd4_8_cluster_markers.Rdata')
all.clus.markers <- unique(clus.markers[, 1])
```

#### 1.2 Load Cancer Markers
```r
load('data_HD/cancer_markers.Rdata')
```

#### 1.3 Load GSEA Data
```r
c5 <- msigdbr(species = 'Homo sapiens', category = 'C5')
c5.go <- c5[which(c5$gs_subcat != 'HPO'), ]
```

## Step 2: Prepare Spatial Data

```r
tcell_markers=c('CD3D','CD3E','CD3G','CD4','CD8A','CD8B2')
ex_markers=c('PDCD1','LAG3','HAVCR2','TIGIT','ENTPD1','CTLA4','TOX')


dat1=data_preparation_HD(folder='/data_HD/square_016um',  
                         filenm = "filtered_feature_bc_matrix.h5", 
                         res=NA,
                         loc_file=NA,
                         ex.list=ex_markers,
                         cell=tcell_markers,
                         cell_cut=1,
                         n.vgg=3000)


sp.counts.r=dat1[[1]]
sp.counts.norm=dat1[[2]]
sp.counts.tpm=dat1[[3]]
vgg=dat1[[4]]
all.gg=colnames(sp.counts.r)
meta_dat=dat1[[6]] 
loc.raw=dat1[[5]][rownames(meta_dat),c('x','y')]
meta_dat=data.frame(cbind(meta_dat,loc.raw))

pt=meta_dat$ex.score
names(pt)=rownames(meta_dat)
rm(dat1)
```

## Step 3: Assign small clusters for pixels adjacent to each other with similar exhaustion scores


#### 3.1 Get unit distance between adjacent pixels

spot.sdist=as.matrix(dist(loc.raw[which(loc.raw[,'x']< min(loc.raw[,'x'])+1000),],diag = T)) 
r.unit.raw=min(spot.sdist[which(spot.sdist!=0)])*1.5
rm(spot.sdist)
print(r.unit.raw)

#### 3.2 Assign sub-clusters based on spatial distance and similarity in exhaustion scores

nbr_radius=r.unit.raw*5 #radius to define neighbors

subclusters=assign_subcluster_greedy(meta_dat,cut.t=nbr_radius,pt_thresh = 0.01)
sub.clusters=subclusters[[2]]
print(length(sub.clusters)) 

#### 3.3 Get meta data for the subclusters (nodes)
meta.cts=meta_subclusters(sub.clusters,loc.raw,pt=pt,meta_dat)
loc.cts=meta.cts[,c('x','y')]
pt.cts=meta.cts$pt
names(sub.clusters)=rownames(meta.cts)


#### 3.4 Get neighbors for each node
coords <- as.matrix(meta.cts[, c("x","y")])
rownames(coords) <- rownames(meta.cts)
pt.cts=meta.cts$pt
names(pt.cts)=rownames(meta.cts)

nbr <- frNN(coords, eps = nbr_radius, sort = FALSE)
nbr.list=nbr$id 
nbr.list=lapply(nbr.list,function(x) meta.cts[x,])

nbr.large <- frNN(coords, eps = nbr_radius*5, sort = FALSE)
nbr.list.large=nbr.large$id 
nbr.list.large=lapply(nbr.list.large,function(x) rownames(meta.cts[x,]))


## Step 4: Identify Potential Trails and Filter Trails by Correlation

#### 4.1 get potential trails
```r
potential.paths=get_all_paths_HD(loc.cts=loc.cts, #rescaled XY coordinate of node data
                                 pt.cts=pt.cts, # Exhaustion score of nodes
                                 nbr.list=nbr.list,
                                 cut.length = 800,
                                 min.length=7)


```

#### 4.2 Select trails with large correlation between consecutive pairs of nodes

```r

#Take a look at Seurat clusters and define tumor clusters
tmp.mean=plot_seurat_clusters(loc=loc.raw,meta_dat=meta_dat,sp.counts.norm=sp.counts.norm,sp.counts.r=sp.counts.r,cancer.markers=cancer.markers$lung)
tumor.clus=c(0,1,3) #pick tumor clusters based on expression of cancer markers

gpt.norm=path_selection(paths=potential.paths,
                            meta_dat=meta_dat,
                            meta.cts=meta.cts,
                            loc.cts=loc.cts,
                            pt.cts=pt.cts,
                            sub.clusters=sub.clusters,
                            sp.counts.norm=sp.counts.norm,
                            tumor.clus=tumor.clus,
                            cor.cut=0.001,
                            cell=tcell_markers,
                            ex.markers=ex_markers,
                            top_percent = 50,
                            p_threshold=0.001)

print(length(gpt.norm))

#check length of trails
ll=unlist(lapply(gpt.norm,length))
print(ll)

#check seurat clusters for each trail
path.clus=unlist(lapply(gpt.norm,function(x)check_path_cluster(x,meta.cts=meta.cts)))
table(path.clus)

#get spot level data for each trail
gpt.norm.spot=lapply(gpt.norm,function(x) do.call(c,sub.clusters[x]))

```

## Step 5: Visualize Trails
```r
plot_all_shortest_paths(all.paths=gpt.norm,loc=loc.cts,length.cut = 4)
points(loc.cts[which(rownames(loc.cts) %in% rownames(meta.cts[which(meta.cts$seurat_clusters %in% tumor.clus),])),],col=alpha('red',0.1),cex=0.2)
```

## Step 6: Systematic analysis (optional)


#### 6.1 Get control trail sets
```r
perm.pathset.norm.all=get_control_paths(paths=gpt.norm,nbr.list=nbr.list,
                                        meta.cts=meta.cts,s.range1=max.size,s.range2=min.size,
                                        n.paths=30000,
                                        loc.cts=loc.cts,
                                        n.set=5000,
                                        path.size=path_size)
```

#### 6.2 Calculate the mean expression of each gene in each control set and the real migration trail set
```r

#10000 control sets
controlset.mean.allgg=controlset.mean.fun(genes=all.gg,sp.dd=sp.counts.norm)
gpt.mat.norm=get_gpt_gg1(dd=sp.counts.norm,pathlist=gpt.norm,gg.markers=all.gg) 
pathset.mean.all=apply(gpt.mat.norm,2,mean) #real pathset: 1*2957

gg.all=vocano_plot1(pathset.mean=pathset.mean.all,
                    controlset.means=controlset.mean.allgg,plot=F,
                    pvalue.cut=0.05,fc.cut.up=1.5,fc.cut.down=1/1.5)
```


#### 6.3 GSEA
```r
gg.all.enrich=enricher(gene=gg.all[[1]],TERM2GENE = c5.go[,c('gs_name','gene_symbol')])
enrich_barplot(gg.all.enrich,1:40)
```



