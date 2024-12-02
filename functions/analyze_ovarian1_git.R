

source("functions/get_paths_git.R")
source("functions/path_downstream_git.R")
source("functions/packaged_functions_git.R")
source("functions/figures_tables_git.R")



#############################Get the trails############################################################################################

#1. CD4, CD8 T cell stage markers
load('data/cd4_8_cluster_markers.Rdata') 
all.clus.markers=unique(clus.markers[,1])


#2. Markers for different cancers
load('data/cancer_markers.Rdata') 

#3. Data for GSEA
c5=msigdbr(species = 'Homo sapiens',category = 'C5')
c5.go=c5[which(c5$gs_subcat!='HPO'),]



#4. Prepare spatial data
dat1=data_preparation(folder="data/visium_ovarian",  
                      filenm = "Visium_FFPE_Human_Ovarian_Cancer_filtered_feature_bc_matrix.h5", 
                      res=0.2,
                      loc_file='data/visium_ovarian/spatial/tissue_positions_list.csv',
                      ex.list=c('PDCD1','GZMB','LAG3','HAVCR2','TIGIT','CD37','CTLA4'),
                      cell=c('CD3D','CD3E','CD3G','CD4','CD8A','CD8B'),   #Use 'CD8B2' for Spatial Gene Expression dataset analyzed using Space Ranger 2.0.0
                      cell_cut=5,
                      cut.t=1e-10,
                      n.vgg=3000,
                      clustering=T)


sp.counts.r=dat1[[1]]
sp.counts.norm=dat1[[2]]
vgg=dat1[[3]]
loc.raw=dat1[[4]]
meta_dat=dat1[[5]]
keep.ids=dat1[[6]]
meta.cts=dat1[[7]]
loc.cts=dat1[[8]]
pt.cts=dat1[[9]]
r.unit=dat1[[10]]
sub.clusters=dat1[[11]]
sp.counts.tpm=dat1[[12]]


#5. Take a look at Seurat clusters
#Define your own tumor region of clusters that have high expression of cancer markers

tmp.mean=plot_seurat_clusters(loc=loc.raw,meta_dat=meta_dat,sp.counts.norm=sp.counts.norm,sp.counts.r=sp.counts.r,cancer.markers=cancer.markers$ovarian)
tumor.clus=order(tmp.mean[[2]][,1],decreasing=T)[1:3]-1 #potential clusters that are tumor


#6. Take a look at the tumor region
#Red spots are tumor region
#Green spots are spots with T cells
#This method only works on sample with high-level of T cell infiltration

loc=loc.raw
plot_tumor(loc=loc.raw,meta_dat=meta_dat,tumor_cluster=tumor.clus)
points(loc[keep.ids,'x'],loc[keep.ids,'y'],col=alpha('green',0.4),pch=19,cex=.5)

#7. Check exhaustion score-need to have a wide range
plot_pt(loc=loc.raw,meta_dat=meta_dat)
boxplot_pt(loc=loc.raw,meta_dat=meta_dat)



#8. Get potential trails
Sys.time()
potential.paths=get_all_paths(loc.cts=loc.cts, #rescaled XY coordinate of node data
                              pt.cts=pt.cts, # Exhaustion score of nodes
                              cut.xydist=r.unit*1.3, # maximum XY distance of adjacent node on a path that is allowed
                              cut.sdist.n=6, # (n-1) thousandth of the pairwise xy distance-threshold for nodes to 3-D line 
                              pt.dist.cuts.n=11, # (n-1) percent of the pairwise exhaustion score differences-threshold for maximum allowed drawback along the trail
                              min.length=5,
                              r.unit=r.unit,
                              nr.graph=1.3)

Sys.time() #It took about 40 mins




#9. Select trails with large correlation between consecutive pairs of nodes
gpt.norm=path_selection(paths=potential.paths,
                        meta_dat=meta_dat,
                        meta.cts=meta.cts,
                        loc.cts=loc.cts,
                        pt.cts=pt.cts,
                        r.unit=r.unit,
                        tcell_spots=keep.ids,
                        sub.clusters=sub.clusters,
                        sp.counts.norm=sp.counts.norm,
                        tumor.clus=tumor.clus,
                        cor.cut=0.2,
                        cell=c('CD3D','CD3E','CD3G','CD4','CD8A','CD8B'), #Use 'CD8B2' for Spatial Gene Expression dataset analyzed using Space Ranger 2.0.0
                        ex.markers=c('PDCD1','GZMB','LAG3','HAVCR2','TIGIT','CD37','CTLA4'))





#10. Take a look of the trails


par(mfrow=c(1,1))
plot_all_shortest_paths(all.paths=gpt.norm,loc=loc.cts,length.cut = 5)

### If a reasonable number of trails can be identified using your sample, then continue with the systematic analysis













#############################Systematic analysis############################################################################################


### Check which clsuters these trails belong to
path.clus=unlist(lapply(gpt.norm,function(x)check_path_cluster(x,meta.cts=meta.cts)))
table(path.clus)

### 1. Get 10000 controlsets
perm.pathset.norm.all=get_control_paths(paths=gpt.norm,meta.cts=meta.cts,nr=1.3,s.range1=7*r.unit,s.range2=3*r.unit,n.paths=90000,loc.cts=loc.cts,r.unit=r.unit,n.set=10000)

### 2. Calculate path means of each gene for the 10000 control set

#Path mean of the 10000 control sets-28*2957 of length 10000
all.means.vgg=c()
for (i in 1:length(perm.pathset.norm.all)){
  pathset=perm.pathset.norm.all[[i]]
  pathset.spots=lapply(pathset,function(x) unlist(sub.clusters[x]))
  pathset.genes=lapply(pathset.spots,function(x) sp.counts.tpm[x,vgg])
  pathset.mean=do.call(rbind,lapply(pathset.genes,function(x) apply(x,2,mean)))
  all.means.vgg=c(all.means.vgg,list(pathset.mean))
}

#Mean for each gene for each set
controlset.means=do.call(rbind,lapply(all.means.vgg,function(x) apply(x, 2, mean))) #10000 control sets: 10000*2957


#Mean for each gene of each path (path1-path28) among 10000 random pathsets, 28*2957
all.path.mean=c() 
for (i in 1:length(gpt.norm)) {
  path.mat=do.call(rbind,lapply(all.means.vgg,function(x)x[i,])) #Gene expression of Path i in the 10000 control sets, 10000*2957
  path.mean=apply(path.mat,2,mean) #Mean expression of 2957 genes of path i
  all.path.mean=rbind(all.path.mean,path.mean)
}

rm(all.means.vgg)

### 3. Path mean of the real set-28*2957
gpt.mat.norm=get_gpt_gg1(dd=sp.counts.tpm,pathlist=gpt.norm,gg.markers=vgg) 
pathset.mean.all=apply(gpt.mat.norm,2,mean) #real pathset: 1*2957

### 4. Upregulated genes along the migration trails and GSEA
gg.all=vocano_plot1(pathset.mean=pathset.mean.all,controlset.means=controlset.means,plot=TRUE)
gg.all.enrich=enricher(gene=gg.all[[1]],TERM2GENE = c5.go[,c('gs_name','gene_symbol')])
enrich_barplot(gg.all.enrich,1:50)





### 5. Phenotype analysis

#Download single cell data for T cells from lung cancer patients (GSE99254)
# Prepare the data with Tcells_lung_cancer.R and get lung_seurat

norm.rna.mat=GetAssayData(lung_seurat, slot = "data")
clusters=Idents(lung_seurat)
table(clusters) #0:12
n.cluster=length(table(clusters))-1

### Tsne plot
DimPlot(lung_seurat, reduction = "tsne")
cluster_enrich(gg.up.pathset=gg.all[[1]])




### 6. Divide the trails into 3 groups by direction related to HNR

direction=rep(0,length(gpt.norm))
direction[c(1,11)]='out'
direction[c(2,4,8)]='out'
direction[3]='towards'
direction[5]='out'
direction[6]='out'
direction[7]='out'
direction[9]='out'
direction[c(10,14)]='none'
direction[12]='none'
direction[c(13,24)]='towards'
direction[15]='out'
direction[16]='out'
direction[17]='out'
direction[18]='out'
direction[19]='out'
direction[c(20,23)]='towards'
direction[c(21,25)]='none'
direction[c(22)]='none'
direction[c(26,27)]='out'
direction[c(28)]='towards'

names(direction)=paste0('Trail',1:length(direction))


out.mean=apply(gpt.mat.norm[direction=='out',],2,mean)
towards.mean=apply(gpt.mat.norm[direction=='towards',],2,mean)
unrelated.mean=apply(gpt.mat.norm[direction=='none',],2,mean)

gg.out=vocano_plot1(pathset.mean=out.mean,controlset.means=controlset.means,plot=TRUE)
gg.towards=vocano_plot1(pathset.mean=towards.mean,controlset.means=controlset.means,plot=TRUE)
gg.unrelated=vocano_plot1(pathset.mean=unrelated.mean,controlset.means=controlset.means,plot=TRUE)

gg.out.enrich=enricher(gene=gg.out[[1]],TERM2GENE = c5.go[,c('gs_name','gene_symbol')])
gg.towards.enrich=enricher(gene=gg.towards[[1]],TERM2GENE = c5.go[,c('gs_name','gene_symbol')])
gg.unrelated.enrich=enricher(gene=gg.unrelated[[1]],TERM2GENE = c5.go[,c('gs_name','gene_symbol')])

enrich_barplot(gg.out.enrich,1:50)
enrich_barplot(gg.towards.enrich,1:50)
enrich_barplot(gg.unrelated.enrich,1:50)

cluster_enrich(gg.up.pathset=gg.out[[1]])
cluster_enrich(gg.up.pathset=gg.towards[[1]])
cluster_enrich(gg.up.pathset=gg.unrelated[[1]])
