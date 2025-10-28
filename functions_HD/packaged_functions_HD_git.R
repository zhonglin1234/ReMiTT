
options(future.globals.maxSize = 5 * 1024^3)  # 5 GB

require(Seurat)
require(arrow)
require(ggplot2)
require(dplyr)
require(monocle)
require(igraph)
library(plotly)
library(RColorBrewer)
require(openxlsx)
require(LearnGeom)
require(StereoMorph)
require(openxlsx)
require(msigdbr)
library(dbscan)
library(Matrix)

###Colors
colfunc<-colorRampPalette(c("royalblue","springgreen","yellow",'red'))


###1. Load and prepare data

read.spatial.dat<-function(folder,filenm,resolut=0.3,loc_file,n.vgg=3000){ # data can be spot or HD
  spdat <- Load10X_Spatial(
    data.dir      =  folder,
    filename      = filenm,  # use filtered counts
    assay         = "Spatial",
    slice         = "HD_16um",
    filter.matrix = TRUE
  )
  spdat$total_UMI <- Matrix::colSums(spdat@assays$Spatial@layers$counts)
  spdat <- subset(spdat, subset = total_UMI >= 50)
  spdat <- NormalizeData(spdat, normalization.method = "LogNormalize", scale.factor = 1e4)
  spdat <- FindVariableFeatures(spdat, selection.method = "vst", nfeatures = 3000)
  
  
  all.genes=rownames(spdat@assays$Spatial$counts)
  vgg <- VariableFeatures(spdat)
  genes=unique(c(vgg,"PDCD1","LAG3","HAVCR2","TIGIT","ENTPD1","CTLA4","TOX",'CD3D','CD3E','CD3G','CD4','CD8A','CD8B2',
                 'CXCL9','CXCL10','CXCL11','CCL4','CCL5','ITGA8','ITGA11','TJP3', 'LGALS3','MMP1','MMP2'))
  genes=intersect(all.genes,genes)
  
  tcells=apply(spdat@assays$Spatial$counts[c('CD3D','CD3E','CD3G','CD4','CD8A','CD8B2'),],2,function(x) max(x,na.rm=T))
  spdat@assays$Spatial$count=spdat@assays$Spatial$counts[genes,names(tcells)[which(tcells>0)]]
  spdat@assays$Spatial$data=spdat@assays$Spatial$data[genes,names(tcells)[which(tcells>0)]]
  
  #get clusters
  projection=read.csv(paste0(folder,'/analysis/umap/gene_expression_2_components/projection.csv'))
  set.seed(123)  
  k <- 5         
  km <- kmeans(projection[, c("UMAP.1","UMAP.2")], centers = k)
  projection$cluster <- km$cluster-1
  rownames(projection)=projection$Barcode
  
  meta_dat=spdat@meta.data
  meta_dat=meta_dat[names(tcells)[which(tcells>0)],]
  projection=projection[rownames(meta_dat),]
  meta_dat$seurat_clusters=projection$cluster
  
  loc.raw=GetTissueCoordinates(spdat)
  loc.raw=loc.raw[which(tcells>0),]
  
  #plot(loc.raw$x,loc.raw$y)
  return(c(list(as.matrix(t(spdat@assays$Spatial$count))),
           list(as.matrix(t(spdat@assays$Spatial$data))),
           list(vgg),
           list(meta_dat),
           list(loc.raw))
  )
}






### sum of cell markers and sum of exhaustion markers
add_sum_to_meta<-function(ex.list,cell,sp.counts.r,meta_dat){
  ex.sum=rowSums(sp.counts.r[,ex.list])
  cell.sum=rowSums(sp.counts.r[,cell])
  
  print(paste('Sum of ',paste(cell,collapse = ',')))
  print(quantile(cell.sum,probs=seq(0,1,length=21)))
  
  print(paste('Sum of ',paste(ex.list,collapse = ',')))
  print(quantile(round(ex.sum),probs=seq(0,1,length=201)))
  
  names(cell.sum)=rownames(sp.counts.r)
  names(ex.sum)=rownames(sp.counts.r)
  
  meta_dat$cell.sum=cell.sum[rownames(meta_dat)]
  meta_dat$ex.sum=ex.sum[rownames(meta_dat)]
  
  return(meta_dat)
}


## Exhaustion score
add_pt_to_meta<-function(keep.ids,meta_dat){
  meta_dat_cell=meta_dat[keep.ids,]
  tmp=lm(ex.sum~cell.sum,data=meta_dat_cell)
  ex.score=tmp$residuals
  ex.score=ex.score+abs(min(ex.score))
  meta_dat$ex.score=NA
  meta_dat[keep.ids,'ex.score']=ex.score
  print(head(meta_dat))
  return(meta_dat)
}


# Assign subclusters

# meta_dat: data.frame with columns x, y, ex.score (and optionally an ID column)
# cut.t    : spatial radius (same units as x,y) to gather candidate neighbors
# pt_thresh: threshold on |ex.score_i - ex.score_j| relative to the SEED point
# id_col   : optional column with unique spot IDs; if NULL, uses rownames
# min_size : minimum cluster size to keep; smaller groups become singletons (or NA if drop_singletons=TRUE)
# order_by : how to order seeds ("index" | "degree_desc")
assign_subcluster_greedy <- function(meta_dat,
                                     cut.t,
                                     pt_thresh,
                                     id_col = NULL,
                                     min_size = 1,
                                     drop_singletons = FALSE,
                                     order_by = c("index", "degree_desc")) {
  stopifnot(all(c("x","y","ex.score") %in% colnames(meta_dat)))
  order_by <- match.arg(order_by)
  
  # IDs
  if (is.null(id_col)) {
    if (is.null(rownames(meta_dat))) stop("Provide id_col or set rownames(meta_dat) to spot IDs.")
    ids <- rownames(meta_dat)
  } else {
    ids <- as.character(meta_dat[[id_col]])
  }
  n <- nrow(meta_dat)
  
  # coords & pt
  coords <- as.matrix(meta_dat[, c("x","y")])
  rownames(coords) <- ids
  pt <- setNames(as.numeric(meta_dat$ex.score), ids)
  
  # radius neighbors (includes self)
  nbr <- frNN(coords, eps = cut.t, sort = FALSE)
  
  # optional: process high-degree seeds first (often larger, more stable clusters)
  seed_order <- seq_len(n)
  if (order_by == "degree_desc") {
    deg <- vapply(nbr$id, length, integer(1))
    seed_order <- order(deg, decreasing = TRUE)
  }
  
  cluster_id <- rep(NA_integer_, n)
  names(cluster_id) <- ids
  visited <- rep(FALSE, n)
  cid <- 0L
  
  for (k in seed_order) {
    if (visited[k]) next
    
    seed_id <- ids[k]
    nb_idx  <- nbr$id[[k]]
    nb_ids  <- ids[nb_idx]
    
    # keep neighbors whose pt is close to the SEED's pt
    keep_ids <- nb_ids[ abs(pt[seed_id] - pt[nb_ids]) < pt_thresh ]
    
    # ensure at least the seed
    if (length(keep_ids) == 0) keep_ids <- seed_id
    
    # enforce min_size (optionally drop singletons)
    if (length(keep_ids) < min_size) {
      if (drop_singletons) {
        visited[match(seed_id, ids)] <- TRUE
        next
      } else {
        keep_ids <- seed_id
      }
    }
    
    cid <- cid + 1L
    cluster_id[keep_ids] <- cid
    visited[match(keep_ids, ids)] <- TRUE
  }
  
  # build cluster list
  clusters <- split(names(cluster_id)[!is.na(cluster_id)], cluster_id[!is.na(cluster_id)])
  
  # small report
  size_tbl <- table(vapply(clusters, length, integer(1)))
  message("Formed ", length(clusters), " clusters; size distribution:\n", paste(capture.output(print(size_tbl)), collapse = "\n"))
  
  list(
    cluster_id = data.frame(spot_id = ids, cluster_id = cluster_id, stringsAsFactors = FALSE),
    clusters   = clusters
  )
}




#make meta data for the data by subclusters
meta_subclusters<-function(sub.clusters,loc,pt,meta_dat){
  loc.cts=do.call(rbind,lapply(sub.clusters,function(x) colSums(loc[x,])/length(x)))
  rownames(loc.cts)=paste0('cts',1:nrow(loc.cts))
  colnames(loc.cts)=c('x','y')
  pt.cts=unlist(lapply(sub.clusters,function(x) mean(pt[x])))
  names(pt.cts)=rownames(loc.cts)
  
  tmp.cluster=lapply(sub.clusters,function(x) meta_dat[x,'seurat_clusters'])
  tmp.cluster=lapply(tmp.cluster,function(x) prop.table(table(x)))
  tmp.cluster=unlist(lapply(tmp.cluster,function(x) names(x)[which(x==max(x))][1]))
  cts.meta=data.frame(loc.cts,pt.cts,tmp.cluster)
  colnames(cts.meta)=c('x','y','pt','seurat_clusters')
  return(cts.meta)
}



data_preparation_HD<-function(folder,
                              filenm,
                              loc_file,
                              res,
                              ex.list,
                              cell,
                              cell_cut, # cutoff of sum of T cells markers to define T cell spot
                              n.vgg) { #cutoff of exhaustion score difference to cluster spots into one single nodes
  
  ###Load data
  seurat_items=read.spatial.dat(folder=folder,  filenm = filenm, resolut=res,
                                loc_file=loc_file,n.vgg=n.vgg)
  
  sp.counts.r=seurat_items[[1]]
  sp.counts.norm=seurat_items[[2]]
  meta_dat=seurat_items[[4]]
  loc.raw=seurat_items[[5]]
  vgg=seurat_items[[3]]
  
  ex.list=intersect(ex.list,colnames(sp.counts.r))
  print(ex.list)
  cell=intersect(cell,colnames(sp.counts.r))
  print(cell)
  
  sp.counts.tpm=apply(sp.counts.r,2,function(x) as.numeric(x)/meta_dat[rownames(sp.counts.r),'nCount_Spatial']*1e6)
  rownames(sp.counts.tpm)=rownames(sp.counts.r)
  
  ### Add exhaustion score to meta data
  #Add sum of cell markers and sum of exhaustion scores
  meta_dat=add_sum_to_meta(ex.list=ex.list,cell=cell,sp.counts.r=sp.counts.r,meta_dat=meta_dat)
  meta_dat=add_pt_to_meta(keep.ids=rownames(meta_dat),meta_dat=meta_dat)
  pt=meta_dat$ex.score
  
  output=c(list(sp.counts.r),list(sp.counts.norm),list(sp.counts.tpm),list(vgg),list(loc.raw),list(meta_dat))
  
  names(output)=c('raw count data','normalized count data','TPM data','Top variable genes','XY coordinate of the raw data',
                  'meta data of the spot data')
  return(output)
}





###2.  Plot Seurat clusters
plot_seurat_clusters<-function(loc,meta_dat,sp.counts.norm,sp.counts.r,cancer.markers){
  loc=loc[rownames(meta_dat),]
  loc[,'x']=as.numeric(loc[,'x'])
  loc[,'y']=as.numeric(loc[,'y'])
  tmp.clusters=meta_dat$seurat_clusters
  number.clusters=length(table(tmp.clusters))
  col.list=c('blue','green','yellow','orange','brown','red')
  plot(loc,col='white',xlim=c(min(loc[,'x'])*0.9,max(loc[,'x'])*1.2),main='Seurat clusters',cex=0.1)
  for (i in 1:number.clusters) {
    points(loc[tmp.clusters==(i-1),'x'],loc[tmp.clusters==(i-1),'y'],col=col.list[i],pch=19,cex=.1)
  }
  legend(x=max(loc[,'x'])*1.05,y=median(loc[,'y']),legend=0:(number.clusters-1),col=col.list[1:number.clusters],pch=19,cex=0.7,bty='n')
  
  tmp=clustermeans(count_dat=sp.counts.norm,meta_dat=meta_dat,cancer.markers=cancer.markers)
  print('Median expression of each cancer markers in normalized data:')
  print(tmp[[1]])
  print('Sum of median expression of each cancer markers in each cluster with normalized data:')
  print(rowSums(tmp[[1]]))
  
  tmp=clustermeans(count_dat=sp.counts.r,meta_dat=meta_dat,cancer.markers=cancer.markers)
  print('Median expression of each cancer markers in raw count data:')
  print(tmp[[1]])
  print('Sum of median expression of each cancer markers in each cluster with raw count data:')
  print(rowSums(tmp[[1]]))
  return(tmp)
}



###3. Plot tumor area
plot_tumor<-function(loc,meta_dat,tumor_cluster){
  loc=loc[rownames(meta_dat),]
  tmp.clusters=meta_dat$seurat_clusters
  tumor.ind=rownames(meta_dat[which(meta_dat$seurat_clusters %in% tumor_cluster),])
  plot(loc,col='grey',xlim=c(min(loc[,'x'])*0.9,max(loc[,'x'])*1.2),main='Tumor region (red)')
  points(loc[tumor.ind,'x'],loc[tumor.ind,'y'],col=alpha('red',0.4),pch=19)
}


### 4. Plot T cell spots
plot_Tcell_spots<-function(keep.ids,loc){
  print(paste('Number of T cell spots:',length(keep.ids),'out of',nrow(loc),'spots'))
  plot(loc,col='grey',xlim=c(min(loc[,'x'])*0.9,max(loc[,'x'])*1.2),main='Tcell spots (pink)')
  points(loc[keep.ids,'x'],loc[keep.ids,'y'],col=alpha('pink',0.4),pch=19,cex=.5)
}



### 5. Plot exhaustion score
plot_pt<-function(loc,meta_dat){
  meta_dat_cell=meta_dat[which(!is.na(meta_dat$ex.score)),]
  loc=loc[rownames(meta_dat),]
  loc=loc[rownames(meta_dat_cell),]
  ex.score=meta_dat_cell$ex.score
  print(quantile(ex.score,probs=seq(0,1,length=21)))
  ggplot(data.frame(loc), aes(x=loc[,1], y=loc[,2], color = ex.score)) +
    geom_point(size=0.1) +
    scale_color_gradient(low = "blue", high = "red") +  labs(title = "Exhaustion score")+
    labs(x='x',y='y')
}


### 6. Boxplot of exhaustion score in each cluster
boxplot_pt<-function(loc,meta_dat){
  meta_dat_cell=meta_dat[which(!is.na(meta_dat$ex.score)),]
  loc=loc[rownames(meta_dat),]
  loc=loc[rownames(meta_dat_cell),]
  ex.score=meta_dat_cell$ex.score
  boxplot(ex.score~meta_dat_cell$seurat_clusters,xlab='Seurat clusters',ylab='Exhaustion score',main='Exhaustion score by Seurat clusters')
}


### 7. To find tumor clusters
clustermeans<-function(count_dat,meta_dat,cancer.markers){
  tmp.markers=cancer.markers[which(cancer.markers %in% colnames(count_dat))]
  tmp.dat=count_dat[,tmp.markers]
  meta_dat=meta_dat[rownames(tmp.dat),]
  
  #print("Median")
  med=apply(tmp.dat,2,function(x) aggregate(x~meta_dat$seurat_clusters,data=tmp.dat,median)[,2])
  #print(med)
  
  #print("Mean")
  mm=apply(tmp.dat,2,function(x) aggregate(x~meta_dat$seurat_clusters,data=tmp.dat,mean)[,2])
  #print(mm)
  return(c(list(med),list(mm)))
}









##################################################################################################################################
#Get paths by minimum spanning tree
###################################################################################################################################


# Get edge weight (difference in Exhaustion score)
edge_wt<-function(spot,nbr.list,pt) {
  nbs=rownames(nbr.list[spot][[1]])
  if (length(nbs)>0) {
    tmp.edges=as.matrix(cbind(rep(spot,length(nbs)),nbs))
    drct=pt[nbs]-pt[spot]
    wt1=abs(drct) #difference in pt
    dist.mat=as.matrix(dist(loc.cts[c(spot,nbs),]))
    wt2=dist.mat[spot,nbs] #location distance
    
    #swap if pt1 < pt2
    swch=which(drct<0)
    if (length(swch)>0) {
      col1= tmp.edges[swch,2]
      col2= tmp.edges[swch,1]
      tmp.edges[swch,1]=col1
      tmp.edges[swch,2]=col2
    }
    out=cbind(tmp.edges,wt1,wt2)
    rownames(out)=paste(spot,1:nrow(out),sep='_')
    colnames(out)=c('edge0','edge1','dist_pt','dist_s')
    return(out)
  }
}



# Make the graph
make_graph<-function(spots,nbr.list,pt,ws,wp){
  g.dat=do.call(rbind,sapply(spots,function(x) edge_wt(spot=x,nbr.list,pt=pt)))
  g.dat=unique(g.dat)
  rownames(g.dat)=paste(g.dat[,'edge0'],g.dat[,'edge1'],sep='_')
  
  edges <- data.frame(from = g.dat[,'edge0'],to = g.dat[,'edge1'])
  weights <- as.numeric(g.dat[,'dist_s'])*ws+as.numeric(g.dat[,'dist_pt'])*wp   #Weight is a combination of space distance and pseudotime distance
  
  graph <- graph_from_data_frame(edges, directed = TRUE)
  E(graph)$weight <- weights
  return(c(list(graph),list(g.dat)))
}



# Get MSP 
get_msp1<-function(start.point,
                   spots,
                   pt,
                   graph,
                   cut.length, # mininum distance between start and end point
                   min.length){ # minimum number of nodes on the trail
  if (!start.point%in%names(V(graph))) {
    print('Start point not in the graph')
    return(NULL)
  }
  end.points=spots[which(pt[spots]-pt[start.point]>quantile(pt,probs = seq(0,1,length=11))[9])] #has to have bigger pt than start point for at least 80% 
  end.points=end.points[which(end.points %in% names(V(graph)))] #has to be in graph
  space.distance=sqrt((loc.cts[end.points,1]-loc.cts[start.point,][1])^2+(loc.cts[end.points,2]-loc.cts[start.point,][2])^2)
  end.points=end.points[which(space.distance>=cut.length)] #cannot be too close to start point
  
  if (length(end.points)==0) {return(NULL)}
  
  path.list=c()
  for (end.point in end.points) {
    msp <- shortest_paths(graph, from = start.point, to = end.point, weights = E(graph)$weight)
    vv=unlist(msp$vpath)
    if (length(vv)>=min.length) {path.list=c(path.list,list(names(vv)))}
  }
  return(path.list)
}





get_all_paths_HD <-function(loc.cts, #rescaled XY coordinate of node data
                            pt.cts, # Exhaustion score of nodes
                            nbr.list,
                            cut.length=800,
                            min.length=8){ # Edge can exist between two nodes with at most nr.graph*r.unit distance
  loc.cts=loc.cts[names(pt.cts),]
  keep.cts=names(pt.cts)
  
  ### Set up thresholds for start and end points
  
  cut.low=quantile(pt.cts,probs=seq(0,1,length=5))[2]
  cut.high=quantile(pt.cts,probs=seq(0,1,length=5))[4]
  
  print(paste('low threshold of exhaustion score for start points:',cut.low))
  print(paste('high threshold of exhaustion score for end points:',cut.high))
  
  start.pts=names(pt.cts)[which(pt.cts<=cut.low)] #First 1/4 most naive nodes as start nodes
  end.pts=names(pt.cts)[which(pt.cts>=cut.high)] #last 1/4 most exhaustive nodes as end nodes
  print(paste('Number of start nodes:',length(start.pts))) 
  print(paste('Number of end nodes:',length(end.pts))) 
  
  #Plot start and end nodes:
  par(mfrow=c(1,1))
  plot(loc.cts[start.pts,],main='Start(black) and end nodes (red)')
  points(loc.cts[end.pts,'x'],loc.cts[end.pts,'y'],col='red')
  
  
  ###Get MS path
  
  #Make the graph with distance being weight
  tmp.graph=make_graph(spots=keep.cts,nbr.list=nbr.list,pt=pt.cts,ws=1,wp=0)
  graph2=tmp.graph[[1]]
  g.dat2=tmp.graph[[2]]
  
  #Make the graph with difference in exhaustion score as weight
  tmp.graph=make_graph(spots=keep.cts,nbr.list=nbr.list,pt=pt.cts,ws=0,wp=1)
  graph3=tmp.graph[[1]]
  g.dat3=tmp.graph[[2]]
  
  #Get paths
  paths_msp_spacew=c()
  paths_msp_exscorew=c()
  n=0
  for (i in start.pts) {
    n=n+1
    #print(n)
    tmp.out2=get_msp1(start.point=i,spots=keep.cts,pt=pt.cts,
                      graph=graph2,
                      cut.length=cut.length,
                      min.length=min.length)
    if (length(tmp.out2)>0) {paths_msp_spacew=c(paths_msp_spacew,list(tmp.out2))}
    
    tmp.out3=get_msp1(start.point=i,spots=keep.cts,pt=pt.cts,
                      graph=graph3,
                      cut.length=cut.length,
                      min.length=min.length)
    if (length(tmp.out3)>0) {paths_msp_exscorew=c(paths_msp_exscorew,list(tmp.out3))}
  }
  
  ll=unlist(lapply(paths_msp_spacew,function(x) length(x[[1]])))
  print('The lengh of paths from shortest path, space as weight:')
  print(table(ll))
  
  ll=unlist(lapply(paths_msp_exscorew,function(x) length(x[[1]])))
  print('The lengh of paths from shortest path, exhaustion score as weight:')
  print(table(ll))
  
  output=c(list(paths_msp_spacew),list(paths_msp_exscorew))
  names(output)=c('paths from shortest path with distance being weight','paths from shortest path with difference in exhaustion score being weight')
  return(output)
}








##################################################################################################################################
#Select paths based on correlation of t cell stage markers between consecutive spots along the trail.
###################################################################################################################################


# select genes for comparing correlation (between consecutive nodes)
cell.markers3<-function(meta_dat,cell,dat,p_threshold,tumor.clus){
  sc=meta_dat[rownames(dat),'seurat_clusters']
  dat=dat[which(sc %in% tumor.clus),]
  markerlevel=rowSums(dat[,cell])
  names(markerlevel)=rownames(dat)
  
  tmp.cor=apply(dat,2,function(x) cor.test(x,markerlevel,method='s')[c(3,4)])
  pt=unlist(lapply(tmp.cor,function(x)x[[1]]))
  tmp.corr=unlist(lapply(tmp.cor,function(x)x[[2]]))
  genes=colnames(dat)[which(pt<p_threshold)]
  names(tmp.corr)=colnames(dat)
  tmp.corr=tmp.corr[genes]
  tmp.corr=tmp.corr[order(tmp.corr,decreasing = T)]
  return(tmp.corr)
}


edge_wt_perm<-function(spot,nbs.list) {
  nbs=rownames(nbr.list[spot][[1]])
  if (length(nbs)>0) {
    tmp.edges=as.matrix(cbind(rep(spot,length(nbs)),nbs))
    dist.mat=as.matrix(dist(loc.cts[c(spot,nbs),]))
    wt=dist.mat[spot,nbs] #location distance
    out=cbind(tmp.edges,wt)
    rownames(out)=paste(spot,1:nrow(out),sep='_')
    colnames(out)=c('edge0','edge1','dist_s')
    return(out)
  }
}

make_graph_perm<-function(spots,nbs.list){ #Undirected graph with weight being spatial distance
  g.dat=do.call(rbind,sapply(spots,function(x) edge_wt_perm(spot=x,nbs.list=nbs.list)))
  g.dat=unique(g.dat)
  rownames(g.dat)=paste(g.dat[,'edge0'],g.dat[,'edge1'],sep='_')
  
  edges <- data.frame(from = g.dat[,'edge0'],to = g.dat[,'edge1'])
  weights <- as.numeric(g.dat[,'dist_s'])  
  
  graph <- graph_from_data_frame(edges, directed = FALSE)
  E(graph)$weight <- weights
  return(c(list(graph),list(g.dat)))
}



#Get alternative paths

get_short_paths<-function(start.p,end.p,
                          tmp.path,
                          n.path,
                          g.dat.perm,
                          nbr.list.large
){
  nb.pps=unique(do.call(c,nbr.list.large[tmp.path]))
  path.list=list(tmp.path)
  
  for (i in 1:n.path) {
    rm.nodes=unique(do.call(c,path.list))
    rm.nodes=rm.nodes[which(!rm.nodes %in% c(start.p,end.p))]
    #print(i)
    #print(rm.nodes)
    
    tmp.gdat=data.frame(g.dat.perm[which(g.dat.perm[,'edge0'] %in% nb.pps | g.dat.perm[,'edge1'] %in% nb.pps),])
    rm.rows=which(tmp.gdat$edge0 %in% rm.nodes | tmp.gdat$edge1 %in% rm.nodes)
    if (length(rm.rows)>0) {tmp.gdat=tmp.gdat[-rm.rows,]}
    
    tmp.edges <- data.frame(from = tmp.gdat[,'edge0'],to = tmp.gdat[,'edge1'])
    tmp.weights <- as.numeric(tmp.gdat[,'dist_s'])  
    tmp.graph <- graph_from_data_frame(tmp.edges, directed = FALSE)
    E(tmp.graph)$weight <- tmp.weights
    #print(dim(tmp.edges))
    if ((start.p %in% tmp.edges[,1]|start.p %in% tmp.edges[,2]) & (end.p %in% tmp.edges[,1]|end.p %in% tmp.edges[,2])) {
      msp <- shortest_paths(tmp.graph, from = start.p, to = end.p, weights = E(tmp.graph)$weight)
      vv=names(unlist(msp$vpath))
      if (length(vv)==0) {break}
      if (length(vv)>2)
        path.list=c(path.list,list(vv))
    }
    else {break}
  }
  path.list=path.list[-1] #remove the gradient path
  return(path.list)
}


get_other_paths2<-function(tmp.path,n.path,g.dat.perm,nbr.list.large){
  realstart=tmp.path[1]
  realend=tmp.path[length(tmp.path)]
  realpath=tmp.path
  
  path.list1=get_short_paths(start.p=realstart,end.p=realend,tmp.path=tmp.path,n.path=n.path,g.dat.perm=g.dat.perm,nbr.list.large=nbr.list.large)
  if (length(path.list1)==0) return(NULL)
  
  tmp.path=path.list1[[1]]
  path.list2=get_short_paths(start.p=tmp.path[2],end.p=tmp.path[length(tmp.path)-1],tmp.path=tmp.path,n.path=n.path,g.dat.perm=g.dat.perm,nbr.list.large=nbr.list.large)
  path.list2=lapply(path.list2,function(x) c(realstart,x,realend))
  rm.ind=unlist(lapply(path.list2,function(x) length(intersect(x,realpath[2:(length(realpath)-1)]))))
  path.list2=path.list2[-which(rm.ind>0)]
  
  
  if (length(path.list1)==1) return(c(path.list1,path.list2)) #there could be duplicates, or path that has nodes in real path
  
  tmp.path=path.list1[[2]]
  path.list3=get_short_paths(start.p=tmp.path[2],end.p=tmp.path[length(tmp.path)-1],tmp.path=tmp.path,n.path=n.path,g.dat.perm=g.dat.perm,nbr.list.large=nbr.list.large)
  path.list3=lapply(path.list3,function(x) c(realstart,x,realend))
  rm.ind=unlist(lapply(path.list3,function(x) length(intersect(x,realpath[2:(length(realpath)-1)]))))
  path.list3=path.list3[-which(rm.ind>0)]
  
  if (length(path.list1)==2) return(c(path.list1,path.list2,path.list3)) #there could be duplicates, or path that has nodes in real path
  
  tmp.path=path.list1[[3]]
  path.list4=get_short_paths(start.p=tmp.path[2],end.p=tmp.path[length(tmp.path)-1],tmp.path=tmp.path,n.path=n.path,g.dat.perm=g.dat.perm,nbr.list.large=nbr.list.large)
  path.list4=lapply(path.list4,function(x) c(realstart,x,realend))
  rm.ind=unlist(lapply(path.list4,function(x) length(intersect(x,realpath[2:(length(realpath)-1)]))))
  path.list4=path.list4[-which(rm.ind>0)]
  
  return(unique(c(path.list1,path.list2,path.list3,path.list4))) #there could be path that has nodes in real path
}




### 5. calculate consecutive pairwise correlation
pairwaise_correlation<-function(tmp.path,cor.mat,sub.clusters){
  #Get spots
  tmp.spots=do.call(c,sub.clusters[tmp.path])
  #Get consecutive pairwise correlation
  tmp.cor=cor.mat[tmp.spots,tmp.spots]
  pair.cor=apply(cbind(1:(length(tmp.path)-1),2:length(tmp.path)),1,function(x) tmp.cor[x[1],x[2]])
  return(c(list(pair.cor),list(mean(pair.cor))))
}


# (1) Get paths mean pairwise correlation of a list of path
get_mean_cor<-function(path.list,cor.mat,sub.clusters){
  path.cor.list=lapply(path.list,function(x) pairwaise_correlation(tmp.path=x,cor.mat=cor.mat,sub.clusters=sub.clusters))
  mean.cor.list=unlist(lapply(path.cor.list,function(x) x[[2]]))
  names(mean.cor.list)=names(path.list)
  mean.cor.list=mean.cor.list[order(mean.cor.list,decreasing = T)]
  return(mean.cor.list)
}

# (2) Get paths with largest correlation compared to its alternative paths

get_large_cor<-function(dat,path.list,alter.paths,pna,top_percent=25,sub.clusters,meta.cts){
  cor.mat2=cor(t(dat)) #Pearson's correlation between spots
  diag(cor.mat2)=NA
  
  gpt.cor2=get_mean_cor(path.list=path.list,cor.mat=cor.mat2,sub.clusters=sub.clusters)
  
  if (length(pna)>0) {
    print('pna correlations')
    print(gpt.cor2[pna])
  }
  
  alter.cor2=lapply(alter.paths,function(x) get_mean_cor(path.list=x,cor.mat=cor.mat2,sub.clusters=sub.clusters))
  alter.max2=unlist(lapply(alter.cor2,max))
  names(alter.max2)=names(alter.paths)
  
  if (length(pna)>0) {
    gpt.cor2a=gpt.cor2[which(!names(gpt.cor2) %in% names(pna))]
  }
  
  if (length(pna)==0) {
    gpt.cor2a=gpt.cor2
  }
  
  alter.max2=alter.max2[names(gpt.cor2a)]
  cut.trds=quantile(do.call(c,alter.cor2),seq(0,1,0.01),na.rm = T)
  print(cut.trds)
  
  gpt1=names(gpt.cor2a)[which(gpt.cor2a>=alter.max2)]
  print('Number of migration paths with largest mean pairwise correlation:')
  print(length(gpt1))
  #print('The mean pairwise correlations:')
  #gpt.cor2a[gpt1]
  
  
  gpt2=intersect(gpt1,names(gpt.cor2)[which(gpt.cor2>cut.trds[(100-top_percent)+1])]) # >75% percentile of all the adjcent alternative paths 
  print(paste0('Number of migration paths with mean correlation >',100-top_percent,' percentile of all the alternative paths:'))
  print(length(gpt2))
  
  print('Path mean correlations:')
  print(gpt.cor2[gpt2])
  
  print('Maximum alternative mean correlations')
  print(alter.max2[gpt2])
  
  print("Number of alternative paths for each migration path")
  n.alter=unlist(lapply(alter.paths[gpt2],length))
  print(table(n.alter))
  
  print('Migration paths that only have one alternative path:')
  tmp.pp=gpt2[which(n.alter==1)]
  print(tmp.pp)
  gpt2=gpt2[which(!gpt2 %in% tmp.pp)] #Temporarily remove them
  
  tmp.keep=tmp.pp[which(gpt.cor2[tmp.pp]>cut.trds[81])]
  print('The ones that has mean cor > 80th percentile of all alternative paths, keep them')
  print(tmp.keep)
  if (length(tmp.keep>0)) {gpt2=c(gpt2,tmp.keep)}
  
  if (length(pna)>0) {
    print('For those migration paths that have no alternative path, if its mean cor > 90th percentile of all alternative paths, keep them')
    tmp.keep.pna=pna[which(gpt.cor2[pna]>cut.trds[91])]
    print(tmp.keep.pna)
    if (length(tmp.keep.pna>0)) {gpt2=c(gpt2,tmp.keep.pna)}
  }
  gpt2.list=remove_overlaps(path.list[gpt2])
  print(length(gpt2.list))
  print(table(unlist(lapply(gpt2.list,length))))
  return(gpt2.list)
}

#Remove duplicated paths
n.overlap<-function(v1,v.list){
  unlist(lapply(v.list,function(x) length(intersect(v1,x))))
}

remove_overlaps<-function(path.list){
  rem.list=c()
  for (i in 1:length(path.list)) {
    path.ll=length(path.list[[i]])
    tmp.no=n.overlap(path.list[[i]],path.list[-i])
    if (max(tmp.no==path.ll)) {rem.list=c(rem.list,i)}
  }
  if (length(rem.list)>0) {path.list=path.list[-rem.list]}
  return(path.list)
}


path_selection<-function(paths,
                         meta_dat,
                         sp.counts.norm,
                         tumor.clus,
                         meta.cts,
                         loc.cts,
                         pt.cts,
                         r.unit,
                         tcell_spots,
                         sub.clusters,
                         cor.cut=0.001,
                         cell,
                         ex.markers,
                         top_percent=25,
                         p_threshold = 0.001
){
  
  loc.cts=loc.cts[names(pt.cts),]
  meta.cts=meta.cts[names(pt.cts),]
  keep.cts=names(pt.cts)
  #Space distance matrix among nodes
  #cts.sdist=as.matrix(dist(loc.cts))
  
  ### Load the potential paths
  paths_msp_spacew=paths[[1]]
  paths_msp_exscorew=paths[[2]]
  
  path.list=c(lapply(paths_msp_spacew,function(x) x[[1]]),lapply(paths_msp_exscorew,function(x) x[[1]]))
  print(length(path.list)) 
  names(path.list)=paste0('path',1:length(path.list))
  
  ll=unlist(lapply(path.list,length))
  print(table(ll))
  
  ### Get T cell stage related gene markers for correlation filtering
  cor.info=cell.markers3(meta_dat=meta_dat,cell=cell,dat=sp.counts.norm,p_threshold = p_threshold,tumor.clus = tumor.clus)
  markers=names(cor.info)[which(cor.info>cor.cut)]
  gg=unique(c(intersect(clus.markers[,1],markers),ex.markers))
  print(length(gg))
  print(gg)
  
  ### Get undirected graph, only nodes with cd3+cd4+cd8>=threshold were included
  tmp.graph=make_graph_perm(spots=rownames(meta.cts),nbs.list=nbs.list)
  graph.perm.all=tmp.graph[[1]]
  g.dat.perm.all=tmp.graph[[2]]
  
  ### Get alternative paths connecting start and end node of each path in adjacent areas
  alter.paths=c()
  for (i in 1:length(path.list)) {
    tmp.path=path.list[[i]]
    tmp.alter=get_other_paths2(tmp.path=tmp.path,n.path=6,g.dat.perm = g.dat.perm.all,nbr.list.large = nbr.list.large )
    alter.paths=c(alter.paths,list(tmp.alter))
  }
  names(alter.paths)=names(path.list)
  ll.alter=unlist(lapply(alter.paths,length))
  print(table(ll.alter))
  
  #paths with no alternative paths
  pna=which(ll.alter==0) 
  print(length(pna))
  if (length(pna)>0) {
    alter.paths=alter.paths[-pna]
  }
  print(length(alter.paths))
  ll.alter=unlist(lapply(alter.paths,length))
  
  path.spots=unique(do.call(c,sub.clusters[unique(do.call(c,path.list))]))
  alter.cts=unique(do.call(c,do.call(c,alter.paths)))
  alter.spots=unique(do.call(c,sub.clusters[alter.cts]))
  
  gpt.norm=get_large_cor(dat=sp.counts.norm[unique(c(path.spots,alter.spots)),gg],
                         path.list=path.list,
                         alter.paths=alter.paths,
                         pna=pna,
                         top_percent=top_percent,
                         sub.clusters = sub.clusters,
                         meta.cts=meta.cts)
  
  return(gpt.norm)
}










##################################################################################################################################
#Plot paths
###################################################################################################################################

# Plot a path on xy-plane
plot_2d_path<-function(tmp.path,loc=loc.cts) {
  plot(loc,cex=0.5)
  points(loc[tmp.path,'x'],loc[tmp.path,'y'],col=colfunc(length(tmp.path)),pch=19,cex=0.5)
}


# Plot output paths
plot_all_shortest_paths<-function(all.paths,loc,length.cut){
  plot(loc,col='lightgrey',cex=.7)
  for (k in 1:length(all.paths)) {
    tmp.path=all.paths[[k]]
    start.point=tmp.path[1]
    end.point=tmp.path[length(tmp.path)]
    if (length(tmp.path)>length.cut) {
      points(loc[start.point,'x'],loc[start.point,'y'],pch=20,col='green',xlab='x',ylab='y',cex=0.7)
      points(loc[end.point,'x'],loc[end.point,'y'],pch=20,col='chocolate4',xlab='x',ylab='y',cex=0.7)
      for (i in 1:(length(tmp.path)-1)) {
        x1=loc[tmp.path[i],'x']
        y1=loc[tmp.path[i],'y']
        
        x2=loc[tmp.path[i+1],'x']
        y2=loc[tmp.path[i+1],'y']
        arrows(x1, y1, x2, y2, length = 0.05, angle = 15, col = "blue")
      } 
    }
  }
}



###3. Plot tumor area
plot_tumor<-function(loc,meta_dat,tumor_cluster){
  loc=loc[rownames(meta_dat),]
  tmp.clusters=meta_dat$seurat_clusters
  tumor.ind=rownames(meta_dat[which(meta_dat$seurat_clusters %in% tumor_cluster),])
  plot(loc,col='grey',xlim=c(min(loc[,'x'])*0.9,max(loc[,'x'])*1.2),main='Tumor region (red)')
  points(loc[tumor.ind,'x'],loc[tumor.ind,'y'],col=alpha('red',0.4),pch=19)
}


### 4. Plot T cell spots
plot_Tcell_spots<-function(keep.ids,loc){
  print(paste('Number of T cell spots:',length(keep.ids),'out of',nrow(loc),'spots'))
  plot(loc,col='grey',xlim=c(min(loc[,'x'])*0.9,max(loc[,'x'])*1.2),main='Tcell spots (pink)')
  points(loc[keep.ids,'x'],loc[keep.ids,'y'],col=alpha('pink',0.4),pch=19,cex=.5)
}



### 5. Plot exhaustion score
plot_pt<-function(loc,meta_dat){
  meta_dat_cell=meta_dat[which(!is.na(meta_dat$ex.score)),]
  loc=loc[rownames(meta_dat),]
  loc=loc[rownames(meta_dat_cell),]
  ex.score=meta_dat_cell$ex.score
  print(quantile(ex.score,probs=seq(0,1,length=21)))
  ggplot(data.frame(loc), aes(x=loc[,1], y=loc[,2], color = ex.score)) +
    geom_point(size=0.1) +
    scale_color_gradient(low = "blue", high = "red") +  labs(title = "Exhaustion score")+
    labs(x='x',y='y')
}


### 6. Boxplot of exhaustion score in each cluster
boxplot_pt<-function(loc,meta_dat){
  meta_dat_cell=meta_dat[which(!is.na(meta_dat$ex.score)),]
  loc=loc[rownames(meta_dat),]
  loc=loc[rownames(meta_dat_cell),]
  ex.score=meta_dat_cell$ex.score
  boxplot(ex.score~meta_dat_cell$seurat_clusters,xlab='Seurat clusters',ylab='Exhaustion score',main='Exhaustion score by Seurat clusters')
}


### 7. To find tumor clusters
clustermeans<-function(count_dat,meta_dat,cancer.markers){
  tmp.markers=cancer.markers[which(cancer.markers %in% colnames(count_dat))]
  tmp.dat=count_dat[,tmp.markers]
  meta_dat=meta_dat[rownames(tmp.dat),]
  
  #print("Median")
  med=apply(tmp.dat,2,function(x) aggregate(x~meta_dat$seurat_clusters,data=tmp.dat,median)[,2])
  #print(med)
  
  #print("Mean")
  mm=apply(tmp.dat,2,function(x) aggregate(x~meta_dat$seurat_clusters,data=tmp.dat,mean)[,2])
  #print(mm)
  return(c(list(med),list(mm)))
}




