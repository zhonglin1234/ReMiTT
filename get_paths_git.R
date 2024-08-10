require(Seurat)
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


colfunc<-colorRampPalette(c("royalblue","springgreen","yellow",'red'))
#plot(rep(1,50),col=(colfunc(50)), pch=19,cex=2)


####################################Data preparation###########################################################################


##### Load spatial data to Seurat
read.spatial.dat<-function(folder,filenm,resolut=0.3,loc_file,n.vgg=3000){
  spdat <- Load10X_Spatial(folder,
                           filename = filenm,
                           assay = "Spatial", 
                           slice = "slice1",
                           filter.matrix = TRUE,
                           to.upper = FALSE,
                           image = NULL)
  
  spdat=SCTransform(spdat,assay = 'Spatial',variable.features.n = n.vgg) 
  spdat=RunPCA(spdat)
  spdat=FindNeighbors(spdat,dims=1:30)
  spdat=FindClusters(spdat,resolution=resolut)
  spdat=RunUMAP(spdat,dims=1:30)
  DimPlot(spdat,reduction='umap',label = TRUE,repel=TRUE)
  
  sp.counts.r=t(as.matrix(spdat@assays$Spatial$counts))
  sp.counts.norm=t(as.matrix(as.matrix(spdat@assays$SCT@data))) #Normalized data
  sp.counts.r=sp.counts.r[which(rownames(sp.counts.r) %in% rownames(sp.counts.norm)),
                          which(colnames(sp.counts.r) %in% colnames(sp.counts.norm))]
  meta_dat=spdat@meta.data
  
  sp.counts.norm=sp.counts.norm[rownames(meta_dat),]
  sp.counts.r=sp.counts.r[rownames(meta_dat),]
  sp.counts.scale=t(spdat@assays$SCT$scale.data)
  sp.counts.scale=sp.counts.scale[rownames(meta_dat),]
  
  #print(sp.counts.r[1:5,1:5])
  row.sum=rowSums(sp.counts.r)
  col.sum=colSums(sp.counts.r)
  print(paste("Minimum number of scripts per row:",min(row.sum)))
  print(paste("Minimum number of scripts per colum:",min(col.sum)))
  
  ###Location
  tmp.loc=read.csv(loc_file,sep=',',header=F)
  #dim(tmp.loc)
  tmp.loc=tmp.loc[which(tmp.loc[,1] %in% rownames(sp.counts.r)),]
  rownames(tmp.loc)=tmp.loc[,1]
  tmp.loc=tmp.loc[,5:6]
  colnames(tmp.loc)=c('x','y')
  loc.raw=data.frame(tmp.loc)
  
  loc.raw$x=as.numeric(loc.raw$x)
  loc.raw$y=as.numeric(loc.raw$y)
  
  loc.raw=loc.raw[rownames(meta_dat),]
  
  out=c(list(sp.counts.r),list(sp.counts.norm),list(sp.counts.scale),list(meta_dat),list(loc.raw))
  
  names(out)=c('sp.counts.r','sp.counts.norm','sp.counts.scale','meta_dat','loc.raw')
  return(out)
}


#####Take a look at the data

### Plot Seurat clusters
plot_seurat_clusters<-function(loc,meta_dat,sp.counts.norm,sp.counts.r,cancer.markers){
  loc=loc[rownames(meta_dat),]
  loc[,'x']=as.numeric(loc[,'x'])
  loc[,'y']=as.numeric(loc[,'y'])
  tmp.clusters=meta_dat$seurat_clusters
  number.clusters=length(table(tmp.clusters))
  col.list=c('darkblue','blue','lightblue','cyan','green','yellow','orange','chocolate','brown','red')
  plot(loc,col='white',xlim=c(min(loc[,'x'])*0.9,max(loc[,'x'])*1.2),main='Seurat clusters')
  for (i in 1:number.clusters) {
    points(loc[tmp.clusters==(i-1),'x'],loc[tmp.clusters==(i-1),'y'],col=col.list[i],pch=19,cex=.5)
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


### Plot tumor area
plot_tumor<-function(loc,meta_dat,tumor_cluster){
  loc=loc[rownames(meta_dat),]
  tmp.clusters=meta_dat$seurat_clusters
  tumor.ind=rownames(meta_dat[which(meta_dat$seurat_clusters %in% tumor_cluster),])
  plot(loc,col='grey',xlim=c(min(loc[,'x'])*0.9,max(loc[,'x'])*1.2),main='Tumor region (red)')
  points(loc[tumor.ind,'x'],loc[tumor.ind,'y'],col=alpha('red',0.4),pch=19)
}


### Plot T cell spots
plot_Tcell_spots<-function(keep.ids,loc){
  print(paste('Number of T cell spots:',length(keep.ids),'out of',nrow(loc),'spots'))
  plot(loc,col='grey',xlim=c(min(loc[,'x'])*0.9,max(loc[,'x'])*1.2),main='Tcell spots (pink)')
  points(loc[keep.ids,'x'],loc[keep.ids,'y'],col=alpha('pink',0.4),pch=19,cex=.5)
}

### Plot exhaustion score
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

boxplot_pt<-function(loc,meta_dat){
  meta_dat_cell=meta_dat[which(!is.na(meta_dat$ex.score)),]
  loc=loc[rownames(meta_dat),]
  loc=loc[rownames(meta_dat_cell),]
  ex.score=meta_dat_cell$ex.score
  boxplot(ex.score~meta_dat_cell$seurat_clusters,xlab='Seurat clusters',ylab='Exhaustion score',main='Exhaustion score by Seurat clusters')
}


### Find tumor clusters
clustermeans<-function(count_dat,meta_dat,cancer.markers){
  tmp.markers=cancer.markers[which(cancer.markers %in% colnames(count_dat))]
  tmp.dat=count_dat[,tmp.markers]
  meta_dat=meta_dat[rownames(tmp.dat),]
  
  #print("Median")
  med=apply(tmp.dat,2,function(x) aggregate(x~meta_dat$seurat_clusters,data=meta_dat,median)[,2])
  #print(med)
  
  #print("Mean")
  mm=apply(tmp.dat,2,function(x) aggregate(x~meta_dat$seurat_clusters,data=meta_dat,mean)[,2])
  #print(mm)
  return(c(list(med),list(mm)))
}

### Assign Exhaustion score (PT) to each spot
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


###Assign subclusters and create meta data for subclusters

assign_subcluster<-function(cut.t,spots,spot.sdist,nr,r.unit,spot.pt.dist){ 
  #cut.t is the threshold of pseudotime distance between spots that can form a subcluster,nr is the threshold of number of space distance unit
  tmp.list=c()
  for (i in spots) {
    nbs=names(which(spot.sdist[i,]<nr*r.unit))
    if (length(nbs)==0) {print('error')}
    tmp.clus=nbs[which(spot.pt.dist[i,nbs]<cut.t)]
    tmp.list=c(tmp.list,list(tmp.clus))
  }
  tmp.list=unique(tmp.list)
  
  ll=unlist(lapply(tmp.list,length))
  if (max(ll)<=2) {
    return(tmp.list)
  }
  if (max(ll)>2) {
    final.list=c()
    for (j in 2:(max(ll)-1)){
      s.clus=tmp.list[ll==j]
      b.clus=tmp.list[ll>j]
      max.overlap=unlist(lapply(s.clus,function(x) max(n.overlap(v1=x,v.list=b.clus))))
      if (max(max.overlap)==j) {
        s.clus=s.clus[-which(max.overlap==j)]
      }
      final.list=c(final.list,s.clus)
    }
    final.list=c(tmp.list[ll==1],final.list,tmp.list[ll==max(ll)])
    print(sum(unlist(lapply(final.list,length))))
    print(table(unlist(lapply(final.list,length))))
    return(final.list)
  }
}
meta_subclusters<-function(sub.clusters,loc,pt,meta_dat,keep.ids){
  meta_dat_cell=meta_dat[keep.ids,]
  
  loc.cts=do.call(rbind,lapply(sub.clusters,function(x) colSums(loc[x,])/length(x)))
  rownames(loc.cts)=paste0('cts',1:nrow(loc.cts))
  colnames(loc.cts)=c('x','y')
  pt.cts=unlist(lapply(sub.clusters,function(x) mean(pt[x])))
  names(pt.cts)=rownames(loc.cts)
  
  tmp.cluster=lapply(sub.clusters,function(x) meta_dat_cell[x,'seurat_clusters'])
  tmp.cluster=lapply(tmp.cluster,function(x) prop.table(table(x)))
  tmp.cluster=unlist(lapply(tmp.cluster,function(x) names(x)[which(x==max(x))][1]))
  cts.meta=data.frame(loc.cts,pt.cts,tmp.cluster)
  colnames(cts.meta)=c('x','y','pt','seurat_clusters')
  return(cts.meta)
}




###################Path identification#################################################################################################################################################

##### 3-D line related functions

### Rescale location coordinates based on range of exhaustion score 
rescale_loc<-function(keep.ids=keep.ids,loc,pt,aa=1){ #aa act as a weight of physical distance
  loc=data.frame(loc[keep.ids,])
  loc$x=as.numeric(loc$x)
  loc$y=as.numeric(loc$y)
  
  loc$x=loc$x-min(loc$x)
  loc$y=loc$y-min(loc$y)
  loc$x=loc$x*max(pt)*aa/max(loc$x)
  loc$y=loc$y*max(pt)*aa/max(loc$y)
  plot(loc)
  return(loc)
}


### For a node, get the nodes within certain distance from it, on the xy-plane
get_local_cts<-function(stp,cts.sdist,s.range){
  local.cts=names(which(cts.sdist[stp,]<=s.range))
  #print(length(local.cts))
  
  plot(loc.cts,col='grey')
  points(loc.cts[local.cts,'x'],loc.cts[local.cts,'y'],col='lightblue',pch=19,cex=0.5)
  points(loc.cts[stp,'x'],loc.cts[stp,'y'],col='green',pch=19,cex=0.7)
  return(local.cts)
}

### Get point to line distance
p_to_l_dist<-function(p,p.start,p.end){ #m0 is the point out of line, m1, m2 are points on the line
  projected.p=orthogonalProjectionToLine(p=p, l1=list(l1=p.start, l2=p.end))
  return(dist(rbind(p,projected.p)))
}


### Get distance to line for all nodes in one area
get_dists<-function(stp,end.pts,cts.sdist,s.range1,s.range2,loc,pt.cts){
  local.cts=names(which(cts.sdist[stp,]<=s.range1))
  local.pt=pt.cts[local.cts]
  #print(length(local.cts))
  
  tmp.end.pts=local.cts[which(local.cts %in% end.pts)]
  tmp.end.pts=names(which(cts.sdist[stp,tmp.end.pts]>=s.range2))
  #print(length(tmp.end.pts))
  
  start=as.numeric(c(loc[stp,],pt.cts[stp]))
  if (length(tmp.end.pts)==0) {return(NULL)}
  
  #Sys.time()
  dds.mat=matrix(NA,nrow=length(local.cts),ncol=length(tmp.end.pts))
  rownames(dds.mat)=local.cts
  colnames(dds.mat)=tmp.end.pts
  
  for (edp in tmp.end.pts) {
    end=as.numeric(c(loc[edp,],pt.cts[edp]))
    dds.mat[,edp]=sapply(local.cts,function(x) p_to_l_dist(p=as.numeric(c(loc[x,],pt.cts[x])),p.start=start,p.end=end))
  }
  
  #Sys.time()
  colnames(dds.mat)=paste(stp,tmp.end.pts,sep='_')
  #print(dds.mat[1:5,1:5])
  return(dds.mat)
}

###Get info of each path

get_path_info1<-function(path,loc=loc.cts,pt.cts,cts.sdist){
  tmp.pt=pt.cts[path]
  s.dists=cts.sdist[path[1:(length(path)-1)],path[2:length(path)]]
  s.dists=sapply(1:(length(path)-1),function(x) s.dists[x,x])
  pt.dists=tmp.pt[2:length(path)]-tmp.pt[1:(length(path)-1)]
  
  dist.mat=cbind(s.dists,pt.dists)
  colnames(dist.mat)=c('Space_dist','PT_dist')
  rownames(dist.mat)=paste(path[1:(length(path)-1)],path[2:length(path)],sep='_')
  
  out=c(list(dist.mat),list(pt.cts[path]))
  names(out)=c('Distance between cts','PT')
  return(out)
} #not plot the path


### Remove trails with gap between consecutive pairs
filter_by_gap<-function(dds.mat,cut.sdist,cut.xydist,cut.pt.dist,loc,min.length,remove_names=T,
                        pt.cts,cts.sdist){ #to use get_path_info2 to plot path, set remove_names to be F
  if (is.null(dds.mat)) {return(NULL)}
  local.cts=rownames(dds.mat)
  tmp=apply(dds.mat,2,function(x) local.cts[which(x<=cut.sdist)])
  ll=unlist(lapply(tmp,length))
  
  if (max(ll)<5) {return(NULL)}
  
  dds.mat=dds.mat[,ll>=5]
  tmp=tmp[which(ll>=5)] #number of nodes must>5
  ll=unlist(lapply(tmp,length))
  
  path.list=c()
  for (j in 1:length(tmp)) {
    start=unlist(strsplit(names(tmp)[j],'\\_'))[1]
    end=unlist(strsplit(names(tmp)[j],'\\_'))[2]
    path=tmp[[j]]
    
    #Order nodes on the path along the line between start and endpoint
    start.loc=as.numeric(loc[start,])
    end.loc=as.numeric(loc[end,])
    
    Line <- CreateLinePoints(start.loc, end.loc)
    #Draw(Line, "red")
    projection=t(sapply(path,function(x) orthogonalProjectionToLine(p=as.numeric(loc[x,]), l1=start.loc, l2=end.loc)))
    path=path[order(projection[,1])]
    if (pt.cts[path[1]] > pt.cts[path[length(path)]]) {path=path[length(path):1]}
    
    tmp.info=get_path_info1(path=path,loc=loc.cts, pt.cts =pt.cts, cts.sdist=cts.sdist)[[1]]
    
    #Remove those nodes with wrong z gaps
    while (min(tmp.info[,'PT_dist'])<cut.pt.dist*(-1)) {    
      rm1=which(tmp.info[,'PT_dist'] < cut.pt.dist*(-1))
      if (length(rm1)>0) {path=path[-rm1]}
      if (length(path)<min.length) {break}
      tmp.info=get_path_info1(path=path,loc=loc.cts, pt.cts =pt.cts, cts.sdist=cts.sdist)[[1]]
    }
    
    #Remove paths with big space gaps
    #print(max(tmp.info[,'Space_dist'])<cut.xydist)
    
    if (max(tmp.info[,'Space_dist'])<cut.xydist & length(path)>=min.length) {
      names(path)=paste(start,end,sep='_')
      path.list=c(path.list,list(path))
    }
  }
  if (remove_names==T & length(path.list)>0) {
    for (i in 1:length(path.list)) {
      names(path.list[[i]])=NULL
    }
  }
  return(path.list)
}

### Get the trails

#Get distance to 3D line for each start points 
get_all_dists<-function(start.pts,end.pts,cts.sdist,s.range1,s.range2,loc,pt.cts){
  dds.list=c()
  for (k in start.pts) {
    print(k)
    tmp.dds=get_dists(stp=k,end.pts=end.pts,cts.sdist=cts.sdist,s.range1=s.range1,s.range2=s.range2,loc=loc,pt.cts=pt.cts)
    dds.list=c(dds.list,list(tmp.dds))
  }
  names(dds.list)=start.pts
  return(dds.list)
}

# Get potential trails
paths_from_3d_line<-function(dds.list,cut.sdist,cut.pt.dist,cut.xydist,loc,min.length,pt.cts,cts.sdist){
  all.paths=c()
  for (i in 1:length(dds.list)) {
    print(i)
    dds.mat=dds.list[[i]]
    all.paths=c(all.paths,filter_by_gap(dds.mat=dds.mat,cut.sdist=cut.sdist,cut.xydist=cut.xydist,
                                        cut.pt.dist=cut.pt.dist,loc=loc,
                                        min.length=min.length,
                                        pt.cts=pt.cts,cts.sdist=cts.sdist))
  }
  return(all.paths)
}



##### Minimum spanning path method


### Make the graph with all T cell nodes

# Get neighbors of a spot, not including the spot
get_nbs<-function(spot,dist.mat,nr,r.unit){
  tmp=names(which(dist.mat[spot,]<=nr*r.unit))
  return(tmp[-which(tmp==spot)])
}

# Get neighbors of a spot, including the spot
get_nbs1<-function(spot,dist.mat,nr,r.unit){
  tmp=names(which(dist.mat[spot,]<=nr*r.unit))
  return(tmp)
}

# Get edge weight (difference in Exhaustion score)
edge_wt<-function(spot,dist.mat,pt,nr,r.unit) {
  nbs=get_nbs(spot=spot,dist.mat=dist.mat,nr=nr,r.unit=r.unit)
  if (length(nbs)>0) {
    tmp.edges=as.matrix(cbind(rep(spot,length(nbs)),nbs))
    drct=pt[nbs]-pt[spot]
    wt1=abs(drct)
    wt2=dist.mat[spot,nbs]
    
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
make_graph<-function(spots,cts.sdist,nr,pt,ws,wp,r.unit){
  g.dat=do.call(rbind,sapply(spots,function(x) edge_wt(spot=x,dist.mat=cts.sdist[spots,spots],nr=nr,pt=pt,r.unit=r.unit)))
  g.dat=unique(g.dat)
  rownames(g.dat)=paste(g.dat[,'edge0'],g.dat[,'edge1'],sep='_')
  
  edges <- data.frame(from = g.dat[,'edge0'],to = g.dat[,'edge1'])
  weights <- as.numeric(g.dat[,'dist_s'])*ws+as.numeric(g.dat[,'dist_pt'])*wp   #Weight is a combination of space distance and pseudotime distance
  
  graph <- graph_from_data_frame(edges, directed = TRUE)
  E(graph)$weight <- weights
  return(c(list(graph),list(g.dat)))
}


### Get MSP 
get_msp1<-function(start.point,spots,dist.mat.s,pt,graph,cut.length,min.length){
  if (!start.point%in%names(V(graph))) {
    print('Start point not in the graph')
    return(NULL)
  }
  end.points=spots[which(pt[spots]-pt[start.point]>quantile(pt,probs = seq(0,1,length=10))[3])] #has to have bigger pt than start point for at least 20% 
  end.points=end.points[which(end.points %in% names(V(graph)))] #has to be in graph
  end.points=end.points[which(dist.mat.s[start.point,end.points]>cut.length)] #cannot be too close to start point
  
  if (length(end.points)==0) {return(NULL)}
  
  path.list=c()
  for (end.point in end.points) {
    msp <- shortest_paths(graph, from = start.point, to = end.point, weights = E(graph)$weight)
    vv=unlist(msp$vpath)
    if (length(vv)>=min.length) {path.list=c(path.list,list(names(vv)))}
  }
  return(path.list)
}


### Make one undirected graph for getting control trails, assign weight to be euclidean distance 
edge_wt_perm<-function(spot,dist.mat,nr,r.unit) {
  nbs=get_nbs(spot=spot,dist.mat=dist.mat,nr=nr,r.unit=r.unit)
  if (length(nbs)>0) {
    tmp.edges=as.matrix(cbind(rep(spot,length(nbs)),nbs))
    wt=dist.mat[spot,nbs]
    out=cbind(tmp.edges,wt)
    rownames(out)=paste(spot,1:nrow(out),sep='_')
    colnames(out)=c('edge0','edge1','dist_s')
    return(out)
  }
}

### Make the graph
make_graph_perm<-function(spots,cts.sdist,nr,r.unit){ #Undirected graph with weight being spatial distance
  g.dat=do.call(rbind,sapply(spots,function(x) edge_wt_perm(spot=x,dist.mat=cts.sdist[spots,spots],nr=nr,r.unit=r.unit)))
  g.dat=unique(g.dat)
  rownames(g.dat)=paste(g.dat[,'edge0'],g.dat[,'edge1'],sep='_')
  
  edges <- data.frame(from = g.dat[,'edge0'],to = g.dat[,'edge1'])
  weights <- as.numeric(g.dat[,'dist_s'])  
  
  graph <- graph_from_data_frame(edges, directed = FALSE)
  E(graph)$weight <- weights
  return(c(list(graph),list(g.dat)))
}


### Remove duplicated trails

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




##############Path validation##################################################################################################################################################

##### Plot the trails

### Plot trails
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



### plot alternative paths
plot_alternative_paths<-function(tmp.path,loc=loc.cts,cts.sdist){
  alter.paths=get_other_paths2(tmp.path=tmp.path,n.path=100,cts.sdist=cts.sdist,g.dat.perm = g.dat.perm)
  par(mfrow=c(1,1))
  plot(loc,col='lightgrey')
  points(loc[tmp.path,'x'],loc[tmp.path,'y'],pch=20,col=colfunc(length(tmp.path)),xlab='x',ylab='y',cex=.5)
  for (i in 1:length(alter.paths)) {
    vv=alter.paths[[i]]
    print(vv)
    lines(loc[vv,'x'],loc[vv,'y'],col="pink",lty=1)
  }
}



##### Select genes for comparing pairwise correlation between trail vs. alternative routes

### Correlation with the expression of certain markers in tumor region
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


### Get alternative paths around gradient paths and select the path with largest correlation among its alternative paths

get_short_paths<-function(start.p,end.p,tmp.path,n.path,g.dat.perm,cts.sdist,r.unit){
  nb.pps=c()
  for (pp in tmp.path) {
    nb.pps=c(nb.pps,get_nbs(spot=pp,dist.mat=cts.sdist,nr=10,r.unit=r.unit))
  }
  nb.pps=unique(nb.pps)
  
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

get_other_paths2<-function(tmp.path,n.path,g.dat.perm,cts.sdist){
  realstart=tmp.path[1]
  realend=tmp.path[length(tmp.path)]
  realpath=tmp.path
  
  path.list1=get_short_paths(start.p=realstart,end.p=realend,tmp.path=tmp.path,n.path=n.path,g.dat.perm=g.dat.perm,cts.sdist=cts.sdist,r.unit)
  if (length(path.list1)==0) return(NULL)
  
  tmp.path=path.list1[[1]]
  path.list2=get_short_paths(start.p=tmp.path[2],end.p=tmp.path[length(tmp.path)-1],tmp.path=tmp.path,n.path=n.path,g.dat.perm=g.dat.perm,cts.sdist=cts.sdist,r.unit)
  path.list2=lapply(path.list2,function(x) c(realstart,x,realend))
  rm.ind=unlist(lapply(path.list2,function(x) length(intersect(x,realpath[2:(length(realpath)-1)]))))
  path.list2=path.list2[-which(rm.ind>0)]
  
  
  if (length(path.list1)==1) return(c(path.list1,path.list2)) #there could be duplicates, or path that has nodes in real path
  
  tmp.path=path.list1[[2]]
  path.list3=get_short_paths(start.p=tmp.path[2],end.p=tmp.path[length(tmp.path)-1],tmp.path=tmp.path,n.path=n.path,g.dat.perm=g.dat.perm,cts.sdist=cts.sdist,r.unit)
  path.list3=lapply(path.list3,function(x) c(realstart,x,realend))
  rm.ind=unlist(lapply(path.list3,function(x) length(intersect(x,realpath[2:(length(realpath)-1)]))))
  path.list3=path.list3[-which(rm.ind>0)]
  
  if (length(path.list1)==2) return(c(path.list1,path.list2,path.list3)) #there could be duplicates, or path that has nodes in real path
  
  tmp.path=path.list1[[3]]
  path.list4=get_short_paths(start.p=tmp.path[2],end.p=tmp.path[length(tmp.path)-1],tmp.path=tmp.path,n.path=n.path,g.dat.perm=g.dat.perm,cts.sdist=cts.sdist,r.unit)
  path.list4=lapply(path.list4,function(x) c(realstart,x,realend))
  rm.ind=unlist(lapply(path.list4,function(x) length(intersect(x,realpath[2:(length(realpath)-1)]))))
  path.list4=path.list4[-which(rm.ind>0)]
  
  return(unique(c(path.list1,path.list2,path.list3,path.list4))) #there could be path that has nodes in real path
}



### Calculate consecutive pairwise correlation
pairwaise_correlation<-function(tmp.path,cor.mat,sub.clusters){
  #Get spots
  tmp.spots=do.call(c,sub.clusters[tmp.path])
  #Get consecutive pairwise correlation
  tmp.cor=cor.mat[tmp.spots,tmp.spots]
  pair.cor=apply(cbind(1:(length(tmp.path)-1),2:length(tmp.path)),1,function(x) tmp.cor[x[1],x[2]])
  return(c(list(pair.cor),list(mean(pair.cor))))
}


###Get mean pairwise correlation of a list of trails
get_mean_cor<-function(path.list,cor.mat,sub.clusters){
  path.cor.list=lapply(path.list,function(x) pairwaise_correlation(tmp.path=x,cor.mat=cor.mat,sub.clusters=sub.clusters))
  mean.cor.list=unlist(lapply(path.cor.list,function(x) x[[2]]))
  names(mean.cor.list)=names(path.list)
  mean.cor.list=mean.cor.list[order(mean.cor.list,decreasing = T)]
  return(mean.cor.list)
}

###Get trails with largest correlation compared to its alternative routes
get_large_cor<-function(dat,path.list,alter.paths,pna,top_percent=25,sub.clusters,meta.cts){
  cor.mat2=cor(t(dat)) #Pearson's correlation between spots
  diag(cor.mat2)=NA
  #cut.trds=quantile(cor.mat2,seq(0,1,0.05),na.rm = T)
  #print(cut.trds)
  
  gpt.cor2=get_mean_cor(path.list=path.list,cor.mat=cor.mat2,sub.clusters=sub.clusters)
  
  if (length(pna)>0) {
    print('pna correlations')
    print(gpt.cor2[pna])
  }
  
  
  #Correlation with pairwise correlation and UMI
  #plot(unlist(lapply(path.list,get_path_umi)),gpt.cor2)
  #print(cor(unlist(lapply(path.list,get_path_umi)),gpt.cor2))
  
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
  
  gpt1=names(gpt.cor2a)[which(gpt.cor2a>alter.max2)]
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




##### Get matched control set

### Get lots of random control paths with different lengths

get_permutated_paths<-function(meta.cts,cts.sdist,s.range1=7*r.unit,s.range2=3*r.unit,graph.perm=graph.perm.tumor,n.paths=30000){
  start.end.pairs=c()
  tmp.pp=rownames(meta.cts)[which(rownames(meta.cts)%in% names(V(graph.perm)))]
  #print(length(tmp.pp))
  tmp.sdist=cts.sdist[tmp.pp,tmp.pp]
  for (i in 1:ncol(tmp.sdist)){
    tmp.dd=tmp.sdist[,i]
    keep.dd=which(tmp.dd>=s.range2 & tmp.dd<=s.range1)
    if (length(keep.dd)<2) {next}
    tmp.pairs=cbind(rownames(tmp.sdist)[keep.dd],rep(colnames(tmp.sdist)[i],length(keep.dd)))
    tmp.pairs=t(apply(tmp.pairs,1,function(x)x[order(x)]))
    start.end.pairs=rbind(start.end.pairs,tmp.pairs)
  }
  
  start.end.pairs=unique(start.end.pairs)
  print(length(start.end.pairs))
  
  tmp.pairs=start.end.pairs[sample(1:nrow(start.end.pairs),n.paths),]
  
  perm.paths.list=c()
  for (i in 1:n.paths) {
    msp=shortest_paths(graph.perm, from = tmp.pairs[i,1], to = tmp.pairs[i,2])
    vv=names(unlist(msp$vpath))
    perm.paths.list=c(perm.paths.list,list(vv))
  }
  
  ll=unlist(lapply(perm.paths.list,function(x)length(x)))
  print(table(ll))
  return(perm.paths.list[which(ll>0)])
}


### Check which cluster each control trail belongs to

check_path_cluster<-function(tmp.path,meta.cts){
  s.clus=meta.cts[tmp.path,'seurat_clusters']
  tmp.tab=table(s.clus)
  return(names(tmp.tab)[which(tmp.tab==max(tmp.tab))[1]])
}


### Get a controlset based on length and cluster of each migration trail
perm.paths <-function(perm.paths.list,real.path.list,meta.cts){
  out=c()
  for (i in 1:length(real.path.list)) {
    tmp.ll=real.ll[i]
    tmp.clus=path.clus[i]
    tmp.pool=perm.paths.list[which(perm.ll==tmp.ll & perm.path.clus==tmp.clus)]
    out=c(out,tmp.pool[sample(1:length(tmp.pool),1)])
  }
  return(out)
}

get_control_paths<-function(paths,meta.cts,nr,s.range1,s.range2,n.paths,loc.cts,r.unit,n.set){
  
  ###Get alternative paths of all kinds of lengths
  cts.sdist=as.matrix(dist(loc.cts))
  paths.nodes=unique(do.call(c,paths))
  path.clus=unlist(lapply(paths,function(x)check_path_cluster(x,meta.cts=meta.cts)))
  all.nodes=rownames(meta.cts)[which(meta.cts$seurat_clusters %in% unique(path.clus))]
  other.nodes=all.nodes[which(!all.nodes %in% paths.nodes)]
  tmp.graph=make_graph_perm(spots=other.nodes,cts.sdist=cts.sdist,nr=nr,r.unit=r.unit)
  perm.paths.list=get_permutated_paths(meta.cts=meta.cts,s.range1=s.range1,s.range2=s.range2,
                                       graph.perm=tmp.graph[[1]],n.paths=n.paths,cts.sdist=cts.sdist)
  
  
  #Length
  real.ll=unlist(lapply(paths,length))
  #Cluster
  path.clus=unlist(lapply(paths,function(x)check_path_cluster(x,meta.cts=meta.cts)))
  
  #Length
  perm.ll=unlist(lapply(perm.paths.list,length))
  #Cluster
  perm.path.clus=unlist(lapply(perm.paths.list,function(x)check_path_cluster(x,meta.cts=meta.cts)))
  
  ll.list=names(table(real.ll))
  clus.list=names(table(path.clus))
  
  combs=expand.grid(ll.list,clus.list)
  
  perm.pathset.norm=c()
  k=0 
  while(k<n.set) {
    k=k+1
    pathset=as.list(rep(NA, length(paths)))
    for (i in 1:nrow(combs)) {
      ll=as.numeric(as.character(combs[i,1]))
      clus=as.numeric(as.character(combs[i,2]))
      path.ids=which(real.ll==ll & path.clus==clus)
      n.path=length(path.ids)
      if (n.path==0) {next}
      tmp.pool=perm.paths.list[which(perm.ll==ll & perm.path.clus==clus)]
      pathset[path.ids]=tmp.pool[sample(1:length(tmp.pool),n.path)]
    }
    perm.pathset.norm=c(perm.pathset.norm,list(pathset))
  }
  return(perm.pathset.norm)
}


