##### Prepare data

data_preparation<-function(folder,
                           filenm,
                           loc_file,
                           res,
                           ex.list,
                           cell,
                           cell_cut, # cutoff of sum of T cells markers to define T cell spot
                           cut.t=1e-10,
                           n.vgg ) { #cutoff of exhaustion score difference to cluster spots into one single nodes
  
  ###Load data
  seurat_items=read.spatial.dat(folder=folder,  filenm = filenm, resolut=res,
                                loc_file=loc_file,n.vgg=n.vgg)
  sp.counts.r=seurat_items[[1]]
  #print(sp.counts.r[1:5,1:5])
  sp.counts.norm=seurat_items[[2]]
  meta_dat=seurat_items[[4]]
  loc.raw=seurat_items[[5]]
  vgg=colnames(seurat_items[[3]]) #sp.counts.scale
  
  
  sp.counts.tpm=apply(sp.counts.r,2,function(x) as.numeric(x)/meta_dat[rownames(sp.counts.r),'nCount_Spatial']*1e6)
  rownames(sp.counts.tpm)=rownames(sp.counts.r)
  
  ### Add exhaustion score to meta data
  
  #Add sum of cell markers and sum of exhaustion scores
  meta_dat=add_sum_to_meta(ex.list=ex.list,cell=cell,sp.counts.r=sp.counts.r,meta_dat=meta_dat)
  
  #Define T cell spots
  keep.ids=rownames(meta_dat[which(meta_dat$cell.sum>=cell_cut),]) 
  print(paste("Number of T cell spots are:",length(keep.ids)))
  
  vgg.mean=apply(sp.counts.r[keep.ids,vgg],2,mean)
  quantile(vgg.mean,seq(0,1,0.1))
  vgg=vgg[which(vgg.mean>0.1)]
  
  #Add exhaustion score
  meta_dat=add_pt_to_meta(keep.ids=keep.ids,meta_dat=meta_dat)
  meta_dat_cell=meta_dat[keep.ids,] #meta data of the T cell spots
  
  pt=meta_dat_cell$ex.score #Exhaustion score 
  names(pt)=rownames(meta_dat_cell)
  
  ### Rescale xy coordinates
  loc.rsc=rescale_loc(loc=loc.raw,pt=pt,aa=1)
  
  #Coordinates of T cell spots
  loc1=loc.rsc[keep.ids,]
  
  #Exhaustion score distance matrix among spots
  spot.pt.dist=as.matrix(dist(pt))
  
  rownames(spot.pt.dist)=rownames(loc1)
  colnames(spot.pt.dist)=rownames(loc1)
  
  #unit spatial distance
  spot.sdist=as.matrix(dist(loc1,diag = T))
  r.unit=min(spot.sdist[which(spot.sdist!=0)])*1.015
  
  ### Make subclusters and get their info
  sub.clusters=assign_subcluster(cut.t=cut.t,nr=1.1,spots=keep.ids,spot.sdist=spot.sdist,r.unit=r.unit,spot.pt.dist=spot.pt.dist)
  ll=unlist(lapply(sub.clusters,length))
  print('Number of spots in each node:')
  print(table(ll))
  
  #Get meta data for subclusters
  meta.cts=meta_subclusters(sub.clusters=sub.clusters,loc=loc1,pt=pt,meta_dat=meta_dat,keep.ids=keep.ids)
  names(sub.clusters)=rownames(meta.cts)
  
  #Get location of subclusters
  loc.cts=meta.cts[,c('x','y')]
  pt.cts=meta.cts$pt
  names(pt.cts)=rownames(meta.cts)
  
  
  output=c(list(sp.counts.r),list(sp.counts.norm),list(vgg),list(loc.raw),list(meta_dat),
           list(keep.ids),
           list(meta.cts),list(loc.cts),list(pt.cts),list(r.unit),
           list(sub.clusters),list(sp.counts.tpm))
  
  names(output)=c('raw count data','normalized count data','Top variable genes','XY coordinate of the raw data',
                  'meta data of the spot data','Barcodes of T cell spots',
                  'meta data of the node data',
                  'Rescaled XY coordinate of the node data','Exhaustion scores of the nodes','Distance between adjacent nodes',
                  'Spots in nodes','TPM data')
  return(output)
}




##### Get migration trails
get_all_paths <-function(loc.cts, #rescaled XY coordinate of node data
                         pt.cts, # Exhaustion score of nodes
                         cut.xydist, # maximum XY distance of adjacent node on a path that is allowed
                         cut.sdist.n, # (n-1) thousandth of the pairwise xy distance-threshold for nodes to 3-D line 
                         pt.dist.cuts.n, # (n-1) percent of the pairwise exhaustion score differences-threshold for maximum allowed drawback along the trail
                         min.length, # require a path to at least have 5 nodes
                         r.unit,
                         nr.graph=1.3){ # Edge can exist between two nodes with at most nr.graph*r.unit distance
  
  
  loc.cts=loc.cts[names(pt.cts),]
  keep.cts=names(pt.cts)
  #Space distance matrix among nodes
  cts.sdist=as.matrix(dist(loc.cts))
  #Exhaustion score distance matrix among spots
  cts.pt.dist=as.matrix(dist(pt.cts))
  
  
  ### Set up thresholds
  sdist.cuts=quantile(cts.sdist,seq(0,1,length=1001)) # This is for the threshold of 3-D distance of nodes to straight line
  pt.dist.cuts=quantile(cts.pt.dist,seq(0,1,length=101)) # On a path, PT should be increasing, although there can be small zigzags. This is the threshold of the size of allowable zigzag
  
  cut.sdist=sdist.cuts[cut.sdist.n] #threshold of 3-D distance of nodes to straight line
  cut.pt.dist=pt.dist.cuts[pt.dist.cuts.n] #threshold of the size of allowable PT zigzag
  cut.xydist=cut.xydist #There cannot be big spatial gaps between any consecutive pairs along the path. This is the threshold of allowable xy-plane gap
  
  start.pts=names(pt.cts)[which(pt.cts<=quantile(pt.cts,probs=seq(0,1,length=4))[2])] #First 1/3 most naive nodes as start nodes
  end.pts=names(pt.cts)[which(pt.cts>=quantile(pt.cts,probs=seq(0,1,length=4))[3])] #last 1/3 most exhaustive nodes as end nodes
  print(paste('Number of start nodes:',length(start.pts))) 
  print(paste('Number of end nodes:',length(end.pts))) 
  
  #Plot start and end nodes:
  par(mfrow=c(1,1))
  plot(loc.cts[start.pts,],main='Start(black) and end nodes (red)')
  points(loc.cts[end.pts,'x'],loc.cts[end.pts,'y'],col='red')
  
  
  ### Get 3_D paths
  dds.list=get_all_dists(start.pts=start.pts,end.pts=end.pts,cts.sdist=cts.sdist,s.range1=r.unit*10,s.range2=r.unit*3,
                         loc=loc.cts,pt.cts=pt.cts)
  
  paths_3d=paths_from_3d_line(dds.list=dds.list,
                              cut.sdist=cut.sdist,cut.pt.dist=cut.pt.dist,
                              cut.xydist=cut.xydist,loc=loc.cts,min.length=min.length,
                              pt.cts=pt.cts,cts.sdist=cts.sdist)
  
  #Take a look
  ll=unlist(lapply(paths_3d,length))
  print('The lengh of paths from 3-D line method:')
  print(table(ll))
  
  ###Get MS path
  
  #Make the graph with distance being weight
  tmp.graph=make_graph(spots=keep.cts,cts.sdist=cts.sdist,nr=nr.graph,pt=pt.cts,ws=1,wp=0,r.unit=r.unit)
  graph2=tmp.graph[[1]]
  g.dat2=tmp.graph[[2]]
  
  #Make the graph with difference in exhaustion score as weight
  tmp.graph=make_graph(spots=keep.cts,cts.sdist=cts.sdist,nr=nr.graph,pt=pt.cts,ws=0,wp=1,r.unit=r.unit)
  graph3=tmp.graph[[1]]
  g.dat3=tmp.graph[[2]]
  
  #Get paths
  paths_msp_spacew=c()
  paths_msp_exscorew=c()
  for (i in start.pts) {
    tmp.out2=get_msp1(start.point=i,spots=keep.cts,
                      dist.mat.s=cts.sdist[keep.cts,keep.cts],pt=pt.cts,
                      graph=graph2,cut.length = 4*r.unit,min.length=min.length)
    if (length(tmp.out2)>0) {paths_msp_spacew=c(paths_msp_spacew,list(tmp.out2))}
    tmp.out3=get_msp1(start.point=i,spots=keep.cts,
                      dist.mat.s=cts.sdist[keep.cts,keep.cts],pt=pt.cts,graph=graph3,
                      cut.length = 4*r.unit,min.length=min.length)
    if (length(tmp.out3)>0) {paths_msp_exscorew=c(paths_msp_exscorew,list(tmp.out3))}
  }
  
  ll=unlist(lapply(paths_msp_spacew,length))
  print('The lengh of paths from shortest path, space as weight:')
  print(table(ll))
  
  output=c(list(paths_3d),list(paths_msp_spacew),list(paths_msp_exscorew))
  names(output)=c('paths from 3D','paths from shortest path with distance being weight','paths from shortest path with difference in exhaustion score being weight')
  return(output)
}




##### Select trails based on pairwise correlation of T cell stage-related genes
path_selection<-function(paths,meta_dat,sp.counts.norm,tumor.clus,
                         meta.cts,
                         loc.cts,
                         pt.cts,
                         r.unit,
                         tcell_spots,
                         sub.clusters,
                         cor.cut=0.2,
                         cell,
                         ex.markers){
  
  loc.cts=loc.cts[names(pt.cts),]
  meta.cts=meta.cts[names(pt.cts),]
  keep.cts=names(pt.cts)
  #Space distance matrix among nodes
  cts.sdist=as.matrix(dist(loc.cts))
  
  ### Load the potential paths
  paths_3d=paths[[1]]
  paths_msp_spacew=paths[[2]]
  paths_msp_exscorew=paths[[3]]
  
  path.list=unique(c(paths_3d,do.call(c,lapply(paths_msp_spacew,function(x) x)),do.call(c,lapply(paths_msp_exscorew,function(x) x))))
  print(length(path.list)) 
  names(path.list)=paste0('path',1:length(path.list))
  
  ll=unlist(lapply(path.list,length))
  print(table(ll))
  
  ### Get T cell stage related markers for correlation filtering
  cor.info=cell.markers3(meta_dat=meta_dat,cell=cell,dat=sp.counts.norm,p_threshold = 0.001,tumor.clus = tumor.clus)
  markers=names(cor.info)[which(cor.info>cor.cut)]
  gg=unique(c(intersect(clus.markers[,1],markers),ex.markers))
  print(length(gg))
  print(gg)
  
  ### Get undirected graph, only nodes with cd3+cd4+cd8>=5 were included
  tmp.graph=make_graph_perm(spots=rownames(meta.cts),
                            cts.sdist=cts.sdist,
                            nr=1.3,
                            r.unit=r.unit)
  graph.perm.all=tmp.graph[[1]]
  g.dat.perm.all=tmp.graph[[2]]
  
  ### Get alternative paths connecting start and end node of each path in adjcent areas
  alter.paths=c()
  for (i in 1:length(path.list)) {
    tmp.path=path.list[[i]]
    tmp.alter=get_other_paths2(tmp.path=tmp.path,n.path=6,g.dat.perm = g.dat.perm.all,cts.sdist=cts.sdist)
    alter.paths=c(alter.paths,list(tmp.alter))
  }
  names(alter.paths)=names(path.list)
  ll.alter=unlist(lapply(alter.paths,length))
  print(table(ll.alter))
  
  
  #paths with no alternative paths
  pna=which(ll.alter==0) 
  print(length(pna))
  if (length(pna)>0) {
    plot_all_shortest_paths(all.paths=path.list[pna],loc=loc.cts,length.cut = 0)
    alter.paths=alter.paths[-pna]
  }
  print(length(alter.paths))
  ll.alter=unlist(lapply(alter.paths,length))
  table(ll.alter)
  
  gpt.norm=get_large_cor(dat=sp.counts.norm[tcell_spots,gg],path.list=path.list,alter.paths=alter.paths,pna=pna,top_percent=25,
                         sub.clusters = sub.clusters,meta.cts=meta.cts)
  
  return(gpt.norm)
}



##### Select trails based on pairwise correlation of exhaustion markers
path_selection1<-function(paths,meta_dat,sp.counts.norm,tumor.clus,
                          meta.cts,
                          loc.cts,
                          pt.cts,
                          r.unit,
                          tcell_spots,
                          sub.clusters,
                          cor.cut=0.2,
                          cell,
                          ex.markers){
  
  loc.cts=loc.cts[names(pt.cts),]
  meta.cts=meta.cts[names(pt.cts),]
  keep.cts=names(pt.cts)
  #Space distance matrix among nodes
  cts.sdist=as.matrix(dist(loc.cts))
  
  ### Load the potential paths
  paths_3d=paths[[1]]
  paths_msp_spacew=paths[[2]]
  paths_msp_exscorew=paths[[3]]
  
  path.list=unique(c(paths_3d,do.call(c,lapply(paths_msp_spacew,function(x) x)),do.call(c,lapply(paths_msp_exscorew,function(x) x))))
  print(length(path.list)) 
  names(path.list)=paste0('path',1:length(path.list))
  
  ll=unlist(lapply(path.list,length))
  print(table(ll))
  
  ### Get T cell stage related markers for correlation filtering
  cor.info=cell.markers3(meta_dat=meta_dat,cell=cell,dat=sp.counts.norm,p_threshold = 0.001,tumor.clus = tumor.clus)
  gg=ex.markers
  
  ### Get undirected graph, only nodes with cd3+cd4+cd8>=5 were included
  tmp.graph=make_graph_perm(spots=rownames(meta.cts),
                            cts.sdist=cts.sdist,
                            nr=1.3,
                            r.unit=r.unit)
  graph.perm.all=tmp.graph[[1]]
  g.dat.perm.all=tmp.graph[[2]]
  
  ### Get alternative paths connecting start and end node of each path in adjcent areas
  alter.paths=c()
  for (i in 1:length(path.list)) {
    tmp.path=path.list[[i]]
    tmp.alter=get_other_paths2(tmp.path=tmp.path,n.path=6,g.dat.perm = g.dat.perm.all,cts.sdist=cts.sdist)
    alter.paths=c(alter.paths,list(tmp.alter))
  }
  names(alter.paths)=names(path.list)
  ll.alter=unlist(lapply(alter.paths,length))
  print(table(ll.alter))
  
  
  #paths with no alternative paths
  pna=which(ll.alter==0) 
  print(length(pna))
  if (length(pna)>0) {
    plot_all_shortest_paths(all.paths=path.list[pna],loc=loc.cts,length.cut = 0)
    alter.paths=alter.paths[-pna]
  }
  print(length(alter.paths))
  ll.alter=unlist(lapply(alter.paths,length))
  table(ll.alter)
  
  gpt.norm=get_large_cor(dat=sp.counts.norm[tcell_spots,gg],path.list=path.list,alter.paths=alter.paths,pna=pna,top_percent=25,
                         sub.clusters = sub.clusters,meta.cts=meta.cts)
  
  return(gpt.norm)
}





