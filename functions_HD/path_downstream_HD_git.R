##############Universal functions for migration path finding project##########################
#Date started: 03/14/2023

library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
require(lme4)
require(multcomp)
require(ggrepel)
require(ggplot2)
library(ggrepel)
library(msigdbr)
library(clusterProfiler)
library(StereoMorph)
require(Seurat)
require(dplyr)
require(igraph)
library(plotly)
library(RColorBrewer)
require(openxlsx)
require(LearnGeom)
require(StereoMorph)
library(gplots)
require(ppcor)



################Path validation ############################################################################################################################################




######  Get control sets for systematic analysis   ###### 

### 1. Get random control path sets

# (1) Get lots of random control paths with different lengths

# Requires: dbscan, igraph
# - meta.cts: data.frame with rownames = spot IDs (used to pick valid nodes)
# - loc.cts:  data.frame/matrix with rows named by spot IDs and (x,y) coords (or more dims)
# - s.range1: upper distance bound
# - s.range2: lower distance bound
# - graph.perm: igraph with vertex names matching spot IDs
# - n.paths: number of random start–end pairs to try
get_permutated_paths <- function(meta.cts,
                                 loc.cts,
                                 s.range1 = 7 * r.unit.raw,
                                 s.range2 = 3 * r.unit.raw,
                                 graph.perm = graph.perm.tumor,
                                 n.paths = 30000,
                                 seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1) Keep only spots present in the graph
  tmp.pp <- intersect(rownames(meta.cts), igraph::V(graph.perm)$name)

  # 2) Coordinates for those spots (keeps order aligned with tmp.pp)
  loc.pp <- loc.cts[tmp.pp,]
  
  # 3) Radius neighbors up to s.range1 (much cheaper than a full distance matrix)
  #    frNN returns, for each point i, a vector of neighbor indices and distances within eps
  fr <- dbscan::frNN(as.matrix(loc.pp), eps = s.range1, sort = TRUE)
  
  # 4) Build candidate start–end pairs with distance in [s.range2, s.range1]
  #    Avoid self-pairs and de-duplicate by sorting each pair's names
  ids <- names(tmp.pp) <- tmp.pp
  
  make_pairs <- function(i) {
    nbr_idx <- fr$id[[i]]
    if (length(nbr_idx) == 0) return(NULL)
    nbr_dist <- fr$dist[[i]]
    keep <- which(nbr_dist >= s.range2 & nbr_dist <= s.range1 & nbr_idx != i)
    if (!length(keep)) return(NULL)
    a <- tmp.pp[i]
    b <- tmp.pp[nbr_idx[keep]]
    # sort each pair so (a,b) == (b,a)
    p <- cbind(pmin(a, b), pmax(a, b))
    unique(p)
  }
  
  if (length(tmp.pp)>8000) {
    pair_list <- lapply(seq_along(sample(tmp.pp,8000)), make_pairs)
  }
  else {pair_list <- lapply(seq_along(tmp.pp), make_pairs)}
  
  if (length(pair_list) == 0) {
    message("No candidate pairs found in the requested distance range.")
    return(list())
  }
  start.end.pairs <- unique(do.call(rbind, pair_list))
  colnames(start.end.pairs) <- c("from", "to")
  
  message(sprintf("Candidate pairs found: %d", nrow(start.end.pairs)))
  if (nrow(start.end.pairs) == 0) return(list())
  
  # 5) Sample pairs (without replacement) up to requested n.paths
  m <- min(n.paths, nrow(start.end.pairs))
  samp_idx <- sample.int(nrow(start.end.pairs), m, replace = FALSE)
  tmp.pairs <- start.end.pairs[samp_idx, , drop = FALSE]
  
  # 6) Compute shortest paths for sampled pairs
  #    Keep only non-empty paths
  perm.paths.list=c()
  for (i in seq_len(m)) {
    sp <- igraph::shortest_paths(graph.perm, from = tmp.pairs[i, 1], to = tmp.pairs[i, 2])
    vv <- names(unlist(sp$vpath, use.names = TRUE))
    if (length(vv)>5) {perm.paths.list=c(perm.paths.list,list(vv))}
  }
  
  ll <- vapply(perm.paths.list, length, integer(1))
  message("Path length table:")
  print(table(ll))
  
  return(perm.paths.list)
}



# (2) Check which cluster path is from

check_path_cluster<-function(tmp.path,meta.cts){
  s.clus=meta.cts[tmp.path,'seurat_clusters']
  tmp.tab=table(s.clus)
  return(names(tmp.tab)[which(tmp.tab==max(tmp.tab))[1]])
}


# (3) Get a random path set based on length of each migration paths
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

# (4) Get all the sets
get_control_paths<-function(paths, # final trails
                            nbr.list, 
                            meta.cts,
                            s.range1,s.range2, # largest size and smallest size (distance) of permutated trails
                            n.paths, #total number of permutated trails to make
                            loc.cts,
                            n.set, #total number of control trail sets
                            path.size){ # size of each trail
  
  ###Get alternative paths of all kinds of lengths
  paths.nodes=unique(do.call(c,paths))
  path.clus=unlist(lapply(paths,function(x)check_path_cluster(x,meta.cts=meta.cts)))
  all.nodes=rownames(meta.cts)[which(meta.cts$seurat_clusters %in% unique(path.clus))]
  other.nodes=all.nodes[which(!all.nodes %in% paths.nodes)] #nodes that not on the path but in the same clusters
  
  tmp.graph=make_graph_perm(spots=other.nodes,nbs.list=nbs.list)
  perm.paths.list=get_permutated_paths(meta.cts=meta.cts,s.range1=s.range1,s.range2=s.range2,loc.cts = loc.cts,
                                       graph.perm=tmp.graph[[1]],n.paths=n.paths)
  
  
  #Length of real trails
  real.size=unlist(path.size)
  #Cluster of real trails
  path.clus=unlist(lapply(paths,function(x)check_path_cluster(x,meta.cts=meta.cts)))
  
  #Length of permurated trails
  str.end=lapply(perm.paths.list, function(x) c(x[1], x[length(x)]))
  perm.size=lapply(str.end,function(x) sqrt((loc.cts[x[1],1]-loc.cts[x[2],1])^2+(loc.cts[x[1],2]-loc.cts[x[2],2])^2))
  perm.size=unlist(perm.size)
  print(quantile(perm.size))
  
  #Cluster of permutated trails
  perm.path.clus=unlist(lapply(perm.paths.list,function(x)check_path_cluster(x,meta.cts=meta.cts)))
  
  perm.path.pool=c()
  for (i in 1:length(paths)) {
    tmp.clus=path.clus[i]
    tmp.size=real.size[i]
    #print(tmp.size)
    tmp.pool=perm.paths.list[which(perm.path.clus==tmp.clus)]
    pool.size=perm.size[which(perm.path.clus==tmp.clus)]
    tmp.pool=tmp.pool[which(pool.size>=tmp.size*0.8 & pool.size<=tmp.size*1.2)]
    print(length(tmp.pool))
    matched.pool=tmp.pool[sample(1:length(tmp.pool),size=n.set,replace = TRUE)]
    perm.path.pool=c(perm.path.pool,list(matched.pool))
  }
  perm.pathset.norm=c()
  for (k in 1:n.set) {
    tmp.set=lapply(perm.path.pool,function(x) x[[k]])
    perm.pathset.norm=c(perm.pathset.norm,list(tmp.set))
  }
  return(perm.pathset.norm)
}






####################################################
#Genes upregulated and downregulated along the path
###################################################



controlset.mean.fun<-function(genes,sp.dd){
  controlset.means=c()
  for (i in 1:length(perm.pathset.norm.all)){
    print(i)
    pathset=perm.pathset.norm.all[[i]]
    pathset.spots=lapply(pathset,function(x) unlist(sub.clusters[x]))
    pathset.genes=lapply(pathset.spots,function(x) sp.dd[x,genes])
    
    pathset.mean=do.call(rbind,lapply(pathset.genes,function(x) apply(x,2,mean)))
    set.mean=apply(pathset.mean, 2, mean)
    controlset.means=rbind(controlset.means,set.mean)
  }
  return(controlset.means)
}








### Get mean UMI of a path
get_path_umi<-function(tmp.path){
  tmp.spots=do.call(c,sub.clusters[tmp.path])
  mean.umi=mean(meta_dat[tmp.spots,'nCount_Spatial'])
  return(mean.umi)
}


### Get mean expression of a gene of some spots
get_spots_gene_mean<-function(dd,tmp.spots,gene){
  gene.mean=mean(dd[tmp.spots,gene])
  return(gene.mean)
}

### Get mean expression of a gene along a path
get_gene_mean<-function(dd,tmp.path,gene){
  tmp.spots=do.call(c,sub.clusters[tmp.path])
  mean(dd[tmp.spots,gene])
}

### Get mean expression of some genes of  one path set
get_gpt_gg1<-function(dd,pathlist,gg.markers) {
  tmp.list=lapply(pathlist,function(path) sapply(gg.markers,function(x) get_gene_mean(dd=dd,tmp.path=path,gene=x)))
  do.call(rbind,tmp.list)
}

### Get mean expression of some genes of many random control path sets, to get empirical distribution

get_perm_gg_mat1<-function(dd,perm.pathset.list,gg.markers,start.index,end.index) {
  perm.pathset.list=perm.pathset.list[start.index:end.index]
  lapply(perm.pathset.list,function(x) get_gpt_gg1(dd=dd,pathlist=x,gg.markers=gg.markers))
}

### show gene expression on a trail
show_dat<-function(tmp.path,gg){
  tmp.spots=do.call(c,sub.clusters[tmp.path])
  tmp.dat=sp.counts.r[tmp.spots,gg]
  
  tmp.meta=meta_dat_cell[tmp.spots,]
  tmp.perc=tmp.dat/tmp.meta[,'nCount_Spatial']
  print(tmp.dat)
  print(tmp.perc)
  #cell.sd=sd(tmp.meta[,'cell.sum'])
  #cell.mean=mean(tmp.meta[,'cell.sum'])
  return(c(cell.sd,cell.mean))
}


### Empirical distribution of mean gene expression of the random control path sets
emp_dist<-function(perm.mat.list,gg){
  pool=unlist(lapply(perm.mat.list,function(x) mean(x[,gg])))
  #hist(pool)
  return(pool)
}


### Get empirical P value
get_pvalue<-function(gg.mat,perm.mat.list){
  p.list=c()
  for (i in 1:ncol(gg.mat)) {
    mm.gpt=mean(gg.mat[,i])
    mm.perm=emp_dist(perm.mat.list,colnames(gg.mat)[i])
    pavlue_1sided=length(which(mm.perm>=mm.gpt))/length(perm.mat.list)
    p.list=c(p.list,pavlue_1sided)
  }
  names(p.list)=colnames(gg.mat)
  p.list=p.list[order(p.list)]
  return(p.list)
}

### Get FC
get_fold_change<-function(gg.mat,perm.mat.list){
  fc.list=c()
  for (i in 1:ncol(gg.mat)) {
    mm.gpt=mean(gg.mat[,i])
    mm.perm=emp_dist(perm.mat.list,colnames(gg.mat)[i])
    fc=mm.gpt/mean(mm.perm)
    fc.list=c(fc.list,fc)
  }
  names(fc.list)=colnames(gg.mat)
  return(fc.list)
}




### Compare path set mean to mean of all keep.ids
compare_to_all_keepids<-function(gg.markers){
  out=c()
  for (gg in gg.markers) {
    gpt.dat=unlist(lapply(gpt.norm,function(x) get_gene_mean(dd=sp.counts.tpm,tmp.path=x,gg)))
    all.spots.dat=sp.counts.tpm[keep.ids,gg]
    mms=c(median(all.spots.dat),median(gpt.dat))
    tmp.ttest=t.test(gpt.dat,sp.counts.tpm[keep.ids,gg])
    mms=tmp.ttest$estimate
    pv=tmp.ttest$p.value
    out=rbind(out,c(mms,pv))
  }
  rownames(out)=gg.markers
  colnames(out)=c('Mean TPM on real paths','Mean TPM in all Tcell spots','P values')
  return(out)
}


### Plot migration path mean compared with empirical path sets distribution

plot_emp<-function(gg,ypos,legend=F) {
  par(mar=c(4,4,1,1))
  gpt.dat=unlist(lapply(gpt.norm,function(x) get_gene_mean(dd=sp.counts.tpm,tmp.path=x,gg)))
  gpt.mean=mean(gpt.dat)
  emp.means=unlist(lapply(random.pathset.TMPmean.norm,function(x) mean(x[,gg])))
  pvalue1=length(which(gpt.mean<=emp.means))/length(emp.means)
  pvalue2=length(which(gpt.mean>=emp.means))/length(emp.means)
  pvalue1=round(pvalue1,5)
  pvalue2=round(pvalue2,5)
  pvalue=min(pvalue1,pvalue2)
  hist(emp.means,breaks=100,xlim=c(min(gpt.mean*0.9,min(emp.means)),max(emp.means*1.1,gpt.mean)),main='Distribution of mean TPM among random path sets',xlab=paste0('TPM of ',gg),ylab='Frequency')
  abline(v=gpt.mean,col='red')
  if (legend==T) {legend(x=max(emp.means)*0.8,y=ypos+20,legend=': Mean TMP along real paths',col='red',lty=1,cex=0.8,bty = 'n')}
  text(max(gpt.mean,max(emp.means))*0.95,ypos,labels=paste('P value:',pvalue),cex=0.8)
}




###Get gene trend along a trail

# Boxplot by groups
group_boxplot<-function(dat,gg){
  p= ggplot(dat, aes(group, eval(parse(text=gg)), fill = eval(parse(text=gg)))) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    guides(fill = "none") +
    labs(x = "Path", y = gg)
  p
}


# Paneled scatter plot
fig1<-function(gg,path.list){
  par(mfrow=c(5,7),mar=c(2,2,1,1))
  for (i in 1:length(gpt.norm)) {
    tmp.spots=do.call(c,sub.clusters[path.list[[i]]])
    x=1:length(tmp.spots)
    y=sp.counts.tpm[tmp.spots,gg]
    plot(y=y,x=x)
  }
}


# Get long-format data
long_dat<-function(path.list,trt,dd,gg.markers,method='median'){ #trt=1, if gpt paths; trt=0, if random paths
  dat=c()
  for (i in 1:length(path.list)){
    path.spots=sub.clusters[path.list[[i]]]
    x=1:length(path.spots)
    len.list=lapply(path.spots,function(x) length(x))
    y.list=lapply(path.spots,function(x) dd[x,gg.markers])
    y.list <- lapply(1:length(y.list), function(x) {
      obj=y.list[[x]]
      if (len.list[x]==1) c(x,meta_dat[path.spots[[x]],'ex.score'],obj) else cbind(rep(x,len.list[x]),meta_dat[path.spots[[x]],'ex.score'],obj)
    })
    
    y.list.median <- lapply(1:length(y.list), function(x) {
      obj=y.list[[x]]
      if (len.list[x]==1) obj else apply(obj,2,median)
    })
    
    y.list.mean <- lapply(1:length(y.list), function(x) {
      obj=y.list[[x]]
      if (len.list[x]==1) obj else apply(obj,2,mean)
    })
    
    if (method=='mean')
    y=do.call(rbind,y.list.mean)
    
    if (method=='median')
      y=do.call(rbind,y.list.median)

    colnames(y)[1]='ex_score'
    
    if (trt==1) {
      path_label=paste0('Trail',i)
      path_id=rep(path_label,length(path.spots))
    }
    
    if (trt==0) {
      path_label=paste0('random_',i)
      path_id=rep(path_label,length(path.spots))
    }
    
    tmp.dat=data.frame(path_id,x,y)
    
    dat=rbind(dat,tmp.dat)
  }
  trt=rep(trt,nrow(dat))
  dat=data.frame(dat,trt)
  dat$x=as.factor(dat$x)
  return(dat)
}



# Mixed linear model
mlm_pvalue<-function(dat,gg.markers) {
  p.list=c()
  beta.list=c()
  gg.list=c()
  for (i in gg.markers) {
    if (length(which(dat[,i]==0))/nrow(dat) > 0.7) {next}
    else{
      print(i)
      gg.list=c(gg.list,i)
      m1=lmer(eval(parse(text=i)) ~ as.numeric(ex_score) + as.numeric(x)  +  (1 | path_id), data =dat)
      tmp=summary(m1)
      print(tmp)
      beta=tmp$coefficients[2,1]
      p.value=pf(anova(m1)$`F value`, 1, (nrow(dat)-1), lower.tail = FALSE)[2]
      print(p.value)
      p.list=c(p.list,p.value)
      beta.list=c(beta.list,beta)
    }
  }
  names(p.list)=gg.list
  names(beta.list)=gg.list
  p.list.adj=p.adjust(p.list,'BY')
  p.list.adj=p.list.adj[order(p.list.adj)]
  beta.list=beta.list[names(p.list.adj)]
  return(c(list(beta.list),list(p.list.adj)))
}




mlm_pvalue1<-function(dat,gg.markers) {
  p.list=c()
  beta.list=c()
  gg.list=c()
  for (i in gg.markers) {
    if (length(which(dat[,i]==0))/nrow(dat) > 0.7) {next}
    else{
      print(i)
      gg.list=c(gg.list,i)
      m1=lmer(eval(parse(text=i)) ~ as.numeric(ex_score) +  (1 | path_id), data =dat)
      residuals_m1 <- residuals(m1)
      m2=lmer(residuals_m1 ~  as.numeric(x)  +  (1 | path_id), data =dat)
      
      
      tmp=summary(m2)
      #print(tmp)
      beta=tmp$coefficients[2,1]
      p.value=pf(anova(m2)$`F value`, 1, (nrow(dat)-1), lower.tail = FALSE)
      print(p.value)
      p.list=c(p.list,p.value)
      beta.list=c(beta.list,beta)
    }
  }
  names(p.list)=gg.list
  names(beta.list)=gg.list
  p.list.adj=p.adjust(p.list,'BY')
  p.list.adj=p.list.adj[order(p.list.adj)]
  beta.list=beta.list[names(p.list.adj)]
  return(c(list(beta.list),list(p.list.adj)))
}



# Linear model
lm_pvalue<-function(dat,gg.markers) {
  p.list=c()
  r.list=c()
  beta.list=c()
  for (i in gg.markers) {
    print(i)
    m1=lm(eval(parse(text=i)) ~ as.numeric(x), data =dat)
    tmp=summary(m1)
    p.value=tmp$coefficients[2,4]
    r.square=tmp$r.squared
    beta=tmp$coefficients[2,1]
    beta.list=c(beta.list,beta)
    p.list=c(p.list,p.value)
    r.list=c(r.list,r.square)
  }
  gg.markers=gsub('\\.','\\-',gg.markers)
  names(p.list)=gg.markers
  names(r.list)=gg.markers
  names(beta.list)=gg.markers
  p.list.adj=p.adjust(p.list,'BY')
  p.list.adj=p.list.adj[order(p.list.adj)]
  beta.list=beta.list[names(p.list.adj)]
  r.list=r.list[names(p.list.adj)]
  return(c(list(beta.list),list(p.list.adj),list(r.list)))
}


# Report genes with trend
trend_gg<-function(longdat,vgg.rename=vgg.rename){
  p.norm.mlm=mlm_pvalue(dat=longdat,gg.markers=vgg.rename)
  trend.norm.adjpvalues=p.norm.mlm[[2]]
  gg.trend.norm=names(trend.norm.adjpvalues)[which(trend.norm.adjpvalues<0.05)]
  
  beta.norm=p.norm.mlm[[1]][gg.trend.norm]
  trend.norm.up=gg.trend.norm[which(beta.norm>0)]
  trend.norm.down=gg.trend.norm[which(beta.norm<0)]
  return(list(trend.norm.up,trend.norm.down))
}



### Write barcodes to a file to import into Loupe

#For one path
get_rm_barcodes<-function(tmp.path){
  tmp.spots=do.call(c,sub.clusters[tmp.path])
  tmp.others=rownames(meta_dat)[which(!rownames(meta_dat) %in% tmp.spots)]
  write.table(tmp.others,file='barcodes.csv',sep=',',row.names = F,col.names = F)
}


#For a path list
get_rm_barcodes<-function(path.list,fname){
  tmp.keep=c()
  for (i in 1:length(path.list)) {
    tmp.path=path.list[[i]]
    tmp.spots=do.call(c,sub.clusters[tmp.path])
    tmp.keep=c(tmp.keep,tmp.spots)
  }
  tmp.keep=unique(tmp.keep)
  tmp.others=rownames(loc.raw)[which(!rownames(loc.raw) %in% tmp.keep)]
  write.table(tmp.others,file=fname,sep=',',row.names = F,col.names = F)
}


get_barcodes_loc<-function(tmp.path){
  tmp.spots=do.call(c,sub.clusters[tmp.path])
  out=loc.raw[tmp.spots,]
  out=data.frame(rownames(out),out[,2],out[,1])
  colnames(out)=c('Barcode','X Coordinate','Y Coordinate')
  write.csv(out,file='barcodes.csv',row.names = F)
}


###Check trails that have a big overlap

check_trails_overlap<-function(gpt.norm){
  ov.mat=c()
  for (i in 1:length(gpt.norm)) {
    tmp.path=gpt.norm[[i]]
    ovlp=c()
    for (j in 1:length(gpt.norm)) {
      tmp.p=gpt.norm[[j]]
      ovlp=c(ovlp,length(intersect(tmp.path,tmp.p)))
    }
    ov.mat=rbind(ov.mat,ovlp)
  }
  
  diag(ov.mat)=0
  ov.mat[upper.tri(ov.mat)]=0
  
  for (i in 1:28) {
    for (j in 1:28) {
      if (ov.mat[i,j]>2) {
        print(c(i,j))
      }
    }
  }
}

### Plot trails
par(mfrow=c(1,1))

###Figure 2

plot_trails<-function(paths){
  ex.score=round(meta.cts[unlist(paths),'pt'])
  colors=colfunc(round(max(ex.score) -  min(ex.score))+1)
  
  plot(loc.cts,pch=19,col='grey',xlab='',ylab='')
  points(meta.cts[which(meta.cts$seurat_clusters %in% tumor.clus),c('x','y')],col=alpha('red',0.4))
  for (i in 1:length(paths)) {
    ex.score=meta.cts[unlist(paths[[i]]),'pt']-4
    cols=colors[round(ex.score)]
    points(loc.cts[paths[[i]],],col=cols,pch=19,cex=0.8)
    lines(loc.cts[paths[[i]],],pch=19,cex=0.6,lwd=2,col=alpha('black',0.5))
  }
}




vocano_plot1<-function(pathset.mean,controlset.means,plot,pvalue.cut=0.1,fc.cut.up=1,fc.cut.down=1) {
  null.means=apply(controlset.means,2,mean) #Mean of the mean of 10000 controlsets
  fc=pathset.mean/null.means #fold change of pathset mean
  
  #Empirical pvalues
  pathset.raw.pvalues=sapply(1:length(pathset.mean),function(x) length(which(controlset.means[,x]>pathset.mean[x])))/nrow(controlset.means)
  pathset.raw.pvalues[which(fc<=1)]=1-pathset.raw.pvalues[which(fc<=1)]
  names(pathset.raw.pvalues)=names(pathset.mean)
  
  #Adjusted p values
  pathset.pvalues=p.adjust(pathset.raw.pvalues,method='BH')
  names(pathset.pvalues)=names(pathset.raw.pvalues)
  names(fc)=names(pathset.raw.pvalues)
  
  output.pathset=cbind(pathset.pvalues,fc,pathset.raw.pvalues)
  colnames(output.pathset)=c('pvalue','fc','raw.pvalue')
  
  output.up=output.pathset[which(output.pathset[,'fc']>fc.cut.up & output.pathset[,'pvalue']<pvalue.cut),]
  output.down=output.pathset[which(output.pathset[,'fc']<fc.cut.down & output.pathset[,'pvalue']<pvalue.cut),]
  
  gg.up.pathset=rownames(output.up)
  gg.down.pathset=rownames(output.down)
  
  
  
  #Data for vacano plot
  
  pvalues=output.pathset[,'pvalue']
  fc=output.pathset[,'fc']
  
  pvalues[which(pvalues==0)]=1e-3
  df=data.frame(-log(pvalues,10),log(fc,2))
  colnames(df)=c('Pvalue','fc')
  sig.gg=rep(1,nrow(df))
  sig.gg[which(df$Pvalue>=-log10(0.05) | df$fc>log2(1.5)|df$fc<log2(1/1.5))]=3
  sig.gg[which((df$Pvalue>=-log10(0.05) & df$fc>log2(1.5)) | (df$Pvalue>=-log10(0.05) & df$fc<log2(1/1.5)) )]=2
  table(sig.gg)
  df=data.frame(df,sig.gg)
  colnames(df)[3]='sig_genes'
  df$genes=rownames(df)
  df$genes[which(df$sig_genes!=2)]=''
  
  if (plot==TRUE) {
    p=ggplot(data = df, aes(x = fc, y = Pvalue, label = genes)) +
      geom_vline(xintercept = c(0, log2(1.5), log2(1/1.5)), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = 'dashed') + 
      geom_point(size = 0.3, col = df$sig_genes) +
      coord_cartesian(ylim = c(0, 4), xlim = c(-2, 2)) + # setting limits
      labs(color = 'Significance', 
           x = expression("log"[2]*"FC"), 
           y = expression("-log"[10]*"p-value")) + 
      scale_x_continuous(breaks = seq(-10, 10, 2)) + # customizing x-axis breaks
      ggtitle('') + # Plot title 
      geom_text_repel(max.overlaps = 100, size = 3) + # To show all labels
      theme(
        axis.title.x = element_text(size = 15), # Enlarging x-axis label
        axis.title.y = element_text(size = 15),  # Enlarging y-axis label
        axis.text.x = element_text(size = 15),  # Enlarging x-axis tick marks
        axis.text.y = element_text(size = 15)   # Enlarging y-axis tick marks
      )
    print(p)
  }

  return(list(gg.up.pathset,gg.down.pathset,output.pathset))
}


enrich_barplot<-function(gg.up.pathset.enrich,range){
  enrich.out=gg.up.pathset.enrich@result[range,]
  enrich.out[,'ID'] = gsub('_', ' ', enrich.out[,'ID'])
  enrich.out[,'ID']=tolower(enrich.out[,'ID'])
  df = data.frame(enrich.out[range, c('Count', 'ID', 'p.adjust')])
  colnames(df) = c('counts', 'group', 'p_values')
  df$group=tolower(df$group)
  df[,'group'] = gsub('_', ' ', df[,'group'])
  
  df$group <- reorder(df$group, -df$p_values)
  
  p=ggplot(df, aes(x = counts, y = group, fill = p_values)) +
    geom_bar(stat = "identity") + # Create the bars
    scale_fill_gradient(low = "red", high = "blue", 
                        guide = guide_colorbar(reverse = TRUE)) + # Reverse the color gradient in the legend
    labs(
      x = "Count",
      y = "",
      fill = "p-value"
    ) + 
    theme_minimal() + # Use a minimal theme
    theme(
      axis.text.x = element_text(size = 14), # Enlarge x-axis text
      axis.text.y = element_text(size = 13), # Enlarge y-axis text
      axis.title.x = element_text(size = 16), # Enlarge x-axis label
      axis.title.y = element_text(size = 16), # Enlarge y-axis label
      legend.text = element_text(size = 12),  # Enlarge legend bar labels
      legend.title = element_text(size = 14)  # Enlarge legend title
    )
  print(p)
  return(df)
}



cluster_enrich<-function(gg.up.pathset){
  ggs=gg.up.pathset[which(gg.up.pathset %in% rownames(norm.rna.mat))]
  norm.rna.mat=norm.rna.mat[ggs,]
  norm.rna.mat=t(as.matrix(norm.rna.mat))
  
  ggs.mean.mat=c()
  for (i in ggs) {
    ggs.mean.mat=rbind(ggs.mean.mat,aggregate(norm.rna.mat[,i] ~ clusters, data=norm.rna.mat, FUN=mean)[,2])
  }
  
  rownames(ggs.mean.mat)=ggs
  colnames(ggs.mean.mat)=paste0('Cluster',0:n.cluster)
  n.cluster=length(unique(clusters))-1
  
  
  ##let x denote gene and y denote cluster
  tmp.loc=expand.grid(1:length(ggs),1:(n.cluster+1))
  colnames(tmp.loc)=c('gene','cluster')
  value=ggs.mean.mat[as.matrix(tmp.loc)]
  
  rank.mat=t(apply(ggs.mean.mat,1,function(x) rank(as.numeric(x))))
  rank.mat[which(rank.mat<7)]=7
  rownames(rank.mat)=ggs
  colnames(rank.mat)=colnames(ggs.mean.mat)
  rank.value=rank.mat[as.matrix(tmp.loc)]
  
  
  rank.mean=apply(rank.mat,2,mean)
  rank.mean=rank.mean+rnorm(n=length(rank.mean),mean=0,sd=0.0001)
  reorder=rank(rank.mean)
  names(reorder)=paste0('Cluster',0:12)
  
  #Order columns by average rank
  rank.mat.new=rank.mat[,names(reorder)[order(reorder)]]
  
  table(rank.mat.new[,ncol(rank.mat.new)])
  table(rank.mat.new[,ncol(rank.mat.new)-1])
  table(rank.mat.new[,1])
  
  ggs.mean.mat.new=ggs.mean.mat[,order(rank.mean)]
  rank.value.new=rank.mat.new[as.matrix(tmp.loc)]
  
  #Get the most enriched cluster's fc
  enrich.clus=which(rank.mean==max(rank.mean))-1
  fc.list=c()
  for (i in ggs) {
    tmp.dd1=norm.rna.mat[which(clusters==enrich.clus),i]
    tmp.dd2=norm.rna.mat[which(clusters!=enrich.clus),i]
    fc.list=c(fc.list,mean(tmp.dd1)/mean(tmp.dd2))
  }
  names(fc.list)=ggs
  ggs.new=names(fc.list[order(fc.list)])
  fc.list[ggs.new]
  
  
  #rank.mat.new=rank.mat.new[ggs.new,]
  rank.value.new=rank.mat.new[as.matrix(tmp.loc)]
  
  tmp.loc$cluster=factor(tmp.loc$cluster)
  tmp.loc$gene=factor(tmp.loc$gene)
  #Heatmap
  p <- ggplot(tmp.loc, aes(y = cluster, x = gene, size = rank.value.new, fill = rank.value.new)) +
    geom_point(shape = 21, color = "black", alpha = 0.9) +  # Add a border to the dots
    scale_size_continuous(range = c(.1, 1.5), breaks = seq(0, 1, length = 6)) +  # Customize size range and breaks
    scale_fill_gradient(low = "blue", high = "red",labels = function(x) ifelse(x == 7, "<=7", x)) +
    labs(title = "",
         y = "Cluster",
         x = "Gene",
         size = "Rank",
         fill = "Rank") +
    theme_minimal() +
    theme(#axis.text.x = element_text(angle = 90, hjust = 1,size=8.5),  # Rotate x-axis labels
      axis.text.x = element_blank(), 
      axis.text.y = element_text(size=13),
      axis.title.x = element_text(size = 16),  # Enlarge x-axis label
      axis.title.y = element_text(size = 16),  # Enlarge y-axis label
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)) +
    scale_y_discrete(breaks=1:13,labels = colnames(rank.mat.new)) +
    scale_x_discrete(breaks=1:length(ggs.new),labels=ggs.new)
  # Customize x-axis ticks (if needed)
  # Print the customized plot
  print(p)
  print(fc.list[ggs.new])
}




#vocano plot of single cell data
sc_vacanoplot<-function(clus.markers){
  pvalues=clus.markers$p_val_adj
  pvalues[pvalues<=1e-100]=1e-100
  fc=clus.markers$avg_log2FC
  
  df=data.frame(-log(pvalues,10),fc)
  colnames(df)=c('Pvalue','fc')
  sig.gg=rep(1,nrow(df))
  sig.gg[which(df$Pvalue>=90 | df$fc>log2(10)|df$fc<log2(1/50))]=3
  sig.gg[which((df$Pvalue>=90 & df$fc>log2(4)) | (df$Pvalue>=90 & df$fc<log2(1/50)) )]=2
  table(sig.gg)
  df=data.frame(df,sig.gg)
  colnames(df)[3]='sig_genes'
  rownames(df)=rownames(clus.markers)
  df$genes=rownames(df)
  df$genes[which(df$sig_genes!=2)]=''
  
  
  ggplot(data = df, aes(x = fc, y = Pvalue,label=genes)) +
    geom_vline(xintercept = c(0, log2(10),log2(1/50)), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = 'dashed') + 
    geom_point(size = 0.3,col=df$sig_genes) +
    coord_cartesian(ylim = c(0, 210), xlim = c(-6, 15)) + # since some genes can have minuslog10padj of inf, we set these limits
    labs(color = 'Significance', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
    scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
    ggtitle('Mean expression on potential migration trail vs. matched random trails') + # Plot title 
    geom_text_repel(max.overlaps = 100,size=2.5) # To show all labels 
}


par_cor<-function(path.list,dd,gg.markers){
  cor.mat=c()
  for (i in 1:length(path.list)){
    path=path.list[[i]]
    tmp.spots=do.call(c,sub.clusters[path]) #This is by spots
    x=1:length(tmp.spots)
    z=meta_dat[tmp.spots,'ex.score']
    cor.list=c()
    for (j in gg.markers) {
      y=dd[tmp.spots,j]
      sample_data <- data.frame(x,y,z)
      print(sample_data)
      print(pcor( sample_data,method='spearman' ))
      print('==============================')
      tmp.cor=pcor( sample_data ,method='spearman')$estimate [1,2]
      cor.list=c(cor.list,tmp.cor)
    }
    cor.mat=rbind(cor.mat,cor.list)
  }
  colnames(cor.mat)=gg.markers
  rownames(cor.mat)=paste('trail',1:length(path.list),sep='_')
  return(cor.mat)
}


#Spagetti_plot in panels
spagetti_plot<-function(gg,dat,pathset){
  for (i in 1:length(pathset)) {
    path=pathset[[i]]
    tmp.spots=unlist(sub.clusters[path])
    tmp.dd=dat[tmp.spots,gg]
    x=1:length(tmp.dd)
    m1=lm(tmp.dd~x)
    rr=summary(m1)$r.squared
    beta=summary(m1)$coef[2,1]
    plot(y=tmp.dd,x=1:length(tmp.dd),type='l')
  }
}

