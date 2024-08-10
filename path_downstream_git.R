library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
require(lme4)
require(multcomp)
require(ggrepel)
require(ggplot2)


################Path validation ############################################################################################################################################

#####Descriptive analysis

### Get mean UMI of a path
get_path_umi<-function(tmp.path){
  tmp.spots=do.call(c,sub.clusters[tmp.path])
  mean.umi=mean(meta_dat[tmp.spots,'nCount_Spatial'])
  return(mean.umi)
}

### Get mean expression of a gene of a few spots
get_spots_gene_mean<-function(dd,tmp.spots,gene){
  gene.mean=mean(dd[tmp.spots,gene])
  return(gene.mean)
}

### Get mean expression of a gene along a path
get_gene_mean<-function(dd,tmp.path,gene){
  tmp.spots=do.call(c,sub.clusters[tmp.path])
  mean(dd[tmp.spots,gene])
}

### Get mean expression of some genes of one path set
get_gpt_gg1<-function(dd,pathlist,gg.markers) {
  tmp.list=lapply(pathlist,function(path) sapply(gg.markers,function(x) get_gene_mean(dd=dd,tmp.path=path,gene=x)))
  do.call(rbind,tmp.list)
}

### Get mean expression of some genes of many random control path sets, to get empirical distribution

get_perm_gg_mat1<-function(dd,perm.pathset.list,gg.markers,start.index,end.index) {
  perm.pathset.list=perm.pathset.list[start.index:end.index]
  lapply(perm.pathset.list,function(x) get_gpt_gg1(dd=dd,pathlist=x,gg.markers=gg.markers))
}

### Show gene expression of a path
show_dat<-function(tmp.path,gg){
  tmp.spots=do.call(c,sub.clusters[tmp.path])
  tmp.dat=sp.counts.r[tmp.spots,gg]
  
  tmp.meta=meta_dat_cell[tmp.spots,]
  tmp.perc=tmp.dat/tmp.meta[,'nCount_Spatial']
  print(tmp.dat)
  print(tmp.perc)
  cell.sd=sd(tmp.meta[,'cell.sum'])
  cell.mean=mean(tmp.meta[,'cell.sum'])
  return(c(cell.sd,cell.mean))
}


### Empirical distribution of mean gene expression of the random control path sets
emp_dist<-function(perm.mat.list,gg){
  pool=unlist(lapply(perm.mat.list,function(x) mean(x[,gg])))
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
compare_to_all_keepids<-function(paths=gpt.norm,gg.markers){
  out=c()
  for (gg in gg.markers) {
    gpt.dat=unlist(lapply(paths,function(x) get_gene_mean(dd=sp.counts.tpm,tmp.path=x,gg)))
    all.spots.dat=sp.counts.tpm[keep.ids,gg]
    mms=c(median(all.spots.dat),median(gpt.dat))
    tmp.ttest=t.test(gpt.dat,sp.counts.tpm[keep.ids,gg])
    mms=tmp.ttest$estimate
    pv=tmp.ttest$p.value
    out=rbind(out,c(mms,pv))
  }
  rownames(out)=gg.markers
  colnames(out)=c('Mean TPM on trails','Mean TPM in all Tcell spots','P values')
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



##### Gene trend along the path


### Boxplot by groups
group_boxplot<-function(dat,gg){
  p= ggplot(dat, aes(group, eval(parse(text=gg)), fill = eval(parse(text=gg)))) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    guides(fill = "none") +
    labs(x = "Path", y = gg)
  p
}


### Paneled scatter plot
fig1<-function(gg,path.list){
  par(mfrow=c(5,7),mar=c(2,2,1,1))
  for (i in 1:length(gpt.norm)) {
    tmp.spots=do.call(c,sub.clusters[path.list[[i]]])
    x=1:length(tmp.spots)
    y=sp.counts.tpm[tmp.spots,gg]
    plot(y=y,x=x)
  }
}

### Get long-format data
long_dat<-function(path.list,trt,dd,gg.markers){ #trt=1, if gpt paths; trt=0, if random paths
  dat=c()
  for (i in 1:length(path.list)){
    path=path.list[[i]]
    tmp.spots=do.call(c,sub.clusters[path]) #This is by spots
    x=1:length(tmp.spots)
    y=dd[tmp.spots,gg.markers]
    
    if (trt==1) {
      path_label=paste0('Trail',i)
      path_id=rep(path_label,length(tmp.spots))
    }
    
    if (trt==0) {
      path_label=paste0('random_',i)
      path_id=rep(path_label,length(tmp.spots))
    }
    tmp.dat=data.frame(path_id,x,y)
    
    dat=rbind(dat,tmp.dat)
  }
  trt=rep(trt,nrow(dat))
  dat=data.frame(dat,trt)
  dat$x=as.factor(dat$x)
  return(dat)
}

### Mixed linear model to test trend
mlm_pvalue<-function(dat,gg.markers) {
  p.list=c()
  beta.list=c()
  gg.list=c()
  for (i in gg.markers) {
    print(i)
    if (length(which(dat[,i]==0))/nrow(dat) > 0.7) {next}
    gg.list=c(gg.list,i)
    m1=lmer(eval(parse(text=i)) ~ as.numeric(x) + (1 | path_id), data =dat)
    tmp=summary(m1)
    beta=tmp$coefficients[2,1]
    p.value=pf(anova(m1)$`F value`, 1, (nrow(dat)-1), lower.tail = FALSE)
    p.list=c(p.list,p.value)
    beta.list=c(beta.list,beta)
  }
  names(p.list)=gg.list
  names(beta.list)=gg.list
  p.list.adj=p.adjust(p.list,'BY')
  p.list.adj=p.list.adj[order(p.list.adj)]
  beta.list=beta.list[names(p.list.adj)]
  return(c(list(beta.list),list(p.list.adj)))
}

### Linear model to test trend
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


### Genes with trend

trend_gg<-function(longdat,vgg.rename=vgg.rename){
  p.norm.mlm=mlm_pvalue(dat=longdat,gg.markers=vgg.rename)
  trend.norm.adjpvalues=p.norm.mlm[[2]]
  gg.trend.norm=names(trend.norm.adjpvalues)[which(trend.norm.adjpvalues<0.05)]
  
  beta.norm=p.norm.mlm[[1]][gg.trend.norm]
  trend.norm.up=gg.trend.norm[which(beta.norm>0)]
  trend.norm.down=gg.trend.norm[which(beta.norm<0)]
  return(list(trend.norm.up,trend.norm.down))
}



### Write barcodes of trails to a file to import into Loupe
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