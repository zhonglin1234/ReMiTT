library(ggrepel)
library(msigdbr)
library(clusterProfiler)
library(StereoMorph)
require(Seurat)
require(ggplot2)
require(igraph)
library(plotly)
library(RColorBrewer)
require(openxlsx)
require(LearnGeom)
#require(lme4)
library(gplots)



##### Vocano plot
vocano_plot1<-function(pathset.mean,controlset.means,plot) {
  null.means=apply(controlset.means,2,mean) #Mean of the mean of 10000 controlsets
  fc=pathset.mean/null.means #fold change of pathset mean
  
  #Empirical pvalues
  pathset.raw.pvalues=sapply(1:length(pathset.mean),function(x) length(which(controlset.means[,x]>pathset.mean[x])))/nrow(controlset.means)
  pathset.raw.pvalues[which(fc<=1)]=1-pathset.raw.pvalues[which(fc<=1)]
  names(pathset.raw.pvalues)=names(pathset.mean)
  
  #Adjusted p values
  pathset.pvalues=p.adjust(pathset.raw.pvalues,method='BH')
  names(pathset.pvalues)=vgg
  names(fc)=vgg
  
  output.pathset=cbind(pathset.pvalues,fc)
  colnames(output.pathset)=c('pvalue','fc')
  
  output.up=output.pathset[which(output.pathset[,'fc']>1 & output.pathset[,'pvalue']<0.05),]
  output.down=output.pathset[which(output.pathset[,'fc']<1 & output.pathset[,'pvalue']<0.05),]
  
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


##### Barplot of GSEA

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






##### Bubble plot to check T clusters enriched with trail markers
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
  print('FC of mean in the most enriched cluster compared to the mean of other clusters')
  print(fc.list[ggs.new])
}




