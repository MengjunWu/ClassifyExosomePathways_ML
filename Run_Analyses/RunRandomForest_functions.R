suppressMessages(library(caret))
suppressMessages(library(optparse))
suppressMessages(library(matrixStats))
suppressMessages(library(doParallel))
suppressMessages(library(ROCR))

library(grid)
neatload <- function(dataname){
  temp.space <- new.env()
  bar <- load(dataname, temp.space)
  the.object <- get(bar, temp.space)
  rm(temp.space)
  the.object
}

CalcMetrics <- function(prediction, truth){
  confusionTable <- table(prediction, truth)
  res <- c()
  prec <- posPredValue(confusionTable)
  rec <- sensitivity(confusionTable)
  F1score <- (2 * prec * rec) / (prec + rec)
  res <- data.frame(precision=prec,
                    recall=rec,
                    F1.score=F1score)
  res
}

calc_average_metric <- function(metrics.files){
  testres.mean <- c()
  testres.sd <- c()
  testres <- c()
  for (i in 1: length(metrics.files)){
    res.tmp <- read.table(metrics.files[i], header = T)
    testres <- rbind(testres, res.tmp)
  }
  test.metrics.mean <- colMeans(as.matrix(testres))
  test.metrics.sd<- colSds(as.matrix(testres))
  test.metrics.res <- data.frame(test.metrics.mean)
  test.metrics.res$SD <- test.metrics.sd
  colnames(test.metrics.res) <- c("mean", "sd")
  
  test.metrics.res
}

calc_average_importance <- function(importance.files){
  importance.average <- c()
  res1 <- read.table(importance.files[1], header = T)
  order.kmer <- row.names(res1)
  for(i in 1:length(importance.files)){
    res.tmp <- read.table(importance.files[i], header = T)
    res.tmp <- res.tmp[order.kmer, ]
    res.tmp2 <- res.tmp[,1]
    importance.average <- cbind(importance.average, res.tmp2)
  }
  
  importance.res <- rowMeans(importance.average)
  res <- data.frame(importance.score=importance.res)
  row.names(res) <- order.kmer
  res
}

calc_average_importance_multi <- function(importance.files){
  importance.average1 <- c()
  importance.average2 <- c()
  importance.average3 <- c()
  res1 <- read.table(importance.files[1], header = T)
  order.kmer <- row.names(res1)
  for(i in 1:length(importance.files)){
    res.tmp <- read.table(importance.files[i], header = T)
    res.tmp <- res.tmp[order.kmer, ]
    importance.average1 <- cbind(importance.average1, res.tmp[,1])
    importance.average2 <- cbind(importance.average2, res.tmp[,2])
    importance.average3 <- cbind(importance.average3, res.tmp[,3])
  }
  
  importance.res1 <- rowMeans(importance.average1)
  importance.res2 <- rowMeans(importance.average2)
  importance.res3 <- rowMeans(importance.average3)
  res <- data.frame(insensitive.score=importance.res1,
                    z1.score=importance.res2,
                    z8.score=importance.res3)
  row.names(res) <- order.kmer
  res
}


rankbyAverage <-function(data){
  data$average <- rowMeans(data)
  data <- data[order(data$average, decreasing = T), ]
  data <- data[,!colnames(data) %in% "average"]
  data
}
rankbyZ1 <-function(data){
  data <- data[order(data$Z1.importance, decreasing = T), ]
  data
}
rankbyZ8 <-function(data){
  data <- data[order(data$Z8.importance, decreasing = T), ]
  data
}

plotImportance_class14 <- function(datalist, class, zscore, color="grey42"){
  i <- class
  for (j in names(datalist)){
    data <- datalist[[j]]
    data.tmp <- data[data$group==i,]
    if(class=="class1.TSS"){
      data.tmp$Features <- unlist(lapply(data.tmp$feature, function(x)strsplit(x, split="TSS_")[[1]][2]))
    }
    if(class=="class4"){
      data.tmp$Features <- unlist(lapply(data.tmp$feature, function(x)strsplit(x, split="TES_")[[1]][2]))
    }
    #data.tmp$Features <- data.tmp$feature
    data.tmp <- data.tmp[, c("importance.score", "Features")]
    data.tmp <- data.tmp[data.tmp$importance.score > zscore, ]
    data.tmp <- data.tmp[order(data.tmp$importance.score, decreasing = T), ][1:min(nrow(data.tmp), 15),]
    data.tmp <- data.tmp[order(data.tmp$importance.score, decreasing = F), ]
    data.tmp$Features <- factor(data.tmp$Features, levels = data.tmp$Features)
    p <- ggplot(data=data.tmp, 
                aes(x=Features,y=importance.score)) + 
      geom_bar(position="dodge",stat="identity", width = 0.7, fill=color) + 
      coord_flip() +  
      sensitivity.genomicfeature.mytheme2 + ylab("importance score") + 
      xlab("Features") 
    pdf(paste0(outplot, j, "_",i, "_importance.score.pdf"),width=6, height=3)
    print(p)
    dev.off()
    print(p)
  }
}


plotImportance_class_all <- function(datalist, class,zscore){
  i <- class
  for (j in names(datalist)){
    data <- datalist[[j]]
    data.tmp <- data[data$group==i,]
    data.tmp$Features <- data.tmp$feature
    data.tmp <- data.tmp[, c("importance.score", "Features")]
    data.tmp <- data.tmp[data.tmp$importance.score >zscore, ]
    data.tmp <- data.tmp[order(data.tmp$importance.score, decreasing = T), ][1:min(nrow(data.tmp),15),]
    data.tmp <- data.tmp[order(data.tmp$importance.score, decreasing = F), ]
    data.tmp$Features <- factor(data.tmp$Features, levels = data.tmp$Features)
    p <- ggplot(data=data.tmp, 
                aes(x=Features,y=importance.score)) + geom_bar(position="dodge",stat="identity", width = 0.7, fill="grey42") + 
      coord_flip() +  
      mytheme + ylab("importance score") + 
      xlab("Features") 
    pdf(paste0(outplot, j, "_",i, "_importance.score.pdf"),width=7, height=5)
    print(p)
    dev.off()
  }
}
plotImportance_class23 <- function(datalist, class,zscore, color="grey42"){
  i <- class
  for (j in names(datalist)){
    print(j)
    data <- datalist[[j]]
    data.tmp <- data[data$group==i,]
    if(class=="class2.TSS"){
      data.tmp$Features <- unlist(lapply(data.tmp$feature, function(x)strsplit(x, split="class2.TSS_")[[1]][2]))
    }
    if(class=="class3"){
      data.tmp$Features <- unlist(lapply(data.tmp$feature, function(x)strsplit(x, split="class3.TES_")[[1]][2]))
    }
    print(data.tmp[1:3,])
    data.tmp$Motif_ID <- data.tmp$Features
    data.tmp2 <- merge(data.tmp, RBP.motifID_name, by="Motif_ID")
    if(!dim(data.tmp2)[1] == 0){
      data.tmp2$Features <- paste0(data.tmp2$RBP_name,",",data.tmp2$Motif_ID)
      data.tmp3 <- data.tmp[!data.tmp$Features %in% data.tmp2$Motif_ID,][,c("importance.score", "Features")]
      data.tmp4 <- data.tmp2[,c("importance.score", "Features")]
      data.tmp <- rbind(data.tmp3, data.tmp4)
      
    }else{
      data.tmp <- data.tmp[,c("importance.score", "Features")]
    }
    print(nrow(data.tmp))
    data.tmp <- data.tmp[data.tmp$importance.score >zscore, ]
    data.tmp <- data.tmp[order(data.tmp$importance.score, decreasing = T), ][1:min(nrow(data.tmp), 15),]
    data.tmp <- data.tmp[order(data.tmp$importance.score, decreasing = F), ]
    data.tmp$Features <- factor(data.tmp$Features, levels = data.tmp$Features)
    p <- ggplot(data=data.tmp, 
                aes(x=Features,y=importance.score)) + 
      geom_bar(position="dodge",stat="identity", width = 0.7, fill=color) + 
      coord_flip() +  
      sensitivity.genomicfeature.mytheme2  + ylab("importance score") + 
      xlab("Features") 
    pdf(paste0(outplot, j, "_",i, "_importance.score_top15.pdf"),width=7, height=4)
    print(p)
    dev.off()
    print(p)
  }
}


plotImportance_class23_2 <- function(datalist, class,zscore){
  i <- class
  for (j in names(datalist)){
    print(j)
    data <- datalist[[j]]
    data.tmp <- data[data$group==i,]
    if(class=="class2.TSS"){
      data.tmp$Features <- unlist(lapply(data.tmp$feature, function(x)strsplit(x, split="class2.TSS_")[[1]][2]))
    }
    if(class=="class3"){
      data.tmp$Features <- unlist(lapply(data.tmp$feature, function(x)strsplit(x, split="class3.TES_")[[1]][2]))
    }
    data.tmp$Motif_ID <- data.tmp$Features
    data.tmp2 <- merge(data.tmp, RBP.motifID_name, by="Motif_ID")
    if(!nrow(data.tmp2)==0){
      #print("what's wrong")
      #print(dim(data.tmp2))
      data.tmp2$Features <- paste0(data.tmp2$RBP_name,",",data.tmp2$Motif_ID)
      data.tmp3 <- data.tmp[!data.tmp$Features %in% data.tmp2$Motif_ID,][,c("importance.score", "Features")]
      data.tmp4 <- data.tmp2[,c("importance.score", "Features")]
      data.tmp <- rbind(data.tmp3, data.tmp4)
     
      
    }else{
      data.tmp <- data.tmp[,c("importance.score", "Features")]
    }
    data.tmp <- data.tmp[data.tmp$importance.score >zscore, ]
    featurefile1 <- as.vector(data.tmp$Features)
    #print(data.tmp[1:20,])
    featuresname1 <- unlist(lapply(featurefile1, function(x)strsplit(x, split="_")[[1]][1]))
    featuresname1 <- unlist(lapply(featuresname1, function(x)strsplit(x, split=",")[[1]][1]))
    data.tmp$featurename <- featuresname1
    #print(data.tmp[1:20,])
    write.table(data.tmp, file=paste0(outfile, j, "_",i, "_zscore",zscore, "_all.txt"), sep = "\t", row.names = F,
                col.names = F, quote = F)
    write.table(data.tmp[,c("featurename")], file=paste0(outfile, j, "_",i, "_zscore",zscore, "_names.txt"), sep = "\t", row.names = F,
                col.names = F, quote = F)
    
  }
}

##merge four classes
plotfeatureCor <- function(data, datagroup, cor.data){
  data.sig <- data[data$importance.score > 1,]
  data.sig <- data.sig[data.sig$group %in% featureclasses, ]
  print(nrow(data.sig))
  feature.cor.data <- cor.data[unique(data.sig$feature) , ]
  feature.cor.data <- feature.cor.data[, colnames(feature.cor.data) %in% unique(data.sig$feature)]
  feature.cor.data.dist <- as.dist(1-abs(feature.cor.data))
  feature.cor.data.dist.tree <- hclust(feature.cor.data.dist, method="ward.D2")
  
  
  feature.cor.data <- feature.cor.data[row.names(feature.cor.data), row.names(feature.cor.data)]
  
  # generate the color scheme to use
  clusterlabel <- c()
  rn <- c()
  for(i in featureclasses){
    data.tmp <- data.sig[data.sig$group == i,]
    clusterlabel <- c(clusterlabel, rep(i, nrow(data.tmp)))
    rn <- c(rn, data.tmp$feature)
  }
  
  colorsplot <- list(group=c(class1.TSS="dodgerblue3",class2.TSS="indianred2",class3="seagreen3",class4="goldenrod2"))
  cluster.annotation <- data.frame(group=factor(clusterlabel, labels = c("class1.TSS","class2.TSS","class3","class4")))
  row.names(cluster.annotation) <- rn
  color.scheme <- rev(heat.colors(11))
  p <- pheatmap(abs(feature.cor.data), 
                color =  color.scheme, border_color = NA,annotation_row = cluster.annotation,annotation_colors =colorsplot[1],
                fontsize_col=8, clustering_method = "ward.D2", fontsize_row=6,cluster_cols =T, cluster_rows = T, cellwidth =5, cellheight = 5, 
                show_rownames = T, show_colnames = F, na_col = "black")
  pdf(file=paste0(outplot,datagroup, "_featurespace_correlation_heatmap.pdf"), width=7, height=7)
  print(p)  
  dev.off()
} 

############individual correlation 

plotfeatureCor_individualClasses <- function(data, datagroup, cor.data, featureclass,zscore){
  data <- data[data$importance.score > zscore,]
  for(i in featureclass){
    data.sig <- data[data$group %in% i, ]
    data.sig <- data.sig[order(data.sig$importance.score, decreasing = T), ][1:min(nrow(data.sig), 15),]
    data.sig <- data.sig[order(data.sig$importance.score, decreasing = F), ]
    print(data.sig)
    print(nrow(data.sig))
    feature.cor.data <- cor.data[unique(data.sig$feature) , ]
    feature.cor.data <- feature.cor.data[, colnames(feature.cor.data) %in% unique(data.sig$feature)]
    feature.cor.data.dist <- as.dist(1-abs(feature.cor.data))
    feature.cor.data.dist.tree <- hclust(feature.cor.data.dist, method="ward.D2")
    
    feature.cor.data <- feature.cor.data[row.names(feature.cor.data), row.names(feature.cor.data)]
    
    
    # generate the color scheme to use
    
    
    color.scheme <- rev(heat.colors(11))
    p <- pheatmap(abs(feature.cor.data), 
                  color =  color.scheme, border_color = NA,
                  fontsize_col=8, clustering_method = "ward.D2", fontsize_row=6,cluster_cols =T, cluster_rows = T, cellwidth =5, cellheight = 5, 
                  show_rownames = T, show_colnames = F, na_col = "black")
    pdf(file=paste0(outplot,datagroup,"_",i,"_featurespace_correlation_heatmap.pdf"), width=7, height=7)
    print(p)  
    dev.off()
  }
  
  
} 

multiplot <- function(..., plotlist = NULL, cols = cols, layout = NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout),
                                               ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

getFeatureTable <- function(featurelist){
  rowname.ref <- row.names(class1_TSS_nascent)
  featuredata <- c()
  fnames.tmp <- names(featurelist)
  fnames <- c()
  for(i in 1: length(featurelist)){
    print(fnames.tmp[i])
    prefix.tmp <- strsplit(fnames.tmp[i],split="_")[[1]][1]
    print(prefix.tmp)
    tmp <- featurelist[[i]]
    row.names(tmp) <- unlist(lapply(row.names(tmp), function(x)strsplit(x, split="\\(")[[1]][1]))
    cn.tmp <- colnames(tmp)
    cn.tmp <- paste0(prefix.tmp, "_",cn.tmp)
    fnames <- c(fnames, cn.tmp)
    tmp <- tmp[rowname.ref, ]
    featuredata <- cbind(featuredata, as.matrix(tmp))
  }
  print(length(fnames))
  print(dim(featuredata))
  colnames(featuredata) <- fnames
  featuredata <- data.frame(featuredata)
  featuredata
}

expressionmatch_subsampling <- function(feature.table, small.group="positive", larger.group= "negative"){
  feature_subset <- feature.table[,c("expressionlevel", "class")]
  feature_subset$name <- row.names(feature_subset)
  group1 <- feature_subset[feature_subset$class==small.group,]$expressionlevel
  group2 <- feature_subset[feature_subset$class==larger.group,]$expressionlevel
  breaks <- seq(-1, ceiling(max(group1))+1000, by=(ceiling(max(group1))+1000)/1000)
  
  ##for subsampling
  group1.data <- feature_subset[feature_subset$class==small.group,]
  group1.data <- dplyr::mutate(group1.data, 
                               expressionRange = cut(group1.data$expressionlevel, breaks = breaks))
  
  group2.data <- feature_subset[feature_subset$class==larger.group,]
  group2.data <- dplyr::mutate(group2.data, expressionRange = cut(group2.data$expressionlevel, 
                                                                  breaks = breaks))
  ###for calculating frequencies
  group1_cut <- table(cut(group1, breaks = breaks))
  group2_cut <- table(cut(group2, breaks = breaks))
  overR <- unique(names(group1_cut[group1_cut > group2_cut]))
  underR <- unique(names(group1_cut[group1_cut <= group2_cut]))
  ##sum over
  diff.over <- sum(group1_cut[overR]) - sum(group2_cut[overR])
  ##sum of 
  under.total <- sum(group1_cut[underR]) + diff.over
  
  group1.under.percen <- group1_cut[underR]/sum(group1_cut[underR])
  group2.sample <- round(under.total * group1.under.percen)
  
  tmp <- unique(group2.data$expressionRange) 
  tmp <- tmp[!is.na(tmp)]
  group2.sample <- group2.sample[names(group2.sample) %in% tmp]
  group2.sample <- group2.sample[group2.sample>0]
  res.down <- c()
  for(i in names(group2.sample)){
    group2.data.tmp <- group2.data[group2.data$expressionRange %in% i, ]
    group2.data.tmp.sub <- dplyr::sample_n(group2.data.tmp,group2.sample[i])
    res.down <- rbind(res.down, group2.data.tmp.sub)
  }
  res <- rbind(res.down, group2.data[group2.data$expressionRange %in% overR, ])
  res2 <- feature.table[row.names(feature.table) %in% res$name,]
  res2
}

expressionmatch_subsampling2 <- function(feature.table, small.group="positive", larger.group= "negative"){
  feature_subset <- feature.table[,c("expressionlevel", "class")]
  feature_subset$name <- row.names(feature_subset)
  group1 <- feature_subset[feature_subset$class==small.group,]$expressionlevel
  group2 <- feature_subset[feature_subset$class==larger.group,]$expressionlevel
  breaks <- seq(-1, ceiling(max(group1))+1000, by=(ceiling(max(group1))+1000)/2000)
  ##for subsampling
  group1.data <- feature_subset[feature_subset$class==small.group,]
  group1.data <- dplyr::mutate(group1.data, 
                               expressionRange = cut(group1.data$expressionlevel, breaks = breaks))
  
  group2.data <- feature_subset[feature_subset$class==larger.group,]
  group2.data <- dplyr::mutate(group2.data, expressionRange = cut(group2.data$expressionlevel, 
                                                                  breaks = breaks))
  ###for calculating frequencies
  group1_cut <- table(cut(group1, breaks = breaks))
  group2_cut <- table(cut(group2, breaks = breaks))
  overR <- unique(names(group1_cut[group1_cut > group2_cut]))
  underR <- unique(names(group1_cut[group1_cut <= group2_cut]))
  #######differences
  overR <- unique(names(group1_cut[group1_cut > group2_cut]))
  underR <- unique(names(group1_cut[group1_cut <= group2_cut]))
  
  
  ###sample proportionally in each categories
  group1.sample.overR <- group2_cut[overR]
  group2.sample.underR <- group1_cut[underR]
  
  ##result for group1
  
  tmp1 <- unique(group1.data$expressionRange) 
  tmp1 <- tmp1[!is.na(tmp1)]
  group1.sample.overR <- group1.sample.overR[names(group1.sample.overR) %in% tmp1]
  #group1.sample.overR <- group1.sample.overR[group1.sample.overR>0]
  res.sample1 <- c()
  for(i in names(group1.sample.overR)){
    group1.data.tmp <- group1.data[group1.data$expressionRange %in% i, ]
    group1.data.tmp.sub <- dplyr::sample_n(group1.data.tmp, group1.sample.overR[i])
    res.sample1 <- rbind(res.sample1, group1.data.tmp.sub)
  }
  res1 <- rbind(res.sample1, group1.data[group1.data$expressionRange %in% underR, ])
  
  ##result for group1
  
  tmp2 <- unique(group2.data$expressionRange) 
  tmp2 <- tmp2[!is.na(tmp2)]
  group2.sample.underR <- group2.sample.underR[names(group2.sample.underR) %in% tmp2]
  group2.sample.underR <- group2.sample.underR[group2.sample.underR>0]
  res.sample2 <- c()
  for(i in names(group2.sample.underR)){
    group2.data.tmp <- group2.data[group2.data$expressionRange %in% i, ]
    group2.data.tmp.sub <- dplyr::sample_n(group2.data.tmp, group2.sample.underR[i])
    res.sample2 <- rbind(res.sample2, group2.data.tmp.sub)
  }
  res2 <- rbind(res.sample2, group2.data[group2.data$expressionRange %in% overR, ])
  
  res1_1 <- feature.table[row.names(feature.table) %in% res1$name,]
  res2_1 <- feature.table[row.names(feature.table) %in% res2$name,]
  reslist <- list()
  reslist[["positive"]] <- res1_1
  reslist[["negative"]] <- res2_1
  reslist
}


expressionmatch_subsampling3 <- function(feature.table, group1.tag, group2.tag, group3.tag){
  feature_subset <- feature.table[,c("expressionlevel", "class")]
  feature_subset$name <- row.names(feature_subset)
  group1 <- feature_subset[feature_subset$class==group1.tag,]$expressionlevel
  group2 <- feature_subset[feature_subset$class==group2.tag,]$expressionlevel
  group3 <- feature_subset[feature_subset$class==group3.tag,]$expressionlevel
  
  min.tmp <- max(c(min(group1),min(group2),min(group3)))
  max.tmp <- min(c(max(group1),max(group2),max(group3)))
  
  breaks <- seq(min.tmp-1, ceiling(max.tmp) + 1000, by=(ceiling(max.tmp)+1000)/2000)
  
  ##for subsampling
  group1.data <- feature_subset[feature_subset$class==group1.tag,]
  group1.data <- dplyr::mutate(group1.data, 
                               expressionRange = cut(group1.data$expressionlevel, breaks = breaks))
  
  group2.data <- feature_subset[feature_subset$class==group2.tag,]
  group2.data <- dplyr::mutate(group2.data, expressionRange = cut(group2.data$expressionlevel, 
                                                                  breaks = breaks))
  
  group3.data <- feature_subset[feature_subset$class==group3.tag,]
  group3.data <- dplyr::mutate(group3.data, 
                               expressionRange = cut(group3.data$expressionlevel, 
                                                     breaks = breaks))
  
  data.list <- list()
  data.list[[group1.tag]] <- group1.data
  data.list[[group2.tag]] <- group2.data
  data.list[[group3.tag]] <- group3.data
  
  ###for calculating frequencies
  group1_cut <- table(cut(group1, breaks = breaks))
  group2_cut <- table(cut(group2, breaks = breaks))
  group3_cut <- table(cut(group3, breaks = breaks))
  
  sample.size.data <- data.frame(cbind(group1_cut, group2_cut, group3_cut))
  sample.size <- do.call(pmin, sample.size.data)
  names(sample.size) <- row.names(sample.size.data)
  res <- c()
  for (i in c(group1.tag, group2.tag, group3.tag)){
    tmp <- data.list[[i]]
    tmp1 <- unique(tmp$expressionRange) 
    tmp1 <- tmp1[!is.na(tmp1)]
    res.sample.tmp <- c()
    for(j in names(sample.size)){
      tmp2 <- tmp[tmp$expressionRange %in% j, ]
      tmp2.sub <- dplyr::sample_n(tmp2, sample.size[j])
      res.sample.tmp <- rbind(res.sample.tmp, tmp2.sub)
    }
    res <- rbind(res, res.sample.tmp)
    res
  }
  res
}



sensitivity.genomicfeature.mytheme2 <- theme_bw() + theme(plot.title = element_text(size = 13),
                                                          axis.title = element_text(size = 13),
                                                          axis.text = element_text(size = 13),
                                                          panel.grid.major = element_blank(), 
                                                          panel.grid.minor = element_blank(),
                                                          legend.title = element_text( size=20),
                                                          legend.text = element_text(size = 18))

mytheme <- theme_bw() + theme(plot.title = element_text(size = 14),
                              axis.title = element_text(size = 12),
                              axis.text = element_text(size = 12),
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = "black"),
                              legend.title = element_text( size=14),
                              legend.text = element_text(size = 14))