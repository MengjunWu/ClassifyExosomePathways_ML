library("optparse")
source("/binf-isilon/sandelin/people/mengjun//Exosome_ML/Run_Analyses/RunRandomForest_functions.R")

library(optparse)
library(matrixStats)
option_list = list(
  make_option("--filedir", type="character", default=NULL, 
              help="directory for the binary exosome targets pair", metavar="character"),
  make_option("--featureclass", type="character", default=NULL, 
              help="feature class", metavar="character"),
  make_option("--zscore", type="numeric", default=NULL, 
              help="zscore cutoff", metavar="number")
); 


opt = parse_args(OptionParser(option_list=option_list))
filedir <- opt$filedir
featureclass <- opt$featureclass
zscore <- opt$zscore

filedir.iterative <- paste0(filedir, "iterative_zscore", zscore, "_", featureclass)
outfiledir <- paste0(filedir, "iterative_zscore", zscore, "_summary_output/")
dir.create(outfiledir)
setwd(filedir.iterative)

##take all the iterative files
iterative.files <- list.files(pattern = "runiterative")
filepattern <- unlist(lapply(iterative.files, function(x)strsplit(x, split="_")[[1]][1]))
filepattern <- as.numeric(unique(unlist(lapply(filepattern, function(x)strsplit(x, split="tive")[[1]][2]))))
filepattern <- filepattern[order(filepattern, decreasing = F)]
iterative.pattern <- paste0("runiterative", filepattern)

iterative.files.importance <- iterative.files[grep(pattern = "importance", x=iterative.files)]
iterative.files.metrics <- iterative.files[grep(pattern = "metrics", x=iterative.files)]


###generate files#####

###calculate the metrics 
res.metric <- c()
for (i in 1: length(iterative.pattern)){
  tmp.name <- iterative.pattern[i]
  #print(tmp.name)
  tmp.name2 <- paste0("^", tmp.name,"_")
  tmp.files <- iterative.files.metrics[grep(pattern = tmp.name2, x=iterative.files.metrics)]
  #print(tmp.files)
  metrics.tmp <- calc_average_metric(tmp.files)
  metrics.tmp$metric <- row.names(metrics.tmp)
  metrics.tmp$group <- tmp.name
  res.metric <- rbind(res.metric, metrics.tmp)
}

####calculate the importance 
importance.n <- c()
for (i in 1: length(iterative.pattern)){
  tmp.name <- iterative.pattern[i]
  #print(tmp.name)
  tmp.name2 <- paste0("^", tmp.name,"_")
  tmp.files <- iterative.files.importance[grep(pattern = tmp.name2, x=iterative.files.importance)]
  importance.tmp <- calc_average_importance(tmp.files)
  n.tmp <- nrow(importance.tmp)
  importance.n <- c(importance.n, n.tmp)
}
res.importance <- data.frame(group=iterative.pattern, feature.number=importance.n)
res.metric.importance <- merge(res.metric, res.importance)
res.metric.importance$group <- factor(res.metric.importance$group, levels = iterative.pattern)
res.metric.importance.AUC <- res.metric.importance[res.metric.importance$metric=="AUC", ]

res.metric.importance.F1score <- res.metric.importance[res.metric.importance$metric=="F1.score", ]


if(featureclass %in% c("class2", "class3", "TSS", "TES","all")){
  factor <- 250
}
if(featureclass %in% c("class1", "class4")){
  factor <- 50
}

pdf(file=paste0(outfiledir,featureclass, "_iteration_AUC.pdf"), width = 9, height = 4)
ggplot(res.metric.importance.AUC, aes(x=group,y=mean)) + 
  geom_bar( stat="identity", width=0.8) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                width=.2, position=position_dodge(0.05)) + 
  geom_text(aes(label=round(mean,digits = 2)), vjust=3, color="white", size=2.5) + 
  geom_line(aes(y=feature.number/factor, group=1), size=1) + scale_y_continuous(
    # Features of the first axis
    name = "Average AUC",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*factor, name="Number of features")
  ) + mytheme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file=paste0(outfiledir,featureclass, "_iteration_F1.score.pdf"), width = 9, height = 4)
ggplot(res.metric.importance.F1score, aes(x=group,y=mean)) + 
  geom_bar( stat="identity", width=0.8) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                width=.2, position=position_dodge(0.05)) + 
  geom_text(aes(label=round(mean,digits = 2)), vjust=3, color="white", size=2.5) + 
  geom_line(aes(y=feature.number/factor, group=1), size=1) + scale_y_continuous(
    # Features of the first axis
    name = "Average F1 score",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*factor, name="Number of features")
  ) + mytheme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

###write the metric file
res.metric.iterative5 <- res.metric.importance[res.metric.importance$group == "runiterative5",] 
res.metric.iterative5$group <- featureclass 

write.table(iterative5.importance, file=paste0(outfiledir,featureclass, "_metric.subfeature.txt"),
            row.names = F, col.names = T, sep = "\t", quote=F)

###write the subfeature file of iterative 5
tmp.last <- paste0("^", iterative.pattern[6], "_")
iterative5.importance <- iterative.files.importance[grep(pattern = tmp.last, x=iterative.files.importance)]
iterative5.importance <- calc_average_importance(iterative5.importance)
iterative5.importance$features <- row.names(iterative5.importance)
iterative5.importance$group <- featureclass 
write.table(iterative5.importance, file=paste0(outfiledir,featureclass, "_importance.subfeature.txt"),
            row.names = F, col.names = T, sep = "\t", quote=F)