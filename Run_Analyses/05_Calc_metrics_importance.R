library("optparse")
source("/binf-isilon/sandelin/people/mengjun//Exosome_ML/Run_Analyses/RunRandomForest_functions.R")

option_list = list(
  make_option("--filedir", type="character", default=NULL, 
              help="directory for the binary exosome targets pair", metavar="character")
); 


opt = parse_args(OptionParser(option_list=option_list))
filedir <- opt$filedir

outfiledir <- paste0(filedir, "summary_output/")
dir.create(outfiledir)

setwd(filedir)

metrics.files <- list.files(pattern = "_metrics.txt")

files.pattern <- unique(unlist(lapply(metrics.files, function(x)strsplit(x, split="_rep")[[1]][1])))

group.plot <- c("class1", "class2",
                "TSS", "class3", "class4",
                "TES", "all")

###generate files#####
res.metric <- c()
for (i in 1: length(files.pattern)){
  tmp.name <- files.pattern[i]
  print(tmp.name)
  tmp.name2 <- paste0("^", tmp.name, "_")
  tmp.files <- metrics.files[grep(pattern = tmp.name2, x=metrics.files)]
  print(tmp.files)
  metrics.tmp <- calc_average_metric(tmp.files)
  importance.tmp <- calc_average_importance(tmp.files)
  metrics.tmp$metric <- row.names(metrics.tmp)
  metrics.tmp$group <- tmp.name
  res.metric <- rbind(res.metric, metrics.tmp)
}

res.metric <- res.metric[res.metric$group %in% group.plot, ]
res.metric$group <- factor(res.metric$group, 
                           levels = group.plot)
###plot upset plot

barplot_metrics <- function(data, metric, outdir){
  data.metric <- data[data$metric == metric, ]
  
  p <- ggplot(data.metric, aes(x=group, y=mean)) +
    geom_bar(stat="identity", color="gray23", width=0.6,) + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                  width=.2, position=position_dodge(0.05)) + 
    mytheme + coord_cartesian(ylim=c(0.5,1)) + ylab(paste0("Average ", metric))
  
  pdf(file=paste0(outdir,metric, "_plot.pdf"), width = 5, height = 3)
  print(p)
  dev.off()
}

##plot precision
barplot_metrics(res.metric, metric="precision", outdir=outfiledir)
##plot recall
barplot_metrics(res.metric, metric="recall", outdir=outfiledir)
##plot F1 score
barplot_metrics(res.metric, metric="F1.score", outdir=outfiledir)
##plot AUC
barplot_metrics(res.metric, metric="AUC", outdir=outfiledir)



write.table(res.metric, file=paste0(outfiledir,"metricsall.txt"), row.names = T, 
            col.names = T, sep = "\t", quote=F)


####calculate the average importance 
importance.files <- list.files(pattern = "_importance.txt")
importance.files.pattern <- unique(unlist(lapply(importance.files, function(x)strsplit(x, split="_rep")[[1]][1])))


###generate files#####
res.importance <- c()
for (i in 1: length(importance.files.pattern)){
  tmp.name <- importance.files.pattern[i]
  tmp.name2 <- paste0("^", tmp.name, "_")
  tmp.files <- importance.files[grep(pattern = tmp.name2, x=importance.files)]
  print(tmp.files)
  importance.tmp <- calc_average_importance(tmp.files)
  importance.tmp$feature <- row.names(importance.tmp)
  importance.tmp$group <- tmp.name
  res.importance <- rbind(res.importance, importance.tmp)
}

write.table(res.importance, file=paste0(outfiledir,"importanceall.txt"), row.names = T, 
            col.names = T, sep = "\t", quote=F)


