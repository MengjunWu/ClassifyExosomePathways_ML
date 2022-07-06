library("optparse")
source("/binf-isilon/sandelin/people/mengjun//Exosome_ML/Run_Analyses/RunRandomForest_functions.R")
option_list = list(
  make_option("--filedir", type="character", default=NULL, 
              help="directory for the binary exosome targets pair", metavar="character"),
  make_option("--iterative", type="numeric", default=NULL, 
              help="number of interative round", metavar="number"),
  make_option("--zscore", type="numeric", default=NULL, 
              help="zscore cutoff", metavar="number")
  
); 

opt = parse_args(OptionParser(option_list=option_list))
filedir <- opt$filedir
iterative <- opt$iterative
zscore <- opt$zscore
zscore <- as.numeric(zscore)

setwd(filedir)
importance.files <- list.files(pattern = paste0("runiterative", iterative, "_"))
importance.files <- importance.files[grep(pattern = "importance", x=importance.files)]
importance.tmp <- calc_average_importance(importance.files)

importance.tmp$features <- row.names(importance.tmp)
importance.tmp <- importance.tmp[order(importance.tmp$importance.score, decreasing = F),]
print("number of features before filter:")
print(nrow(importance.tmp))

print("zscore: ")
print(zscore)

importance.tmp.sub <- importance.tmp[importance.tmp$importance.score > zscore, ]
print("number of features after filter:")
print(nrow(importance.tmp.sub))

if(nrow(importance.tmp.sub)==0){
  write.table(importance.tmp, file=paste0("all_zscoreShouldbelessthan0.5_after",iterative2, "iters.txt"),
              row.names = T, col.names = T, sep = "\t", quote=F)
  stop(paste0("no feature with importance score larger than: ", zscore))
}else{
  iterative2 <- as.numeric(iterative)+1
  write.table(importance.tmp.sub, file=paste0("selectedfeatures_iterative", iterative2, "zscore",zscore,"_importance.txt"), 
              row.names = T, col.names = T, sep = "\t", quote=F)
}

