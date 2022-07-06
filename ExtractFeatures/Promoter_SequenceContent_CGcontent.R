library("optparse")
library("rtracklayer")
library("dplyr")
library("reshape2")
library("Biostrings")
source("/binf-isilon/sandelin/people/mengjun/Exosome_ML/ExtractFeatures/FeatureExtraction_function.R")

option_list = list(  
  make_option("--fastafiledir", type="character", default=NULL, 
              help="directory for fastafile", metavar="character"),
  make_option("--slidingWindow", type="integer", default=25, 
              help="sliding window to calculate CG content", metavar="number"),
  make_option("--windowcg", type="integer", default=25, 
              help="windowsize to calculate CG content", metavar="number"),
  make_option("--outdir", type="character", default=NULL, 
              help="directory for the output file", metavar="character")
); 

opt = parse_args(OptionParser(option_list=option_list))
fastafiledir <- opt$fastafiledir
SlidingWindow <- opt$slidingWindow
windowcg <- opt$windowcg
outdir <- opt$outdir

###
SequenceSet <- readDNAStringSet(fastafiledir)
print("calculating the CG content new...")
GC_content <- CalculateCGFromSequenceSet(SequenceSet, letters=c("G", "C"), window=SlidingWindow)
###center 
center.position <- 1101
start.position <- center.position - windowcg
end.position <- center.position + windowcg - 1
GC_content.window <- GC_content[, start.position: end.position]
GC_content.window <- data.frame(rowSums(GC_content.window))
colnames(GC_content.window) <- "CG_content"
row.names(GC_content.window) <- unlist(lapply(row.names(GC_content.window), function(x)strsplit(x, split="\\:")[[1]][1]))
save(GC_content.window, file=paste0(outdir, "GC_content.window.rda"))
###window for GC spread
##fixed calculation for 500 and 1000

end.position.spread1 <- center.position + 1000 - 1
end.position.spread2 <- center.position + 500 - 1

GC_content.spread1 <- GC_content[, center.position: end.position.spread1]
GC_content.spread2 <- GC_content[, center.position: end.position.spread2]

GC_content.spread1 <- GC_content.spread1 - 0.5
GC_content.spread1[GC_content.spread1 < 0] <- 0

GC_content.spread2 <- GC_content.spread2 - 0.5
GC_content.spread2[GC_content.spread2 < 0] <- 0

##caculate GC spread with different parameters
  
GC_content.spread1000_quantile0.5 <- GC_content.spread(GC_content.spread1, quantile = 0.5, 
                                                       decay = 0.5, width = 1000, 
                                                       avg_signal = 0.1)

GC_content.spread1000_quantile0.75 <- GC_content.spread(GC_content.spread1, quantile = 0.75, 
                                                       decay = 0.5, width = 1000, 
                                                       avg_signal = 0.1)

GC_content.spread500_quantile0.5 <- GC_content.spread(GC_content.spread2, quantile = 0.5, 
                                                       decay = 0.5, width = 500, 
                                                       avg_signal = 0.1)

GC_content.spread500_quantile0.75 <- GC_content.spread(GC_content.spread2, quantile = 0.75, 
                                                        decay = 0.5, width = 500, 
                                                        avg_signal = 0.1)

GC_content_spread_res <- cbind(GC_quantile0.5_spread1000 = GC_content.spread1000_quantile0.5, 
                               GC_quantile0.75_spread1000 = GC_content.spread1000_quantile0.75,
                               GC_quantile0.5_spread500 = GC_content.spread500_quantile0.5, 
                               GC_quantile0.75_spread500 = GC_content.spread500_quantile0.75)

row.names(GC_content_spread_res) <- row.names(GC_content)
row.names(GC_content_spread_res) <- unlist(lapply(row.names(GC_content_spread_res), function(x)strsplit(x, split="\\:")[[1]][1]))
print(GC_content_spread_res[1:3,])
save(GC_content_spread_res, file=paste0(outdir, "GC_content_spread_res.rda"))
