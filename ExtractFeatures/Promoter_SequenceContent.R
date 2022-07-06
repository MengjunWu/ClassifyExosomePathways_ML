library("optparse")
library("rtracklayer")
library("dplyr")
library("reshape2")
library("Biostrings")
source("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/scripts/FeatureExtraction_function.R")
option_list = list(  
  make_option("--fastafiledir", type="character", default=NULL, 
              help="directory for fastafile", metavar="character"),
  make_option("--SlidingWindow", type="integer", default=NULL, 
              help="sliding window to calculate CG content", metavar="number"),
  make_option("--windowcg", type="integer", default=NULL, 
              help="window size of cg content", metavar="number"),
  make_option("--windowcgspread", type="integer", default=NULL, 
              help="window size of cg spread", metavar="number"),
  make_option("--outdir", type="character", default=NULL, 
              help="directory for the output file", metavar="character")
); 

opt = parse_args(OptionParser(option_list=option_list))
fastafiledir <- opt$fastafiledir
SlidingWindow <- opt$SlidingWindow
windowcg <- opt$windowcg
windowcgspread <- opt$windowcgspread
outdir <- opt$outdir
###
SequenceSet <- readDNAStringSet(fastafiledir)
GC_content <- CalculateGCFromSequenceSet(SequenceSet, letters=c("G", "C"), window=SlidingWindow)

###center 
center.position <- 1101
start.position <- center.position - windowcg
end.position <- center.position + windowcg - 1
GC_content.window <- GC_content[, start.position: end.position]
GC_content.window <- data.frame(rowSums(GC_content.window))
colnames(GC_content.window) <- "CG_content"

CalculateCGFromSequenceSet <- function(SequenceSet, letters=c("G", "C"), window){
  tmp <- SequenceSet
  tmp.no <- length(tmp)
  tmp.names <- names(tmp)
  res <- c()
  seqname <- c()
  for(i in 1: tmp.no){
    tmp.seq <- rowSums(letterFrequencyInSlidingView(tmp[[i]], view.width=window, letters=letters))/window
    res <- rbind(res, tmp.seq)
    seqname <- c(seqname, tmp.names[i])
  }
  row.names(res) <- seqname
  res
}


