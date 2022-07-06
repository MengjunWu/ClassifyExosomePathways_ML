library("optparse")
library("rtracklayer")
library("dplyr")
library("reshape2")
library("Biostrings")
source("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/scripts/ExtractFeatures/FeatureExtraction_function.R")
option_list = list(  
  make_option("--fastafiledir", type="character", default=NULL, 
              help="directory for fastafile", metavar="character"),
  make_option("--slidingWindow", type="integer", default=25, 
              help="sliding window to calculate CG content", metavar="number"),
  make_option("--outdir", type="character", default=NULL, 
              help="directory for the output file", metavar="character")
); 

opt = parse_args(OptionParser(option_list=option_list))
fastafiledir <- opt$fastafiledir
SlidingWindow <- opt$slidingWindow
outdir <- opt$outdir
###

SequenceSet <- readDNAStringSet(fastafiledir)
print("calculating the CG content...")
GC_content <- CalculateCGFromSequenceSet(SequenceSet, letters=c("G", "C"), window=SlidingWindow)

#save(GC_content, file=paste0(outdir, "GC_content.rda"))

###estimate the background from 500 bp and 1000 bp downstream

GC_content <- GC_content - 0.5
GC_content[GC_content < 0] <- 0

print(GC_content[1:10, 2000:2050])

center.position <- 1101
window.500 <- 500 
window.1000 <- 1000 
signal500 <- get_avgCG_signal(GC_content, start.postion = center.position, window=window.500)
signal1000 <- get_avgCG_signal(GC_content, start.postion = center.position, window=window.1000)

res <- data.frame(signal500=signal500, signal1000=signal1000)
print(res)
write.table(res, file=paste0(outdir, "CG_background_Estimate_result.txt"))

