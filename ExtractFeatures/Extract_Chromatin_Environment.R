library("optparse")
library("tidyverse")
library("rtracklayer")
library("dplyr")
library("reshape2")
source("/binf-isilon/sandelin/people/mengjun/Exosome_ML/ExtractFeatures/FeatureExtraction_function.R")

option_list = list(  
  make_option("--bedfile", type="character", default=NULL, 
              help="filename of the bed files", metavar="character"),
  make_option("--bigwigdir", type="character", default=NULL, 
              help="directory for all the bigwig files", metavar="character"),
  make_option("--windowsize", type="integer", default=NULL, 
              help="size of flank window around the end, this is symmetric", metavar="number"),
  make_option("--outdir", type="character", default=NULL, 
              help="directory for the output file", metavar="character")
); 

opt = parse_args(OptionParser(option_list=option_list))
bedfile <- opt$bedfile
bigwigdir <- opt$bigwigdir
windowsize <- opt$windowsize
outdir <- opt$outdir

bwfiles <- list.files(bigwigdir)
bw_fnames_base <- dir(bigwigdir) %>% keep(grepl('hg38.bw', .)) %>% sub('_hg38.bw', '', .)

head(bw_fnames_base)

##pay attention, to use import, it will bed file is 0-based while GRange is 1-based
##so when using import, the grange will increase 1 base of the bed file

bedfile.grange.original <- rtracklayer::import(bedfile)

##the bed file by default has extended 1100bp window from 5end to both sides
##so we extract the 5end first from the file
bedfile.grange.center <- narrow(bedfile.grange.original, start=1101, end=-1100)
bedfile.grange <- promoters(bedfile.grange.center, upstream=windowsize, downstream=windowsize)
print("window size:")
print(unique(width(bedfile.grange)))

ChIPseq_result <- lapply(bw_fnames_base, function(x)get_data_for_ChIPseq(bw_base=x, ranges=bedfile.grange)) %>% bind_rows
ChIPseq_result$gene <- factor(ChIPseq_result$gene, levels = unique(ChIPseq_result$gene))
ChIPseq_result$ChIPseq <- factor(ChIPseq_result$ChIPseq, levels = unique(ChIPseq_result$ChIPseq))
ChIPseq_result <- acast(ChIPseq_result, gene~ChIPseq, value.var="value")

print(ChIPseq_result[1:3,1:10])
save(ChIPseq_result, file=paste0(outdir, "ChIPseq_result.rda"))

