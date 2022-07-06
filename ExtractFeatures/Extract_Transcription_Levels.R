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
  make_option("--windowup", type="integer", default=NULL, 
              help="window size upstream ", metavar="number"),
  make_option("--windowdn", type="integer", default=NULL, 
              help="window size downstream ", metavar="number"),
  make_option("--outdir", type="character", default=NULL, 
              help="directory for the output file", metavar="character")
);

opt = parse_args(OptionParser(option_list=option_list))
bedfile <- opt$bedfile
bigwigdir <- opt$bigwigdir
windowup <- opt$windowup
windowdn <- opt$windowdn
outdir <- opt$outdir

bw_fnames_base <- dir(bigwigdir) %>% keep(grepl('plus.bw', .)) %>% sub('_plus.bw', '', .)

head(bw_fnames_base)

bedfile.grange.original <- rtracklayer::import(bedfile)

##the bed file by default has extended 1kb window from 5end to both sides
##so we extract the 5end first from the file
bedfile.grange.center <- narrow(bedfile.grange.original, start=1101, end=-1100)
bedfile.grange <- promoters(bedfile.grange.center, upstream=windowup, downstream=windowdn)
print("window size:")
print(unique(width(bedfile.grange)))

bedfile.grange_plus <- bedfile.grange[strand(bedfile.grange)=='+',]
bedfile.grange_minus <- bedfile.grange[strand(bedfile.grange)=='-',]


Nascentseq_result <- lapply(bw_fnames_base, function(x)get_data_for_strandedSeq(bw_base = x, 
                                                                                ranges.plus = bedfile.grange_plus,
                                                                                ranges.minus = bedfile.grange_minus)) %>% bind_rows


Nascentseq_result$gene <- factor(Nascentseq_result$gene, levels = unique(Nascentseq_result$gene))
Nascentseq_result$Nascentseq <- factor(Nascentseq_result$Nascentseq, levels = unique(Nascentseq_result$Nascentseq))
Nascentseq_result <- acast(Nascentseq_result, gene~Nascentseq, value.var="value")

print(Nascentseq_result[1:3,])

save(Nascentseq_result, file=paste0(outdir, "Nascentseq_result.rda"))

