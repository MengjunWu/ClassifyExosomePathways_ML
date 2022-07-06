library("optparse")
library("rtracklayer")
library("dplyr")
library("reshape2")
library("Biostrings")
library("seqPattern")
source("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/scripts/ExtractFeatures/FeatureExtraction_function.R")

option_list = list(  
  make_option("--datadir", type="character", default=NULL, 
              help="file with overlap information of ends with CLIP binding proteins", metavar="character"),
  make_option("--overlapfile", type="character", default=NULL, 
              help="file with overlap information of ends with CLIP binding proteins", metavar="character"),
  make_option("--outdir", type="character", default=NULL, 
              help="directory for the output file", metavar="character")
); 

opt = parse_args(OptionParser(option_list=option_list))
datadir <- opt$datadir
overlapfile <- opt$overlapfile
outdir <- opt$outdir

load("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/RBP_id_unique.rda")
RBP_id_unique_dataframe <- data.frame(RBP_id = RBP_id_unique)
RBP_id_unique_dataframe$RBP_id <- as.vector(RBP_id_unique_dataframe$RBP_id)

RBP.overlap <- data.table::fread(paste0(datadir, overlapfile))
RBP.overlap.count <- RBP.overlap[, c("V4", "V10")]
colnames(RBP.overlap.count) <- c("gene_id", "RBP_id")
RBP.overlap.count$unique.id <- paste0(RBP.overlap.count$gene_id, ".", RBP.overlap.count$RBP_id)
RBP.overlap.count.id <- plyr::count(RBP.overlap.count, vars="unique.id")
RBP.overlap.count <- merge(RBP.overlap.count, RBP.overlap.count.id)
RBP.overlap.count <- RBP.overlap.count[!duplicated(RBP.overlap.count),]
res <- c()
geneids <- c()
for(i in unique(RBP.overlap.count$gene_id)){
  geneids <- c(geneids, i)
  tmp <- RBP.overlap.count[RBP.overlap.count$gene_id == i, ]
  tmp2 <- full_join(tmp, RBP_id_unique_dataframe)
  tmp2 <- tmp2[,c("RBP_id", "freq")]
  tmp2[is.na(tmp2)] <- 0
  tmp2 <- data.frame(tmp2)
  row.names(tmp2) <- tmp2$RBP_id
  tmp2 <- tmp2[RBP_id_unique, ]
  res <- rbind(res, tmp2$freq)
}
colnames(res) <- RBP_id_unique
row.names(res) <- geneids
classdata <- read.table(paste0(datadir, "classid.txt"), header = T)
noCLIPdatagene <-  classdata$name[!classdata$name %in% row.names(res)]
noCLIPdatagene.data <- data.frame(matrix(0, nrow=length(noCLIPdatagene), ncol = ncol(res)))
row.names(noCLIPdatagene.data) <- noCLIPdatagene
colnames(noCLIPdatagene.data) <- colnames(res)

res <- rbind(res, noCLIPdatagene.data)
save(res, file=paste0(outdir, "RBP_binding_CLIP.rda"))
