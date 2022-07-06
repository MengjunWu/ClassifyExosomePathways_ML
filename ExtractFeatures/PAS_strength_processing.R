library("optparse")
option_list = list(  
  make_option("--filedir", type="character", default=NULL, 
              help="directory for fastafile", metavar="character"),
  make_option("--outdir", type="character", default=NULL, 
              help="directory for the output file", metavar="character")
); 

opt = parse_args(OptionParser(option_list=option_list))
filedir <- opt$filedir
outdir <- opt$outdir

PAS.score <- read.csv(filedir)
PAS.score.name <- unlist(lapply(as.vector(PAS.score$X), function(x)strsplit(x, split="\\:")[[1]][1]))

PAS.score <- data.frame(PAS_score=PAS.score$PAS_score)
row.names(PAS.score) <- PAS.score.name

save(PAS.score, file=paste0(outdir, "PAS.score.rda"))

load("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/PAS.score.rda")

