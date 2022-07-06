library("optparse")
library("rtracklayer")
library("dplyr")
library("reshape2")
library("Biostrings")
library("seqPattern")
source("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/scripts/ExtractFeatures/FeatureExtraction_function.R")

option_list = list(  
  make_option("--fastafiledir", type="character", default=NULL, 
              help="directory for fastafile", metavar="character"),
  make_option("--outdir", type="character", default=NULL, 
              help="directory for the output file", metavar="character")
); 

opt = parse_args(OptionParser(option_list=option_list))
fastafiledir <- opt$fastafiledir
outdir <- opt$outdir

SequenceSet <- readDNAStringSet(fastafiledir)

##For PAS relative location is fixed###
##PAS, -5 to -50

center.position <- 1101
start.PAS <- center.position - 50
end.PAS <- center.position - 1
SequenceSet.PAS <- DNAStringSet(x=SequenceSet, start=start.PAS, end=end.PAS)

load("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/pas.rda")

PASmotif <- list(PAS=pas, RBP.ID=c("PAS"))
motif.PAS.freq <- data.frame(motif_pwm_hits(SequenceSet.PAS, motif.pwm.list = PASmotif))
row.names(motif.PAS.freq) <- unlist(lapply(row.names(motif.PAS.freq), function(x)strsplit(x, split="\\:")[[1]][1]))
save(motif.PAS.freq, file=paste0(outdir, "motif.PAS.TESupstream.rda"))

####plot PAS location#####
pdf(paste0(outdir, "TES_pas_motif.pdf"),width=10, height=6)
plotMotifOccurrenceAverage(SequenceSet.PAS, motifPWM= pas, smoothingWindow = 3)
dev.off()
