library("optparse")
library("rtracklayer")
library("dplyr")
library("reshape2")
library("Biostrings")
library("seqPattern")
source("/binf-isilon/sandelin/people/mengjun/Exosome_ML/ExtractFeatures/FeatureExtraction_function.R")

option_list = list(  
  make_option("--fastafiledir", type="character", default=NULL, 
              help="directory for fastafile", metavar="character"),
  make_option("--onesidewindow", type="integer", default=NULL, 
              help="The flank window from the end, it is one sided", metavar="number"),
  make_option("--twosidewindow", type="integer", default=NULL, 
              help="The flank window around the end, it is two sided", metavar="number"),
  make_option("--outdir", type="character", default=NULL, 
              help="directory for the output file", metavar="character")
); 

opt = parse_args(OptionParser(option_list=option_list))
fastafiledir <- opt$fastafiledir
onesidewindow <- opt$onesidewindow
twosidewindow <- opt$twosidewindow
outdir <- opt$outdir

SequenceSet <- readDNAStringSet(fastafiledir)
center.position <- 1101
start.position.oneside <- center.position - onesidewindow + 1
SequenceSet.oneside <- DNAStringSet(x=SequenceSet, start=start.position.oneside, end=center.position)


load("/binf-isilon/sandelin/people/mengjun/Exosome_ML/Input_feature_data/RBP.pwm.rda")
load("/binf-isilon/sandelin/people/mengjun/Exosome_ML/Input_feature_data/pas.rda")
load("/binf-isilon/sandelin/people/mengjun/Exosome_ML/Input_feature_data/ss5.rda")


U1PAS <- list(PAS=pas, ss5=ss5, RBP.ID=c("PAS", "ss5"))

print("scanning the motif...")
motif.RBP.freq <- data.frame(motif_pwm_hits(SequenceSet.oneside, motif.pwm.list = RBP.pwm))
motif.U1PAS.freq <- data.frame(motif_pwm_hits(SequenceSet.oneside, motif.pwm.list = U1PAS))
motif.U1PAS.freq <- motif.U1PAS.freq[row.names(motif.RBP.freq), ]

motif.RBP.U1PAS.freq <- cbind(motif.RBP.freq, motif.U1PAS.freq)
row.names(motif.RBP.U1PAS.freq) <- unlist(lapply(row.names(motif.RBP.U1PAS.freq), function(x)strsplit(x, split="\\:")[[1]][1]))

save(motif.RBP.U1PAS.freq, file=paste0(outdir, "TES_motif_RBP_U1PAS_freq_oneside.rda"))

