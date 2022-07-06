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
  make_option("--onesidewindow", type="integer", default=NULL, 
              help="The flank window from the end, it is one sided", metavar="number"),
  make_option("--outdir", type="character", default=NULL, 
              help="directory for the output file", metavar="character")
); 

opt = parse_args(OptionParser(option_list=option_list))
fastafiledir <- opt$fastafiledir
##for TSS one only care about downstream, so the window the distance from the start site
window <- opt$onesidewindow
outdir <- opt$outdir
#fastafiledir <- "/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivityNew/RunPrediction/Z8sensitiveVSallinsensitive_denovo/TSS.window1100.fa"
SequenceSet <- readDNAStringSet(fastafiledir)
center.position <- 1101
end.position <- center.position + window -1
SequenceSet <- DNAStringSet(x=SequenceSet, start=center.position, end=end.position)

load("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/RBP.pwm.rda")
load("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/pas.rda")
load("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/ss5.rda")
RBP.pwm[["PAS"]] <-  pas
RBP.pwm[["ss5"]] <- ss5
RBP.pwm$RBP.ID <- c(RBP.pwm$RBP.ID, "PAS", "ss5")
print("scanning the motif...")

motif.RBP.U1PAS.freq <- data.frame(motif_pwm_hits(SequenceSet, motif.pwm.list = RBP.pwm))
row.names(motif.RBP.U1PAS.freq) <- unlist(lapply(row.names(motif.RBP.U1PAS.freq), function(x)strsplit(x, split="\\:")[[1]][1]))

motif.RBP.U1PAS.freq$names <- row.names(motif.RBP.U1PAS.freq)
motif.RBP.U1PAS.freq <- merge(motif.RBP.U1PAS.freq, trial)

save(motif.RBP.U1PAS.freq, file=paste0(outdir, "TSS_motif_RBP_U1PAS_freq.rda"))


