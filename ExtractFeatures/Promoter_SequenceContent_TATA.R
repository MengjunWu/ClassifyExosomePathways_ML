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

#fastafiledir="/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/RunPrediction/Z8sensitiveVsZ8insensitive_denovo/TSS.window1200.fa"
SequenceSet <- readDNAStringSet(fastafiledir)
print(unique(width(SequenceSet)))
##For TATA and InR, the window size is fixed###
##TATA, -10 to -50
##InR, +1 to -10

center.position <- 1101
end.TATA <- center.position 
start.TATA<- center.position - 50
SequenceSet.TATA <- DNAStringSet(x=SequenceSet, start=start.TATA, end=end.TATA)
print("what's wrong")
end.Inr<- center.position + 10 
start.Inr <- center.position - 15
SequenceSet.Inr <- DNAStringSet(x=SequenceSet, start=start.Inr, end=end.Inr)



load("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/TATA_box.rda")
load("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/INR.rda")
TATA.motif.list <- list(RBP.ID=c("TATA_box"), TATA_box=TATA_box)
Inr.motif.list <- list(RBP.ID=c("INR"), INR=INR)

motif.TATA.freq <- data.frame(motif_pwm_hits(SequenceSet.TATA, motif.pwm.list = TATA.motif.list))
motif.Inr.freq <- data.frame(motif_pwm_hits(SequenceSet.Inr, motif.pwm.list = Inr.motif.list))
res <- cbind(motif.TATA.freq, motif.Inr.freq)
row.names(res) <- unlist(lapply(row.names(res), function(x)strsplit(x, split="\\:")[[1]][1]))
print(res[1:3,])
save(res, file=paste0(outdir, "Prompter_TATA_INR.rda"))

##plot out the feature also to do sanity check###

pdf(paste0(outdir, "Prompter_TATA_motif.pdf"),width=10, height=6)
plotMotifOccurrenceAverage(SequenceSet.TATA, motifPWM= TATA_box, smoothingWindow = 3)
dev.off()

pdf(paste0(outdir, "Prompter_Inr_motif.pdf"),width=10, height=6)
plotMotifOccurrenceAverage(SequenceSet.Inr, motifPWM= INR, smoothingWindow = 3)
dev.off()
