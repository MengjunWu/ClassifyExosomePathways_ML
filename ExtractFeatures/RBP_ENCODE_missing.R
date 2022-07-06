RBPs_bindSites <- data.table::fread("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivityNew/Input_feature_preparation/RBP_bindingSite_CLIP_sorted.bed")
RBP_ENCODE <- data.table::fread("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/data/experiment_report_2022_5_30_14h_38m.tsv")
RBPs_bindSites.cell <- sapply(RBPs_bindSites$V4, function(x)strsplit(x, split = "_")[[1]][3])

RBPs_bindSites.K562 <- RBPs_bindSites[RBPs_bindSites$cell == "K562",]
RBPs_bindSites.K562.protein <- sapply(RBPs_bindSites.K562$V4, function(x)strsplit(x, split = "_")[[1]][1])
RBPs_bindSites.K562$RBP <- RBPs_bindSites.K562.protein
RBPs_bindSites.K562.protein <- unique(RBPs_bindSites.K562.protein)

save(RBPs_bindSites.K562, file="/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/data/RBP/RBPs_bindSites.K562.rda")

RBP_ENCODE_protein <- RBP_ENCODE$`Target of assay`
RBP_ENCODE_protein[!RBP_ENCODE_protein %in% RBPs_bindSites.K562.protein]
RBPs_bindSites.K562.protein[!RBPs_bindSites.K562.protein %in%RBP_ENCODE_protein]

RPB.bw <- list.files(path="/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/data/RBP/hg38_bw/",
                    pattern="eCLIP")

RPB.bw.protein <- sapply(RPB.bw, function(x)strsplit(x, split = "_")[[1]][1])
length(unique(RPB.bw.protein))
RBP.bw.K562 <- unique(RPB.bw.protein[RPB.bw.protein%in%RBPs_bindSites.K562.protein])

###summarize the signal
bigwigdir <- "/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/data/RBP/hg38_bw/"
bw_fnames_base <- dir(bigwigdir) %>% keep(grepl('eCLIP_plus.bw', .)) %>% sub('_eCLIP_plus.bw', '', .)


library(GenomicRanges)
res <- 
  lapply(split(RBPs_bindSites.K562, RBPs_bindSites.K562$RBP), function(i){
    GRanges(seqnames = i$V1,
            ranges = IRanges(start = i$V2,
                             end = i$V3,
                             names = i$names),
            strand=i$V6)
  })
j=0
get_data_for_strandedSeq <- function(bw_base, ranges.plus, ranges.minus){
  print(paste0(' doing: ', bw_base))
  plus_bw <- paste0(bigwigdir,bw_base,'_eCLIP_plus.bw')
  minus_bw <- paste0(bigwigdir,bw_base,'_eCLIP_minus.bw')
  
  sumsignal <- bind_rows(get_bw_sum(plus_bw, ranges.plus), 
                         get_bw_sum(minus_bw, ranges.minus)) %>% 
    mutate(RBP = bw_base)
  sumsignal
}
get_bw_sum <- function(bw, ranges){
  bw_file <- rtracklayer::BigWigFile(bw)
  seqlevels(ranges, pruning.mode="coarse") <- seqnames(seqinfo(bw_file))
  values <- unlist(lapply(rtracklayer::import(bw_file, 
                                              which=ranges, 
                                              as = 'NumericList'),sum)) 
  data.frame(peak=names(ranges),
             value=values)
}

RBP_signal <- c()
res.name2 <- names(res)[!names(res) %in% c("DGCR8","HNRNPU")]
for(i in res.name2){
  K562.tmp <- res[[i]]
  bedfile.grange_plus <- K562.tmp[strand(K562.tmp)=='+',]
  bedfile.grange_minus <- K562.tmp[strand(K562.tmp)=='-',]
  bw.signal.tmp <- get_data_for_strandedSeq(bw_base = i, 
                           ranges.plus = bedfile.grange_plus, 
                           ranges.minus = bedfile.grange_minus)
  RBP_signal <- rbind(RBP_signal, bw.signal.tmp)
  
}
RBP_signal_K562 <- RBP_signal
save(RBP_signal_K562, file="/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/data/RBP/RBP_signal_K562.rda")
RBPs_bindSites.K562.bed <- merge(RBPs_bindSites.K562, RBP_signal_K562)
RBPs_bindSites.K562.bed <- RBPs_bindSites.K562.bed[!duplicated(RBPs_bindSites.K562.bed),]
RBPs_bindSites.K562.bed <- RBPs_bindSites.K562.bed[,c("V1","V2","V3","peak","value", "V6")]
write.table(RBPs_bindSites.K562.bed, file="/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/data/RBP/RBPs_bindSites.K562.bed", sep = "\t",
            quote = T, col.names = F, row.names = F)

###major isoform 
library(GenomicFeatures)
gtf.anno <- "gencode.v29.annotation.gtf"
library(rtracklayer)
chr.inf <- SeqinfoForUCSCGenome("hg38")
TxDB.tmp <- GenomicFeatures::makeTxDbFromGFF(gtf.anno, chrominfo=chr.inf)

threeUTR <- unlist(threeUTRsByTranscript(TxDB.tmp, use.names=T))
threeUTR$transcript <- names(threeUTR)
threeUTR <- data.frame(threeUTR)

fiveUTR <-  unlist(fiveUTRsByTranscript(TxDB.tmp, use.names=T))              
fiveUTR$transcript <- names(fiveUTR)
fiveUTR <- data.frame(fiveUTR)

exon <- unlist(exonsBy(TxDB.tmp, "tx", use.names=T))
exon$transcript <- names(exon)
exon <- data.frame(exon)

transcripts <- unlist(transcriptsBy(TxDB.tmp, "gene"))
transcripts$gene <- names(transcripts)
transcripts <- data.frame(transcripts)
transcripts$transcript <- transcripts$tx_name

transcript_length <- dplyr::group_by(exon, transcript)
transcript_length <- data.frame(dplyr::mutate(transcript_length, exonic.transcripts.length = rep(sum(width), length(width))))

fiveUTRlength <- dplyr::group_by(fiveUTR, transcript)
fiveUTRlength <- data.frame(dplyr::mutate(fiveUTRlength, fiveUTR.length = rep(sum(width), length(width))))

threeUTRlength <- dplyr::group_by(threeUTR, transcript)
threeUTRlength <- data.frame(dplyr::mutate(threeUTRlength, threeUTR.length = rep(sum(width), length(width))))

geneLength <- dplyr::group_by(transcripts, gene)
geneLength <- data.frame(dplyr::mutate(geneLength, gene.length = rep(max(width), length(width))))

#####################
transcript_length.merge <- transcript_length[,c("transcript", "exonic.transcripts.length")]
transcript_length.merge <- transcript_length.merge[!duplicated(transcript_length.merge),]

fiveUTRlength.merge <- fiveUTRlength[,c("transcript", "fiveUTR.length")]
fiveUTRlength.merge <- fiveUTRlength.merge[!duplicated(fiveUTRlength.merge),]

threeUTRlength.merge <- threeUTRlength[,c("transcript", "threeUTR.length")]
threeUTRlength.merge <- threeUTRlength.merge[!duplicated(threeUTRlength.merge),]

geneLength.merge <- geneLength[,c("gene", "transcript","gene.length")]

res <- left_join(transcript_length.merge, fiveUTRlength.merge)
res <- left_join(res, threeUTRlength.merge)
res <- left_join(res, geneLength.merge)

###gene annotation
##transcripts length -- exons
##5'UTR length -- transcript
##3'UTR length -- transcript
##gene length
load("ENCODE_K652_RNAseq.rda")
