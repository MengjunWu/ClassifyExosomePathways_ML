filedir <- "/binf-isilon/sandelin/people/mengjun/Exosome_ML/data/PAS_score_denovo_TES205.csv"
outdir <- "/binf-isilon/sandelin/people/mengjun/Exosome_ML/Input_feature_data/"

PAS.score <- read.csv(filedir)
PAS.score.name <- unlist(lapply(as.vector(PAS.score$X), function(x)strsplit(x, split="\\:")[[1]][1]))

PAS.score <- data.frame(PAS_score=PAS.score$PAS_score)
row.names(PAS.score) <- PAS.score.name

save(PAS.score, file=paste0(outdir, "PAS.score.rda"))

