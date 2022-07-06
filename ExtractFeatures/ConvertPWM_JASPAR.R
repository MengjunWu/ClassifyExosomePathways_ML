library(TFBSTools)
library(JASPAR2016)
db <- file.path(system.file("extdata", package="JASPAR2016"),
                "JASPAR2016.sqlite")
TATA.pfm <- getMatrixByName(db, name="TATA-Box")
INR.pfm <- getMatrixByName(db, name="INR")

TATA_box <- toPWM(TATA.pfm)@profileMatrix
INR <- toPWM(INR.pfm)@profileMatrix

save(TATA_box, file="/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/TATA_box.rda")
save(INR, file="/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/INR.rda")


ss5 <- as.matrix(data.table::fread("/binf-isilon/sandelin/people/mengjun/ASAP_example/jasparMotif/ss5/ss5_readableform"))
row.names(ss5) <- c("A", "C", "G", "T")
ss5 <- PWM(ss5)
save(ss5, file="/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/ss5.rda")


pas <- as.matrix(data.table::fread("/binf-isilon/sandelin/people/mengjun/ASAP_example/jasparMotif/awtaaa/pas_readableform"))
row.names(pas) <- c("A", "C", "G", "T")
pas <- PWM(pas)
save(pas, file="/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivity/Input_feature_preparation/pas.rda")
