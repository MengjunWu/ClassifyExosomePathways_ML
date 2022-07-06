library(TFBSTools)
library(JASPAR2016)
db <- file.path(system.file("extdata", package="JASPAR2016"),
                "JASPAR2016.sqlite")
TATA.pfm <- getMatrixByName(db, name="TATA-Box")
INR.pfm <- getMatrixByName(db, name="INR")

TATA_box <- toPWM(TATA.pfm)@profileMatrix
INR <- toPWM(INR.pfm)@profileMatrix

save(TATA_box, file="/binf-isilon/sandelin/people/mengjun/Exosome_ML/Input_feature_data/TATA_box.rda")
save(INR, file="/binf-isilon/sandelin/people/mengjun/Exosome_ML/Input_feature_data/INR.rda")


ss5 <- as.matrix(data.table::fread("/binf-isilon/sandelin/people/mengjun/Exosome_ML/data/ss5_readableform"))
row.names(ss5) <- c("A", "C", "G", "T")
ss5 <- PWM(ss5)
save(ss5, file="/binf-isilon/sandelin/people/mengjun/Exosome_ML/Input_feature_data/Input_feature_preparation/ss5.rda")


pas <- as.matrix(data.table::fread("/binf-isilon/sandelin/people/mengjun/Exosome_ML/data/pas_readableform"))
row.names(pas) <- c("A", "C", "G", "T")
pas <- PWM(pas)
save(pas, file="/binf-isilon/sandelin/people/mengjun/Exosome_ML/Input_feature_data/pas.rda")
