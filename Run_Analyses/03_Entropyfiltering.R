library(entropy)
library(matrixStats)
source("/binf-isilon/sandelin/people/mengjun/Exosome_ML/ExtractFeatures/FeatureExtraction_function.R")

####################get full feature table for all classes########################
filedir.Z8 <- "/binf-isilon/sandelin/people/mengjun/Exosome_ML/Z8sensitiveVSinsensitive_all/"
filedir.Z1 <- "/binf-isilon/sandelin/people/mengjun/Exosome_ML/Z1sensitiveVSinsensitive_all/"
filedir.Z1Z8<- "/binf-isilon/sandelin/people/mengjun/Exosome_ML/Z1sensitiveVSZ8sensitive_all/"


feature.Z8 <- getfulltable(filedir.Z8)
feature.Z1 <- getfulltable(filedir.Z1)
feature.Z1Z8 <- getfulltable(filedir.Z1Z8)

z8.sensitive <- feature.Z8$positive
z8.sensitive$class <- "Z8"
allinsensitive <- feature.Z8$negative
allinsensitive$class <- "insensitive"
z1.sensitive <- feature.Z1$positive
z1.sensitive$class <- "Z1"

features_fulltable <- rbind(z8.sensitive,
                            z1.sensitive,
                            allinsensitive)


outputdir <- "/binf-isilon/sandelin/people/mengjun/Exosome_ML/Input_feature_data/"
save(features_fulltable, file=paste0(outputdir, "features_fulltable_GC1000.rda"))

##############################filter by entropy#######################################################
entropyvalue <- c()
for (i in 1: ncol(features_fulltable)){
  a <- table(features_fulltable[,i])
  entropyvalue <- c(entropyvalue, entropy(a, unit="log2"))
  
}

res <- data.frame(feature=colnames(features_fulltable), 
                  entropyv=entropyvalue)
res2 <- res[res$entropyv >= 0.8,]
save(res2$feature, file=paste0(outputdir, "EntropyFilteredFeature_GC1000.rda"))

