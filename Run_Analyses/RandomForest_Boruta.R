library("optparse")
library(Boruta)
source("/binf-isilon/sandelin/people/mengjun//Exosome_ML/Run_Analyses/RunRandomForest_functions.R")

option_list = list(  
  make_option("--datadir", type="character", default=NULL, 
              help="directory for the binary exosome targets pair", metavar="character"),
  make_option("--featureclass", type="character", default=NULL, 
              help="which combination of features to use", metavar="character")
); 

opt = parse_args(OptionParser(option_list=option_list))
datadir <- opt$datadir
featureclass <- opt$featureclass


outdir <- paste0(datadir, "RandomForest_Boruta/")
dir.create(outdir)

classdata <- read.table(paste0(datadir, "classid.txt"), header = T)

##import all featuredata

class1_TSS_chrEnv <- neatload(paste0(datadir, "features/TSS_ChIPseq_result.rda"))
class1_TSS_nascent <- neatload(paste0(datadir, "features/Nascentseq_result.rda"))
class1_TSS_GCcontent <- neatload(paste0(datadir, "features/GC_content.window.rda"))
class1_TSS_GCspread <- neatload(paste0(datadir, "features/GC_content_spread_res.rda"))
class1_TSS_GCspread0.75 <- data.frame(GC_quantile0.75_spread1000=class1_TSS_GCspread[,"GC_quantile0.75_spread1000"])

class1_TSS_TATAInr <- neatload(paste0(datadir, "features/Prompter_TATA_INR.rda"))

class2_TSS_u1pasRBP <-  neatload(paste0(datadir, "features/TSS_motif_RBP_U1PAS_freq.rda"))
class2_TSS_RBPclip <- neatload(paste0(datadir, "features/TSS_overlap500_RBP_binding_CLIP.rda"))

class3_TES_u1pasRBP <- neatload(paste0(datadir, "features/TES_motif_RBP_U1PAS_freq_oneside.rda"))
class3_TES_RBPclip <- neatload(paste0(datadir, "features/TES_upstream500_RBP_binding_CLIP.rda"))
class3_pas_upstream <- neatload(paste0(datadir, "features/motif.PAS.TESupstream.rda"))
class3_pas_strength <- neatload("/binf-isilon/sandelin/people/mengjun//Exosome_ML/Input_feature_data/PAS.score.rda")
class4_TES_chrEnv <- data.frame(neatload(paste0(datadir, "features/TES_ChIPseq_result.rda"))) 

featurelist <- list(class1.TSS_chrEnv=class1_TSS_chrEnv,
                    class1.TSS_nascent=class1_TSS_nascent, 
                    class1.TSS_GCcontent=class1_TSS_GCcontent, 
                    class1.TSS_GCspread=class1_TSS_GCspread0.75, 
                    class1.TSS_TATAInr=class1_TSS_TATAInr,
                    class2.TSS_u1pasRBP=class2_TSS_u1pasRBP,
                    class2.TSS_RBPclip=class2_TSS_RBPclip,
                    class3.TES_u1pasRBP=class3_TES_u1pasRBP, 
                    class3.TES_RBPclip=class3_TES_RBPclip, 
                    class3.TES_pas_upstream=class3_pas_upstream, 
                    class3.TES_pas_strength=class3_pas_strength,
                    class4.TES_chrEnv=class4_TES_chrEnv)

getFeatureTable <- function(featurelist){
  rowname.ref <- row.names(class1_TSS_nascent)
  featuredata <- c()
  fnames.tmp <- names(featurelist)
  fnames <- c()
  for(i in 1: length(featurelist)){
    print(fnames.tmp[i])
    prefix.tmp <- strsplit(fnames.tmp[i],split="_")[[1]][1]
    print(prefix.tmp)
    tmp <- featurelist[[i]]
    row.names(tmp) <- unlist(lapply(row.names(tmp), function(x)strsplit(x, split="\\(")[[1]][1]))
    cn.tmp <- colnames(tmp)
    cn.tmp <- paste0(prefix.tmp, "_",cn.tmp)
    fnames <- c(fnames, cn.tmp)
    tmp <- tmp[rowname.ref, ]
    featuredata <- cbind(featuredata, as.matrix(tmp))
  }
  print(length(fnames))
  print(fnames[1:3])
  print(dim(featuredata))
  colnames(featuredata) <- fnames
  featuredata <- data.frame(featuredata)
  featuredata
}

####extract features table
entropyfeatures<- neatload("/binf-isilon/sandelin/people/mengjun/Exosome_SLICCAGE_3end/Determinants_ExosomeSensitivityNew/Input_feature_preparation/EntropyFilteredFeature_GC1000.rda")

featureTable <- getFeatureTable(featurelist)
print("feature table rows")
print(nrow(featureTable))
featureTable <- featureTable[classdata$name, ]
print(nrow(featureTable))
expressionlevel <- featureTable$class1.TSS_hg38_HeLa_NETseq_pooled
featureTable <- featureTable[,colnames(featureTable) %in% entropyfeatures]

if(featureclass != "all"){
  featureTable <- featureTable[, grep(pattern = featureclass, x=colnames(featureTable))]
}
featureTable$expressionlevel <- expressionlevel

print("number of features:")
print(ncol(featureTable))

featureTable$name <- row.names(featureTable)
featureTable <- merge(featureTable, classdata, by="name")
row.names(featureTable) <- featureTable$name
featureTable <- featureTable[, !names(featureTable)%in% c("name")]

####set up random forest#####

featureTable2 <- featureTable

featureTable.pos <- featureTable2[featureTable2$class == "positive", ]
print("positive cases...")
print(nrow(featureTable.pos))

featureTable.neg <- featureTable2[featureTable2$class == "negative", ]
print("negative cases...")
print(nrow(featureTable.neg))

res <- expressionmatch_subsampling2(featureTable2)
featureTable.pos <- res[["positive"]]
featureTable.pos$class <- TRUE
featureTable.neg <- res[["negative"]]
featureTable.neg$class <- FALSE

featureTable <- rbind(featureTable.pos, featureTable.neg)

featureTable <- featureTable[,colnames(featureTable) !="expressionlevel"]

Boruta.imp <- Boruta(class~.,data=featureTable, getImp=getImpRfZ, maxRuns=300)
Boruta.imp.table <- attStats(Boruta.imp)

write.table(Boruta.imp.table, 
            file=paste0(outdir,featureclass, "_Boruta.imp.table.txt"), 
            quote=F, sep = "\t",row.names = T, col.names = T)







