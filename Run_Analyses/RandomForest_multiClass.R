library("optparse")
library(pROC)
source("/binf-isilon/sandelin/people/mengjun//Exosome_ML/Run_Analyses/RunRandomForest_functions.R")

option_list = list(  
  make_option("--datadir", type="character", default=NULL, 
              help="directory for all exosome targets", metavar="character"),
  make_option("--featureclass", type="character", default=NULL, 
              help="which combination of features to use", metavar="character")
); 

opt = parse_args(OptionParser(option_list=option_list))
datadir <- opt$datadir
featureclass <- opt$featureclass

outdir <-  paste0(datadir, "RandomForest_importantFeatures_stableFeatures/")
dir.create(outdir)

classdata <- read.table(paste0(datadir, "classid.txt"), header = T)

##import all featuredata

class1_TSS_chrEnv <- neatload(paste0(datadir, "features/TSS_ChIPseq_result.rda"))
class1_TSS_nascent <- neatload(paste0(datadir, "features/Nascentseq_result.rda"))
class1_TSS_GCcontent <- neatload(paste0(datadir, "features/GC_content.window.rda"))
class1_TSS_GCspread <- neatload(paste0(datadir, "features/GC_content_spread_res.rda"))

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
                    class1.TSS_GCspread=class1_TSS_GCspread, 
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
    #print(fnames.tmp[i])
    prefix.tmp <- strsplit(fnames.tmp[i],split="_")[[1]][1]
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
featureTable <- getFeatureTable(featurelist)
expressionlevel <- featureTable$class1.TSS_hg38_HeLa_NETseq_pooled

##read in subfeatures
subfeaturedir <- "/binf-isilon/sandelin/people/mengjun//Exosome_ML/outfile/"
class1.TSS <- read.table(paste0(subfeaturedir, "class1_iterativeSelectedFeatures.txt"), header=F, stringsAsFactors = F)$V1
class2.TSS <- read.table(paste0(subfeaturedir, "class2_iterativeSelectedFeatures.txt"), header=F, stringsAsFactors = F)$V1
class3 <- read.table(paste0(subfeaturedir, "class3_iterativeSelectedFeatures.txt"), header=F, stringsAsFactors = F)$V1
class4 <- read.table(paste0(subfeaturedir, "class4_iterativeSelectedFeatures.txt"), header=F, stringsAsFactors = F)$V1

subfeatureList <- list()
subfeatureList[["class1.TSS"]] <- class1.TSS
subfeatureList[["class2.TSS"]] <- class2.TSS
subfeatureList[["class3"]] <- class3
subfeatureList[["class4"]] <- class4

subfeatureList[["TSS"]] <- c(class1.TSS, class2.TSS)
subfeatureList[["TES"]] <- c(class3, class4)
subfeatureList[["all"]] <- c(class1.TSS, class2.TSS, class3, class4)

subfeatures <- subfeatureList[[featureclass]]

print("number of subfeatures:")
print(length(subfeatures))

featureTable <- featureTable[, subfeatures]
print("number of features:")
print(ncol(featureTable))

featureTable$expressionlevel <- expressionlevel

featureTable$name <- row.names(featureTable)
featureTable <- merge(featureTable, classdata, by="name")
row.names(featureTable) <- featureTable$name
featureTable <- featureTable[, !names(featureTable)%in% c("name")]

####set up random forest#####
rep.split <- 10
ctrl <- trainControl(method = "repeatedcv", 
                     number = 5, 
                     repeats = 5, 
                     verboseIter = FALSE)

featureTable2 <- featureTable
for (i in 1: rep.split){
  res <- expressionmatch_subsampling3(feature.table = featureTable2, group1.tag = "z1",
                                      group2.tag = "z8",
                                      group3.tag = "insensitive")
  
  featureTable.z1 <- res[res$class == "z1", ]
  featureTable.z1 <- featureTable2[row.names(featureTable2) %in% featureTable.z1$name,]
  print("z1 cases...")
  print(nrow(featureTable.z1))
  
  featureTable.z8 <- res[res$class == "z8", ]
  featureTable.z8 <- featureTable2[row.names(featureTable2) %in% featureTable.z8$name,]
  print("z8 cases...")
  print(nrow(featureTable.z8))
  
  featureTable.insensitive <- res[res$class == "insensitive", ]
  featureTable.insensitive <- featureTable2[row.names(featureTable2) %in% featureTable.insensitive$name,]
  print("insensitive cases...")
  print(nrow(featureTable.insensitive))
  
  featureTable <- rbind(featureTable.z1, featureTable.z8, featureTable.insensitive)
  featureTable <- featureTable[,colnames(featureTable) !="expressionlevel"]
  
  index <- createDataPartition(featureTable$class, p = 0.7, list = FALSE)
  train_data <- featureTable[index, ]
  test_data  <- featureTable[-index, ]
  
  cores <- 20
  cl <- makePSOCKcluster(cores)
  registerDoParallel(cl)
  print("runing random forest...")
  model <- caret::train(class ~ .,
                        data = train_data,
                        method = "parRF",
                        preProcess = c("scale", "center"),
                        trControl = ctrl,
                        metric= "Accuracy",
                        importance = TRUE)
  ###calculate metrics
  print("calculating metrics...")
  predicted.classes <- predict(model, test_data)
  observed.classes <- test_data$class
  
  predicted.classes <- factor(predicted.classes, levels = c("z1","z8", "insensitive"))
  observed.classes <- factor(observed.classes, levels = c("z1","z8", "insensitive"))
  
  confusion.matrix.res <- caret::confusionMatrix(data=predicted.classes, reference=observed.classes)
  print("writing confusion Matrix...")
  write.csv(confusion.matrix.res$table,
            file=paste0(outdir, featureclass, "_rep",i, "_confusionMatrix.csv"))
  
  confusem.res <- confusion.matrix.res$overall
  confusem.Accuracy <- confusem.res["Accuracy"]
  
  AUC_prob <- multiclass.roc(observed.classes, predict(model, test_data,type="prob"))
  AUC_prob <- as.numeric(AUC_prob$auc)
  
  metric.res <- data.frame(Accurary=confusem.Accuracy,
                    AUC_prob=AUC_prob)

  
  print("calculating importance...")
  write.table(metric.res, 
              file=paste0(outdir, featureclass, "_rep",i, "_metrics.txt"), 
              quote=F, sep = "\t",row.names = F, col.names = T)
  
  ####write feature importance
  importance.score <- varImp(model, scale = F)$importance
  
  write.table(importance.score, 
              file=paste0(outdir, featureclass, "_rep",i, "_importance.txt"), 
              quote=F, sep = "\t",row.names = T, col.names = T)
}




