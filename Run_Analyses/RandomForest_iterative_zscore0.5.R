library("optparse")
source("/binf-isilon/sandelin/people/mengjun//Exosome_ML/Run_Analyses/RunRandomForest_functions.R")

option_list = list(  
  make_option("--datadir", type="character", default=NULL, 
              help="directory for the binary exosome targets pair", metavar="character"),
  make_option("--iterativefiledir", type="character", default=NULL, 
              help="directory for iteratively selected feature file", metavar="character"),
  make_option("--iterative", type="numeric", default=NULL, 
              help="iterative round", metavar="number"),
  make_option("--outdir", type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option("--zscore", type="numeric", default=NULL, 
              help="zscore cutoff", metavar="number")
  
); 

opt = parse_args(OptionParser(option_list=option_list))
datadir <- opt$datadir
iterativefiledir <- opt$iterativefiledir
iterative <- opt$iterative
outdir <- opt$outdir
zscore <- opt$zscore
zscore <- as.numeric(zscore)

iterativefile <- paste0(iterativefiledir,"selectedfeatures_iterative", iterative, "zscore",zscore,"_importance.txt")
if(!file.exists(iterativefile)){
  stop("iterativefeature does not exist...")
}else{
  print(iterativefile)
  subfeatures <- read.table(iterativefile, header = T)
  subfeatures <- as.vector(subfeatures$features)
  classdata <- read.table(paste0(datadir, "classid.txt"), header = T)
  
  ##import all featuredata
  if(length(subfeatures) <= 1){
    stop("too few features to run the model")
  }else{
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
    ####extract features table
    entropyfeatures<- neatload("/binf-isilon/sandelin/people/mengjun//Exosome_ML/Input_feature_data/EntropyFilteredFeature_GC1000.rda")
    featureTable <- getFeatureTable(featurelist)
    featureTable <- featureTable[classdata$name, ]
    expressionlevel <- featureTable$class1.TSS_hg38_HeLa_NETseq_pooled
    
    featureTable <- featureTable[,colnames(featureTable) %in% entropyfeatures]
    
    featureTable <- featureTable[, subfeatures]
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
      featureTable.pos <- featureTable2[featureTable2$class == "positive", ]
      print("positive cases...")
      print(nrow(featureTable.pos))
      
      featureTable.neg <- featureTable2[featureTable2$class == "negative", ]
      print("negative cases...")
      print(nrow(featureTable.neg))
      
      res <- expressionmatch_subsampling2(featureTable2)
      featureTable.pos <- res[["positive"]]
      featureTable.neg <- res[["negative"]]
      featureTable <- rbind(featureTable.pos, featureTable.neg)
      featureTable <- featureTable[,colnames(featureTable) !="expressionlevel"]
      
      print("positive cases after sampling...")
      print(nrow(featureTable.pos))
      
      print("negative cases after sampling...")
      print(nrow(featureTable.neg))
      
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
                            importance = TRUE)
      ###calculate metrics
      print("calculating metrics...")
      predicted.classes <- predict(model, test_data)
      observed.classes <- test_data$class
      predicted.classes <- factor(predicted.classes, levels = c("positive","negative"))
      observed.classes <- factor(observed.classes, levels = c("positive","negative"))
      
      metric.res <- CalcMetrics(prediction = predicted.classes, truth = observed.classes)
      
      prediction_auc <- predict(model, test_data,type="prob")[,2]
      pred <- prediction(prediction_auc,observed.classes)
      pred.auc <- performance(pred, measure = "auc")@y.values[[1]]
      metric.res$AUC <- pred.auc
      
      print("calculating importance...")
      write.table(metric.res, 
                  file=paste0(outdir, "runiterative", iterative, "_rep",i, "_metrics.txt"), 
                  quote=F, sep = "\t",row.names = F, col.names = T)
      
      ####write feature importance
      importance.score <- varImp(model, scale = F)$importance
      
      write.table(importance.score, 
                  file=paste0(outdir, "runiterative", iterative, "_rep",i, "_importance.txt"), 
                  quote=F, sep = "\t",row.names = T, col.names = T)
    }
  }
  
}

