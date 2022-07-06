library("optparse")

option_list = list(  
  make_option("--outputdir", type="character", default=NULL, 
              help="directory for the output file", metavar="character"),
  make_option("--targetfirst", type="character", default=NULL, 
              help="one of the exosome target categories: PAXT, NEXT, or insensitive", metavar="character"),
  make_option("--targetsecond", type="character", default=NULL, 
              help="another the exosome target category from PAXT, NEXT, or insensitive", metavar="character")
); 
source("/binf-isilon/sandelin/people/mengjun//Exosome_ML/Run_Analyses/ExosomeTargets_SamplePreparation_function.R")

opt = parse_args(OptionParser(option_list=option_list))
outputdir <- opt$outputdir
targetfirst <- opt$targetsecond
targetsecond <- opt$targetsecond
targetPair <- c(targetfirst, targetsecond)


###############################read in the exosome targets files######################################

denovo_hela_dir <- "/binf-isilon/sandelin/people/mengjun/Exosome_ML/data/"
denova_hela_anno <- data.frame(fread(paste0(denovo_hela_dir, "SLAannopublished_MeolaSillaNEXTPAXTtype.tsv")))
load(paste0(denovo_hela_dir, "NEXT_PAXT_nonNEXTPAXT_targets.rda"))
denova_hela_anno2 <- merge(denova_hela_anno, PAXT_NEXT_target)
denova_hela_anno.PAXT <- denova_hela_anno2[denova_hela_anno2$NEXTPAXT_class_new =="PAXT",]
denova_hela_anno.NEXT <- denova_hela_anno2[denova_hela_anno2$NEXTPAXT_class_new =="NEXT",]
denova_hela_anno.insensitive <- denova_hela_anno2[denova_hela_anno2$NEXTPAXT_class_new =="NS",]
Z1Z8insensitive.RNA <- rbind(denova_hela_anno.PAXT, denova_hela_anno.NEXT, denova_hela_anno.insensitive)

#####################prepare files for feature extraction in binary classification#############################
#####process the binary comparison information from the input
binaryComparisons1 <- c("NEXT", "insensitive")
binaryComparisons2 <- c("PAXT", "insensitive")
binaryComparisons3 <- c("NEXT", "PAXT")
index1 <- sum(targetPair %in% binaryComparisons1)
index2 <- sum(targetPair %in% binaryComparisons2)
index3 <- sum(targetPair %in% binaryComparisons3)

index.pair <- which(c(index1, index2, index3) == 2)

if(length(index.pair) == 0){
  stop("could not find corresponding binary comparisons, 
       please check the two input target categories")
}

if(index.pair==1){
  filename <- "Z8sensitiveVSinsensitive"
}

if(index.pair==2){
  filename <- "Z1sensitiveVSinsensitive"
}

if(index.pair==3){
  filename <- "Z1sensitiveVSZ8sensitive"
}

#####generate files for subsequent feature extraction

##all targets
outputdir <- paste0(outputdir, filename, "_all/")
positiveclass <- list(denova_hela_anno.NEXT, denova_hela_anno.PAXT, 
                      denova_hela_anno.PAXT)
negativeclass <- list(denova_hela_anno.insensitive, denova_hela_anno.insensitive, 
                      denova_hela_anno.NEXT)

PrepareInputData_symetricWindow(positiveclass = positiveclass[index.pair],
                                negativeclass = negativeclass[index.pair],
                                output = outputdir)

PrepareInputData_TSS_downstreamWindow(positiveclass = positiveclass[index.pair],
                                      negativeclass = negativeclass[index.pair],
                                      output = outputdir)

PrepareInputData_TES_upstreamWindow(positiveclass = positiveclass[index.pair],
                                    negativeclass = negativeclass[index.pair],
                                    output = outputdir)

##biotypes
biotypes <- c("protein_coding", "lncRNA", "PROMPT", "enhancer")

for(i in biotypes){
  Z1Z8insensitive.RNA1 <- Z1Z8insensitive.RNA[!Z1Z8insensitive.RNA$group %in% c("Z1","Z8"),]
  Z1Z8insensitive.RNA2 <- Z1Z8insensitive.RNA[Z1Z8insensitive.RNA$group %in% c("Z1","Z8"),]
  Z1Z8insensitive.subset.tmp <- data.frame(Z1Z8insensitive.RNA2[Z1Z8insensitive.RNA2$biotype == i,])
  Z1Z8insensitive.subset <- rbind(Z1Z8insensitive.subset.tmp, Z1Z8insensitive.RNA1)
  
  outputdir.biotype <- paste0(outputdir, filename, "_", i, "/")
  PrepareInputData_symetricWindow(positiveclass = positiveclass[index.pair],
                                  negativeclass = negativeclass[index.pair],
                                  subset = Z1Z8insensitive.subset$name,
                                  output = outputdir.biotype)
  
  PrepareInputData_TSS_downstreamWindow(positiveclass = positiveclass[index.pair],
                                        negativeclass = negativeclass[index.pair],
                                        output = outputdir.biotype)
  
  PrepareInputData_TES_upstreamWindow(positiveclass = positiveclass[index.pair],
                                      negativeclass = negativeclass[index.pair],
                                      output = outputdir.biotype)
  
}


