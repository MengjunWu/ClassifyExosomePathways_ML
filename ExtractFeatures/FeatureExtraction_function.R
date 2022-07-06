

get_bw_sum <- function(bw, ranges){
  bw_file <- rtracklayer::BigWigFile(bw)
  seqlevels(ranges, pruning.mode="coarse") <- seqnames(seqinfo(bw_file))
  values <- unlist(lapply(rtracklayer::import(bw_file, 
                                              which=ranges, 
                                              as = 'NumericList'),sum)) 
  data.frame(gene=ranges$name,
             value=values)
}

get_data_for_ChIPseq <- function(bw_base, ranges){
  print(paste0(' doing: ', bw_base))
  bwfile <- paste0(bigwigdir,bw_base, "_hg38.bw")
  sumsignal<- get_bw_sum(bwfile, ranges) %>%  
    mutate(ChIPseq = bw_base)
  sumsignal
}

get_data_for_strandedSeq <- function(bw_base, ranges.plus, ranges.minus){
  print(paste0(' doing: ', bw_base))
  plus_bw <- paste0(bigwigdir,bw_base,'_plus.bw')
  minus_bw <- paste0(bigwigdir,bw_base,'_minus.bw')
  
  sumsignal <- bind_rows(get_bw_sum(plus_bw, ranges.plus), 
                         get_bw_sum(minus_bw, ranges.minus)) %>% 
    mutate(Nascentseq = bw_base)
  sumsignal
}


CalculateCGFromSequenceSet <- function(SequenceSet, letters=c("G", "C"), window){
  tmp <- SequenceSet
  tmp.no <- length(tmp)
  tmp.names <- names(tmp)
  res <- c()
  seqname <- c()
  for(i in 1: tmp.no){
    tmp.seq <- rowSums(letterFrequencyInSlidingView(tmp[[i]], view.width=window, letters=letters))/window
    res <- rbind(res, tmp.seq)
    seqname <- c(seqname, tmp.names[i])
  }
  row.names(res) <- seqname
  res
}



CalculateGCspread <- function(data, GC.background){
  name <- c
  N = rowsum(data)
  m = 0
  row.n <- nrow(data)
  col.n <- ncol(data)
  fix.term <- c(col.n:1)
  data.tmp <- data - GC.background
  data.tmp[data.tmp < 0] <- 0
  res <- c()
  for (i in 1:row.n){
    N.tmp <- N[i]
    v.tmp <- fix.term * data.tmp[i]/N.tmp
    res <- c(res, v.tmp)
  }
  
}


motif_pwm_hits <- function(SequenceSet, motif.pwm.list, min.score="90%"){
  motif.ids <- motif.pwm.list$RBP.ID
  seq.res <- c()
  for(i in 1: length(SequenceSet)){
    seq.tmp <- SequenceSet[[i]]
    seq.res.tmp <- c()
    for(j in motif.ids){
      hits.count <- countPWM(motif.pwm.list[[j]], 
                             seq.tmp, min.score = min.score)
      print(hits.count)
      seq.res.tmp <- c(seq.res.tmp, hits.count)
    }
    seq.res <- rbind(seq.res, seq.res.tmp)
  }
  colnames(seq.res) <- motif.ids
  row.names(seq.res) <- names(SequenceSet)
  seq.res
}

neatload <- function(dataname){
  temp.space <- new.env()
  bar <- load(dataname, temp.space)
  the.object <- get(bar, temp.space)
  rm(temp.space)
  the.object
}

##estimate CG background
get_avgCG_signal <- function(data, start.postion, window){
  end.position <- start.postion + window -1
  data <- data[, start.postion: end.position]
  print(dim(data))
  data.mean <- mean(data)
  data.mean
}

##estimate GC spread by quantile, 
##and assign different weight to GC in different position (decay with the position)
##the weight is modelled as decay, the distribution: math.pow(2,-x**2/decay**2)
##setting bp limit (width)
##position correction by "avg signal" 
GC_content.spread <- function(data, quantile, decay, width, avg_signal){
  row.n <- nrow(data)
  col.n <- min(width,ncol(data))
  pos <- rep(0, row.n)
  if (decay < 1) {
    decay <- round(col.n * decay)
  }
  scale_decay_position <-  scale_position(1:col.n, decay=decay)
  for(i in 1: row.n){
    total.scale <- sum(data[i, 1:col.n]*scale_decay_position)
    target.scale <-  total.scale * quantile
    rt <- 0
    for (j in 1: col.n){
      dr.tmp <- data[i,j] * scale_decay_position[j]
      rt <- rt + as.numeric(dr.tmp)
      if(rt > target.scale){
        #print (j)
        pos[i] <- j
        #calculate the non-scaled signal
        signal_found <- mean(data[i, 1: j])
        if(signal_found < avg_signal){
          pos[i] <- round(j*(signal_found/avg_signal))
        }
        break
      }
    }
  }
  pos
}

scale_position <- function(vector, decay){
  res <- c()
  for(i in vector){
    dr.tmp <- decay_rate(i, decay = decay)
    res <- c(res, dr.tmp)
  }
  res
}
decay_rate <- function(x, decay){
  dr <- 2^(-x^2/decay^2)
  dr
}

CalculateCGFromSequenceSet <- function(SequenceSet, letters=c("G", "C"), window){
  tmp <- SequenceSet
  tmp.no <- length(tmp)
  tmp.names <- names(tmp)
  res <- c()
  seqname <- c()
  for(i in 1: tmp.no){
    tmp.seq <- rowSums(letterFrequencyInSlidingView(tmp[[i]], view.width=window, letters=letters))/window
    res <- rbind(res, tmp.seq)
    seqname <- c(seqname, tmp.names[i])
  }
  row.names(res) <- seqname
  res
}

getfulltable <- function(datadir){
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
  class4_TES_chrEnv <- class4_TES_chrEnv[, !names(class4_TES_chrEnv)%in% c("GSM733759_Hela_Pol2b", "GSM733785_Hela_CTCF", "GSM1003520_Hela_Ezh2")]
  
  
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
  featureTable <- getFeatureTable(featurelist)
  
  classdata <- read.table(paste0(datadir, "classid.txt"), header = T)
  featureTable$name <- row.names(featureTable)
  featureTable <- merge(featureTable, classdata, by="name")
  row.names(featureTable) <- featureTable$name
  featureTable <- featureTable[, !names(featureTable)%in% c("name")]
  featureTable.pos <- featureTable[featureTable$class=="positive", ]
  featureTable.neg <- featureTable[featureTable$class=="negative", ]
  res <- list()
  featureTable.pos <- featureTable.pos[, !names(featureTable)%in% c("class")]
  featureTable.neg <- featureTable.neg[, !names(featureTable)%in% c("class")]
  res[["positive"]] <- featureTable.pos
  res[["negative"]] <- featureTable.neg
  res
}




getFeatureTable <- function(featurelist){
  nt <- names(featurelist)
  nascentname <- nt[grep(pattern = "nascent", nt)]
  rowname.ref <- row.names(featurelist[[nascentname]])
  featuredata <- c()
  fnames.tmp <- names(featurelist)
  fnames <- c()
  for(i in 1: length(featurelist)){
    #print(fnames.tmp[i])
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
