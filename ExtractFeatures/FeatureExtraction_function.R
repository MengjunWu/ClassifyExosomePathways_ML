

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

