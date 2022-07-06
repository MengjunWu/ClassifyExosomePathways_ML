PrepareInputData_symetricWindow <- function(positiveclass, negativeclass, subset=NULL, output, window=1100){
  positiveclass$class <- "positive"
  negativeclass$class <- "negative"
  tmpfile <- rbind(positiveclass, negativeclass)
  if(!is.null(subset)){
    tmpfile <- tmpfile[tmpfile$name %in% subset, ]
  }
  classid <- tmpfile[,c("name", "class")]
  tmpfile <- data.frame(chr=tmpfile$seqnames, start=tmpfile$start, end=tmpfile$end,
                        id=tmpfile$name, score=tmpfile$RNA_size, strand=tmpfile$strand)
  tmpfile.plus <- tmpfile[tmpfile$strand=="+", ]
  tmpfile.minus <- tmpfile[tmpfile$strand=="-", ]
  tmpfile.plus.TES <- tmpfile.plus
  tmpfile.plus.TSS <- tmpfile.plus
  tmpfile.minus.TES <- tmpfile.minus
  tmpfile.minus.TSS <- tmpfile.minus
  
  tmpfile.plus.TES$start <- tmpfile.plus$end -window
  tmpfile.plus.TES$end <- tmpfile.plus$end + window -1
  
  
  tmpfile.plus.TSS$start <- tmpfile.plus$start - window
  tmpfile.plus.TSS$end <- tmpfile.plus$start + window -1
  
  tmpfile.minus.TES$start <- tmpfile.minus$start -window
  tmpfile.minus.TES$end <- tmpfile.minus$start + window -1
  
  tmpfile.minus.TSS$start <- tmpfile.minus$end - window
  tmpfile.minus.TSS$end <- tmpfile.minus$end + window -1
  
  tmpfile.TSS <- rbind(tmpfile.plus.TSS, tmpfile.minus.TSS)
  tmpfile.TES <- rbind(tmpfile.plus.TES, tmpfile.minus.TES)
  tmpfile.TSS <- tmpfile.TSS[order(tmpfile.TSS$start),]
  tmpfile.TSS <- tmpfile.TSS[order(tmpfile.TSS$chr),]
  tmpfile.TES <- tmpfile.TES[order(tmpfile.TES$start),]
  tmpfile.TES <- tmpfile.TES[order(tmpfile.TES$chr),]
  
  write.table(classid, file=paste0(output, "classid.txt"), col.names = T,
              row.names = F, quote = F, sep = "\t")
  write.table(tmpfile.TSS ,file=paste0(output, "TSS.window", window, ".bed"), col.names = F,
              row.names = F, quote = F, sep = "\t")
  write.table(tmpfile.TES ,file=paste0(output, "TES.window", window, ".bed"), col.names = F,
              row.names = F, quote = F, sep = "\t")
  
}

PrepareInputData_TSS_downstreamWindow <- function(positiveclass, negativeclass, subset=NULL, 
                                                  output, downstreamwindow=500){
  positiveclass$class <- "positive"
  negativeclass$class <- "negative"
  tmpfile <- rbind(positiveclass, negativeclass)
  if(!is.null(subset)){
    tmpfile <- tmpfile[tmpfile$name %in% subset, ]
  }
  classid <- tmpfile[,c("name", "class")]
  tmpfile <- data.frame(chr=tmpfile$seqnames, start=tmpfile$start, end=tmpfile$end,
                        id=tmpfile$name, score=tmpfile$RNA_size, strand=tmpfile$strand)
  tmpfile.plus <- tmpfile[tmpfile$strand=="+", ]
  tmpfile.minus <- tmpfile[tmpfile$strand=="-", ]
  tmpfile.plus.TSS <- tmpfile.plus
  tmpfile.minus.TSS <- tmpfile.minus
  
  tmpfile.plus.TSS$end <- tmpfile.plus$start + downstreamwindow 
  
  
  tmpfile.minus.TSS$start <- tmpfile.minus$end - downstreamwindow
  
  tmpfile.TSS <- rbind(tmpfile.plus.TSS, tmpfile.minus.TSS)
  tmpfile.TSS <- tmpfile.TSS[order(tmpfile.TSS$start),]
  tmpfile.TSS <- tmpfile.TSS[order(tmpfile.TSS$chr),]
  
  write.table(tmpfile.TSS ,file=paste0(output, "TSS.downstream", downstreamwindow, ".bed"), col.names = F,
              row.names = F, quote = F, sep = "\t")
  
}


PrepareInputData_TES_upstreamWindow <- function(positiveclass, negativeclass, subset=NULL, 
                                                output, upstreamwindow=500){
  positiveclass$class <- "positive"
  negativeclass$class <- "negative"
  tmpfile <- rbind(positiveclass, negativeclass)
  if(!is.null(subset)){
    tmpfile <- tmpfile[tmpfile$name %in% subset, ]
  }
  classid <- tmpfile[,c("name", "class")]
  tmpfile <- data.frame(chr=tmpfile$seqnames, start=tmpfile$start, end=tmpfile$end,
                        id=tmpfile$name, score=tmpfile$RNA_size, strand=tmpfile$strand)
  tmpfile.plus <- tmpfile[tmpfile$strand=="+", ]
  tmpfile.minus <- tmpfile[tmpfile$strand=="-", ]
  tmpfile.plus.TSS <- tmpfile.plus
  tmpfile.minus.TSS <- tmpfile.minus
  
  tmpfile.plus.TSS$start <- tmpfile.plus$end - upstreamwindow
  
  tmpfile.minus.TSS$end <- tmpfile.minus$start + upstreamwindow 
  tmpfile.TSS <- rbind(tmpfile.plus.TSS, tmpfile.minus.TSS)
  tmpfile.TSS <- tmpfile.TSS[order(tmpfile.TSS$start),]
  tmpfile.TSS <- tmpfile.TSS[order(tmpfile.TSS$chr),]
  
  write.table(tmpfile.TSS ,file=paste0(output, "TES.upstream", upstreamwindow, ".bed"), col.names = F,
              row.names = F, quote = F, sep = "\t")
  
}
