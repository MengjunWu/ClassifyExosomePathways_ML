z1sins <- "/binf-isilon/sandelin/people/mengjun/Exosome_ML/Z1sensitiveVSinsensitive_all/"
z8sins <- "/binf-isilon/sandelin/people/mengjun/Exosome_ML/Z8sensitiveVSinsensitive_all/"
z1z8s <- "/binf-isilon/sandelin/people/mengjun/Exosome_ML/Z1sensitiveVSZ8sensitive_all/"

z1sensitiveVsInsensitive <- read.table(paste0(z1sins, "classid.txt"), header=T)
z8sensitiveVsInsensitive <- read.table(paste0(z8sins, "classid.txt"), header=T)
z1sensitiveVsz8sensitive <- read.table(paste0(z1z8s, "classid.txt"), header=T)

z1sensitiveVsInsensitive.pos <- z1sensitiveVsInsensitive[z1sensitiveVsInsensitive$class=="positive",]
z1sensitiveVsInsensitive.neg <- z1sensitiveVsInsensitive[z1sensitiveVsInsensitive$class=="negative",]

z8sensitiveVsInsensitive.pos <- z8sensitiveVsInsensitive[z8sensitiveVsInsensitive$class=="positive",]
z8sensitiveVsInsensitive.neg <- z8sensitiveVsInsensitive[z8sensitiveVsInsensitive$class=="negative",]

z1sensitiveVsz8sensitive.pos <- z1sensitiveVsz8sensitive[z1sensitiveVsz8sensitive$class=="positive",]
z1sensitiveVsz8sensitive.neg <- z1sensitiveVsz8sensitive[z1sensitiveVsz8sensitive$class=="negative",]

z1sensitiveVsz8sensitive.pos$class <- "z1"
z1sensitiveVsz8sensitive.neg$class <- "z8"
z1sensitiveVsInsensitive.neg$class <- "insensitive"

z1.z8.insensitive <- rbind(z1sensitiveVsz8sensitive.pos, z1sensitiveVsz8sensitive.neg,
                           z1sensitiveVsInsensitive.neg)
newrowname <- as.vector(z1.z8.insensitive$name)
write.table(z1.z8.insensitive, file="/binf-isilon/sandelin/people/mengjun/Exosome_ML/Z1sensitiveVSZ8sensitiveVSallinsensitive/classid.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)
####features
datadir1 <- z1sins
datadir2 <- z8sins 
datadir3 <- z1z8s 
output <- "/binf-isilon/sandelin/people/mengjun/Exosome_ML/Z1sensitiveVSZ8sensitiveVSallinsensitive/features/"
filename <- list.files(path=paste0(z8sins, "features/"),
                       pattern="*.rda")

neatload <- function(dataname){
  temp.space <- new.env()
  bar <- load(dataname, temp.space)
  the.object <- get(bar, temp.space)
  rm(temp.space)
  the.object
}
extractFeatures <- function(filename){
  for(i in filename){
    res1 <- neatload(paste0(datadir1, "features/", i))
    res2 <- neatload(paste0(datadir2, "features/", i))
    res3 <- neatload(paste0(datadir3, "features/", i))
    res <- rbind(res1, res2, res3)
    ct <- colnames(res)
    res <- data.frame(res[newrowname, ])
    colnames(res) <- ct
    row.names(res) <- newrowname
    save(res, file=paste0(output, i))
  }
   
}

extractFeatures(filename)
