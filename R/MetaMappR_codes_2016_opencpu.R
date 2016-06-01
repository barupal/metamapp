
getChemSimNet <- function (cids, cutoff=0.7) {
  library(RCurl)
  #cids <- c(1:5)
  #cutoff <- 0.7
  if (length(cids)>200) {
      hex <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f")
      bin <- c("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111")
      cidlist <- split(cids, ceiling(seq_along(cids)/200))
      subkeys <- do.call(rbind,lapply(cidlist,function(x) { read.csv(paste(c('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',paste(x,collapse=","),'/property/Fingerprint2D/csv'),collapse=""))}))
      m <- t(sapply(subkeys[,2],function(x) { as.integer(strsplit(paste(sapply(strsplit(paste(base64Decode(x,"raw")[5:115],collapse=""),"")[[1]],function(x) {bin[which(hex==x)]}),collapse=""),"")[[1]][1:881] ) }))
      mat <- m%*%t(m)
      len <- length(m[,1])
      s <- mat.or.vec(len,len)
      for (i in 1:len) {
        for (j in 1:len){
          s[i,j] <- mat[i,j]/(mat[i,i]+mat[j,j]-mat[i,j])
        }
      }
      diag(s) <- 0
      s[lower.tri(s)]<-0
      chemsimdf <- do.call(rbind,sapply(1:length(cids), function (k) { if(length(which(s[k,]>cutoff))>0) {cbind(cids[k],"tmsim",cids[which(s[k,]>cutoff)])}} ))
      write.table(chemsimdf,file=paste(c("chemsim_",gsub("[.]","",as.character(cutoff)),".sif"),collapse=""), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)  ## To write the cytoscape network file as an output
      return(chemsimdf)
  } else{

  hex <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f")
  bin <- c("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111")
  subkeys <- read.csv(paste(c('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',paste(cids,collapse=","),'/property/Fingerprint2D/csv'),collapse=""))
  m <- t(sapply(subkeys[,2],function(x) { as.integer(strsplit(paste(sapply(strsplit(paste(base64Decode(x,"raw")[5:115],collapse=""),"")[[1]],function(x) {bin[which(hex==x)]}),collapse=""),"")[[1]][1:881] ) }))
  mat <- m%*%t(m)
  len <- length(m[,1])
  s <- mat.or.vec(len,len)
  for (i in 1:len) {
    for (j in 1:len){
      s[i,j] <- mat[i,j]/(mat[i,i]+mat[j,j]-mat[i,j])
    }
  }
  diag(s) <- 0
  s[lower.tri(s)]<-0
  chemsimdf <- do.call(rbind,sapply(1:length(cids), function (k) { if(length(which(s[k,]>cutoff))>0) {cbind(cids[k],"tmsim",cids[which(s[k,]>cutoff)])}} ))
  write.table(chemsimdf,file=paste(c("chemsim_",gsub("[.]","",as.character(cutoff)),".sif"),collapse=""), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)  ## To write the cytoscape network file as an output
  return(chemsimdf)
}
}
## Tanimoto score calculation was adopted from http://data2quest.blogspot.com/2013/10/fast-tanimoto-similarity-calculation.html

krp <- read.table("KRPlinks.txt", sep="\t")

getKEGGRpairs <- function (cids, keggids, cutoff) {
  krp.1 <- match(krp[,1],keggids)
  krp.2 <- match(krp[,2],keggids)
  krp.cbind <- cbind (krp.1,krp.2)
  krp.net <- subset(krp.cbind, krp.1!="NA" & krp.2!="NA")
  cid.krp.2 <- cids[krp.net[,2]]
  cid.krp.1 <- cids[krp.net[,1]]
  krp.cid.net <- cbind(cid.krp.1,"krp",cid.krp.2)
  chemsim <- getChemSimNet(cids,cutoff)
  krp.cid.net <- rbind(krp.cid.net,chemsim)
  write.table(krp.cid.net,file=paste(c("chemsim_krp_",gsub("[.]","",as.character(cutoff)),".sif"),collapse=""), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)  ## To write the cytoscape network file as an output
}


getNodeAttributes <- function (stat_file,samplecount) {
  names <- strsplit(strsplit(stat_file,"\n")[[1]][1],"\t")[[1]][-1]
  df1 <- do.call(rbind,lapply(strsplit(stat_file,"\n")[[1]][-1],function(x){as.integer(strsplit(x,"\t")[[1]])}))
  pvals <- sapply(1:length(df1[,1]),function(i){ t.test(df1[i,2:(samplecount+1)],df1[i,(samplecount+2):(samplecount+1+samplecount)],alternative = c("two.sided", "less", "greater"), var.equal=TRUE)$p.value  })
  cids <- df1[,1]
  statres <- do.call(rbind,lapply(1:length(pvals),function(z){ if(pvals[z]<0.05) {fc <- median(df1[z,2:(samplecount+1)])/median(df1[z,(samplecount+2):(samplecount+1+samplecount)]) ; if (fc>1.0){return( c(cids[z],pvals[z],round(fc,2),"Up")   )}else{  return( c(cids[z],pvals[z],round(1/fc,2),"Down")) }     } else {return(c(cids[z],pvals[z],1,"No Change"))} }))
  colnames(statres) <- c("CID","p.value","fold-change","direction")
  write.table(statres,file=paste(c("ttest_attr_",gsub("[.]| ","",names[1]),"_Vs_",gsub("[.]| ","",names[1+samplecount]),".txt"),collapse=""), quote=FALSE,sep="\t",col.names=T,row.names=FALSE)  ## To write the cytoscape network file as an output
}

hello <- function() {

}


