# @author Dinesh Barupal dinkumar@ucdavis.edu
# @version 1.0 July 2016
# @use call "getMetaMapp" function via a Ajax call from JS.

getChemSimNet <- function (cids, cutoff=0.7) {
  #cids <- c(1:5)
  #cutoff <- 0.7
  if (length(cids)>200) {
      hex <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f")
      bin <- c("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111")
      cidlist <- split(cids, ceiling(seq_along(cids)/200))
      subkeys <- do.call(rbind,lapply(cidlist,function(x) { read.csv(paste(c('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',paste(x,collapse=","),'/property/Fingerprint2D/csv'),collapse=""))}))
      m <- t(sapply(subkeys[,2],function(x) { as.integer(strsplit(paste(sapply(strsplit(paste(RCurl::base64Decode(x,"raw")[5:115],collapse=""),"")[[1]],function(x) {bin[which(hex==x)]}),collapse=""),"")[[1]][1:881] ) }))
      mat <- m%*%t(m)
      len <- length(m[,1])
      s <- mat.or.vec(len,len)
      for (i in 1:len) {
        for (j in 1:len){
          s[i,j] <- mat[i,j]/(mat[i,i]+mat[j,j]-mat[i,j])
        }
      }
      diag(s) <- 0
      dfmax <- cbind(cids,"tmsim",cids[sapply(1:length(cids), function (k) {which.max(s[k,]) })])
      s[lower.tri(s)]<-0
      chemsimdf <- do.call(rbind,sapply(1:length(cids), function (k) { if(length(which(s[k,]>cutoff))>0) {cbind(cids[k],"tmsim",cids[which(s[k,]>cutoff)])}} ))
      chemsimdf <- rbind(chemsimdf, cbind(cids,"tmsim",""), dfmax )

      ## Duplicated edges removal
      chemsimdf <- chemsimdf[-which(chemsimdf[,3]==""),]
      chemsimdf <- rbind(chemsimdf,cbind(chemsimdf[,3], chemsimdf[,2], chemsimdf[,1]))
      chemsimdf <-  chemsimdf[!duplicated( chemsimdf),]
      chemsimdf <- chemsimdf[order(chemsimdf[,3],decreasing = F),]
      chemsimdf <- chemsimdf[order(chemsimdf[,1],decreasing = F),]

      pmids_a <- cids

      for (i in 1:length(pmids_a)) {
        sind <- c((max(which(chemsimdf[,1]==pmids_a[i])) +1) :nrow(chemsimdf))
        chemsimdf[,3][which(chemsimdf[,3][sind]==pmids_a[i]) + (sind[1]-1) ] <- "XX"
      }
      chemsimdf <- chemsimdf[-which(chemsimdf[,3]=="XX"),]


      write.table(chemsimdf,file=paste(c("chemsim_",gsub("[.]","",as.character(cutoff)),".sif"),collapse=""), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)  ## To write the cytoscape network file as an output
      return(chemsimdf)
  } else{

  hex <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f")
  bin <- c("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111")
  subkeys <- read.csv(paste(c('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',paste(cids,collapse=","),'/property/Fingerprint2D/csv'),collapse=""))
  m <- t(sapply(subkeys[,2],function(x) { as.integer(strsplit(paste(sapply(strsplit(paste(RCurl::base64Decode(x,"raw")[5:115],collapse=""),"")[[1]],function(x) {bin[which(hex==x)]}),collapse=""),"")[[1]][1:881] ) }))
  mat <- m%*%t(m)
  len <- length(m[,1])
  s <- mat.or.vec(len,len)
  for (i in 1:len) {
    for (j in 1:len){
      s[i,j] <- mat[i,j]/(mat[i,i]+mat[j,j]-mat[i,j])
    }
  }
  diag(s) <- 0
  dfmax <- cbind(cids,"tmsim",cids[sapply(1:length(cids), function (k) {which.max(s[k,]) })])
  s[lower.tri(s)]<-0
  chemsimdf <- do.call(rbind,sapply(1:length(cids), function (k) { if(length(which(s[k,]>cutoff))>0) {cbind(cids[k],"tmsim",cids[which(s[k,]>cutoff)])}} ))
  chemsimdf <- rbind(chemsimdf, cbind(cids,"tmsim",""), dfmax )

  ## Duplicated edges removal
  chemsimdf <- chemsimdf[-which(chemsimdf[,3]==""),]
  chemsimdf <- rbind(chemsimdf,cbind(chemsimdf[,3], chemsimdf[,2], chemsimdf[,1]))
  chemsimdf <-  chemsimdf[!duplicated( chemsimdf),]
  chemsimdf <- chemsimdf[order(chemsimdf[,3],decreasing = F),]
  chemsimdf <- chemsimdf[order(chemsimdf[,1],decreasing = F),]

  pmids_a <- cids

  for (i in 1:length(pmids_a)) {
    sind <- c((max(which(chemsimdf[,1]==pmids_a[i])) +1) :nrow(chemsimdf))
    chemsimdf[,3][which(chemsimdf[,3][sind]==pmids_a[i]) + (sind[1]-1) ] <- "XX"
  }
  chemsimdf <- chemsimdf[-which(chemsimdf[,3]=="XX"),]

  write.table(chemsimdf,file=paste(c("chemsim_",gsub("[.]","",as.character(cutoff)),".sif"),collapse=""), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)  ## To write the cytoscape network file as an output
  return(chemsimdf)
}
}
## Tanimoto score calculation was adopted from http://data2quest.blogspot.com/2013/10/fast-tanimoto-similarity-calculation.html

krp <- read.table("KRPlinks.txt", sep="\t")

getKEGGRpairs <- function (cids, keggids, cutoff=0.7) {
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

getMetaMapp <- function (rsess, pvalind=1, cutoff=0.7) {
  statres <- paste(rsess,"R/.val/csv",sep="")
  df <- read.csv(statres)
  getKEGGRpairs(df$PubChem[is.na(df$PubChem)==FALSE], as.character(df$KEGG[is.na(df$PubChem)==FALSE]), cutoff)
  fcdf <- df[which(is.na(df$PubChem)==FALSE),grep("^Mean",colnames(df))]
  fcdf$meanratio <- fcdf$Mean.of.t1...None / fcdf$Mean.of.t2...Test.Compound
  pvaldf <- df[which(is.na(df$PubChem)==FALSE),grep("p_value",colnames(df))]
  exportdf <- data.frame(Pubchem=df$PubChem[is.na(df$PubChem)==FALSE])
  exportdf$KEGG <- as.character(df$KEGG[is.na(df$PubChem)==FALSE])
  exportdf$fold.change <- rep(1,length(exportdf[,1]))
  exportdf$pvaldirection <- rep("No Change",length(exportdf[,1]))
  exportdf <- cbind(exportdf,df[which(is.na(df$PubChem)==FALSE),grep("p_value",colnames(df))])
  exportdf$CpdName <- as.character(df$BinBase_name[is.na(df$PubChem)==FALSE])
  sigind <- which(exportdf[,grep("p_value",colnames(exportdf))[1]]<0.05)
  for( x in sigind)  {
    if(fcdf$meanratio[x]<1) {
      exportdf$fold.change[x] <- round(1/fcdf$meanratio[x],1)
      exportdf$pvaldirection[x] <- "Down"
      } else {
        exportdf$fold.change[x] <- round(fcdf$meanratio[x],1)
        exportdf$pvaldirection[x] <- "Up"
        }
  }
  #exportdf$meshanno <- sapply(exportdf$Pubchem, function (x) { paste(fromJSON( paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pccompound&db=mesh&id=",x,"&retmode=json",sep="") )$linksets[[1]]$linksetdbs[[1]]$links,collapse=",")  }   )
  write.table( exportdf, file=paste("node_attributes_chemsim_krp_",gsub("[.]","",as.character(cutoff)) ,".tsv", sep="" ), col.names = T, row.names = F, quote = F, sep = "\t" )
}

runMetaMapp <- function(stat_file, cutoff=0.7) {
  cfile <- strsplit(stat_file,"\n")[[1]]
  df1 <- do.call(rbind, lapply(cfile, function (x) { strsplit(x,"\t")[[1]]  } ))
  colnames(df1) <- df1[1,]
  df1 <- df1[-1,]
  getKEGGRpairs(df1[,1][is.na(df1[,1])==FALSE], df1[,2][is.na(df1[,1])==FALSE], cutoff)
  #exportdf <- data.frame(Pubchem_ID=df1[,1][is.na(df1[,1])==FALSE], KEGG_ID=df1[,2][is.na(df1[,1])==FALSE], CompoundName=df1[,3][is.na(df1[,1])==FALSE])
  exportdf <- as.data.frame(setNames(replicate(  length( grep("foldchange",colnames(df1) ) ), rep("No Change",length(df1[,1])), simplify = F), paste(colnames(df1)[grep("foldchange",colnames(df1))],"_direction",sep="") ), stringsAsFactors=FALSE)
  for (k in grep("pvalue",colnames(df1)  )) {
    df1[,k] <- as.numeric(df1[,k])
    sigind <- which(df1[,k]<0.05)
    df1[which(1:length(df1[,1])%in%sigind==FALSE),(k+1)] <- 1.0  ## convert all the non-significant fold changes to 1.00.
    for( x in sigind)  {
      if(df1[x,(k+1)]<1) {
        df1[x,(k+1)] <- round(1/as.numeric(df1[x,(k+1)] ),1)
        exportdf[x,paste(colnames(df1)[k+1],"_direction",sep="")] <- "Down"
      } else {
        df1[x,(k+1)] <- round(as.numeric(df1[x,(k+1)] ),1)
        exportdf[x,paste(colnames(df1)[k+1],"_direction",sep="")] <- "Up"
      }
    }
  }
  exportdf <- cbind(df1, exportdf)
  write.table( exportdf, file=paste("node_attributes_chemsim_krp_",gsub("[.]","",as.character(cutoff)) ,".tsv", sep="" ), col.names = T, row.names = F, quote = F, sep = "\t" )
}





