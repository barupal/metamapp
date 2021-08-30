# @author Dinesh Barupal dinkumar@ucdavis.edu
# @version August 30th 2021, optimized for google colab

inputfile <- "metamapp_input.xlsx"
cutoff = 0.7
download.file("https://raw.githubusercontent.com/barupal/metamapp/master/KRPlinks.txt", destfile = "KRPlinks.txt")
krp <- read.table("KRPlinks.txt", sep="\t")
ndf <- data.frame(readxl::read_xlsx(path = inputfile, sheet = "input"), stringsAsFactors = F)

getChemSimNet <- function (cids, cutoff=0.7,fps) {
  m <- fps
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

  #write.table(chemsimdf,file=paste(c("chemsim_",gsub("[.]","",as.character(cutoff)),".sif"),collapse=""), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)  ## To write the cytoscape network file as an output
  return(chemsimdf)
}
## Tanimoto score calculation was adopted from http://data2quest.blogspot.com/2013/10/fast-tanimoto-similarity-calculation.html

getKEGGRpairs <- function (cids, keggids,fps, cutoff=0.7) {
  krp.1 <- match(krp[,1],keggids)
  krp.2 <- match(krp[,2],keggids)
  krp.cbind <- cbind (krp.1,krp.2)
  krp.net <- subset(krp.cbind, krp.1!="NA" & krp.2!="NA")
  if(nrow(krp.net)==0) { # export only the chemical similarity map
    chemsim <- getChemSimNet(cids,cutoff,fps)
    write.table(chemsim,file=paste(c("chemsim_krp_",gsub("[.]","",as.character(cutoff)),".sif"),collapse=""), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)  ## To write the cytoscape network file as an output
  } else {
    cid.krp.2 <- cids[krp.net[,2]]
    cid.krp.1 <- cids[krp.net[,1]]
    krp.cid.net <- cbind(cid.krp.1,"krp",cid.krp.2)
    chemsim <- getChemSimNet(cids,cutoff,fps)
    krp.cid.net <- rbind(krp.cid.net,chemsim)
    write.table(krp.cid.net,file=paste(c("chemsim_krp_",gsub("[.]","",as.character(cutoff)),".sif"),collapse=""), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)  ## To write the cytoscape network file as an output
  }
}

runMetaMapp <- function(ndf) {
  df1 <- ndf
  if( length(table(is.na(df1$PubChem_ID)))==2 ) { stop("Missing PubChem CIDs") }
  if( length(which(table(df1$PubChem_ID)>1))!=0) { stop(paste0("Please remove the duplicate CIDs  -   ",paste(names(which(table(c(1,df1$PubChem_ID))>1)),collapse=";"),"   in the input file."))}
  if( (length(which(is.na(df1$SMILES)==TRUE)))) { stop("Missing SMILES codes. Please check the input.") }
  ## Calculate the fingerprints
  fps1 <- t(sapply(1:nrow(df1), function(x) {
    as.character(rcdk::get.fingerprint(rcdk::parse.smiles(df1$SMILES[x])[[1]],type="pubchem"))
  }))
  ### If any of the smiles codes are wrong. Break the code here.
  if( (length(which(fps1==0))>0)) { stop("Incorrect SMILES Code provided. Please check the input") }

  df1.bitmat <- do.call(rbind,lapply(fps1,function(x) as.integer(strsplit(x,"")[[1]][1:881])))

  getKEGGRpairs(df1[,1], df1[,2],df1.bitmat, cutoff)
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

runMetaMapp(ndf)
