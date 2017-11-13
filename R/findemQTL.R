# load("~/KIZ/data/HNLong/Rdata/Wht/off_SP_all/exp_meth/raw_meth.RData")
# load("/home/owht/KIZ/data/HNLong/Rdata/fpkm_merge.RData")
# sampleinfo<-read.table("~/KIZ/data/testData/methylation/sampleinfo.csv",header = T)
# fpkm<-fpkm[,colnames(fpkm) %in% as.character(sampleinfo$sub)]
# aaa<-data.frame(aaa)
# testExp<-fpkm[1:200,]
# testMet<-aaa[1:5000,]
# result<-

##Data clean
# CM<-read.csv("/home/owht/KIZ/data/HNLong/FPKM_Value_Cufflinks_58mRNA.csv")
# WN<-read.csv("/home/owht/KIZ/data/HNLong/FPKM_Value_Cufflinks_63mRNA.csv")
# totalFPKM<-merge(CM,WN,by.x = "X",by.y = "X")
#
# load("/home/owht/KIZ/data/HNLong/Rdata/gencode_hg19.RData")
# pcGenes<-unique(gencode.hg19[gencode.hg19$gene_type=="protein_coding",]$gene_name)
# pcGenes<-as.character(pcGenes)
# cdFPKM<-totalFPKM[totalFPKM$X %in% pcGenes,]
# sampleinfo<-read.table("~/KIZ/data/testData/methylation/sampleinfo.csv",header = T)
#
# gencode.hg19<-gencode.hg19[gencode.hg19$gene_type == "protein_coding",]
# chrGene<-data.frame(gene=as.character(gencode.hg19$gene_name),
#                     chr = as.character(gencode.hg19$seqnames),
#                     stringsAsFactors = F)
# chrGene<-unique(chrGene)
# p<-table(chrGene$gene)
# p<-data.frame(p,stringsAsFactors = F)
# chrGene<-chrGene[chrGene$gene %in% p[p$Freq<2,]$Var1,]
# rownames(chrGene)<-chrGene$gene
# cdFPKM$chr<-chrGene[as.character(cdFPKM$X),]$chr
# cdFPKM<-na.omit(cdFPKM)
# write.csv(cdFPKM,file = "/home/owht/KIZ/data/HNLong/Rdata/cdFPKM.csv")

# cdFPKM<-read.csv("/home/owht/KIZ/data/HNLong/Rdata/cdFPKM.csv")
# load("~/KIZ/data/HNLong/Rdata/Wht/off_SP_all/exp_meth/raw_meth.RData")
# test<-subset(aaa,seqnames == "chr5")
# test<-data.frame(test)
# slmLoc<-test[,1:5]
# test<-test[,-(1:5)]
# sub<-colnames(test)
# library(stringr)
# sub<-str_extract_all(sub,pattern = "[0-9]+")
# sub<-do.call(rbind,sub)
# colnames(test)<-paste0("X",sub)
# cdFPKM<-cbind(cdFPKM$X,
#               cdFPKM[,colnames(cdFPKM) %in% paste0("X",sub)],
#               cdFPKM$chr)
#
# test<-test[,colnames(chr5)[2:45]]
#
# chr5<-cdFPKM[cdFPKM$`cdFPKM$chr` == "chr5",]
# write.csv(cdFPKM,file ="/home/owht/KIZ/data/HNLong/Rdata/meth_45_FPKM.csv")
# colnames(chr5)[1]<-"X"

#This is an inherent function to select the most correlated gene
#for each CpG site
testFun<-function(j,met,exp,site,qFlag = F){ #based on CpG site
  slmLoc<-site
  # result<-data.frame(Gene = as.character(exp$X), ##noticing
  #                    p = rep(0,nrow(exp)),
  #                    slope = rep(0,nrow(exp)))
  ptest<-lapply(1:nrow(exp),FUN=corR,
                Met = met[j,],Exp = exp)
  ptest<-do.call(rbind,ptest)
  result<-na.omit(ptest)
  colnames(result)<-c("Gene","p","slope")
  result<-data.frame(result,stringsAsFactors = F)

  if(qFlag==T){
     result$p<-p.adjust(result$p)}

  result<-result[order(result$p,decreasing = F),]
  result<-cbind(site[j,],result[1,])
  return(result);
}

#This is an inherent function to calculate correlation between A gene
#and A CpG site
#for each CpG site
corR<-function(index,Met,Exp){ #based on Expression genes
  # result<-data.frame(Gene = as.character(Exp$X), ##noticing
  #                    p = rep(0,nrow(Exp)),
  #                    slope = rep(0,nrow(Exp)))
  fit<-cor.test(as.numeric(Met),
                as.numeric(Exp[index,2:(ncol(Exp)-1)]),
                method = "pearson")
  res<-c(as.character(Exp[index,]$X),
         as.numeric(fit$p.value),
         as.numeric(fit$estimate))
  return(res);
}


#'@title findemQTL
#'@description  This function find the potential emQTL by correlation test.
#'
#'@export
#'@param met   a data.frame contains methylation percentation
#'@param exp   a data.frame contains gene and fpkm(or some others)
#'             The 1st column is the gene names and the last is chrom
#'@param cors  a numeric value of numbers of cores you wanna use. (default:1)
#'@param siteFile  a dataframe contains the CpG site info
#'@param aFlag logical(default:F). Whether to use FDR correction.
#'@return a dataframe contains site location info and potential regulated gene
#'and p-value(or qvalue in fact) and slope
#'@rdname findemQTL-method
#'@docType methods
#'@examples
#'data("ExampleemQTL")
#'res<-findemQTL(test[1:20,],chr5,cors = 4,siteFile = siteInfo)
#'res
#'
setGeneric("findemQTL",function(met,exp,cors=1,siteFile,qFlag = F){
  standardGeneric("findemQTL")
})

setMethod("findemQTL","data.frame",
          function(met,exp,cors=1,siteFile,qFlag){
  if(cors>1){
    cl<-makeForkCluster(cors)
    a<-pblapply(1:nrow(met),FUN=testFun,
                met = met,exp = exp,site = siteFile,
                qFlag = qFlag,
                cl = cl)
    stopCluster(cl)
  }
  else{
    a<-lapply(1:nrow(met),FUN=testFun,
                met = met,exp = exp,site = siteFile)
  }
  a<-do.call(rbind,a)
  return(a);
})