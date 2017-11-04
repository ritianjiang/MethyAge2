# load("/home/owht/KIZ/data/testData/methylation/testMethyperc.RData")
# test<-data.frame(testAAA)


#
# sampleInfo<-read.csv("/home/owht/KIZ/data/testData/methylation/Longevity_families_CM_mRNA_WN_mRNA_lncRNA.csv")
# sampleInfo$expre<-rep(0,nrow(sampleInfo))
# sampleInfo$sub<-as.character(sampleInfo$sub)
# aM<-t(test[,-(1:5)])
# rownames(aM)<-str_extract_all(rownames(aM),"\\d+")
# sampleInfo$expre<-aM[sampleInfo$sub,1]
# null.model<-lme(expre~gender+site,
#                 random = ~1|sub,
#                 data = sampleInfo,method = "ML")
# full.model<-lme(expre~gender+site+age,
#                 random = ~1|sub,
#                 data = sampleInfo,method = "ML")
# fit<-anova(null.model,full.model)
# result[1,]$p.value = fit$`p-value`[2]
# result[1,]$slope = full.model$coefficients$fixed[4]

#############################S4 Objects######################

#'This is a stupid DataStructure that storage the sample infomation
#'and the Methylation perc data
#'
#'If you want to create one MetnAge, you should use the function
#'named makeMetnAge.
#'@slot MetaInfo The methylation status for each site
#'@slot SampleInfo The sample info always include age/subid/gender ...
#'@slot locInfo The location infomation of each site (chr, bp, strand)
#'@author Wang Haotian (ritianjiang@gmx.com)
#'@export
#'@rdname MetnAge-class
#'@docType class

setClass("MetnAge",slots = list(MetaInfo = "matrix",
                                SampleInfo = "data.frame",
                                locInfo = "data.frame"))
###############S4 objects END##################################

###############Functions Start###########################

#'This function create the MetnAge object
#'
#'If you want to create one MetnAge, you should use the function
#'named makeMetnAge.
#'@export
#'@param metainfo a dataframe contains the Methylation perc, the
#'                dataframe can be got from the \code{methylKit::percMethylation} -> GRanges
#'                The first 5 cols should be chr,seqnames,start,end,width.
#'                In fact the dataframe is construct from GRanges
#'@param sampleinfo a dataframe contain the sampleinfo. It must include ages.
#'
#'@return a Metnage object
#'@rdname makeMethAge-method
#'@docType methods
#'@examples
#'data("ExampleMethnAge")
#'testobj<-makeMetnAge(test,sampleInfo)
#'

setGeneric("makeMetnAge",function(metainfo,sampleinfo){
  standardGeneric("makeMetnAge")
})

setMethod("makeMetnAge",
          signature(metainfo = "data.frame",sampleinfo = "data.frame"),
          function(metainfo,sampleinfo){
              metainfo1<-t(metainfo[,-(1:5)])
              locinfo<-as.data.frame(metainfo[,1:5])
              rownames(metainfo1)<-str_extract_all(rownames(metainfo1),
                                          "\\d+")
              result<-new("MetnAge",MetaInfo = metainfo1,
                           SampleInfo = sampleinfo,
                           locInfo=locinfo)
              return(result)
})



#This is an invisible(unexport) function. To be the unit of
#ageMetCorr
ageM_each<-function(sampleinfo,am,ni){
  sampleinfo$expre<-am[sampleinfo$sub,ni]
  null.model<-lme(expre~gender+site,
                  random = ~1|sub,
                  data = sampleinfo,method = "ML")
  full.model<-lme(expre~gender+site+age,
                  random = ~1|sub,
                  data = sampleinfo,method = "ML")
  fit<-anova(null.model,full.model)
  result<-c(fit$`p-value`[2],full.model$coefficients$fixed[4])
  return(result)
}


#'@title ageMethCorr
#'@description  This function calculate the roles that age played in methylation
#'              status.we used the mixed linear effect model to calculate thep-value and slope of age by nlme. The null model is
#'              \code{Methperc~gender+site+(1|sub)} and the full mode is
#'              \code{Methperc~gender+site+age+(1|sub)}. To get the LRT, we
#'              adapt the "ML"  instead of "REML" to estimate paramaters
#'
#'@export
#'@param objects MetnAge object
#'@param cors a numeric value of numbers of cores you wanna use. (default:1)
#'
#'@return a dataframe contains site location info and p-value and
#'@rdname ageMetCorr-method
#'@docType methods
#'@examples
#'data("ExampleMethnAge")
#'testobj<-makeMetnAge(test,sampleInfo)
#'cls<-ageMetCorr(testobj,cors = 2)
#'
setGeneric("ageMetCorr",function(object,cors){
  standardGeneric("ageMetCorr")
})

setMethod("ageMetCorr","MetnAge",function(object,cors=1){
    slmSampleInfo<-object@SampleInfo
    slmaM<-object@MetaInfo
    slmLoc<-object@locInfo

    if(cors>1){
      cl<-makeForkCluster(cors)
      a<-pblapply(1:ncol(slmaM),FUN=ageM_each,
                  sampleinfo = slmSampleInfo,
                  am = slmaM,cl = cl)
      stopCluster(cl)
    }
    else{a<-lapply(1:ncol(slmaM),FUN=ageM_each,
                   sampleinfo = slmSampleInfo,
                   am = slmaM)}
    a<-do.call(rbind,a)
    result<-cbind(slmLoc,a)
    colnames(result)[6]<-"pValue"
    colnames(result)[7]<-"slope"
    return(result)

})
