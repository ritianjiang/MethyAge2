#'The function annotates the methylDiff/methylBase obj.
#'
#' @param object: A methylDiff/methyBase/GRanges object
#' @param grAnno: A GRanges object
#' @export
#' @return  A dataframe which contains site location and annotation. In the dataframe, there won't be the metainfo like p-values etc.
#' @examples
#' data(ExampleAnnotate)
#' result<-DMSAnnotate(methylDiff.obj, ExampleGencode)
#'
#'
#'

setGeneric("DMSAnnotate",function(object,grAnno){standardGeneric("DMSAnnotate")})

setMethod("DMSAnnotate","methylDiff",
          function(object,grAnno){
            data<-getData(object);
            dataGR<-with(data,GRanges(chr,IRanges(start,end),strand=strand,diff = meth.diff))
            data.anno<-findOverlaps(dataGR,grAnno)
            annoGrange<-data.frame(grAnno[data.anno@to])
            Result<-cbind(dataGR[data.anno@from],annoGrange,stringsAsFactors = F)
            gc()
            return(Result)
          })

setMethod("DMSAnnotate","methylBase",
          function(object,grAnno){
            data<-getData(object);
            dataGR<-with(data,GRanges(chr,IRanges(start,end),strand=strand,diff = meth.diff))
            data.anno<-findOverlaps(dataGR,grAnno)
            annoGrange<-data.frame(grAnno[data.anno@to])
            Result<-cbind(dataGR[data.anno@from],annoGrange,stringsAsFactors = F)
            gc();
            return(Result)
          })

setMethod("DMSAnnotate","GRanges",
          function(object,grAnno){
            data.anno<-findOverlaps(object,grAnno)
            annoGrange<-data.frame(grAnno[data.anno@to])
            Result<-cbind(object[data.anno@from],annoGrange,stringsAsFactors = F)
            gc();
            return(Result)
          })
