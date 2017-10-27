#'The function annotates the methylDiff/methylBase obj.
#'
#' @param object: A methylDiff/methyBase object
#' @param grAnno: A GRanges object
#' @export
#' @return  A dataframe which contains site location and annotation. In the dataframe, there won't be the metainfo like p-values etc.
#' @examples
#' data(ExampleAnnotate)
#' result<-DMSAnnotate()
#'
#'
#'

setGeneric("DMSAnnotate",function(object,grAnno){standardGeneric("DMSAnnotate")})

setMethod("DMSAnnotate","methylDiff",
           function(object,grAnno){
                 data<-getData(object);
                 dataGR<-with(data,GRanges(chr,IRanges(start,end),strand=strand))
                 data.anno<-findOverlaps(dataGR,grAnno)
                 annoGrange<-data.frame(grAnno[data.anno@to])
                 Result<-cbind(dataGR[data.anno@from],annoGrange)
                 return(Result)
           })

setMethod("DMSAnnotate","methylBase",
          function(object,grAnno){
            data<-getData(object);
            dataGR<-with(data,GRanges(chr,IRanges(start,end),strand=strand))
            data.anno<-findOverlaps(dataGR,grAnno)
            annoGrange<-data.frame(grAnno[data.anno@to])
            Result<-cbind(dataGR[data.anno@from],annoGrange)
            return(Result)
          })
