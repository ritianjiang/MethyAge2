#'Plot methylation status across genic regions(or some other features) and their up/down stream.
#'
#' @param object: A GRanges object contains the methylsites
#' @param grAnno: A TxDb object containing the features. Please provide the strand info if you are ploting
#' @param numbin: bin number to split the feature regions we wanna plot
#' @param wingbin: bin number of the up/down stream
#' @param upstream: numeric; upstream bp to plot; default:2000
#' @param downstream: numeric; downstream bp to plot; default:500
#' @param scoreColumnName: character; The name of score in object.
#' @param Strand: boolean; Wheather use the strand info of object. default: F.
#' @export
#' @return  A dataframe which contains site location and annotation. In the dataframe, there won't be the metainfo like p-values etc.
#' @examples
#' data(ExampleAnnotate)
#' result<-DMSAnnotate(methylDiff.obj, ExampleGencode)
#'
#'
#'
plotGenic<-function(object,grAnno,numbin,wingbin,upstream=2000,
                    downstream=500,scoreColumnName,Strand=FALSE){
  cat("Generating Annotation files###");cat("\n")
  genicRegion<-genes(grAnno) %>% as.data.frame()
  genicRegion$gene_id<-"Gene Body"
  if(Strand==F){object$strand<-"*"}

  ups<-genicRegion
  ups[ups$strand =="+",]$end<-ups[ups$strand =="+",]$start-1
  ups[ups$strand =="+",]$start<-ups[ups$strand =="+",]$start-(upstream+1)
  ups[ups$strand =="-",]$start<-ups[ups$strand =="-",]$end+1
  ups[ups$strand =="-",]$end<-ups[ups$strand =="-",]$end+(upstream+1)
  ups$gene_id<-"Up"

  dws<-genicRegion
  dws[dws$strand =="+",]$start<-dws[dws$strand =="+",]$end+1
  dws[dws$strand =="+",]$end<-dws[dws$strand =="+",]$end+(downstream+1)
  dws[dws$strand =="-",]$end<-dws[dws$strand =="-",]$start-1
  dws[dws$strand =="-",]$start<-dws[dws$strand =="-",]$start-(downstream+1)
  dws$gene_id<-"Down"

  gR<-rbind(genicRegion,ups,dws)
  gR<-as(gR,"GRanges")
  cat("Annotating....");cat("\n")
  anno<-DMSAnnotate(object,gR)
  anno1<-anno[,c(2,8,9,11,12)]
  anno1$score<-anno[,scoreColumnName]
  anno1$loc<-(anno1$start - anno1$start.1)/(anno1$end - anno1$start.1)
  anno1[anno1$strand=="-",]$loc<- 1-anno1[anno1$strand == "-",]$loc
  anno1<-anno1[,c(5,6,7)]
  colnames(anno1)[2]<-"score"
  anno1[anno1$gene_id == "Gene Body",]$loc<-anno1[anno1$gene_id == "Gene Body",]$loc %/% (1/numbin)
  anno1[anno1$gene_id != "Gene Body",]$loc<-anno1[anno1$gene_id != "Gene Body",]$loc %/% (1/wingbin)
  anno2<-aggregate(x = anno1$score,by=list(anno1$gene_id,anno1$loc),FUN=mean)
  return(anno2)
}
