#This is a inherent function
corrSiteEach<-function(i,df1,df2,method){
  cls<-cor.test(x = as.numeric(df1[i,]),
                y = as.numeric(df2[i,]),method = method)
  res<-c(rownames(df1[i,]),
         as.numeric(cls$estimate),
         as.numeric(cls$p.value))
  return(res)
}

#'@title findParentage
#'@description  This function find the potential heritable sites between
#'              p and f1 using correlation test.
#'
#'@export
#'@param DF1   a data.frame contains methylation percentation
#'@param DF2   a data.frame contains methylation percentation
#'@param cors  a numeric value of numbers of cores you wanna use. (default:1)
#'@param method whether "spearman" or "pearson
#'@return a dataframe contains pvalue,corr,location
#'@examples
#' data("ExampleParent")
#' aaa<-findParentage(Father,Son)
#' aaa


setGeneric("findParentage",function(DF1,DF2,cors=4,method = "spearman"){
  standardGeneric("findParentage")
})

setMethod("findParentage","data.frame",function(DF1,DF2,cors,method){
  if(nrow(DF1)!=nrow(DF2)){
    print("Error,The 2 df dont have the same length!")
    return(0);
  }
  if(cors>1){
    cl<-makeForkCluster(cors)
    a<-pblapply(1:nrow(DF1),FUN = corrSiteEach,df1 = DF1,df2=DF2,method = method)
    stopCluster(cl)}
  else{
    a<-lapply(1:nrow(DF1),FUN=corrSiteEach,df1=DF1,df2=DF2,method = method)
  }
  a<-do.call(rbind,a)
  colnames(a)<-c("loc","slope","pvalue")
  a<-as.data.frame(aaa)
  a<-na.omit(a)

  return(a)
})
