#'The function get the correlation between methylation and expression
#'
#' @param df1: A dataframe which contains the methylation,1st column is gene Symbol(or refseq name, gene id etc)
#' @param df2: A dataframe which contains the expression,1st column is gene Symbol(or refseq name, gene id etc)
#' @export
#' @return  A dataframe which contains the p-value & corr coefficients
#' @examples
#' data("ExampleCorr")
#' result<-getMethyVsExpression(exMethyMatrix,exExpMatrix)
#'
#'
#'

getMethyVsExpression<-function(df1,df2){
  dfx<-df1[,-1]
  dfx<-dfx[,order(colnames(dfx))]
  df1<-cbind(df1[,1],dfx)
  dfy<-df2[,-1]
  dfy<-dfy[,order(colnames(dfy))]
  df2<-cbind(df2[,1],dfy)

  if(nrow(df1)!=nrow(df2)){print("The gene numbers differs, please check your input")}
  if(!all(df1[,1]==df2[,1])){print("The gene symbols don't match!!")}
  Result<-data.frame(Symbol=df1[,1],
                     Cor=rep(0,nrow(df1)),
                     p.Value = rep(0,nrow(df1)),stringsAsFactors = F)
  for (i in 1:nrow(df1)){
    fit<-cor.test(x = as.numeric(df1[i,-1]),
                  y = as.numeric(df2[i,-1]))
    Result$Cor[i]<-fit$estimate
    Result$p.Value[i]<-fit$p.value
  }
  Result<-Result %$% .[which(p.Value >0),]
  return (Result)
}
