slidewindow<-function(x){
  start = x$Group.2 - 2 #generate the result matrix
  end<-start + 4
  result<-data.frame(chr = unique(x$Group.1),start,end)
  result2<-matrix(data = 0,nrow = nrow(result),ncol = ncol(x)-2) %>% as.data.frame()
  colnames(result2)<-colnames(x)[-c(1,2)]
  result<-cbind(result,result2)

  for(j in 1:nrow(result)){
    a<-x[x$Group.2 %in% seq(result[j,]$start,result[j,]$end),]
    if(nrow(a) == 0){
      result[j,-c(1,2,3)] = NA}
    else{
      result[j,-c(1,2,3)]<-apply(X=a[,-c(1,2)],MARGIN = 2,FUN = mean)}
    #if(i %% 1000 == 0){print(i)}
  }
  result<-na.omit(result)
  return(result)
}

#'The function find the DMVs defined by Xie et al in Cell 2013.
#'To identify each DMVs in a cell type, the genome was first divided in 1kb
#'bins and the DNA methylation level was averaged within each bin. Then a sliding
#'5kb window (with 1kb step) was used to identify regions that have an averaged
#'methylation level less than 0.15 in a 5kb window. Continuous regions resulting
#'from this analysis were then merged to form DMVs.
#'
#' @param methyMatrix: A matrix contains the methylation ratios. Columns are
#'                     samples and rows are sites. Each rowname is like 'chrX.NNNNN'
#'
#' @export
#' @return  A dataframe shows the methylation of 5kb-1kbstep regions of each sample.
#' @examples
#'
#'
#'




findDMVS<-function(methyMatrix){
  chrMatrix<-methyMatrix[methyMatrix %>% rownames
                         %>% str_sub(1,3) == "chr"
                         %>% head,]
  chrMatrix$chr<-(str_split(chrMatrix %>% rownames,
                           pattern = "\\.",n = 3,
                           simplify = T) %>% as.data.frame())$V1 %>% as.character()
  chrMatrix$loc<-(str_split(chrMatrix %>% rownames,
                            pattern = "\\.",n = 3,
                            simplify = T) %>% as.data.frame(stringsAsFactors = F))$V2 %>%as.numeric()
  chrMatrix$loc <- (chrMatrix$loc %/% 1000)+1
  avMatrix<-aggregate(x = chrMatrix[,1:(ncol(chrMatrix)-2)],
                      FUN = mean,
                      by = list(chrMatrix$chr,chrMatrix$loc))
  chrList<-unique(avMatrix$Group.1)
  pp<-data.frame()
  for(i in 1:length(chrList)){
    temp<-avMatrix[avMatrix$Group.1 == chrList[i],]
    avg<-slidewindow(temp)
    pp<-rbind(pp,avg)
    return(pp)
  }
}
