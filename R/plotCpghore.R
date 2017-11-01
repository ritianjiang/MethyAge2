#'The function deals the cpg methylation status
#'
#' @param object: A GRanges object that is the methylation perc file
#' @param grAnno: A GRanges object that is the cpgi&shores file
#' @param plot: ... bool, plot or not (default: F)
#' @param len: ... int, The length of cpgi you wanna normalized (default: 1000)
#' @export
#' @return  A dataframe which contains cpgi & shore methylation status. the
#' @examples
#' data(ExampleCpgPlot)
#' result<-plotCpghore(testMeatrix,testCPG)
#'
#'
#'

setGeneric("plotCpghore",
           function(object,cpgFile,...){
             standardGeneric("plotCpghore")})
setMethod("plotCpghore","GRanges",
          function(object,cpgFile,plot = F,len = 1000){
            MethAnno<-DMSAnnotate(object,cpgFile)
            MethAnno<-unique(MethAnno)

            down<-MethAnno[MethAnno$ann =="down200",]
            up<-MethAnno[MethAnno$ann == "up200",]
            cpg<-MethAnno[MethAnno$ann != "up200" & MethAnno$ann != "down200",]

            down$dist<-down[,ncol(down)-3] - down$start
            down_mean<-aggregate(down[,6:(ncol(down)-7)],by = list(down$dist),FUN = mean)
            down_mean$Group.1<-down_mean$Group.1 + len

            cpg$dist<-((cpg[,ncol(cpg)-3] - cpg$start)/(cpg[,ncol(cpg)-2]))*len
            cpg$dist<-as.integer(cpg$dist)
            cpg_mean<-aggregate(cpg[,6:(ncol(cpg)-7)],by = list(cpg$dist),FUN = mean)

            up$dist<-up[,ncol(up)-4] - up$start
            up_mean<-aggregate(up[,6:(ncol(up)-7)],by = list(up$dist),FUN = mean)


            totalcpg<-rbind(up_mean,cpg_mean)
            totalcpg<-rbind(totalcpg,down_mean)
            totalcpg<-aggregate(totalcpg[,-1],by = list(totalcpg$Group.1),FUN = mean)
            totalcpg_app<-data.frame(loc = totalcpg$Group.1
                                     ,meth = apply(totalcpg,MARGIN = 1,FUN = mean))
            if(plot == T){
              silda<-filter((totalcpg_app[,2])/10,rep(1,10))
              silda<-data.frame(loc = totalcpg_app$loc,
                                meth = silda)
              p<-ggplot(totalcpg_app,aes(loc,meth))

              p<-p+geom_point()+geom_line(data=silda,col = "green1")

              return(list(totalcpg_app,p))}


            return(totalcpg_app)
          })



