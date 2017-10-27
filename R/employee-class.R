#' Am example

setClass("employee",representation(
      name = "character",
      salary = "numeric",
      union = "logical",
      info = "data.frame"
  )
)

setMethod("show","employee",
          function(object){
            inorout<-ifelse(object@union,"is","is not")
            cat(object@name,"has a salary of", object@salary,
                "and",inorout,"in the union","\n")
          })
