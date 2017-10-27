setGeneric("Display",function(object){standardGeneric("Display")})

setMethod("Display","methylBase",
           function(object){
             data<-getData(object)
             return(data)
           })
