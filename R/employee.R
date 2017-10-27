#' Generate a employee object.
#' @param name the name, defalut is "Joe"
#' @param salary the salary, default is 5000
#' @param union  boolean whether the man is in Union, defualt is T
#' @examples
#' Wht<-employee(name="Wht",salary = 6000)
#' Wht
employee<-function(name = "Joe", salary = 5000, union = T,
                   info = data.frame(x=rep(1,10),y=rep(2,10))){
  joe<-new("employee",name = name, salary = salary,
           union = union,info = info)
  return(joe)
}
