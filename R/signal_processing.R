#' area under the curve
#' 
#' take a numeric vector and return the area under the curve (AUC)
#' 
#' @param data : a numeric vector
#' @param pos : a logical value. if TRUE (default) AUC is calculated on positive values, if FALSE on negative values.
#' @return auc : area under the curve
#' @return max/min : maxiaml value (or minimal if pos=FALSE)
#' @export
auc<-function (data,pos=TRUE) {
  x<-seq_along(data)
  if (pos){
    data.pos<-data>0
    y<-data*data.pos
    auc<-sum(diff(x)*(y[-length(y)]+y[-1])/2)
    point<-max(data)
    return (list(auc=auc,max=point))
  }
  else {
    data.neg<-data<0
    y<-data*data.neg
    auc<-sum(diff(x)*(y[-length(y)]+y[-1])/2)
    point<-min(data)
    return (list(auc=auc,min=point))
  }
}

# Interactive area under the curve
# 
# take vector, ask for start and end points then return positiv and negativ area under the curve (AUC)
# 
# @param data : a numeric vector
# @return aucPos, aucNeg : positiv and negativ area under the curve,
# @return Max, Min : maximal and minimal values
# @export
#auc.interact<-function (data,to=0,interact=F,points=NA, ...) {
#   if (interact){
#           abs<-data.frame(x=points)
#   }
#   else {
#         X11()
#   print(plot.ts(data,type="l",...)) # plot data
#   abline(v=to,col = "red")
#   abline(v=to+5*60,col = "lightgray")
#   abline(v=to+10*60,col = "lightgray")
#   abline(v=to+15*60,col = "lightgray")
#   text(10,50,"line 10,15,20 min",col="gray60")
#   abs<-locator(n=2,type="o") #aske for window points
#   }
#   data<-data[abs$x[1]:abs$x[2]] #subset data
#   
#   x<-seq_along(data) #genetate abscisse sequence
#   
#   data.pos<-data>0 # select positive values
#   y.pos<-data*data.pos # select positive values
#   auc.pos<-sum(diff(x)*(y.pos[-length(y.pos)]+y.pos[-1])/2) 
#   max<-max(data)
# 
#   data.neg<-data<0
#   y.neg<-data*data.neg
#   auc.neg<-sum(diff(x)*(y.neg[-length(y.neg)]+y.neg[-1])/2)
#   min<-min(data)
#   
#   return (list(aucPos=auc.pos,Max=max,aucNeg=auc.neg,Min=min))
# 
# }

# interactive polynomial fit
# 
# plot data then aske for inital and end x points then
# fit an nth degree polynomial 
# (y = CnX^n + Cn-1X^n-1 + . . . + C1X + C0)
# to vectors x and y
# @param y : a numeric vector (x and y must have the same length)
# @param order : a numeric value for the polynomial degree. default is one.
# @return Coef : coeficients in ascending order (i.e. C0, C1, C2, ... , Cn )
# @return R2 : Rsquare goodness of fit
# @return model formula
# @export
# polyfit.interact<-function (y,order=1,...){
#   X11()
#   print(plot(y,type="l",...)) # plot data
#   abs<-locator(n=2,type="o") #aske for window points
#   
#   y<-y[abs$x[1]:abs$x[2]] #subset data
#   x<-seq_along(y) #genetate abscisse sequence
#   
#   fit<-lm(as.formula(paste('y~',paste('I(x^',1:order,')',sep='',collapse='+'))))
#   fit.r2<-summary(fit)$r.squared
#   fit.coef<-coef(fit)
#   
#   return(list(model=fit$call, Coef=fit.coef,R2=fit.r2))
# }

#' polynomial fit
#' 
#' Fit a n\eqn{^{th}}{th} degree polynom 
#' \deqn{(y = C_nX^n + C_{n-1}X^{n-1} + \ldots + C_1X + C_0)}{(y = CnX^n + Cn-1X^n-1 + \ldots + C1X + C0)}
#' to vectors x and y
#' @param x : a numeric vector
#' @param y : a numeric vector (x and y must have the same length)
#' @param order : a numeric value for the polynomial degree. default is one.
#' @return model formula
#' @return Coef: coeficients in ascending order (i.e. C0, C1, C2, \ldots , Cn )
#' @return R2: Rsquare
#' @export
polyfit<-function (x,y,order=1){
  
  
  fit<-lm(as.formula(paste('y~',paste('I(x^',1:order,')',sep='',collapse='+'))))
  fit.r2<-summary(fit)$r.squared
  fit.coef<-coef(fit)
  
  return(list(model=fit$call,Coef=fit.coef,R2=fit.r2))
}

#' polynomial evaluation
#' 
#' evaluation a n\eqn{^{th}}{th} degree polynom 
#' \deqn{(y = C_nX^n + C_{n-1}X^{n-1} + \ldots + C_1X + C_0)}{(y = CnX^n + Cn-1X^n-1 + \ldots + C1X + C0)}
#' at given values of x
#' 
#' @param x : a numeric vector
#' @param coef : a n dimention vector corresponding to the polynom coefficient (ascending order, i.e. C0, C1, C2, ... , Cn )
#' @return y
#' @export
polyval<-function (x,coef){
  
  order<-length(coef)
  y<-eval(parse(text=paste("coef[1] + x^",1:(order-1),
                           "*coef[",2:order,"]",sep="",collapse="+")))
  return(y)
}

#' apply a function FUN on a rooling windows of a vector
#' 
#' @param data : a numeric vector
#' @param width : the size of the rolling window
#' @param size : a logical value indicating if the returned vector have the same length as original data. default is TRUE.
#' @param FUN : the function to apply
#' @param ... : additional argument to pass to the FUN
#' @return a vector with FUN result 
#' @export
roll.funct<-function(data,width,FUN,size=T,...) {
  if (size){
    reslt<-NA
    for (i in 1:(trunc(length(data)))) {
      reslt[i]<-FUN(data[i:(i+width),...])
    }
  }      
  
  else {
    j<-1
    reslt<-NA
    for (i in 1:(trunc(length(data)/width))) {
      reslt[j]<-FUN(data[j:(j+width)],...)
      j<-j+width
    }
  }
  return(reslt)
}

#' Remove artifacts from biosensor signal, based on data's standar deviation
#' 
#' @param data : a numeric vector
#' @param z : a numeric value. Number of SD  over which values should be exculded
#' @param width : size of the window used to roll SD over data (see roll.funct)
#' @return a vector with NA remplacing exculded values
#' @export
noise.na<-function(data,z=20,width=30){
  
  reslt<-roll.funct(data,width,sd)
  data.rm<-reslt>median(reslt,na.rm=T)*z
  data2<-data
  data2[data.rm]<-NA
  return(data2)
}