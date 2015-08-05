
#' Sensor calibration
#' 
#' Fit an n\eqn{^{th}}{th} degree polynom
#' \deqn{(y = C_nX^n + C_{n-1}X^{n-1} + \ldots + C_1X + C_0)}{
#' (y = CnX^n + Cn-1X^n-1 + \ldots + C1X + C0)}
#' to vectors x=volt and y=concentration
#' 
#' @param volt : a numeric vector of the biosensor voltage
#' @param mol : a numeric vector of the molecule conentration in the medium (x and y must have the same length)
#' @param order : a numeric value for the polynomial degree. default is one.
#' @return Coef : coeficients in ascending order (i.e. C0, C1, C2, \ldots , Cn )
#' @return R2 : goodness of fit
#' @import ggplot2
#' @export
calibration<-function(volt,mol,order=1) {
  
  if (length(volt)!=length(mol)) stop ("volt and mol must be the same length")
  dev.new()
  fit<-polyfit(volt,mol,order)
  
  graph1 <- ggplot (data=data.frame(volt=volt,mol=mol),aes(x=volt,y=mol)) +
    geom_point(shape=5,size=5) +
    theme_bw(base_size = 20) +
    stat_smooth(method="lm", formula = eval(fit$model[[2]]) ,se=F) +
    ggtitle(paste("Calibration Curve",
                            paste("order= ",order,"R2= ",fit$R2),sep="\n")) +
    ylab("concentration") +
    xlab("volt")
    
  print(graph1)
  
  return(list(Coef=fit$Coef,R2=fit$R2))
  
}


#' Oxygen Tension correction
#' 
#' Biosensor enzymatic reaction, that underpin amperometric measures, has an asymptotic relation to oxygen tention in the medium :
#' \deqn{m(PO_2,P) = P_1 + (P_2 - P_1)\times \exp(-\exp(P_3)PO_2)}{
#' m(PO2,P) = P1 + (P2 - P1)*exp(-exp(P3)PO2)}
#' parameters are different for each enzyme and has been measured in vitro\cr
#' This function correct biosensor signal depending on PO2 conditions during calibration and experimentation
#' 
#' @param x : a numeric vector of the biosensor voltage
#' @param enz : enzyme on the biosensor i.e. "glucose", "lactate", "glutamate, "daao"
#' @param TPO2 : oxygene tension in the medium during the experiment. the default is 30mmHg measured in anesthetized rat brain.
#' @return  volt.O2.cor : a vector of corected x values for TPO2
#' @references balanca et al 2015
#' @export
correction.TPO2<-function(x,enz,TPO2=28) {
  # Deltapercentvolt ~ SSasymp(PO2,a,b,c) => m(PO2,P) = a + (b - a)*exp(-exp(c)PO2)
    if (enz=="glutamate") {
      Asymp<- -0 ; R0<- -74.2 ; lrc<- -2.9 
    }
    if (enz=="glucose") {
      Asymp<-0 ;R0<- -43.9 ;lrc<- -2.6 
    }
    if (enz=="lactate") {
      Asymp<-0 ;R0<--21.8 ;lrc<- -3.17
    }
    if (enz=="daao") {
      Asymp<-0 ;R0<-1 ;lrc<-1 
    }

    delta<- Asymp + (R0-Asymp) * exp(-exp(lrc)*TPO2)
    return(
      volt.O2.cor<-x/(1+(delta/100))
      )    
}

#' Temperature correction
#' 
#' Biosensor enzymatic reaction, that underpin amperometric measures, has a sigmoid relation to temperature :
#' \deqn{m(x,P)=\frac{P_1+(P_2-P_1)}{(1+\exp((P_3-x)/P_4))}}{
#' m(x,P)=P1+(P2-P1)/(1+exp((P3-x)/P4))}
#' parameters are different for each enzyme and has been measured in vitro.\cr
#' This function correct biosensor signal depending on temperature conditions during calibration and experimentation
#' 
#' @param x : a numeric vector of the biosensor voltage
#' @param enz : enzyme on the biosensor i.e. "glucose", "lactate", "glutamate, "daao"
#' @param temp.calib : temperature of the medium where sensor has been calibrated 
#' the default is 25\eqn{^o}{o}C (i.e. room temperature).
#' @param temp.exp : temperature of the medium during experiment
#' the default is 37\eqn{^o}{o}C (i.e. animal central temperature)
#' @return  volt.temp.cor : a vector of corected x values for temperature
#' @references balanca et al 2015
#' @export
correction.Temp<-function(x,enz,temp.calib=25,temp.exp=37) {
  # SSfpl (x, A, B, xmid, scal) m(x,P)=A+(B-A)/(1+exp((xmid-x)/scal))
  if (enz=="glutamate") {
    A<- -104.2 ; B<- 52.6 ; xmid<- 29.8 ; scal<- 13
    delta.exp<-A+(B-A)/(1+exp((xmid-temp.exp)/scal))
    delta.calib<-A+(B-A)/(1+exp((xmid-temp.calib)/scal))
  }
  if (enz=="glucose") {
    A<- -92.2 ; B<- 57.5 ; xmid<- 34.9 ; scal<- 10.5
    delta.exp<-A+(B-A)/(1+exp((xmid-temp.exp)/scal))
    delta.calib<-A+(B-A)/(1+exp((xmid-temp.calib)/scal))
  }
  if (enz=="lactate") {
    inter<- -76 ; slope<- 1.9
    delta.exp <- slope*temp.exp + inter
    delta.calib <- slope*temp.calib + inter 
  }
  if (enz=="daao") {
    A<-1 ; B<-1 ; xmid<-1 ; scal<-1
    delta.exp<-A+(B-A)/(1+exp((xmid-temp.exp)/scal))
    delta.calib<-A+(B-A)/(1+exp((xmid-temp.calib)/scal))
  }

  return(
    volt.temp.cor<-x*(1+(delta.calib/100))/(1+(delta.exp/100))
    )  
  
}



