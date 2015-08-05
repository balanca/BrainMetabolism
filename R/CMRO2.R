CMRO2<-function (CMR, # to be calcualted
                 CBF, # %age from baseline (My data)
                 TPO2, # mmHg (my data)
                 P50=36, # Hb half saturation - mmHg
                 h=2.7, # hill's coefficient
                 Ca=8, # arterial oxygene content - micormol/ml
                 L, # micromol/100g
                 cbfbase) # from lee et al 2001 JCBFM sprague dawley under alphachloralose - ml/100g/min
{
  cbf<-cbfbase*CBF/100 #ml/100g/min * %age/100
  
  P50 * ((2*Ca*cbf/CMR)-1)^(1/h) - (CMR/(2*L)) - TPO2   
}

# function for L calculation
CMRO2.L<-function (L,
                  CMR, TPO2, P50, h, Ca, cbf) # from lee et al 2001 JCBFM sprague dawley under alphachloralose - ml/100g/min
{  
  P50 * ((2*Ca*cbf/CMR)-1)^(1/h) - (CMR/(2*L)) - TPO2   
}

# graph CRMO2
CMRO2graph <- function (LDF,TPO2,CMR) {
  data<-data.frame(LDF=(LDF-median(LDF[2:5]))/median(LDF[2:5]),TPO2=(TPO2-median(TPO2[2:5]))/median(TPO2[2:5]),CMR=(CMR-median(CMR[2:5]))/median(CMR[2:5]))
  graph<- ggplot(data=data,aes(x=seq(1:length(data$TPO2))))+
    geom_line(aes(y=data$TPO2,colour="1"),size=1,alpha=0.6)+theme_bw()+
    geom_line(aes(y=data$LDF,colour="2"),size=1,alpha=0.6)+
    geom_line(aes(y=data$CMR,colour="3"),size=1.5,alpha=0.8)+
    scale_colour_manual(values=c("#CC0000","#009933","#0066CC"),guide_legend(title="Curve"),labels=c("TPO2","LDF","CMR"))+
    labs(title="delta from baseline",x="time",y="delta value :(x-x0)/x0")
  print(graph)
}

#' Calculate the effective diffusion coefficient of oxygen in brain tissue, L.
#' 
#' take CMRO2 and cerebral blood flow (CBF) to calculate the effective diffusion coefficient of oxygen in brain tissue (L)
#' 
#' \deqn{PbtO_{2} = P_{50} \cdotp \sqrt[\scriptstyle h]{\frac{2 \cdotp C_{a} \cdotp CBF}{CMRO_{2}}-1} - \frac{CMRO_{2}}{2 \cdotp L}}{
#' PbtO2 = P50 x (2 * Ca * CBF/CMRO2-1)^(1/h) - CMRO2/(2*L)} 
#' 
#' @param CMR : Cerebral Metabolic Rate of Oxygen, default is 219 micro mol/100
#' @param TPO2 : numeric value of brain oxygnen pressure (mmHg)
#' @param P50 : half-saturation tension of hemoglobine (default is 36 mmHg)
#' @param h : hill's coeficient
#' @param Ca : oxygen arterial concentration (default is 8 micromol/ml)
#' @param cbf : expected value of CBF (default is 53 ml/100g/min wich was used to calculate L
#' @references Gjedde et al JCBFM (2005) 25(9), 1183
#' @references Piilgaard et al JCBFM (2009) 29, 1517
#' @return L : numeric value of the effective diffusion coefficient of oxygen in brain tissue micomol/100g/mmHg
#' @import ggplot2
#' @export
L.calc<- function (CMR=219, TPO2, P50=36, h=2.7, Ca=8, cbf=53 ) 
{ 
  if (tail(curve(CMRO2.L(x, CMR, TPO2, P50, h, Ca, cbf),from=0,to=100,n=101)$y,n=1)=="NaN"){
    maxindex<-which(curve(CMRO2.L(x, CMR, TPO2, P50, h, Ca, cbf),from=0,to=100,n=101)$y=="NaN")[1]-2
  }
  else {maxindex<-100}
  
  if (tail(curve(CMRO2.L(x, CMR, TPO2, P50, h, Ca, cbf),from=0,to=100,n=101)$y,n=1)<0) {
    stop("L > 100 micomol/100g/mmHg")
  }
  
  
  L <- uniroot(CMRO2.L, c(0, maxindex), CMR, TPO2, P50, h, Ca, cbf)$root
 
  return (L)
}

#' CMRO2 calculation
#' 
#' take tissue oxygen pressure (tPO2) and cerebral blood flow (relative value, e.g. laser doopler lowmetry BPU)
#' 
#' \deqn{ PbtO_{2} = P_{50} \cdotp \sqrt[\scriptstyle h]{\frac{2 \cdotp C_{a} \cdotp CBF}{CMRO_{2}}-1} - \frac{CMRO_{2}}{2 \cdotp L} }{
#' PbtO2 = P50 x (2 * Ca * CBF/CMRO2-1)^(1/h) - CMRO2/(2*L)}
#' 
#' @param TPO2 : a vector/numeric value of brain oxygnen pressure (mmHg)
#' @param LDF : a vector/numeric value of LDF (percentage of baseline, i.e: basal value is 100 \% )
#' @param P50 : half-saturation tension of hemoglobine (default is 36 mmHg)
#' @param h : hill's coeficient
#' @param Ca : oxygen arterial concentration (default is 8 micromol/ml)
#' @param L : effective diffusion coefficient of oxygen in brain tissue, default is 4.03 micomol/100g/mmHg, but one should use the L.calc function to calculate it from their data.
#' @param cbfbase : basal expected value of CBF (default is 53 ml/100g/min wich was used to calculate L) \cr
#' LBF and TPO2 must be the same length
#' @references Gjedde et al JCBFM (2005) 25(9), 1183
#' @references Piilgaard et al JCBFM (2009) 29, 1517
#' @return CMRO2 : vector/numeric value of Cerebral Metabolic Rate of Oxygen in micromol/100g/min
#' @import ggplot2
#' @export
CMRO2.calc<- function (LDF,TPO2,
                       P50=36, # mmHg
                       h=2.7, 
                       Ca=8, # micormol/ml
                       L=4.03, # micromol/100g
                       cbfbase=53 #ml/100g/min
                       ) 
                       { 
  
  if (length(LDF) == length(TPO2)) {
    CMR<-NA
    X11()
    for (i in 1:length(LDF)) {
      if (tail(curve(CMRO2(x,LDF[i],TPO2[i],P50,h,Ca,L,cbfbase),from=0,to=1000,n=1001)[[2]],n=1)=="NaN"){
        maxindex<-which(curve(CMRO2(x,LDF[i],TPO2[i],P50,h,Ca,L,cbfbase),from=0,to=1000,n=1001)[[2]]=="NaN")[1]-2
      }
      
      else {maxindex<-1000}
      CMR[i] <- uniroot(CMRO2, c(0, maxindex),LDF[i],TPO2[i],P50,h,Ca,L,cbfbase)$root
    }
    dev.off()
  }
  
  else {
    stop("LDF and TPO2 must have the same length")
  }

  return (CMR)
}
