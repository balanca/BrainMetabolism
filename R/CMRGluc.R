CMRGluc <- function (CMRG, # to be calculated
                     Gbrain, # micromol/mL (from my data)
                     Vd=0.77, # glucose brai space difusion - ml/g
                     Kt=1.4, # apparent maximal transport rate
                     Tmax=1.27, 
                     Gplasma=7.3) #mmol/L
{

  (Vd*( ( ( (Tmax/CMRG)-1)*Gplasma - Kt ) / ((Tmax/CMRG) + 1) ) ) - Gbrain
}


#' CMRGluc calculation
#' 
#' take extracellular glucose concentration and return brain metabolic rate of glucose, using a reversible Michaelis-Menten model equation:
#' 
#' \deqn{ G_{brain}=V_d \frac{\left(\frac{T_{max}}{CMR_{gluc}}-1\right)\times G_{plasma} - K_t}
#' {\frac{T_{max}}{CMR_{gluc}} +1} }{
#' Gbrain=Vd x ( ((Tmax)/(CMRgluc)-1) * Gplasma - Kt) / ((Tmax)/(CMRgluc) +1) }
#' 
#' @param Gbrain : a vector/numeric value of extracellular glucose concentration in mmol/L, time decay
#' @param Vd : glucose brain space diffusion (default is 0.77 ml/g)
#' @param Kt : glucose apparent maximal transport rate (default is 1.4 mmol/L)
#' @param Tmax : the apparent maximal transport rate (default is 1.27 micormol/g/min)
#' @param Gplasma : plasma glucose concentration (default=6mmol/L)
#' @references Morgenthaler et al. Neurochem Int 2006, 48
#' @references Ori et al. Anesthesiology 1986, 65(2)
#' @references Duarte et al. Front in Neuroenrgetics 2009, 1
#' @references Gruetter et al J Neurochem 1998, 70(1)
#' @return CMRGlucose : numeric value of Cerebral metabolic rate of glucose in micromol/g/min
#' @export

CMRGluc.calc<- function (Gbrain, # mmol/L (from my data)
                        Vd=0.77, # glucose brai space difusion - ml/g
                        Kt=1.4, # apparent maximal transport rate
                        Tmax=1.27,
                        Gplasma=7.3) #mmol/L
  { 
 
    CMR<-NA
    X11()
    for (i in 1:length(Gbrain)) {
      if (tail(curve(CMRGluc(x,Gbrain[i],Vd,Kt,Tmax,Gplasma),from=0.0001,to=10,n=10000)[[2]],n=1)=="NaN"){
        maxindex<-which(curve(CMRGluc(x,Gbrain[i],Vd,Kt,Tmax,Gplasma),from=0.0001,to=10,n=10000)[[2]]=="NaN")[1]-2
      }
      else {maxindex<-10000}
      
      CMR[i] <- uniroot(CMRGluc, c(0.0001,maxindex/1000),Gbrain[i],Vd,Kt,Tmax,Gplasma)$root
    }
    dev.off()

  return (CMR)
}
