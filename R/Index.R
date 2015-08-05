#' Oxidative index calculation
#' 
#' take CMRO2 and LDF to give an oxidative index
#' 
#' \deqn{ IO2 = CMRO2/(cbfbasal*LDF) } 
#' 
#' @param CMRO2 : vector/numeric value of Cerebral metabolic rate of oxygen in micromol/100g/min
#' @param cbfbasal : basale expected value of CBF from the litterature (default is 53 ml/100g/min wich was used to calculate L)
#' @param LDF :  a vector/numeric value of LDF (percentage from baseline, baseline is 100 \%)
#' @references Gjedde et al JCBFM (2000) 20(4), 747
#' @return IO2 : oxidative index, reflect the degree of flow metabolism coupling
#' @export
IO2.calc<-function (CMRO2,LDF,cbfbasal=53) {
  IO2<-CMRO2/(cbfbasal*LDF/100)
  return(IO2)
}

#' Oxygene Glucose index (OGI)
#' 
#' take CMRO2 and CMRGlucose to give an OGI
#' 
#' \deqn{ OGI = CMRO2/CMRGlucose } 
#' 
#' @param CMRO2 : vector (numeric value) of Cerebral metabolic rate of oxygen in micromol/100g/min
#' @param CMRGluc :  a vector (numeric value) of CMRGlucose in micromol/100g/min
#' @return OGI
#' @export
OGI.calc<-function (CMRO2,CMRGluc) {
  OGI<-CMRO2/CMRGluc
  return(OGI)
}
