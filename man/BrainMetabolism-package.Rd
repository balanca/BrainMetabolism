\name{BrainMetabolism-package}
\alias{BrainMetabolism-package}
\alias{BrainMetabolism}
\docType{package}
\title{
calculate brain metabolism rate from extracellular conentrations
}
\description{
calculate CMRO2, CMRgluc, mitochondrialPO2 ...\cr
provide functions to process data from biosensors
}
\details{
\tabular{ll}{
Package: \tab BrainMetabolism\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-03\cr
License: \tab MIT\cr
}

}
\author{
Baptiste Balanca \cr
Lyon Neuroscience Research Center\cr
Team TIGER

Maintainer: baptiste balanca <baptiste.balanca@gmail.com>
}
\references{
Piilgaard Lauritzen JCBFM (2009) 29, 1517\cr
Gjedde et al JCBFM (2005) 25(9), 1183\cr
Gjedde et al JCBFM (2000) 20(4), 747\cr
Du et al JCBFM 2012 32(9)\cr
Gruetter et al J Neurochem 1998 70(1)\cr
balanca et al (2016) in preparation\cr
}

\keyword{ package }

\examples{
MyTPO2<-23 #mmHg
MyLDF<-95 #percent
MyL<-L.calc(CMR=219, TPO2=MyTPO2, P50=36, h=2.7, Ca=8, cbf=53 )
MyCMRO2<-CMRO2.calc(LDF=MyLDF, MyTPO2, P50=36, h=2.7, Ca=8, L=MyL,cbfbase=53)

MyGlucBrain<-1 #mM
MyGlucPlasma<-6 #mM
MyCRMgluc<-CMRGluc.calc(MyGlucbrain, Vd=0.77, Kt=13.4, Tmax=1.35, Gplasma=MyGlucPlasma)
}
