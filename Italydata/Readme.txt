
BAPE： To store the parameter estimations in each model.
ind.BAPE： To store the selected indicators in each model.
EM.BIC: To calculate parameters in each model using EM algorithm and calculate the BIC value based on selected variables in each model at last iteration..
Gibbs.sampler4: To select valuable predictors involving non-ignorable missing values using Gibbs sampler. (Note: The selection model criterion is observed BIC.)
Gibbs.sampler:  To select valuable predictors involving non-ignorable missing values using Gibbs sampler. (Note: The selection model criterion is BIC_Q..)

Before running codes, 
1. R package "MASS" should be installed.
2. Set workplace in the document that includes simulation data (step1_datageneration.RData) and the code file (GibsamplerbicselQ_SiruW.R).
3. The order of the first self-defined parameter (predata) in R function Gibbs.sampler and Gibbs.sampler4 is required as follows:
predictors without missing values Z=(z_{1},z_{2},...,z_{q}); predictors with missing values: X=(x_{1},x_{2},...,x_{p}), the response with missing y
predata should be given as [Z,X,y]=[z_{1},z_{2},...,z_{q},x_{1},x_{2},...,x_{p},y].
4. The order of the second self-defined parameter (predata) in R function Gibbs.sampler and Gibbs.sampler4 is required as follows:
variables of indicators for X: R_X=(r_{1},r_{2},...,r_{p}), indicator variable for y: r_y
rsind should be given as [R_X,r_y]=r_{1},r_{2},...,r_{p}, r_y].

1. BIC_Q as a model criterion.
Open a new R script and running codes are shown as follows:
load("step1_datageneration.RData")
source("GibsamplerbicselQ_SiruW_Italy.R")
rs.indicator=cbind(r.sex,r.fecund,r.contra)
Italy.dummy=as.matrix(Italy.dummy)
indicator.italy=list(rep(1,NCOL(Italy.dummy)),rep(1,NCOL(Italy.dummy)-2),rep(1,NCOL(Italy.dummy)-1),rep(1,NCOL(Italy.dummy)+1),rep(1,NCOL(Italy.dummy)+2),rep(1,NCOL(Italy.dummy)+3))
indicator.italy[[3]][c(13:14)+1]=0
indicator.italy[[4]][c(15,16)+1]=0
indicator.italy[[5]][c(18)+1]=0
indicator.italy[[6]][c(12,14,17)+1]=0

library(MASS)#introduce function ginv()
tab.out=Gibbs.sampler(Italy.dummy,rs.indicator,indicator.italy,90,Weigh)


2. BIC_obs as a model criterion.
Open a new R script and running codes are shown as follows:
load("step1_datageneration.RData")
source("Gibsamplerbicselobs_SiruW_Italy.R")
rs.indicator=cbind(r.sex,r.fecund,r.contra)
Italy.dummy=as.matrix(Italy.dummy)
indicator.italy=list(rep(1,NCOL(Italy.dummy)),rep(1,NCOL(Italy.dummy)-2),rep(1,NCOL(Italy.dummy)-1),rep(1,NCOL(Italy.dummy)+1),rep(1,NCOL(Italy.dummy)+2),rep(1,NCOL(Italy.dummy)+3))
indicator.italy[[3]][c(13:14)+1]=0
indicator.italy[[4]][c(15,16)+1]=0
indicator.italy[[5]][c(18)+1]=0
indicator.italy[[6]][c(12,14,17)+1]=0

library(MASS)#introduce function ginv()
tab.out=Gibbs.sampler4(Italy.dummy,rs.indicator,indicator.italy,90,Weigh)
