
BAPE： To store the parameter estimations in each model.
ind.BAPE： To store the selected indicators in each model.
EM.BIC: To calculate parameters in each model using EM algorithm and calculate the BIC value based on selected variables in each model at last iteration..
Gibbs.sampler4: To select valuable predictors involving non-ignorable missing values using Gibbs sampler. (Note: The selection model criterion is observed BIC.)
Gibbs.sampler:  To select valuable predictors involving non-ignorable missing values using Gibbs sampler. (Note: The selection model criterion is BIC_Q..)

Before running codes, 
1. R package "MASS" should be installed.
2. Set workplace in the document that includes simulation data (step1_misgibsimdata2.RData) and the code file (GibsamplerbicselQ_SiruW.R).
3. The order of the first self-defined parameter (predata) in R function Gibbs.sampler and Gibbs.sampler4 is required as follows:
predictors without missing values Z=(z_{1},z_{2},...,z_{q}); predictors with missing values: X=(x_{1},x_{2},...,x_{p}), the response with missing y
predata should be given as [Z,X,y]=[z_{1},z_{2},...,z_{q},x_{1},x_{2},...,x_{p},y].
4. The order of the second self-defined parameter (predata) in R function Gibbs.sampler and Gibbs.sampler4 is required as follows:
variables of indicators for X: R_X=(r_{1},r_{2},...,r_{p}), indicator variable for y: r_y
rsind should be given as [R_X,r_y]=r_{1},r_{2},...,r_{p}, r_y].

1. BIC_Q as a model criterion.
Open a new R script and running codes are shown as follows:
load("step1_misgibsimdata2.RData")
source("GibsamplerbicselQ_SiruW.R")

n.complete=NROW(na.omit(t(cbind(zz,xx))))#18 covariates without missing values
ind.zx=list(rep(1,NCOL(cbind(zz,xx))+1),rep(1,n.complete+1),rep(1,n.complete+2),rep(1,n.complete+4),rep(1,n.complete+5),rep(1,n.complete+6))

tab.out=Gibbs.sampler(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),cbind(rr.x,ss.y),ind.zx,35)

2. BIC_obs as a model criterion.
Open a new R script and running codes are shown as follows:
load("step1_misgibsimdata2.RData")
source("Gibsamplerbicselobs_SiruW.R")

n.complete=NROW(na.omit(t(cbind(zz,xx))))#18 covariates without missing values
ind.zx=list(rep(1,NCOL(cbind(zz,xx))+1),rep(1,n.complete+1),rep(1,n.complete+2),rep(1,n.complete+4),rep(1,n.complete+5),rep(1,n.complete+6))

tab.out=Gibbs.sampler4(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),cbind(rr.x,ss.y),ind.zx,35)
