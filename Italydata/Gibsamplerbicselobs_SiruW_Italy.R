################################################
#Select important variables in GLM involving nonignorable missing values using BIC_obs
#Codes only used for Italy data
#Siru Wang (siruw@student.unimelb.edu.au)
#################################################
BAPE=function(list.bape){
  #FUNCTION: To store the parameter estimations in each model 
  #Argment: 
  #list.bape: A list which contains a set of vectors. The length of each vector is the number of parameters (intercept included) in the associated model.
  #The order of each vector is fixed.
  #The first vector reprensents the model for the response y regressed on predictor matrix (Z,X).
  #The second set of vectors represents a set of models for predictors with missing values X.
  #The third set of vectors represents a set of models for missing indicators of X.
  #The last vector represents the model for the model for missing indicator of y.
  #Example:
  #There are 2 variables without missing values denoted by Z and 2 variables with missing values denoted by X. The response y also contains missing values.
  #list.bape=list(l1,l2,l3,l4,l5,l6).
  #l1=c(1,1,1,1,1); l2=c(1,1,1);l3=c(1,1,1,1);l4=c(1,1,1,1,1,1);l5=c(1,1,1,1,1,1,1);l6=c(1,1,1,1,1,1,1,1)
  
  Beta=rep(0,sum(list.bape[[1]]))
  Eta=rep(0,sum(list.bape[[NROW(list.bape)]]))
  Alpha=lapply(list.bape[2:(1+0.5*(NROW(list.bape)-2))],function(x){y=rep(0,sum(x));return(y)})
  Phi=lapply(list.bape[(2+0.5*(NROW(list.bape)-2)):(NROW(list.bape)-1)],function(x){y=rep(0,sum(x));return(y)})
  return(list(Beta,Alpha,Phi,Eta))
}
ind.BAPE=function(list.bape){
  #FUNCTION: To store the selected indicators in each model.
  #Argment: 
  #list.bape: A list which contains a set of vectors. The length of each vector is the number of parameters (intercept included) in the associated model.
  ind.beta=list.bape[[1]]
  ind.eta=list.bape[[NROW(list.bape)]]
  ind.alpha=lapply(list.bape[2:(1+0.5*(NROW(list.bape)-2))],identity)
  ind.phi=lapply(list.bape[(2+0.5*(NROW(list.bape)-2)):(NROW(list.bape)-1)],identity)
  return(list(ind.beta,ind.alpha,ind.phi,ind.eta))
}



EM.BIC=function(predata,bape,ind.bape,rs.ind,wgh){
  #FUNCTION: To calculate parameters in each model using EM algorithm and calculate the BIC value based on selected variables in each model at last iteration.
  #Argments£º
  #predata: The combined matrix of predictors without missing values Z, predictors with mssing values X and the response y.
  #bape: A list of temporarily stored parameter estimations in each model at last iteration. 
  #indbape: A list of temporarily stored selected indicators in each model at last iteration. 
  #rs.ind: A matrix of missing indicators for predictors with missing values X and the response with missing values y.
  #wgh: Weights for individuals.
  
  
  #first and second derivative functions for Italydata
  #combination for missing data
  Bitmatrix<-function(n){
    set<-0:(2^n-1)
    rst<-matrix(0,ncol=n,nrow=2^n)
    for (i in 1:n){
      rst[,i]=ifelse((set-rowSums(rst*rep(c(2^((n-1):0)),each=2^n)))/(2^(n-i))>=1,1,0)
    }
    rst
  }
  
  
  
  ###################################################################################
  #conditional weighted probability
  con.w=function(predamis.k,bape,ind.bape,rs.k,wgh.k){#(z.intpt,x,y,bape,ind.bape,r,s,miscate)
    
    p.ymis=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]])))
    p.ymis=as.vector(p.ymis)
    denom.k=dbinom(round(predamis.k[,NCOL(predamis.k)]*wgh.k),round(wgh.k),p.ymis)
    
    predamisr.k=matrix(0.0,nrow=NROW(predamis.k),ncol=NCOL(predamis.k)+2,byrow=TRUE);
    predamisr.k[,1:NCOL(predamis.k)]=predamis.k;
    predamisr.k[,-(1:NCOL(predamis.k))]=matrix(rs.k[-NROW(rs.k)],nrow=NROW(predamis.k),ncol=NROW(rs.k)-1,byrow=TRUE)
    for (i.xr in 1:(NROW(rs.k)-1)){
      p.xrmis=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)])%*%as.matrix(bape[[2]][[i.xr]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)])%*%as.matrix(bape[[2]][[i.xr]])))
      p.xrmis=as.vector(p.xrmis)
      denom.k=denom.k*dbinom(round(predamis.k[,NROW(ind.bape[[2]][[i.xr]])+1]*wgh.k),round(wgh.k),p.xrmis)
      
      p.xrmis=exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)])%*%as.matrix(bape[[3]][[i.xr]]))/(1+exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)])%*%as.matrix(bape[[3]][[i.xr]])))
      p.xrmis=as.vector(p.xrmis)
      denom.k=denom.k*dbinom(round(rs.k[i.xr]*wgh.k),round(wgh.k),p.xrmis)
    }
    
    p.smis=exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
    p.smis=as.vector(p.smis)
    denom.k=denom.k*dbinom(round(rs.k[NROW(rs.k)]*wgh.k),round(wgh.k),p.smis)
    mistype.p=denom.k/sum(denom.k)
    
    as.vector(mistype.p)
  }
  
  ###############################################################################
  # create the function to generate covatiates list
  
  # the indicator of which sample
  comp.ind=complete.cases(predata)
  mis.type=sapply(c(1:3),Bitmatrix)
  
  
  Dy=function(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mis.type){
    predata.inpt=matrix(0.0,NROW(predata),NCOL(predata)+1);
    predata.inpt[,1]=rep(1,NROW(predata));predata.inpt[,-1]=predata
    #observed data
    p.y=exp(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]])/(1+exp(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]]))
    p.y=as.vector(p.y)
    
    D1y=colSums(wgh[(comp.ind==TRUE)]*(predata.inpt[comp.ind==TRUE,NROW(ind.bape[[1]])+1]-p.y)*predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])
    D2y.k=t(apply(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,ind.bape[[1]]==1],1,function(x){x%o%x}))*(-wgh[(comp.ind==TRUE)]*p.y*(1-p.y))
    D2y.k=colSums(D2y.k)
    D2y=matrix(D2y.k,sum(ind.bape[[1]]),sum(ind.bape[[1]]))
    
    #unobserved data
    for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
      mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      predamis.k=matrix(predata.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predata.inpt),byrow=TRUE)
     
      predamis.k[is.na(predamis.k)==T]=mistype.k
      con.P=con.w(predamis.k,bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,],wgh[comp.ind==FALSE][i.mis])
      
      p.y=exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]])/(1+exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]]))
      p.y=as.vector(p.y)
      D1y=D1y+colSums(wgh[comp.ind==FALSE][i.mis]*as.vector(predamis.k[,NCOL(predamis.k)]-p.y)*predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]*con.P)
      D2y.k=t(apply(predamis.k[,1:NROW(ind.bape[[1]])][,ind.bape[[1]]==1],1,function(x){x%o%x}))*wgh[comp.ind==FALSE][i.mis]*(-p.y)*(1-p.y)*con.P
      D2y.k=colSums(D2y.k)
      D2y=D2y+matrix(D2y.k,sum(ind.bape[[1]]),sum(ind.bape[[1]]))
    }
    
    list(D1y,D2y)
  }
  
  
  Dx=function(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mistype,i.xr){
    predata.inpt=matrix(0.0,NROW(predata),NCOL(predata)+1);
    predata.inpt[,1]=rep(1,NROW(predata));predata.inpt[,-1]=predata  
    #D1x=list(rep(0,NROW(bape[[2]][[1]])),rep(0,NROW(bape[[2]][[2]])));
    #D2x=list(matrix(0,NROW(bape[[2]][[1]]),NROW(bape[[2]][[1]])),matrix(0,NROW(bape[[2]][[2]]),NROW(bape[[2]][[2]])));
    
    #observed data
    p.x=exp(predata.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1]%*%bape[[2]][[i.xr]])/
      (1+exp(predata.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1]%*%bape[[2]][[i.xr]]))
    p.x=as.vector(p.x) 
    D1x=colSums(wgh[(comp.ind==TRUE)]*(predata.inpt[comp.ind==TRUE,NROW(ind.bape[[2]][[i.xr]])+1]-p.x)*predata.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])
    D2x.k=t(apply(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1],1,function(x){x%o%x}))*(-wgh[(comp.ind==TRUE)]*p.x*(1-p.x))
    D2x.k=colSums(D2x.k)
    D2x=matrix(D2x.k,sum(ind.bape[[2]][[i.xr]]),sum(ind.bape[[2]][[i.xr]])) 
    
    #unobserved data
    for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
      mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      predamis.k=matrix(predata.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predata.inpt),byrow=TRUE)
      predamis.k[is.na(predamis.k)==T]=mistype.k
      con.P=con.w(predamis.k,bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,],wgh[comp.ind==FALSE][i.mis])
      
      p.x=exp(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]%*%bape[[2]][[i.xr]])/(1+exp(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]%*%bape[[2]][[i.xr]]))
      p.x=as.vector(p.x)
      D1x=D1x+colSums(wgh[comp.ind==FALSE][i.mis]*as.vector(predamis.k[,NROW(ind.bape[[2]][[i.xr]])+1]-p.x)*predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]*con.P)
      D2x.k=t(apply(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1],1,function(x){x%o%x}))*(-p.x)*(1-p.x)*con.P*wgh[comp.ind==FALSE][i.mis]
      D2x.k=colSums(D2x.k)
      D2x=D2x+matrix(D2x.k,sum(ind.bape[[2]][[i.xr]]),sum(ind.bape[[2]][[i.xr]]))
    }
    
    list(D1x,D2x)
  }
  
  Dr=function(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mistype,i.xr){
    predatar.inpt=matrix(0.0,NROW(predata),NCOL(predata)+NCOL(rs.ind))
    predatar.inpt[,1]=rep(1,NROW(predata));
    predatar.inpt[,2:(NCOL(predata)+1)]=predata;predatar.inpt[,(NCOL(predatar.inpt)-1):NCOL(predatar.inpt)]=rs.ind[,-NCOL(rs.ind)]
    
    #observed data
    p.r=exp(predatar.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1]%*%bape[[3]][[i.xr]])/
      (1+exp(predatar.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1]%*%bape[[3]][[i.xr]]))
    p.r=as.vector(p.r) 
    D1r=colSums(wgh[comp.ind==TRUE]*(predatar.inpt[comp.ind==TRUE,NROW(ind.bape[[3]][[i.xr]])+1]-p.r)*predatar.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])
    D2r.k=t(apply(predatar.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1],1,function(x){x%o%x}))*(-p.r*(1-p.r)*wgh[comp.ind==TRUE])
    D2r.k=colSums(D2r.k)
    D2r=matrix(D2r.k,sum(ind.bape[[3]][[i.xr]]),sum(ind.bape[[3]][[i.xr]])) 
    
    #unobserved data
    for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
      mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      predamis.k=matrix(predatar.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatar.inpt),byrow=TRUE)
      predamis.k[is.na(predamis.k)==T]=mistype.k
      con.P=con.w(predamis.k[,1:(NCOL(predata)+1)],bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,],wgh[comp.ind==FALSE][i.mis])
      
      p.r=exp(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]%*%bape[[3]][[i.xr]])/(1+exp(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]%*%bape[[3]][[i.xr]]))
      p.r=as.vector(p.r)
      D1r=D1r+colSums(wgh[comp.ind==FALSE][i.mis]*as.vector(predamis.k[,NROW(ind.bape[[3]][[i.xr]])+1]-p.r)*predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]*con.P)
      D2r.k=t(apply(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1],1,function(x){x%o%x}))*(-p.r)*(1-p.r)*con.P*wgh[comp.ind==FALSE][i.mis]
      D2r.k=colSums(D2r.k)
      D2r=D2r+matrix(D2r.k,sum(ind.bape[[3]][[i.xr]]),sum(ind.bape[[3]][[i.xr]]))
    }
    
    list(D1r,D2r)
    
  }
  
  
  Ds=function(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mistype){
    predatars.inpt=matrix(0.0,NROW(predata),NCOL(predata)+NCOL(rs.ind)+1)
    predatars.inpt[,1]=rep(1,NROW(predata));
    predatars.inpt[,2:(NCOL(predata)+1)]=predata;predatars.inpt[,(NCOL(predatars.inpt)-2):(NCOL(predatars.inpt))]=rs.ind
    
    #observed data
    p.s=exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]])/(1+exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]]))
    p.s=as.vector(p.s)
    D1s=colSums(wgh[comp.ind==TRUE]*(predatars.inpt[comp.ind==TRUE,NROW(ind.bape[[4]])+1]-p.s)*predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])
    D2s.k=t(apply(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1],1,function(x){x%o%x}))*(-p.s*(1-p.s)*wgh[comp.ind==TRUE])
    D2s.k=colSums(D2s.k)
    D2s=matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    
    #unobserved data
    for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
      mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      predamis.k=matrix(predatars.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatars.inpt),byrow=TRUE)
      predamis.k[is.na(predamis.k)==T]=mistype.k
      con.P=con.w(predamis.k[,1:(NCOL(predata)+1)],bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,],wgh[comp.ind==FALSE][i.mis])
      
      p.s=exp(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]])/(1+exp(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]]))
      p.s=as.vector(p.s)
      D1s=D1s+colSums(wgh[comp.ind==FALSE][i.mis]*as.vector(predamis.k[,NCOL(predamis.k)]-p.s)*predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]*con.P)
      D2s.k=t(apply(predamis.k[,1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1],1,function(x){x%o%x}))*(-p.s)*(1-p.s)*con.P*wgh[comp.ind==FALSE][i.mis]
      D2s.k=colSums(D2s.k)
      D2s=D2s+matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    }
    list(D1s,D2s)
  }

  
  i.EM=1
  bape0=bape
  while((i.EM<=24)||((i.EM<=26)&&max(abs(unlist(bape)-unlist(bape0)))>=0.1)){
    bape0=bape
    cat("i.EM",i.EM,"\n")
    if(sum(ind.bape[[1]])==1){
      bape[[1]]=glm(predata[,NCOL(predata)]~1,family=binomial(link=logit))$coefficients
    }else{
      D12y=Dy(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mis.type)
      bape[[1]]=bape[[1]]-ginv(D12y[[2]])%*%D12y[[1]]
    }
    
    for (i.mis in 1:(NCOL(rs.ind)-1)){
      if(sum(ind.bape[[2]][[i.mis]])==1){
        bape[[2]][[i.mis]]=glm(predata[,complete.cases(t(predata))==FALSE][,i.mis]~1,family=binomial(link=logit))$coefficients
      }else{
        D12x=Dx(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mis.type,i.mis)
        bape[[2]][[i.mis]]=bape[[2]][[i.mis]]-ginv(D12x[[2]])%*%D12x[[1]]
      }
      
      if(sum(ind.bape[[3]][[i.mis]])==1){
        bape[[3]][[i.mis]]=glm(rs.ind[,i.mis]~1,family=binomial(link=logit))$coefficients
      }else{
        D12r=Dr(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mis.type,i.mis)
        bape[[3]][[i.mis]]=bape[[3]][[i.mis]]-ginv(D12r[[2]])%*%D12r[[1]]}
    }
    
    if(sum(ind.bape[[4]])==1){
      bape[[4]]=glm(rs.ind[,3]~1,family=binomial(link=logit))$coefficients
    }else{D12s=Ds(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mis.type)
    bape[[4]]=bape[[4]]-ginv(D12s[[2]])%*%D12s[[1]]}
    
    i.EM=i.EM+1

  }
  

  
  # definition of BIC
  BIC.def=function(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mis.type){
    p=NROW(unlist(bape))
    #log information of BIC
    log.BIC=function(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mis.type){
      predatars.inpt=matrix(0.0,NROW(predata),NCOL(predata)+NCOL(rs.ind)+1)
      predatars.inpt[,1]=rep(1,NROW(predata));
      predatars.inpt[,2:(NCOL(predata)+1)]=predata;predatars.inpt[,(NCOL(predatars.inpt)-2):(NCOL(predatars.inpt))]=rs.ind
      
      #observed data
      p.y=exp(as.matrix(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]]))/(1+exp(as.matrix(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]])))
      p.y=as.vector(p.y)
      L.obs=dbinom(round(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[1]])+1]*wgh[comp.ind==TRUE]),round(wgh[comp.ind==TRUE]),p.y)
      for (i.xr in 1:2){
        p.x=exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])%*%as.matrix(bape[[2]][[i.xr]]))/
          (1+exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])%*%as.matrix(bape[[2]][[i.xr]])))
        p.x=as.vector(p.x)
        L.obs=L.obs*dbinom(round(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[2]][[i.xr]])+1]*wgh[comp.ind==TRUE]),round(wgh[comp.ind==TRUE]),p.x)
        p.r=exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])%*%as.matrix(bape[[3]][[i.xr]]))/
          (1+exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])%*%as.matrix(bape[[3]][[i.xr]])))
        p.r=as.vector(p.r)
        L.obs=L.obs*dbinom(round(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[3]][[i.xr]])+1]*wgh[comp.ind==TRUE]),round(wgh[comp.ind==TRUE]),p.r)
      }
      p.s=exp(as.matrix(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
      p.s=as.vector(p.s)
      L.obs=L.obs*dbinom(round(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[4]])+1]*wgh[comp.ind==TRUE]),round(wgh[comp.ind==TRUE]),p.s)
      l.obs=sum(log(L.obs))
      
      #unobserved data
      #  for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
      #    mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      #    predamis.k=matrix(predatars.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatars.inpt),byrow=TRUE)
      #    predamis.k[is.na(predamis.k)==T]=mistype.k
      
      #    p.ymis=exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]])/(1+exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]]))
      #    p.ymis=as.vector(p.ymis)
      #    denom.k=dbinom(predamis.k[,NROW(ind.bape[[1]])+1],1,p.ymis)
      
      #    for (i.xr in 1:(NCOL(rs.ind)-1)){
      #      p.xrmis=exp(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]%*%bape[[2]][[i.xr]])/(1+exp(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]%*%bape[[2]][[i.xr]]))
      #      p.xrmis=as.vector(p.xrmis)
      #      denom.k=denom.k*dbinom(predamis.k[,NROW(ind.bape[[2]][[i.xr]])+1],1,p.xrmis)
      
      #     p.xrmis=exp(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]%*%bape[[3]][[i.xr]])/(1+exp(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]%*%bape[[3]][[i.xr]]))
      #      p.xrmis=as.vector(p.xrmis)
      #      denom.k=denom.k*dbinom(predamis.k[,NROW(ind.bape[[3]][[i.xr]])+1],1,p.xrmis)
      #    }
      
      #   p.smis=exp(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]])/(1+exp(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]]))
      #    p.smis=as.vector(p.smis)
      #    denom.k=denom.k*dbinom(predamis.k[,NROW(ind.bape[[4]])+1],1,p.smis)
      
      #   l.obs=l.obs+log(sum(denom.k))
      #  }
      
      return(l.obs)
    }
    def.val=-2*log.BIC(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mis.type)+p*log(NROW(predata[(comp.ind==TRUE),]))
    
    return(def.val)
  }
  cat("i.EM",i.EM,"\n")
  c(BIC.def(predata,bape,ind.bape,rs.ind,wgh,comp.ind,mis.type),i.EM)
}



Gibbs.sampler<-function(predata,rs.ind,ind,iter.Gibs,wgh){
  #FUNCTION: To select valuable predictors involving non-ignorable missing values using Gibbs sampler.
  #Argments:
  #predata: The combined matrix of predictors without missing values Z, predictors with mssing values X and the response y.
  #rs.ind: A matrix of missing indicators for predictors with missing values X and the response with missing values y.
  #ind: A list which contains a set of vectors. The length of each vector is the number of parameters (intercept included) in the associated model.
  #The order of each vector is fixed.
  #The first vector reprensents the model for the response y regressed on predictor matrix (Z,X).
  #The second set of vectors represents a set of models for predictors with missing values X.
  #The third set of vectors represents a set of models for missing indicators of X.
  #The last vector represents the model for the model for missing indicator of y.
  #Example:
  #There are 2 variables without missing values denoted by Z and 2 variables with missing values denoted by X. The response y also contains missing values.
  #list.bape=list(l1,l2,l3,l4,l5,l6).
  #l1=c(1,1,1,1,1); l2=c(1,1,1);l3=c(1,1,1,1);l4=c(1,1,1,1,1,1);l5=c(1,1,1,1,1,1,1);l6=c(1,1,1,1,1,1,1,1).
  #iter.Gibs: The number of iterations.
  #wgh: Weights of individuals. 
  
  FCP<-function(predata,rs.ind,ind,rr,cc,wgh){#n from the last to the second
    
    ind[[rr]][cc]<-1
    fun1.bic=EM.BIC(predata,BAPE(ind),ind.BAPE(ind),rs.ind,wgh)
    cat("fun1",fun1.bic,"\n")
    if(fun1.bic[2]==27){fun1.bic[1]=10e+5}else{fun1.bic[1]=fun1.bic[1]}
    ###################################################################################
    
    ind[[rr]][cc]<-0
    fun2.bic=EM.BIC(predata,BAPE(ind),ind.BAPE(ind),rs.ind,wgh)
    cat("fun2",fun2.bic,"\n")
    if(fun2.bic[2]==27){
      fun2.bic[1]=10e+5
      res<-1/{1+exp(-fun2.bic[1]+fun1.bic[1])}
    }else{res<-1/{1+exp(0.8*(-fun2.bic[1]+fun1.bic[1]))}}
    
    return(c(res,min(fun1.bic[1],fun2.bic[1])))
  }
  
  #tab.matrix=list(matrix(0,NROW(unlist(ind))*iter.Gibs,NROW(ind[[1]])),matrix(0,NROW(unlist(ind))*iter.Gibs,NROW(ind[[2]])),
  #               matrix(0,NROW(unlist(ind))*iter.Gibs,NROW(ind[[3]])),matrix(0,NROW(unlist(ind))*iter.Gibs,NROW(ind[[4]])),
  #                matrix(0,NROW(unlist(ind))*iter.Gibs,NROW(ind[[5]])),matrix(0,NROW(unlist(ind))*iter.Gibs,NROW(ind[[6]])))
  tab.matrix=lapply(ind,function(x,y){matrix(0,(NROW(unlist(y[[1]]))-NROW(y[[1]]))*y[[2]],NROW(x)+1)},list(ind,iter.Gibs))
  
  del.stru0=list(1:16+1,1:14+1,c(1:12,15)+1,c(1:14,17)+1,c(1:17)+1,c(1:11,13,15,16,18,19)+1)
  i.tabmatrix=1
  for(i.Gib in 1:iter.Gibs){#cat("Gib",i.Gib,"\n")
    for(rr in 1:NROW(ind)){
      for (cc in del.stru0[[rr]]){cat("Gib.iteration",c(i.Gib,rr,cc,i.tabmatrix),"\n")
        print(ind)
        fcp=FCP(predata,rs.ind,ind,rr,cc,wgh)
        ind[[rr]][cc]<-rbinom(1,1,fcp[1])
        for (i.tab in 1:NROW(ind)){tab.matrix[[i.tab]][i.tabmatrix,]=c(fcp[2],ind[[i.tab]])}
        i.tabmatrix=i.tabmatrix+1
      } 
    }
  }
  list(tab.matrix,ind)
}