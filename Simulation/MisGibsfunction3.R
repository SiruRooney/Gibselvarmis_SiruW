BAPE=function(list.bape){
  Beta=rep(0,sum(list.bape[[1]]))
  Eta=rep(0,sum(list.bape[[NROW(list.bape)]]))
  Alpha=lapply(list.bape[2:(1+0.5*(NROW(list.bape)-2))],function(x){y=rep(0,sum(x));return(y)})
  Phi=lapply(list.bape[(2+0.5*(NROW(list.bape)-2)):(NROW(list.bape)-1)],function(x){y=rep(0,sum(x));return(y)})
  return(list(Beta,Alpha,Phi,Eta))
}
ind.BAPE=function(list.bape){
  ind.beta=list.bape[[1]]
  ind.eta=list.bape[[NROW(list.bape)]]
  ind.alpha=lapply(list.bape[2:(1+0.5*(NROW(list.bape)-2))],identity)
  ind.phi=lapply(list.bape[(2+0.5*(NROW(list.bape)-2)):(NROW(list.bape)-1)],identity)
  return(list(ind.beta,ind.alpha,ind.phi,ind.eta))
}

#############################################################

EM.BIC=function(predata,bape,ind.bape,rs.ind){
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
  con.w=function(predamis.k,bape,ind.bape,rs.k){#(z.intpt,x,y,bape,ind.bape,r,s,miscate)
   p.ymis=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]])))
   p.ymis=as.vector(p.ymis)
   denom.k=dbinom(predamis.k[,NCOL(predamis.k)],1,p.ymis)
   
   predamisr.k=matrix(0.0,nrow=NROW(predamis.k),ncol=NCOL(predamis.k)+2,byrow=TRUE);
   predamisr.k[,1:NCOL(predamis.k)]=predamis.k;predamisr.k[,-(1:NCOL(predamis.k))]=matrix(rs.k[-NROW(rs.k)],nrow=NROW(predamis.k),ncol=NROW(rs.k)-1,byrow=TRUE)
   for (i.xr in 1:(NROW(rs.k)-1)){
     p.xrmis=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)])%*%as.matrix(bape[[2]][[i.xr]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)])%*%as.matrix(bape[[2]][[i.xr]])))
     p.xrmis=as.vector(p.xrmis)
     denom.k=denom.k*dbinom(predamis.k[,NROW(ind.bape[[2]][[i.xr]])+1],1,p.xrmis)
     
     p.xrmis=exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)])%*%as.matrix(bape[[3]][[i.xr]]))/(1+exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)])%*%as.matrix(bape[[3]][[i.xr]])))
     p.xrmis=as.vector(p.xrmis)
     denom.k=denom.k*dbinom(rs.k[i.xr],1,p.xrmis)
   }
   print
   p.smis=exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
   p.smis=as.vector(p.smis)
   denom.k=denom.k*dbinom(rs.k[NROW(rs.k)],1,p.smis)
   mistype.p=denom.k/sum(denom.k)
   
   as.vector(mistype.p)
  }
  
  ###############################################################################
  # create the function to generate covatiates list
  
  # the indicator of which sample
  comp.ind=complete.cases(predata)
  mis.type=sapply(c(1:3),Bitmatrix)
  
  
  Dy=function(predata,bape,ind.bape,rs.ind,comp.ind,mis.type){
  predata.inpt=matrix(0.0,NROW(predata),NCOL(predata)+1);
  predata.inpt[,1]=rep(1,NROW(predata));predata.inpt[,-1]=predata
  #observed data
  p.y=exp(as.matrix(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]]))/(1+exp(as.matrix(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]])))
  p.y=as.vector(p.y)

  D1y=colSums((predata.inpt[comp.ind==TRUE,NROW(ind.bape[[1]])+1]-p.y)*predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])
  D2y.k=t(apply(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,ind.bape[[1]]==1],1,function(x){x%o%x}))*(-p.y*(1-p.y))
  D2y.k=colSums(D2y.k)
  D2y=matrix(D2y.k,sum(ind.bape[[1]]),sum(ind.bape[[1]]))
 #unobserved data
  for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
   mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
   predamis.k=matrix(predata.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predata.inpt),byrow=TRUE)
   predamis.k[is.na(predamis.k)==T]=mistype.k
   con.P=con.w(predamis.k,bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,])
   
   p.y=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]])))
   p.y=as.vector(p.y)
   D1y=D1y+colSums(as.vector(predamis.k[,NCOL(predamis.k)]-p.y)*predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]*con.P)
   D2y.k=t(apply(predamis.k[,1:NROW(ind.bape[[1]])][,ind.bape[[1]]==1],1,function(x){x%o%x}))*(-p.y)*(1-p.y)*con.P
   D2y.k=colSums(D2y.k)
   D2y=D2y+matrix(D2y.k,sum(ind.bape[[1]]),sum(ind.bape[[1]]))
  }
  list(D1y,D2y)
  }


  Dx=function(predata,bape,ind.bape,rs.ind,comp.ind,mistype,i.xr){
   predata.inpt=matrix(0.0,NROW(predata),NCOL(predata)+1);
   predata.inpt[,1]=rep(1,NROW(predata));predata.inpt[,-1]=predata  
   D1x=list(rep(0,NROW(bape[[2]][[1]])),rep(0,NROW(bape[[2]][[2]])));
   D2x=list(matrix(0,NROW(bape[[2]][[1]]),NROW(bape[[2]][[1]])),matrix(0,NROW(bape[[2]][[2]]),NROW(bape[[2]][[2]])));
 
  #observed data
  p.x=exp(as.matrix(predata.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])%*%as.matrix(bape[[2]][[i.xr]]))/
     (1+exp(as.matrix(predata.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])%*%as.matrix(bape[[2]][[i.xr]])))
  p.x=as.vector(p.x) 
  D1x[[i.xr]]=colSums((predata.inpt[comp.ind==TRUE,NROW(ind.bape[[2]][[i.xr]])+1]-p.x)*predata.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])
  D2x.k=t(apply(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1],1,function(x){x%o%x}))*(-p.x*(1-p.x))
  D2x.k=colSums(D2x.k)
  D2x[[i.xr]]=matrix(D2x.k,sum(ind.bape[[2]][[i.xr]]),sum(ind.bape[[2]][[i.xr]])) 
  
  #unobserved data
  for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
    mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
    predamis.k=matrix(predata.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predata.inpt),byrow=TRUE)
    predamis.k[is.na(predamis.k)==T]=mistype.k
    con.P=con.w(predamis.k,bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,])
    
    p.x=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)])%*%as.matrix(bape[[2]][[i.xr]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)])%*%as.matrix(bape[[2]][[i.xr]])))
    p.x=as.vector(p.x)
    D1x[[i.xr]]=D1x[[i.xr]]+colSums(as.vector(predamis.k[,NROW(ind.bape[[2]][[i.xr]])+1]-p.x)*predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]*con.P)
    D2x.k=t(apply(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1],1,function(x){x%o%x}))*(-p.x)*(1-p.x)*con.P
    D2x.k=colSums(D2x.k)
    D2x[[i.xr]]=D2x[[i.xr]]+matrix(D2x.k,sum(ind.bape[[2]][[i.xr]]),sum(ind.bape[[2]][[i.xr]]))
  }
  
  list(D1x,D2x)
  }
  
  Dr=function(predata,bape,ind.bape,rs.ind,comp.ind,mistype,i.xr){
    predatar.inpt=matrix(0.0,NROW(predata),NCOL(predata)+NCOL(rs.ind))
    predatar.inpt[,1]=rep(1,NROW(predata));
    predatar.inpt[,2:(NCOL(predata)+1)]=predata;predatar.inpt[,(NCOL(predatar.inpt)-1):(NCOL(predatar.inpt))]=rs.ind[,-NCOL(rs.ind)]
    D1r=list(rep(0,NROW(bape[[3]][[1]])),rep(0,NROW(bape[[3]][[2]])));
    D2r=list(matrix(0,NROW(bape[[3]][[1]]),NROW(bape[[3]][[1]])),matrix(0,NROW(bape[[3]][[2]]),NROW(bape[[3]][[2]])));
    
    #observed data
    p.r=exp(as.matrix(predatar.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])%*%as.matrix(bape[[3]][[i.xr]]))/
      (1+exp(as.matrix(predatar.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])%*%as.matrix(bape[[3]][[i.xr]])))
    p.r=as.vector(p.r) 
    D1r[[i.xr]]=colSums((predatar.inpt[comp.ind==TRUE,NROW(ind.bape[[3]][[i.xr]])+1]-p.r)*predatar.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])
    D2r.k=t(apply(predatar.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1],1,function(x){x%o%x}))*(-p.r*(1-p.r))
    D2r.k=colSums(D2r.k)
    D2r[[i.xr]]=matrix(D2r.k,sum(ind.bape[[3]][[i.xr]]),sum(ind.bape[[3]][[i.xr]])) 
    
    #unobserved data
    for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
      mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      predamis.k=matrix(predatar.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatar.inpt),byrow=TRUE)
      predamis.k[is.na(predamis.k)==T]=mistype.k
      con.P=con.w(predamis.k[,1:(NCOL(predata)+1)],bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,])
      
      p.r=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)])%*%as.matrix(bape[[3]][[i.xr]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)])%*%as.matrix(bape[[3]][[i.xr]])))
      p.r=as.vector(p.r)
      D1r[[i.xr]]=D1r[[i.xr]]+colSums(as.vector(predamis.k[,NROW(ind.bape[[3]][[i.xr]])+1]-p.r)*predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]*con.P)
      D2r.k=t(apply(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1],1,function(x){x%o%x}))*(-p.r)*(1-p.r)*con.P
      D2r.k=colSums(D2r.k)
      D2r[[i.xr]]=D2r[[i.xr]]+matrix(D2r.k,sum(ind.bape[[3]][[i.xr]]),sum(ind.bape[[3]][[i.xr]]))
    }
    
    list(D1r,D2r)

  }
  
  Ds=function(predata,bape,ind.bape,rs.ind,comp.ind,mistype){
    predatars.inpt=matrix(0.0,NROW(predata),NCOL(predata)+NCOL(rs.ind)+1)
    predatars.inpt[,1]=rep(1,NROW(predata));
    predatars.inpt[,2:(NCOL(predata)+1)]=predata;predatars.inpt[,(NCOL(predatars.inpt)-2):(NCOL(predatars.inpt))]=rs.ind
    ###################################################################
    #  if(sum(ind.bape[[4]])==1){
    #   D1s=sum((predatars.inpt[comp.ind==TRUE,NROW(ind.bape[[4]])+1]-p.s)*predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])
    #   D2s.k=(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1]^2)*(-p.s*(1-p.s))
    #  D2s.k=sum(D2s.k)
    #  D2s=matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    
    #unobserved data
    #   for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
    #     mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
    #      predamis.k=matrix(predatars.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatars.inpt),byrow=TRUE)
    #      predamis.k[is.na(predamis.k)==T]=mistype.k
    #      con.P=con.w(predamis.k[,1:(NCOL(predata)+1)],bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,])
    
    #     p.s=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
    #   p.s=as.vector(p.s)
    #   D1s=D1s+sum(as.vector(predamis.k[,NCOL(predamis.k)]-p.s)*predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]*con.P)
    
    #  D2s.k=(predamis.k[,1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1]^2)*(-p.s)*(1-p.s)*con.P
    
    #   D2s.k=sum(D2s.k)
    #   D2s=D2s+matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    # }
    # }else{
    #observed data
    p.s=exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]])/(1+exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]]))
    p.s=as.vector(p.s)
    D1s=colSums((predatars.inpt[comp.ind==TRUE,NROW(ind.bape[[4]])+1]-p.s)*predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])
    D2s.k=t(apply(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1],1,function(x){x%o%x}))*(-p.s*(1-p.s))
    D2s.k=colSums(D2s.k)
    D2s=matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    #unobserved data
    for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
      mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      predamis.k=matrix(predatars.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatars.inpt),byrow=TRUE)
      predamis.k[is.na(predamis.k)==T]=mistype.k
      con.P=con.w(predamis.k[,1:(NCOL(predata)+1)],bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,])
      
      p.s=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
      p.s=as.vector(p.s)
      D1s=D1s+colSums(as.vector(predamis.k[,NCOL(predamis.k)]-p.s)*predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]*con.P)
      D2s.k=t(apply(predamis.k[,1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1],1,function(x){x%o%x}))*(-p.s)*(1-p.s)*con.P
      D2s.k=colSums(D2s.k)
      D2s=D2s+matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    }
    list(D1s,D2s)
  }

  
  ################################################################################################################
  
  i.EM=1
  bape0=bape
  while((i.EM<=20)||((i.EM<=25)&&max(abs(unlist(bape)-unlist(bape0)))<=0.2)){
    bape0=bape
    cat("i.EM",i.EM,"\n")
    D12y=Dy(predata,bape,ind.bape,rs.ind,comp.ind,mis.type)
    bape[[1]]=bape[[1]]-solve(D12y[[2]],D12y[[1]])
    for (i.mis in 1:(NCOL(rs.ind)-1)){
      
      if(sum(ind.bape[[2]][[i.mis]])==1){
        bape[[2]][[i.mis]]=glm(predata[,complete.cases(t(predata))==FALSE][,i.mis]~1,family=binomial(link=logit))$coefficients
        
      }else{
        D12x=Dx(predata,bape,ind.bape,rs.ind,comp.ind,mis.type,i.mis)
        bape[[2]][[i.mis]]=bape[[2]][[i.mis]]-solve(D12x[[2]][[i.mis]],D12x[[1]][[i.mis]])
      }

      
      if(sum(ind.bape[[3]][[i.mis]])==1){
        bape[[3]][[i.mis]]=glm(rs.ind[,i.mis]~1,family=binomial(link=logit))$coefficients
      }else{
      D12r=Dr(predata,bape,ind.bape,rs.ind,comp.ind,mis.type,i.mis)
      bape[[3]][[i.mis]]=bape[[3]][[i.mis]]-solve(D12r[[2]][[i.mis]],D12r[[1]][[i.mis]])}
    }
    if(sum(ind.bape[[4]])==1){
      bape[[4]]=glm(rs.ind[,3]~1,family=binomial(link=logit))$coefficients
    }else{D12s=Ds(predata,bape,ind.bape,rs.ind,comp.ind,mis.type)
    bape[[4]]=bape[[4]]-solve(D12s[[2]],D12s[[1]])}
    #print(class(bape[[4]]))
    i.EM=i.EM+1
  }
  #print(bape)
  cat("i.EM",i.EM,"\n")

      # definition of BIC
    BIC.def=function(predata,bape,ind.bape,rs.ind,comp.ind,mis.type){
      p=NROW(unlist(bape))
      #log information of BIC
      log.BIC=function(predata,bape,ind.bape,rs.ind,comp.ind,mis.type){
        predatars.inpt=matrix(0.0,NROW(predata),NCOL(predata)+NCOL(rs.ind)+1)
        predatars.inpt[,1]=rep(1,NROW(predata));
        predatars.inpt[,2:(NCOL(predata)+1)]=predata;predatars.inpt[,(NCOL(predatars.inpt)-2):(NCOL(predatars.inpt))]=rs.ind
        
        #observed data
        p.y=exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]])/(1+exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]]))
        p.y=as.vector(p.y)
        L.obs=dbinom(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[1]])+1],1,p.y)
        for (i.xr in 1:2){
          p.x=exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])%*%as.matrix(bape[[2]][[i.xr]]))/
            (1+exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])%*%as.matrix(bape[[2]][[i.xr]])))
          p.x=as.vector(p.x)
          L.obs=L.obs*dbinom(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[2]][[i.xr]])+1],1,p.x)
          p.r=exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])%*%as.matrix(bape[[3]][[i.xr]]))/
            (1+exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])%*%as.matrix(bape[[3]][[i.xr]])))
          p.r=as.vector(p.r)
          L.obs=L.obs*dbinom(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[3]][[i.xr]])+1],1,p.r)
        }
        p.s=exp(as.matrix(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
        p.s=as.vector(p.s)
        L.obs=L.obs*dbinom(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[4]])+1],1,p.s)
        l.obs=sum(log(L.obs))
        
        #unobserved data
        for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
          mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
          predamis.k=matrix(predatars.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatars.inpt),byrow=TRUE)
          predamis.k[is.na(predamis.k)==T]=mistype.k
          con.P=con.w(predamis.k[,1:(NCOL(predata)+1)],bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,])
          
          p.ymis=exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]])/(1+exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]]))
          p.ymis=as.vector(p.ymis)
          denom.k=dbinom(predamis.k[,NROW(ind.bape[[1]])+1],1,p.ymis)
          
          for (i.xr in 1:(NCOL(rs.ind)-1)){
            p.xrmis=exp(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]%*%bape[[2]][[i.xr]])/(1+exp(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]%*%bape[[2]][[i.xr]]))
            p.xrmis=as.vector(p.xrmis)
            denom.k=denom.k*dbinom(predamis.k[,NROW(ind.bape[[2]][[i.xr]])+1],1,p.xrmis)
            
            p.xrmis=exp(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]%*%bape[[3]][[i.xr]])/(1+exp(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]%*%bape[[3]][[i.xr]]))
            p.xrmis=as.vector(p.xrmis)
            denom.k=denom.k*dbinom(predamis.k[,NROW(ind.bape[[3]][[i.xr]])+1],1,p.xrmis)
          }
          
          p.smis=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
          p.smis=as.vector(p.smis)
          denom.k=denom.k*dbinom(predamis.k[,NROW(ind.bape[[4]])+1],1,p.smis)
          l.obs=l.obs+sum(log(denom.k)*con.P)
        }
        
        return(l.obs)
      }
      
      def.val=-2*log.BIC(predata,bape,ind.bape,rs.ind,comp.ind,mis.type)+p*log(NROW(predata))
      return(def.val)
    }
    BIC.def(predata,bape,ind.bape,rs.ind,comp.ind,mis.type)
 }




