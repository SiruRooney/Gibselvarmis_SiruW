
fisherinfo=function(predata,bape,ind.bape,rs.ind,wgh){
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
    p.ymis=exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]])/(1+exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]]))
    p.ymis=as.vector(p.ymis)
    denom.k=dbinom(round(predamis.k[,NCOL(predamis.k)]*wgh.k),round(wgh.k),p.ymis)
    
    predamisr.k=matrix(0.0,nrow=NROW(predamis.k),ncol=NCOL(predamis.k)+2,byrow=TRUE);
    predamisr.k[,1:NCOL(predamis.k)]=predamis.k;predamisr.k[,-(1:NCOL(predamis.k))]=matrix(rs.k[-NROW(rs.k)],nrow=NROW(predamis.k),ncol=NROW(rs.k)-1,byrow=TRUE)
    for (i.xr in 1:(NROW(rs.k)-1)){
      p.xrmis=exp(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]%*%bape[[2]][[i.xr]])/(1+exp(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]%*%bape[[2]][[i.xr]]))
      p.xrmis=as.vector(p.xrmis)
      denom.k=denom.k*dbinom(round(predamis.k[,NROW(ind.bape[[2]][[i.xr]])+1]*wgh.k),round(wgh.k),p.xrmis)
      
      p.xrmis=exp(predamisr.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]%*%bape[[3]][[i.xr]])/(1+exp(predamisr.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]%*%bape[[3]][[i.xr]]))
      p.xrmis=as.vector(p.xrmis)
      denom.k=denom.k*dbinom(round(rs.k[i.xr]*wgh.k),round(wgh.k),p.xrmis)
    }
    
    p.smis=exp(predamisr.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]])/(1+exp(predamisr.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]]))
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

    predatars.inpt=matrix(0.0,NROW(predata),NCOL(predata)+NCOL(rs.ind)+1)
    predatars.inpt[,1]=rep(1,NROW(predata));
    predatars.inpt[,2:(NCOL(predata)+1)]=predata;predatars.inpt[,(NCOL(predatars.inpt)-2):(NCOL(predatars.inpt))]=rs.ind
    
    #observed data
    p.y=exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]])/(1+exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]]))
    p.y=as.vector(p.y)

    p.x1=exp(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[1]])][,ind.bape[[2]][[1]]==1]%*%bape[[2]][[1]])/
      (1+exp(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[1]])][,ind.bape[[2]][[1]]==1]%*%bape[[2]][[1]]))
    p.x1=as.vector(p.x1)
    p.x2=exp(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[2]])][,ind.bape[[2]][[2]]==1]%*%bape[[2]][[2]])/
      (1+exp(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[2]])][,ind.bape[[2]][[2]]==1]%*%bape[[2]][[2]]))
    p.x2=as.vector(p.x2)
    p.r1=exp(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[1]])][,ind.bape[[3]][[1]]==1]%*%bape[[3]][[1]])/
      (1+exp(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[1]])][,ind.bape[[3]][[1]]==1]%*%bape[[3]][[1]]))
    p.r1=as.vector(p.r1) 
    p.r2=exp(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[2]])][,ind.bape[[3]][[2]]==1]%*%bape[[3]][[2]])/
      (1+exp(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[2]])][,ind.bape[[3]][[2]]==1]%*%bape[[3]][[2]]))
    p.r2=as.vector(p.r2)
    
    p.s=exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]])/(1+exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]]))
    p.s=as.vector(p.s)
    

    y1T=(predatars.inpt[comp.ind==TRUE,NROW(ind.bape[[1]])+1]-p.y)*predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]
    x1T1=(predatars.inpt[comp.ind==TRUE,NROW(ind.bape[[2]][[1]])+1]-p.x1)*predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[1]])][,ind.bape[[2]][[1]]==1]
    x1T2=(predatars.inpt[comp.ind==TRUE,NROW(ind.bape[[2]][[2]])+1]-p.x2)*predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[2]])][,ind.bape[[2]][[2]]==1]
    r1T1=(predatars.inpt[comp.ind==TRUE,NROW(ind.bape[[3]][[1]])+1]-p.r1)*predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[1]])][,ind.bape[[3]][[1]]==1]
    r1T2=(predatars.inpt[comp.ind==TRUE,NROW(ind.bape[[3]][[2]])+1]-p.r2)*predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[2]])][,ind.bape[[3]][[2]]==1]
    s1T=(predatars.inpt[comp.ind==TRUE,NROW(ind.bape[[4]])+1]-p.s)*predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]
  
    #D1y=colSums(y1T)#;print(D1y)
    #D1x1=colSums(x1T1)#;print(D1x1)
    #D1x2=colSums(x1T2)#;print(D1x2)
    #D1r1=colSums(r1T1)#;print(D1r1)
    #D1r2=colSums(r1T2)#;print(D1r2)
    #D1s=colSums(s1T)#;print(D1s)
    yxrs.T=wgh[(comp.ind==TRUE)]*cbind(y1T,x1T1,x1T2,r1T1,r1T2,s1T)
    D1.T=colSums(t(apply(yxrs.T,1,function(x){x%o%x})))
    D1=matrix(D1.T,sum(unlist(ind.bape)),sum(unlist(ind.bape)))
    D12.T=colSums(t(apply(yxrs.T,1,function(x){x%o%x})))
    D12=matrix(D12.T,sum(unlist(ind.bape)),sum(unlist(ind.bape)))
 
    
    D2y.k=t(apply(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,ind.bape[[1]]==1],1,function(x){x%o%x}))*(-p.y*(1-p.y)*wgh[(comp.ind==TRUE)])
    D2y.k=colSums(D2y.k)
    D2y=matrix(D2y.k,sum(ind.bape[[1]]),sum(ind.bape[[1]]))
    
    D2x1.k=t(apply(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[2]][[1]])][,ind.bape[[2]][[1]]==1],1,function(x){x%o%x}))*(-p.x1*(1-p.x1)*wgh[(comp.ind==TRUE)])
    D2x1.k=colSums(D2x1.k)
    D2x1=matrix(D2x1.k,sum(ind.bape[[2]][[1]]),sum(ind.bape[[2]][[1]])) 
    
    D2x2.k=t(apply(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[2]][[2]])][,ind.bape[[2]][[2]]==1],1,function(x){x%o%x}))*(-p.x2*(1-p.x2)*wgh[(comp.ind==TRUE)])
    D2x2.k=colSums(D2x2.k)
    D2x2=matrix(D2x2.k,sum(ind.bape[[2]][[2]]),sum(ind.bape[[2]][[2]])) 
    
    D2r1.k=t(apply(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[3]][[1]])][,ind.bape[[3]][[1]]==1],1,function(x){x%o%x}))*(-p.r1*(1-p.r1)*wgh[(comp.ind==TRUE)])
    D2r1.k=colSums(D2r1.k)
    D2r1=matrix(D2r1.k,sum(ind.bape[[3]][[1]]),sum(ind.bape[[3]][[1]])) 
    
    D2r2.k=t(apply(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[3]][[2]])][,ind.bape[[3]][[2]]==1],1,function(x){x%o%x}))*(-p.r2*(1-p.r2)*wgh[(comp.ind==TRUE)])
    D2r2.k=colSums(D2r2.k)
    D2r2=matrix(D2r2.k,sum(ind.bape[[3]][[2]]),sum(ind.bape[[3]][[2]]))
    
    D2s.k=t(apply(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1],1,function(x){x%o%x}))*(-p.s*(1-p.s)*wgh[(comp.ind==TRUE)])
    D2s.k=colSums(D2s.k)
    D2s=matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    #D1=c(D1y,D1x1,D1x2,D1r1,D1r2,D1s)
    
    print(any(D12==D1))
    
    
    #unobserved data
    for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){#NROW(predata[comp.ind==FALSE,])
      mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      predamis.k=matrix(predatars.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatars.inpt),byrow=TRUE)
      predamis.k[is.na(predamis.k)==T]=mistype.k
      con.P=con.w(predamis.k[,1:(NCOL(predata)+1)],bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,],wgh[comp.ind==FALSE][i.mis])
      
      p.y=exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]])/(1+exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]]))
      p.y=as.vector(p.y)
      p.x1=exp(predamis.k[,1:NROW(ind.bape[[2]][[1]])][,(ind.bape[[2]][[1]]==1)]%*%bape[[2]][[1]])/(1+exp(predamis.k[,1:NROW(ind.bape[[2]][[1]])][,(ind.bape[[2]][[1]]==1)]%*%bape[[2]][[1]]))
      p.x1=as.vector(p.x1)
      p.x2=exp(predamis.k[,1:NROW(ind.bape[[2]][[2]])][,(ind.bape[[2]][[2]]==1)]%*%bape[[2]][[2]])/(1+exp(predamis.k[,1:NROW(ind.bape[[2]][[2]])][,(ind.bape[[2]][[2]]==1)]%*%bape[[2]][[2]]))
      p.x2=as.vector(p.x2)
      p.r1=exp(predamis.k[,1:NROW(ind.bape[[3]][[1]])][,(ind.bape[[3]][[1]]==1)]%*%bape[[3]][[1]])/(1+exp(predamis.k[,1:NROW(ind.bape[[3]][[1]])][,(ind.bape[[3]][[1]]==1)]%*%bape[[3]][[1]]))
      p.r1=as.vector(p.r1)
      p.r2=exp(predamis.k[,1:NROW(ind.bape[[3]][[2]])][,(ind.bape[[3]][[2]]==1)]%*%bape[[3]][[2]])/(1+exp(predamis.k[,1:NROW(ind.bape[[3]][[2]])][,(ind.bape[[3]][[2]]==1)]%*%bape[[3]][[2]]))
      p.r2=as.vector(p.r2)
      p.s=exp(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]])/(1+exp(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]]))
      p.s=as.vector(p.s)
      
      y1F=as.vector(predamis.k[,NROW(ind.bape[[1]])+1]-p.y)*predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]
      x1F1=as.vector(predamis.k[,NROW(ind.bape[[2]][[1]])+1]-p.x1)*predamis.k[,1:NROW(ind.bape[[2]][[1]])][,(ind.bape[[2]][[1]]==1)]
      x1F2=as.vector(predamis.k[,NROW(ind.bape[[2]][[2]])+1]-p.x2)*predamis.k[,1:NROW(ind.bape[[2]][[2]])][,(ind.bape[[2]][[2]]==1)]
      r1F1=as.vector(predamis.k[,NROW(ind.bape[[3]][[1]])+1]-p.r1)*predamis.k[,1:NROW(ind.bape[[3]][[1]])][,(ind.bape[[3]][[1]]==1)]
      r1F2=as.vector(predamis.k[,NROW(ind.bape[[3]][[2]])+1]-p.r2)*predamis.k[,1:NROW(ind.bape[[3]][[2]])][,(ind.bape[[3]][[2]]==1)]
      s1F=as.vector(predamis.k[,NROW(ind.bape[[4]])+1]-p.s)*predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]
      #print(con.P)
      #D1y=D1y+colSums(y1F*con.P)
      #D1x1=D1x1+colSums(x1F1*con.P)
      #D1x2=D1x2+colSums(x1F2*con.P)
      #D1r1=D1r1+colSums(r1F1*con.P)
      #D1r2=D1r2+colSums(r1F2*con.P)
      #D1s=D1s+colSums(s1F*con.P)
     
      yxrs.F=wgh[(comp.ind==FALSE)][i.mis]*cbind(y1F,x1F1,x1F2,r1F1,r1F2,s1F)
      D1.F=colSums(yxrs.F*con.P)%o%colSums(yxrs.F*con.P)
      D1=D1+D1.F#matrix(D1.F,sum(unlist(ind.bape)),sum(unlist(ind.bape)))
      D12.F=t(apply(yxrs.F,1,function(x){x%o%x}))*con.P
      D12.F=colSums(D12.F)
      D12=D12+matrix(D12.F,sum(unlist(ind.bape)),sum(unlist(ind.bape)))

      D2y.k=t(apply(predamis.k[,1:NROW(ind.bape[[1]])][,ind.bape[[1]]==1],1,function(x){x%o%x}))*(-wgh[(comp.ind==FALSE)][i.mis]*p.y)*(1-p.y)*con.P
      D2y.k=colSums(D2y.k)
      D2y=D2y+matrix(D2y.k,sum(ind.bape[[1]]),sum(ind.bape[[1]]))
      D2x1.k=t(apply(predamis.k[,1:NROW(ind.bape[[2]][[1]])][,ind.bape[[2]][[1]]==1],1,function(x){x%o%x}))*(-wgh[(comp.ind==FALSE)][i.mis]*p.x1)*(1-p.x1)*con.P
      D2x1.k=colSums(D2x1.k)
      D2x1=D2x1+matrix(D2x1.k,sum(ind.bape[[2]][[1]]),sum(ind.bape[[2]][[1]]))
      D2x2.k=t(apply(predamis.k[,1:NROW(ind.bape[[2]][[2]])][,ind.bape[[2]][[2]]==1],1,function(x){x%o%x}))*(-wgh[(comp.ind==FALSE)][i.mis]*p.x2)*(1-p.x2)*con.P
      D2x2.k=colSums(D2x2.k)
      D2x2=D2x2+matrix(D2x2.k,sum(ind.bape[[2]][[2]]),sum(ind.bape[[2]][[2]]))
      D2r1.k=t(apply(predamis.k[,1:NROW(ind.bape[[3]][[1]])][,ind.bape[[3]][[1]]==1],1,function(x){x%o%x}))*(-wgh[(comp.ind==FALSE)][i.mis]*p.r1)*(1-p.r1)*con.P
      D2r1.k=colSums(D2r1.k)
      D2r1=D2r1+matrix(D2r1.k,sum(ind.bape[[3]][[1]]),sum(ind.bape[[3]][[1]]))
      D2r2.k=t(apply(predamis.k[,1:NROW(ind.bape[[3]][[2]])][,ind.bape[[3]][[2]]==1],1,function(x){x%o%x}))*(-wgh[(comp.ind==FALSE)][i.mis]*p.r2)*(1-p.r2)*con.P
      D2r2.k=colSums(D2r2.k)
      D2r2=D2r2+matrix(D2r2.k,sum(ind.bape[[3]][[2]]),sum(ind.bape[[3]][[2]]))
      D2s.k=t(apply(predamis.k[,1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1],1,function(x){x%o%x}))*(-wgh[(comp.ind==FALSE)][i.mis]*p.s)*(1-p.s)*con.P
      D2s.k=colSums(D2s.k)
      D2s=D2s+matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    }

    #print(sum(unlist(ind.bape)))
         
    D2=matrix(0,sum(unlist(ind.bape)),sum(unlist(ind.bape)))
    D2[1:NROW(D2y),1:NCOL(D2y)]=D2y;
    #print(D2y)
    D2[(NROW(D2y)+1):(NROW(D2y)+NROW(D2x1)),(NCOL(D2y)+1):(NCOL(D2y)+NCOL(D2x1))]=D2x1
    D2[-(1:(NROW(D2y)+NROW(D2x1))),-(1:(NCOL(D2y)+NCOL(D2x1)))][1:NROW(D2x2),1:NCOL(D2x2)]=D2x2
    D2[-(1:(NROW(D2y)+NROW(D2x1)+NROW(D2x2))),-(1:(NCOL(D2y)+NCOL(D2x1)+NCOL(D2x2)))][1:NROW(D2r1),1:NCOL(D2r1)]=D2r1
    D2[-(1:(NROW(D2y)+NROW(D2x1)+NROW(D2x2)+NROW(D2r1))),-(1:(NCOL(D2y)+NCOL(D2x1)+NCOL(D2x2)+NCOL(D2r1)))][1:NROW(D2r2),1:NCOL(D2r2)]=D2r2
    D2[-(1:(NROW(D2y)+NROW(D2x1)+NROW(D2x2)+NROW(D2r1)+NROW(D2r2))),-(1:(NCOL(D2y)+NCOL(D2x1)+NCOL(D2x2)+NCOL(D2r1)+NCOL(D2r2)))]=D2s
    #print(NROW(D2y)+NROW(D2x1)+NROW(D2x2)+NROW(D2r1)+NROW(D2r2)+NROW(D2s))
    #D1=c(D1y,D1x1,D1x2,D1r1,D1r2,D1s)
    #print(D1)
    return(-D2-D12+D1)    
}


