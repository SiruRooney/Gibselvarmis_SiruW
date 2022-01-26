

Gibbs.sampler<-function(predata,rs.ind,ind,iter.Gibs,wgh){
  source("F://Missingdatacodes//Missing_data_italy_weighted//MisGibsfunction3.9.3.R")
  
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