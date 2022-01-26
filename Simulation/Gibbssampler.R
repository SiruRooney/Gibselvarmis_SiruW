#Gibbs sampler for GLM in missing data

Gibbs.sampler<-function(predata,rs.ind,ind,iter.Gibs){
  
  ###########################################################################################
  # the definition of full conditional distribution
  FCP<-function(predata,rs.ind,ind,rr,cc){#n from the last to the second
    source("F://Missingdatacodes//MissingGibbs2//MisGibsfunction3.R")
    ind[[rr]][cc]<-1
    fun1.bic=EM.BIC(predata,BAPE(ind),ind.BAPE(ind),rs.ind)
    cat("fun1",fun1.bic,"\n")
    ###################################################################################
    ind[[rr]][cc]<-0
    fun2.bic=EM.BIC(predata,BAPE(ind),ind.BAPE(ind),rs.ind)
    cat("fun2",fun2.bic,"\n")
    res<-1/{1+exp(0.8*(-fun2.bic+fun1.bic))}
    return(c(res,min(fun1.bic,fun2.bic)))
  }
  
  
  #tab.matrix=list(matrix(0,(NROW(unlist(ind))-NROW(ind))*iter.Gibs,NROW(ind[[1]])),matrix(0,(NROW(unlist(ind))-NROW(ind))*iter.Gibs,NROW(ind[[2]])),
  #               matrix(0,(NROW(unlist(ind))-NROW(ind))*iter.Gibs,NROW(ind[[3]])),matrix(0,(NROW(unlist(ind))-NROW(ind))*iter.Gibs,NROW(ind[[4]])),
  #              matrix(0,(NROW(unlist(ind))-NROW(ind))*iter.Gibs,NROW(ind[[5]])),matrix(0,(NROW(unlist(ind))-NROW(ind))*iter.Gibs,NROW(ind[[6]])))
  
  tab.matrix=lapply(ind,function(x,y){matrix(0,(NROW(unlist(y[[1]]))-NROW(y[[1]]))*y[[2]],NROW(x)+1)},list(ind,iter.Gibs))
  
  i.tabmatrix=1
  for(i.Gib in 1:iter.Gibs){#cat("Gib",i.Gib,"\n")
    
    for(rr in NROW(ind):1){
      for (cc in NROW(ind[[rr]]):2){cat("n",c(i.Gib,rr,cc),"\n")
        print(ind)
        fcp=FCP(predata,rs.ind,ind,rr,cc)
        ind[[rr]][cc]<-rbinom(1,1,fcp[1])
        for (i.tab in 1:NROW(ind)){tab.matrix[[i.tab]][i.tabmatrix,]=c(fcp[2],ind[[i.tab]])}
        i.tabmatrix=i.tabmatrix+1
      } 
    }
  }
  list(tab.matrix)
}
