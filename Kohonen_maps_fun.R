require(HistDAWass)
require(Rcpp)
require(RcppArmadillo)
WH_SOM_Topo_error=function(results2,toroidal=FALSE){
  #topographic error
  #validity indices
  ########################################
  
  # results2=WH_2d_Adaptive_Kohonen_maps_test(Age_Pyramids_2014,
  #                                           net=list(xdim=5,ydim=5,topo=c('hexagonal')),
  #                                           simplify = TRUE,qua = 20,schema = 6,
  #                                           weight.sys = "SUM",theta = 2,
  #                                           niter = 50,TMAX = 2,Tmin = 0.4)
  ########################################
  mymap=results2$solution$MAP
  if (!is.null(results2$solution$D2P)){
    D2P=results2$solution$D2P}
  proto=results2$solution$proto
  ID=results2$solution$IDX
  Nprot=get.MatH.nrows(proto)
  WW=results2$solution$weights.comp
  if(results2$solution$Weight.sys=="SUM"){WW=WW^2}
  N=sum(results2$solution$cardinality)
 # N=get.MatH.nrows(data)
#  V=get.MatH.ncols(data)
  
  #topo.error
  #compute distances between objects and protos
  #get the winning proto for each object
  #list neigbours
  distneu=unit.distances(mymap,toroidal = toroidal)#as.matrix(dist(mymap$pts))
  
  distneu[lower.tri(distneu)] <- 0
  neigb=which((distneu<1.1&distneu>0.1),arr.ind=T)
  #
  winning_proto=results2$solution$IDX
  sec_best=rep(0,N)
  DATA_TO_PRO=matrix(0,N,Nprot)
  for (i in 1:N){
    wp=winning_proto[i]
    #search for closest neuron
    Dtmp=rep(0,Nprot)
    Dtmp=D2P[i,]
    
    Dtmp[wp]=.Machine$double.xmax
    sec_best[i]=which.min(Dtmp)
    #  sec_best[i]=which(rank(Dtmp,ties.method = "first")==2)
    #print(i)
  }
  is_connected=rep(1,N)

  for (i in 1:N){
    prox=distneu[min(winning_proto[i],sec_best[i]),max(winning_proto[i],sec_best[i])]
    # sqrt(sum((mymap$pts[winning_proto[i],]-
    #                mymap$pts[sec_best[i],])^2))
    if(prox<1.1){is_connected[i]=0}
    
  }
  
  topo.error=mean(is_connected)
  return(topo.error)
  # print(paste("Topol error ", topo.error))
  # #if (topo.error>0.5){browser()}
  # #topo.error
  # quantiz.error=results2$solution$Crit
  # #quantization error
}
WH_2d_Adaptive_Kohonen_maps_test=function (x,net=list(xdim=4,ydim=3,topo=c('rectangular')), 
                                           kern.param=2, TMAX=-9999, Tmin=-9999, 
                                           niter=30,repetitions ,
                                           simplify=FALSE,
                                           qua=10,
                                           standardize=FALSE, schema=6,
                                           init.weights='EQUAL',weight.sys='PROD',
                                           theta=2,Wfix=FALSE,verbose=FALSE,atleast=1, 
                                           init="random", coords=0,
                                           ASS_DC=FALSE, bubble=FALSE, TOE=FALSE, toroidal=FALSE){
  tol=1e-10
  if(weight.sys=='PROD') {
    weight_sys=1
    theta=1
  }
  if(weight.sys=='SUM') {weight_sys=2}
  # remember to fix TMAX e Tmin
  #require(class)#for somgrid function
  # Check initial conditions and passed parameters
  
  ind=nrow(x@M)
  vars=ncol(x@M)
  ## we homogeneize data for speeding up the code if required
  tmp=Prepare(x,simplify,qua,standardize)
  MM=tmp$MM
  x=tmp$x
  remove(tmp)
  # end of the preprocessing step
  
  TOTSSQ=0
  SSQ_comp=rep(0,vars*2)
  tmpc=as.matrix(get.MatH.stats(x)$mat)
  
  colm=colMeans(tmpc)
  for (v in 1:vars){
    tmpSSQ=WH.SSQ(x[,v])
    TOTSSQ=TOTSSQ+tmpSSQ
    C_comp=sum((tmpc[,v]-colm[v])^2)
    V_comp=tmpSSQ-C_comp
    SSQ_comp[v*2-1]=C_comp
    SSQ_comp[v*2]=V_comp
  }
  ##------------------------------------------------batchKOHONENMAP-----------------------------------------------------
  solutions=list()
  criteria=numeric(0)
  toes=numeric()
  MAP=somgrid(net$xdim,net$ydim,net$topo)
  k=nrow(MAP$pts)
  if (init=="PCA") repetitions=1
  if (missing(repetitions)){repetitions=5}
  
  restart=1
  repet=0
  TopE=1.1
  while(repet<repetitions){
    repet=repet+1
    if (init=="PCA"){
      proto=x[initPCA2(coords,MAP),]
    }else{
      #random selection of prototypes of neurons
      
      if(k<=ind){
        proto=x[sample(1L:ind,k, replace=FALSE),]
      }else{
        proto=x[sample(1L:ind,k, replace=TRUE),]
      }
    }
    rownames(proto@M)=paste0("Prot. ",1:k)
    nhbrdist=unit.distances(MAP,toroidal = toroidal) #as.matrix(dist(MAP$pts))
   
    # dmax=(max(MAP$pts[,1])-min(MAP$pts[,1]))^2+(max(MAP$pts[,2])-min(MAP$pts[,2]))^2
    # dmax=(max(as.matrix(dist(MAP$pts))))^2
    dmax=max(nhbrdist)
    # browser()
    if (TMAX<0){
      TMAX=sqrt(-((0.5*dmax)^2)/(2*log(0.1)))
    }
    if (Tmin<0){
      Tmin=sqrt(-1/(2*log(0.01)))

    print(paste("Tmax=",TMAX,"Tmin=",Tmin))
      }
    
    
    TT=TMAX
    KT=exp(-(nhbrdist^2)/(2*TT^2))
    
    if (bubble) KT[which(nhbrdist<=TT)]=1
    
    KT[which(KT<1e-50)]=0
    #initialize weights
    if (Wfix){lambdas=matrix(1,2*vars,k)}else{
      if (init.weights=='EQUAL'){
        cat("Weights initialization === EQUAL  \n")
        if (weight.sys=='PROD'){
          lambdas=matrix(1,2*vars,k)}
        else{
          if (schema==1||schema==2||schema==3||schema==4){
            lambdas=matrix(1/vars,2*vars,k)}
          if(schema==5||schema==6){
            lambdas=matrix(1/(2*vars),2*vars,k)
          }
        }
      }
      if (init.weights=='PROP'){
        if (weight.sys=='PROD'){
          # 1=A weight for each variable (default) 
          # 2=A weight for the average and the dispersion component of each variable
          # 3=Same as 1 but a different set of weights for each cluster (indipprod=1)
          # 4=Same as 2 but a different set of weights for each cluster (indipprod=1)
          # 5=Same as 1 but a different set of weights for each cluster (allprod=1)
          # 6=Same as 2 but a different set of weights for each cluster (allprod=1)            
          if (schema==1||schema==3||schema==5){
            tmpvals=SSQ_comp[c(1:v)*2-1]+SSQ_comp[c(1:v)*2]
            mg=(prod(tmpvals))^(1/v)
            tmpw=mg/tmpvals
            elew=rep(tmpw,2,each=2)
            lambdas=matrix(rep(elew,k),nrow = 2*vars,ncol = k,byrow = FALSE)
          }
          if (schema==2||4){
            tmpvalsc= SSQ_comp[c(1:v)*2-1]
            tmpvalsv= SSQ_comp[c(1:v)*2]
            mgc=(prod(tmpvalsc))^(1/v)
            mgv=(prod(tmpvalsv))^(1/v)
            
            tmpwc=mgc/tmpvalsc
            tmpwv=mgv/tmpvalsv
            elew=rep(0,2*v)
            elew[c(1:v)*2-1]=tmpwc
            elew[c(1:v)*2]=tmpwv
            
            lambdas=matrix(rep(elew,k),nrow = 2*vars,ncol = k,byrow = FALSE)
          }
          if (schema==6){
            tmpvals=SSQ_comp
            mg=(prod(tmpvals))^(1/(2*v))
            tmpw=mg/tmpvals
            lambdas=matrix(rep(tmpw,k),nrow = 2*vars,ncol = k,byrow = FALSE)
          }
          
        }else{
          if ((schema==1||schema==2||schema==3||schema==4)){
            lambdas=matrix(1/v,nrow = 2*vars,ncol = k,byrow = FALSE)
          }
          if (schema==5||schema==6){
            lambdas=matrix(2/(2*v),nrow = 2*vars,ncol = k,byrow = FALSE)
          }
          
        }
        
      }
      if ((init.weights!='EQUAL')&(init.weights!='PROP')) {#Random initialization
        cat("Weights initialization === RANDOM  \n")
        m1=matrix(runif((vars*k),0.01,0.99),vars,k)
        m2=matrix(runif((vars*k),0.01,0.99),vars,k)
        m1=m1/matrix(rep(apply(m1,2,sum)),vars,k,byrow = TRUE)
        m2=m2/matrix(rep(apply(m2,2,sum)),vars,k,byrow = TRUE)
        if (weight.sys=='PROD'){
          m1=exp(m1*vars-1)
          m2=exp(m2*vars-1)
          
        }
        
        if (schema==1){
          m1=matrix(m1[,1],nrow = vars,ncol = k)
          m2=m1
          
        }
        if (schema==2){
          m1=matrix(rep(m1[,1],k),vars,k)
          m2=matrix(rep(m2[,1],k),vars,k)
        }
        if (schema==3){m2=m1}
        
        lambdas=matrix(0,2*vars,k)
        colnames(lambdas)=paste('Clust',c(1:k),sep="_")
        
        n1=paste('M_Var.',c(1:vars),sep="_")
        n2=paste('C_Var.',c(1:vars),sep="_")
        nr=list()
        for (nn in 1:vars){
          nr[[nn*2-1]]=n1[[nn]]
          nr[[nn*2]]=n2[[nn]]
        }
        rownames(lambdas)=nr
        lambdas[(c(1:vars)*2-1),]=m1
        lambdas[(c(1:vars)*2),]=m2
      }
    }
    
    
    
    #assign objects to closest neuron
    #compute adaptive distances to prototypes
    distances=array(0,dim=c(vars,k,2))
    diINDtoPROT=array(0,dim=c(ind,vars,k,2))
    # compute distances
    MD=(c_ComputeFast_L2_SQ_WASS_DMAT(MM,proto))$Dist
    
    KT[KT<1e-50]=0
################################
#   KT=KT/matrix(rowSums(KT),nrow(KT),ncol(KT),byrow=FALSE)
################################    
    Dind2Neu=MD%*%KT
    #initialize matrix of memberships
    #IDX=apply(Dind2Neu,1,which.min)
    if (ASS_DC){
    IDX=max.col(-Dind2Neu, ties="first")}else{
    IDX=max.col(-MD, ties="first")}
    memb=matrix(0,ind,k)
    for (i in (1:ind)){
      memb[i,IDX[i]]=1
    }
    cat(paste("---------> repetitions  ",repet,"\n"))
    #browser()
    ## initialize clusters and prototypes
    
    SSQ1=matrix(0,k,vars)
    GenCrit=Inf
    ##  compute initial criterion
    tmpIDX=Table2(IDX, k, FALSE)
    
    SSQ=sum(memb*Dind2Neu)
    GenCrit=SSQ
    OK=1 
    #maxit=200
    
    t=0
    ft=20
    fts=0
    OldCrit=GenCrit+1
    while ((abs(TT-Tmin)>1e-10)||(fts<=ft)){
      if (t<(niter-1)){
        
        ################ASSESSING #############################
        #if ((abs(OldCrit-GenCrit)<1e-20)) t=t+1    ############
        #######################################################
        t=t+1
        #browser()
        TT=TMAX*(Tmin/TMAX)^(t/(niter-1))
        KT=exp(-(nhbrdist^2)/(2*TT^2))
        ################################
        # KT=KT/matrix(rowSums(KT),nrow(KT),ncol(KT),byrow=FALSE)
        ################################    
        
        if (bubble) KT[which(nhbrdist<=TT)]=1
        KT[which(KT<1e-20)]=0
      }
      if((t==(niter-1))&(fts<=ft)){
       # print(fts)
        fts=fts+1
      }
   #   if(t==(niter-1)){browser()}
      #computing prototypes
      CardinalityOfNeuron=colSums(memb)
      
      proto=c_PROTO_KOHONEN(proto,k,ind,MM,vars,KT, IDX)
      
      #########################################################################################
      
      
      #first compute distances using kernel
      #  computing distances
      
      ###########    NEW part
      tmpD=c_ComputeFast_L2_SQ_WASS_DMAT(MM,proto)
      #browser()
      GRKT=KT[IDX,1:k]
      distances=array(0,dim=c(vars,k,2))
     
      for (j in 1:vars){
        distances[j,,1]=distances[j,,1]+colSums(tmpD$DET[[j]]$DM*GRKT)
        distances[j,,2]=distances[j,,2]+colSums(tmpD$DET[[j]]$DV*GRKT)
      }
      
      #diINDtoPROT=array(0,dim=c(ind,vars,k,2))
      tmp2=matrix(0,ind,k)
      if (Wfix==FALSE){#Weights computation
        # S2.2) Weights computation
        #NEW
        if (t>0){oldlambda=lambdas}
       # browser()
        lambdas=c_ADA_F_WHEIGHT(list(as.matrix(distances[,,1]),as.matrix(distances[,,2])),
                                k, vars, ind, schema, weight_sys,  theta, 1)
        
        if (t>0){
          
          for (cluster in 1:k){
            if (weight.sys=="PROD"){
              if ((schema>=1 & schema<=4)&
                  ((abs(prod(lambdas[,cluster])-1))>0.0001)){
                print(cat("new",prod(lambdas[,cluster])))
                lambdas[,cluster]=oldlambda[,cluster]
                print(cat("t",t, "CLU", cluster, "new",prod(lambdas[,cluster])))
                
                # browser()
              }
              if ((schema==5||schema==6)&
                  ((abs(prod(lambdas[,cluster])-1))>0.0001)){
                print(cat("new",prod(lambdas[,cluster])))
                lambdas[,cluster]=oldlambda[,cluster]
                print(cat("t",t, "CLU", cluster, "new",prod(lambdas[,cluster])))
              }
            }
            if (weight.sys=="SUM"){
              if ((schema>=1 & schema<=4)&
                  ((abs(sum(lambdas[,cluster])-2))>0.0001)){
                print(cat("new",sum(lambdas[,cluster])))
                lambdas[,cluster]=oldlambda[,cluster]
                print(cat("t",t, "CLU", cluster, "new",sum(lambdas[,cluster])))
                
                # browser()
              }
              if ((schema==5||schema==6)&
                  ((abs(sum(lambdas[,cluster])-1))>0.0001)){
                print(cat("new",sum(lambdas[,cluster])))
                lambdas[,cluster]=oldlambda[,cluster]
                print(cat("t",t, "CLU", cluster, "new",sum(lambdas[,cluster])))
                
              }
            } 
            
            
          }
        }
        wM=(lambdas[((1:vars)*2-1),])^theta
        wV=(lambdas[(1:vars)*2,])^theta
        
        
        #assign data to neurons
        ####################
        if (schema==1|schema==3){
          ####################
          for (variables in (1:vars)){
            tmpM=matrix(wM[variables,],ind,k,byrow = TRUE)
            # diINDtoPROT[,variables,,1]=tmpD$DET[[variables]]$D*tmpM
            tmp2=tmp2+tmpD$DET[[variables]]$D*tmpM
          }
          ###############
        }else{
          for (variables in (1:vars)){
            tmpM=matrix(wM[variables,],ind,k,byrow = TRUE)
            #  diINDtoPROT[,variables,,1]=tmpD$DET[[variables]]$DM*tmpM
            tmp2=tmp2+tmpD$DET[[variables]]$DM*tmpM
            tmpM=matrix(wV[variables,],ind,k,byrow = TRUE)
            #  diINDtoPROT[,variables,,2]=tmpD$DET[[variables]]$DV*tmpM
            tmp2=tmp2+tmpD$DET[[variables]]$DV*tmpM
          }
          
        }
      }
      if (Wfix==TRUE){
        
        ####################
        if (schema==1|schema==3){#one weight for one variable
          for (variables in (1:vars)){
            #  tmpM=matrix(wM[variables,],ind,k,byrow = TRUE)
            # diINDtoPROT[,variables,,1]=tmpD$DET[[variables]]$D
            tmp2=tmp2+tmpD$DET[[variables]]$D
          }
        }else{
          for (variables in (1:vars)){
            #  diINDtoPROT[,variables,,1]=tmpD$DET[[variables]]$DM
            tmp2=tmp2+tmpD$DET[[variables]]$DM
            # diINDtoPROT[,variables,,2]=tmpD$DET[[variables]]$DV
            tmp2=tmp2+tmpD$DET[[variables]]$DV
          }
        }
      }
      remove(tmpD)
      
      Dind2Neu=tmp2%*%KT#matrix(Inf,ind,k)
      
      #recompute criterion
      if(ASS_DC){
      IDX=max.col(-Dind2Neu, ties="first")}else{
      IDX=max.col(-tmp2, ties="first")}
      #browser()
      old_NUM=tmpIDX
      tmpIDX=Table2(IDX, k, FALSE)
      #  table(IDX)
      if (t>5){
        if (length(tmpIDX)<(length(old_NUM)-10)){
          print(t)
          print(old_NUM)
          print(tmpIDX)
          
        }
        
      }
      OldCrit=GenCrit
      GenCrit=0
      WSSQ=0
      for (individual in 1:ind){
        GenCrit=GenCrit+sum(Dind2Neu[individual, IDX[individual]])
        WSSQ=WSSQ+sum(tmp2[individual, IDX[individual]])
      }
      if(GenCrit-OldCrit>0){print(paste(t,GenCrit))}
      if (abs(OldCrit-GenCrit)<1e-50) {
        fts=ft+1
        
        }
      #print(paste(t,OldCrit-GenCrit))
      #  browser()
      if (verbose){cat(paste(t,GenCrit, "\n", sep="---->"))}
      if (is.na(GenCrit)){
        
        print('A GREAT PROBLEM HERE!')
      }
      if (verbose) {print(print_map(MAP,IDX))}
      
    }
    #crisp assignment
    
    cardinality=Table2(IDX, k, FALSE)
   
    names(IDX)=get.MatH.rownames(x)
    colnames(lambdas)=get.MatH.rownames(proto)
    rownames(lambdas)=paste(c("M","D"),rep(get.MatH.varnames(x),each=2),sep=".")
    
    if (length(cardinality)>atleast){
      if(ASS_DC){
        D2P=Dind2Neu
      }else{
        D2P=tmp2
      }
      tmpsol=list(solution=list(MAP=MAP, IDX=IDX,cardinality=cardinality,proto=proto,
                                Crit=GenCrit,weights.comp=lambdas,
                                Weight.sys=weight.sys,within=WSSQ, D2P=D2P))
      if (TOE){
        toe=WH_SOM_Topo_error(tmpsol,toroidal)
        toes=c(toes,toe)
      }
      solutions=c(solutions,tmpsol)
      criteria=c(criteria,GenCrit)
    }else{
      
      repet=repet-1
    }
    
  }
  
  
  if (TOE){
    best.solution=list(solution=solutions[[which.min(toes)]],
                       repetitions=which.min(toes),
                       quality=1-min(criteria[which.min(toes)])/TOTSSQ,dmax=dmax,TMAX=TMAX,Tmin=Tmin,DATA=x,TOE=min(toes))
  }else{
  best.solution=list(solution=solutions[[which.min(criteria)]],
                     repetitions=which.min(criteria),
                     quality=1-min(criteria)/TOTSSQ,dmax=dmax,TMAX=TMAX,Tmin=Tmin,DATA=x)
  }
  FM=best.solution$solution$MAP
  FIDX=best.solution$solution$IDX
  # if (verbose){print_map(FM,FIDX,gr=TRUE)}
  print(print_map(FM,FIDX,gr=TRUE))
  
  memb=matrix(0,ind,k)
  memb[(best.solution$solution$IDX-1)*ind+c(1:ind)]=1
  resTOTSQ=c_WH_ADPT_KMEANS_TOTALSSQ2(x,memb,m=1,
                                      lambdas=best.solution$solution$weights.comp,
                                      proto=best.solution$solution$proto)
  TTSQ2=sum(resTOTSQ$TSQ2_m+resTOTSQ$TSQ2_v)
  #compute Within
  
  #browser()
  best.solution$quality2=1-best.solution$solution$within/TTSQ2
  #    1-min(criteria)/TTSQ2
  return(best.solution)
}

# Batch SOM of HD -----

WH_2d_Kohonen_maps_test =function (x,net=list(xdim=4,ydim=3,topo=c('rectangular')),
                                   kern.param=2, TMAX=2, Tmin=0.2, 
                                   niter=30,repetitions=5 ,      
                                   simplify=FALSE,
                                   qua=10,
                                   standardize=FALSE, verbose=FALSE, 
                                   init="random", coords=0,ASS_DC=FALSE, bubble=FALSE, TOE=FALSE,toroidal=FALSE){
  SOL=WH_2d_Adaptive_Kohonen_maps_test(x,net, kern.param, TMAX, Tmin, 
                                       niter,repetitions ,
                                       simplify,
                                       qua,
                                       standardize, schema=1,
                                       init.weights='EQUAL',
                                       weight.sys='PROD',theta=2,
                                       Wfix=TRUE,verbose=verbose, init=init, 
                                       coords=coords,ASS_DC=ASS_DC,bubble=bubble,TOE=TOE, toroidal=toroidal)
  
  return(SOL)
}

print_map=function(MAP,IDX,gr=FALSE){
  #pp=Table2(IDX, k, FALSE)
  pp=table(IDX)
  print(pp)
  indici=as.numeric(names(pp))
  mappa=MAP$pts
  countv=rep(0,nrow(mappa))
  for (h in 1:length(indici)){
    countv[indici[h]]=pp[h]
  }
  mappa=cbind(mappa,count=countv)
  mappa[,1]=as.factor(mappa[,1])
  mappa[,2]=as.factor(mappa[,2])
  MAPPA=matrix(0,nrow=max(mappa[,1]),ncol=max(mappa[,2]))
  for (i in 1:nrow(mappa)){
    MAPPA[mappa[i,1],mappa[i,2]]=mappa[i,3]
  }
  if (gr){
    plot(MAP)
    text(MAP$pts[,1],MAP$pts[,2],as.character(c(1:nrow(MAP$pts))),pos=2,cex=0.7)
    text(MAP$pts[which(mappa[,3]>0),1],MAP$pts[which(mappa[,3]>0),2],as.character(mappa[which(mappa[,3]>0),3]),pos=1,cex=0.9,col='BLUE')
  }
  return(MAPPA)
}

Prepare =function (x, simplify, qua, standardize){
  resu=c_Prepare(x,simplify,qua,standardize)
  return(resu)
}

Table2=function(x,N, full=TRUE){
  howm=rep(0,N)
  names(howm)=paste0(c(1:N))
  for (i in 1:N){
    howm[i]=length(which(x==i))
  }
  if (full){
    return(howm)}else{
      return(howm[which(howm>0)])
    }
}

#new distance
newD=function(o1,o2){
  dc=(o1@m-o2@m)^2
  cp=c_dotpW(o1,o2)
  x2=o1@s^2+o1@m^2
  y2=o2@s^2+o2@m^2
  df=x2+y2-2*cp
  dv=df-dc
  return(list(D=df,
              Dm=dc,
              Dv=dv, 
              Ds=(o1@s-o2@s)^2,
              Dsh=df-dc-dv,
              rho=(cp-o1@m*o2@m)/(o1@s*o2@s)))
  
}


initPCA=function(data,MAP,m=0.1,M=0.9){
  grid=MAP$pts
  resu=WH.MultiplePCA(data,c(1:get.MatH.ncols(data)),outl = 0.01)
  PCA_coord=resu$global.pca$ind$coord[,1:2]
  minx=quantile(PCA_coord[,1],m)
  maxx=quantile(PCA_coord[,1],M)
  miny=quantile(PCA_coord[,2],m)
  maxy=quantile(PCA_coord[,2],M)
  grid[,1]=(grid[,1]-min(grid[,1]))/(max(grid[,1])-min(grid[,1]))
  grid[,1]=grid[,1]*(maxx-minx)+minx
  
  grid[,2]=(grid[,2]-min(grid[,2]))/(max(grid[,2])-min(grid[,2]))
  grid[,2]=grid[,2]*(maxy-miny)+miny
  #compute distances
  dista=matrix(0,nrow(grid),nrow(PCA_coord))
  for (i in 1:nrow(grid)){
    dista[i,]=sqrt((PCA_coord[,1]-grid[i,1])^2+(PCA_coord[,2]-grid[i,2])^2)
  }
  return(max.col(-dista))
}

initPCA2=function(PCA_coord,MAP,m=0.1,M=0.9){
  grid=MAP$pts
  #resu=WH.MultiplePCA(data,c(1:get.MatH.ncols(data)),outl = 0.01)
  #PCA_coord=resu$global.pca$ind$coord[,1:2]
  minx=quantile(PCA_coord[,1],m)
  maxx=quantile(PCA_coord[,1],M)
  miny=quantile(PCA_coord[,2],m)
  maxy=quantile(PCA_coord[,2],M)
  grid[,1]=(grid[,1]-min(grid[,1]))/(max(grid[,1])-min(grid[,1]))
  grid[,1]=grid[,1]*(maxx-minx)+minx
  
  grid[,2]=(grid[,2]-min(grid[,2]))/(max(grid[,2])-min(grid[,2]))
  grid[,2]=grid[,2]*(maxy-miny)+miny
  #compute distances
  
 
  dista=matrix(0,nrow(grid),nrow(PCA_coord))
  for (i in 1:nrow(grid)){
    dista[i,]=sqrt((PCA_coord[,1]-grid[i,1])^2+(PCA_coord[,2]-grid[i,2])^2)
  }
  return(max.col(-dista))
}

WH_D_to_prot_ADA=function(DATA,proto,WEIGHTS){
  a=Sys.time()
  tmp=Prepare(DATA,simplify=FALSE,qua=20,standardize=FALSE)
  MM=tmp$MM
  
  remove(tmp)
  tmp2=Prepare(proto,simplify=FALSE,qua=20,standardize=FALSE)
  MM2=tmp2$MM
  remove(tmp2)
  
  
  N=get.MatH.nrows(DATA)
  K=get.MatH.nrows(proto)
  V=get.MatH.ncols(DATA)
  Wm=WEIGHTS[c(1:V)*2-1,]
  Wv=WEIGHTS[c(1:V)*2,]
  DM=DV=D=matrix(0,N,K)
  for (v in 1:V){
    for (k in 1:K){
      rr=nrow(MM[[v]])
      C1=(MM[[v]][1:(rr-1),c(1:N)]+MM[[v]][2:rr,c(1:N)])*0.5
      C2=(MM2[[v]][1:(rr-1),k]+MM2[[v]][2:(rr),k])*0.5
      R1=(-MM[[v]][1:(rr-1),c(1:N)]+MM[[v]][2:(rr),c(1:N)])*0.5
      R2=(-MM2[[v]][1:(rr-1),k]+MM2[[v]][2:(rr),k])*0.5
      p=diff(MM[[v]][,(N+1)])
      p=matrix(p,nrow=1,ncol=(rr-1))
      CC2=matrix(C2,nrow = rr-1,ncol = N,byrow=FALSE)
      RR2=matrix(R2,nrow = rr-1,ncol = N,byrow=FALSE)
      
      m1=p%*%C1
      m2=p%*%CC2
      
      dd=p%*%((C1-CC2)^2+1/3*(R1-RR2)^2)
      dm=(m1-m2)^2
      dv=(dd-dm)
      dmp=Wm[v,k]*dm
      dvp=Wv[v,k]*dv
      ddp=dmp+dvp
      D[,k]=D[,k]+ddp
      DM[,k]=DM[,k]+dmp
      DV[,k]=DV[,k]+dvp
    }
  }
 # print(Sys.time()-a)
   resu=list(D=D,DM=DM,DV=DV)
  return(resu)
}

WH_D_to_prot_DECO=function(DATA,proto){
  tmp=Prepare(DATA,simplify=FALSE,qua=20,standardize=FALSE)
  MM=tmp$MM
  tmpd=c_ComputeFast_L2_SQ_WASS_DMAT(MM,proto)
  
  return(tmpd)
}

scroll_and_print_map=function(MAP,IDX,gr=FALSE,L=0,U=1){
  if(U>MAP$ydim){U=U%%MAP$ydim}
  #scroll step
  #shift di 1 verso sopra
 if(U>0){
  old=c(1:(MAP$xdim*MAP$ydim))
  new=old-U*MAP$xdim
  new[new<1]=MAP$xdim*(MAP$ydim-U)+c(1:(U*MAP$xdim))

  for (i in 1:length(IDX)){
    IDX[i]=new[IDX[i]]
  }
 }
  old=c(1:(MAP$xdim*MAP$ydim))
  new=numeric()
  for (i in 0:(MAP$ydim-1)){
    new=c(new,(c(1:MAP$xdim)-1+L)%%MAP$xdim+1+i*MAP$xdim)
  }
  
  for (i in 1:length(IDX)){
    IDX[i]=new[IDX[i]]
  }
  #print the map (or return new indexes)
  print_map(MAP,IDX,gr=TRUE)
#browser()
  #   #pp=Table2(IDX, k, FALSE)
  # pp=table(IDX)
  # print(pp)
  # indici=as.numeric(names(pp))
  # mappa=MAP$pts
  # countv=rep(0,nrow(mappa))
  # for (h in 1:length(indici)){
  #   countv[indici[h]]=pp[h]
  # }
  # mappa=cbind(mappa,count=countv)
  # mappa[,1]=as.factor(mappa[,1])
  # mappa[,2]=as.factor(mappa[,2])
  # MAPPA=matrix(0,nrow=max(mappa[,1]),ncol=max(mappa[,2]))
  # for (i in 1:nrow(mappa)){
  #   MAPPA[mappa[i,1],mappa[i,2]]=mappa[i,3]
  # }
  # if (gr){
  #   plot(MAP)
  #   text(MAP$pts[,1],MAP$pts[,2],as.character(c(1:nrow(MAP$pts))),pos=2,cex=0.7)
  #   text(MAP$pts[which(mappa[,3]>0),1],MAP$pts[which(mappa[,3]>0),2],as.character(mappa[which(mappa[,3]>0),3]),pos=1,cex=0.9,col='BLUE')
  # }
  # return(MAPPA)
  return(list(MAP=MAP, IDX=IDX))
}


