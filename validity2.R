#validity indices
library(matrixStats)

WH_D_to_prot_ADA_Ai=function(DATA,proto,WEIGHTS,ditopro,ID,WSS_M,WSS_V){
  
  N=get.MatH.nrows(DATA)
  K=get.MatH.nrows(proto)
  memb=matrix(0,N,K)
  memb[(ID-1)*N+c(1:N)]=1
  V=get.MatH.ncols(DATA)
  Wm=WEIGHTS[c(1:V)*2-1,]
  Wv=WEIGHTS[c(1:V)*2,]
  
  Nk=colSums(memb)
  Nk[which(Nk<2)]=0
  D1=rowSums(ditopro$D*memb)*rowSums(memb*matrix(Nk,N,K,byrow=TRUE))
  
  WSSm=matrix(0,N,K)
  WSSv=matrix(0,N,K)
  for (v in 1:V){
    WSSm=WSSm+matrix(WW[(v*2-1),]*WSS_M[,v],N,K,byrow=TRUE)
    WSSv=WSSv+matrix(WW[(v*2),]*WSS_V[,v],N,K,byrow=TRUE)
  }
  WSSm=rowSums(WSSm*memb)
  WSSv=rowSums(WSSv*memb)
  
  Nk[which(Nk<2)]=2
  D1=(D1+WSSm+WSSv)/(rowSums(memb*matrix(Nk,N,K,byrow=TRUE))-1)
  # browser()
  return(D1)
}
WH_D_to_prot_ADA_Ai_CAM=function(DATA,proto,ditopro,ID){
  N=get.MatH.nrows(DATA)
  K=get.MatH.nrows(proto)
  memb=matrix(0,N,K)
  memb[(ID-1)*N+c(1:N)]=1
  D1=rowSums(ditopro$D*memb)
  return(D1)
}
WH_D_to_prot_ADA_Bi=function(DATA,proto,WEIGHTS,ID,WSS_M,WSS_V,met=1){
  pippo=WH_D_to_prot_DECO(DATA,proto)
  
  V=get.MatH.ncols(DATA)
  N=get.MatH.nrows(DATA)
  K=get.MatH.nrows(proto)
  B=memb=matrix(0,N,K)
  memb[(ID-1)*N+c(1:N)]=1
  Nk=colSums(memb)
  memb=memb*.Machine$double.xmax
  B[,which(Nk==0)]=.Machine$double.xmax
  B=B+memb
  ID=as.numeric(ID)
  
  for (k in 1:K){
    sel_i=which(ID==k)
    
    if (length(sel_i)>0){
      iter=c(1:K)
      iter=iter[-k]
      for (k2 in iter){
        if(Nk[k2]>0){
          for (v in 1:V){
            # WSS_M_of_i=WSS_M[,k]
            # WSS_V_of_i=WSS_V[,k]
            
            
            WW_M_of_i=WEIGHTS[(v*2-1),k]
            WW_V_of_i=WEIGHTS[(v*2),k]
            WW_M_not_of_i=WEIGHTS[(v*2-1),k2]
            WW_V_not_of_i=WEIGHTS[(v*2),k2]
            WW_M=(WW_M_of_i+WW_M_not_of_i)/2
            WW_V=(WW_V_of_i+WW_V_not_of_i)/2
            if (met==1){
              WW_M=WW_M_not_of_i
              WW_V=WW_V_not_of_i}
            WSS_M_not_of_i=WW_M*WSS_M[k2,v]/Nk[k2]
            WSS_V_not_of_i=WW_V*WSS_V[k2,v]/Nk[k2]
            # if(1%in%sel_i){
            #   print("1 found!")
            # }
            dm=pippo$DET[[v]]$DM[sel_i,k2]*WW_M+WSS_M_not_of_i
            dv=pippo$DET[[v]]$DV[sel_i,k2]*WW_V+WSS_V_not_of_i
            dd=dm+dv
            B[sel_i,k2]=B[sel_i,k2]+dd
            
          }
        }
      }
    }
  }
  
  return(B)
}
# ########################################
# 
# results2=WH_2d_Adaptive_Kohonen_maps_test(Age_Pyramids_2014,
#                                           net=list(xdim=6,ydim=6,topo=c('hexagonal')),
#                                           simplify = TRUE,qua = 20,schema = 6,
#                                           weight.sys = "SUM",theta = 2,
#                                           niter = 50,TMAX = 2,Tmin = 0.1)
# ########################################
data=results2$DATA
#dataM=data@M

mymap=results2$solution$MAP
#Dneuros=as.matrix(dist(mymap$pts))
Dneuros=unit.distances(mymap,toroidal = toroidal)
proto=results2$solution$proto
ID=results2$solution$IDX
Nprot=get.MatH.nrows(proto)
WW=results2$solution$weights.comp
K=Nprot
if(results2$solution$Weight.sys=="SUM"){WW=WW^2}
N=get.MatH.nrows(data)
V=get.MatH.ncols(data)

WSS_M=matrix(0,K,V)
WSS_V=matrix(0,K,V)
proto2=proto
proto2C=proto
for (k in 1:K){
  for (j in 1:V){
    if (length(which(ID==k))>1){
      tmpdata=data[which(ID==k),j]
      proto2@M[k,j][[1]]=WH.vec.mean(tmpdata)
      tmp=as.numeric(WH.SSQ(data[which(ID==k),j]))
      vals=as.numeric(get.MatH.stats(data[which(ID==k),j])$mat)
      WSS_M[k,j]=sum((vals-mean(vals))^2)
      WSS_V[k,j]=tmp-WSS_M[k,j]
    }else{
      if (length(which(ID==k))==1){
        tmpdata=data[which(ID==k),j]
        proto2@M[k,j][[1]]=tmpdata@M[1,1][[1]]
      }
      WSS_M[k,j]=0
      WSS_V[k,j]=0
    }
  }
}
### End of WSS matrices computation
protoOR=proto
proto=proto2
centers_prot=get.MatH.stats(proto)$mat

Proto_C=Center.cell.MatH(proto)

### End of proto decomposition

###### Compute le ai
a1=rep(0,N)
a2=rep(0,N)
a1_CAMP=rep(0,N)

ditopro_sil=WH_D_to_prot_ADA(data,proto,WW)
ditopro_sil_cam=WH_D_to_prot_ADA(data,protoOR,WW)
############# ai
a1=WH_D_to_prot_ADA_Ai(data,proto,WW,ditopro_sil,ID,WSS_M,WSS_V)
a1_CAMP=WH_D_to_prot_ADA_Ai_CAM(data,protoOR,ditopro_sil_cam,ID)
####################
#compute bi
B=WH_D_to_prot_ADA_Bi(data,proto,WW,ID,WSS_M,WSS_V)
memb=matrix(0,N,K)
memb[(ID-1)*N+c(1:N)]=1


Nk=colSums(memb)
memb=memb*.Machine$double.xmax

B_CAMP=ditopro_sil_cam$D+memb
B_CAMP[,which(Nk==0)]=.Machine$double.xmax

b1=rowMins(B)
b1_CAMP=rowMins(B_CAMP)
silhouette.index=mean((b1-a1)/rowMaxs(cbind(b1,a1)))
silhouette.index_CAMP=mean((b1_CAMP-a1_CAMP)/rowMaxs(cbind(b1_CAMP,a1_CAMP)))


## silhouette not adjacent

b2=rep(0,N)
b2_CAMP=rep(0,N)
SILH=rep(0,N)
SILH2=rep(0,N)
SILH_C=rep(0,N)
SILH_C2=rep(0,N)

for (k in 1:K){
  sel_i=which(ID==k)
  if (length(sel_i)>0){
    iter=c(1:K)
    iter=iter[-k]
    for (k2 in iter){
      if(Nk[k2]>0){
        if (Dneuros[k,k2]<1.1){
          B[sel_i,k2]=.Machine$double.xmax
          B_CAMP[sel_i,k2]=.Machine$double.xmax
          }
      }
    }
  }
  
}
b2=rowMins(B)
b2_CAMP=rowMins(B_CAMP)
silhouette.index2=mean((b2-a1)/rowMaxs(cbind(b2,a1)))
silhouette.index2_CAMP=mean((b2_CAMP-a1_CAMP)/rowMaxs(cbind(b2_CAMP,a1_CAMP)))
print(silhouette.index)
print(silhouette.index_CAMP)
print(silhouette.index2)
print(silhouette.index2_CAMP)
