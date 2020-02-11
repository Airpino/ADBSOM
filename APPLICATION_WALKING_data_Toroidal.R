#init
library(HistDAWass)
library(Rcpp)
library(RcppArmadillo)
library(kohonen)
library(ggplot2)

Rcpp::sourceCpp('CPP.cpp')
source('Kohonen_maps_fun.R')
source('Activity_recognition_data.R')

########THE DATA #########################
Work_data=data
############################################

####### The settings
filename="WALK_DC_PROD_toro"
ASS_DC= TRUE


MAP=c(7,7)
MAP2=c(10,10)

#bubble=TRUE
bubble=FALSE
TMAX=-2#3.71#2
TMAX2=-3#3.71#2
Tmin=-0.3#0.2
Tmin2=0.35#0.2
niter = 200#50
init="random"
rep=5
TOE=FALSE
toroidal=TRUE

MAPS=c(rep(paste0(MAP[1],"x",MAP[2]),5),rep(paste0(MAP2[1],"x",MAP2[2]),5))
MET=rep(c("Cl.+st.","P1","P2","P3","P4"),2)
# application
sils=numeric()
sils2=numeric()
silsC=numeric()
silsC2=numeric()
qua=numeric()
topoer=numeric()
QE=numeric()


Maps=c(rep(1,5),rep(2,5))
Met=rep(c("st",rep("ad",4)),2)
schema=rep(c(0,1,5,3,6),2)
WSYS=rep(c("NN",rep("PROD",4)),2)

###########################################
WALK=list()
tries=0
#init="PCA"
if (init=="PCA") {
  resu=WH.MultiplePCA(Work_data,c(1:get.MatH.ncols(Work_data)),outl = 0.01)
  coords=resu$global.pca$ind$coord[,1:2]
}else{coords=0}
ACT_REC_COORDS=coords


######### the workhorse
for (exec in 1:length(Maps)){
  if (Maps[exec]==1){
    MAP_par=MAP
    TMAX_par=TMAX
    Tmin_par=Tmin
  
    }else{
    MAP_par=MAP2
    TMAX_par=TMAX2
    Tmin_par=Tmin2
  }
  if (Met[exec]=="st"){
    a=Sys.time()
    set.seed(12345)#classic stand
    results2=WH_2d_Kohonen_maps_test(Work_data,
                                     net=list(xdim=MAP_par[1],
                                              ydim=MAP_par[2],topo=c('hexagonal')),
                                     simplify = TRUE,repetitions = rep,qua = 10,
                                     standardize = TRUE,
                                     niter = niter,
                                     TMAX = TMAX_par,
                                     Tmin = Tmin_par,
                                     init=init, 
                                     coords = ACT_REC_COORDS,
                                     ASS_DC=ASS_DC,
                                     bubble = bubble, TOE=TOE, toroidal=toroidal)
    
  }else{
    a=Sys.time()
    set.seed(12345)
    results2=WH_2d_Adaptive_Kohonen_maps_test(Work_data,
                                              net=list(xdim=MAP_par[1],
                                                       ydim=MAP_par[2],topo=c('hexagonal')),
                                              TMAX = TMAX_par,
                                              Tmin = Tmin_par, 
                                              niter = niter,repetitions = rep ,
                                              simplify=TRUE,
                                              qua=20,
                                              schema=schema[exec],
                                              weight.sys=WSYS[exec], 
                                              init=init, 
                                              coords = ACT_REC_COORDS,ASS_DC=ASS_DC,
                                              bubble = bubble, TOE=TOE, toroidal = toroidal)
    
  }
  print(paste(exec,Maps[exec],schema[exec],WSYS[exec],sep="-")) ; 
  WALK[[exec]]=results2
  source('validation.R')
  remove(results2)
  b=Sys.time()
  print(difftime(b,a))
  
}
########################

DFRES_WALK=data.frame(MAP=MAPS,MET=MET,silh=sils,silh2=sils2,sil_CAM=silsC,sil_CAM_2=silsC2,
                          QPI=qua,top_err=topoer)
print(DFRES_WALK,digits = 4)
save.image(paste0(filename,".RData"))

