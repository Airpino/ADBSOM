#gini for class
Gini_av_imp=function(OUTcla,INcla){
pippo=table(OUTcla,INcla)
classI=colSums(pippo)
finC=rowSums(pippo)
#finC=finC/sum(finC)
pippo2=(pippo/matrix(finC,nrow(pippo),ncol(pippo),byrow=FALSE))^2
finP=(1-rowSums(pippo2))/((ncol(pippo)-1)/ncol(pippo))
giniI=sum(finP*finC/sum(finC))
return(giniI)
}

