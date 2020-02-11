#activity recognition
require(HistDAWass)
load("means.RData")
############ SELECTOR FUNCTION
SelACT_and_persons=function(sel_activity,sel_person,activity,person){
  selected1=numeric()
  selected2=numeric()
  for(i in 1:length(activity)){
    els=which(activity==sel_activity[i])
    selected1=c(selected1,els)
  }
  for (j in 1:length(person)){
    els=which(person==as.character(sel_person[j]))
    selected2=c(selected2,els)
  }
  selected=intersect(selected1,selected2)
  selected=sort(selected)
  return(selected)
}
######################################
#do labels
#This is for plotting
ROW_labels=get.MatH.rownames(A);
person=rep("a",length(ROW_labels))
person[grep("p1",ROW_labels)]="1"#"F 21.8bmi 25y 170h 63k "
person[grep("p2",ROW_labels)]="2"#"F 20.6bmi 20y 162h 54k "
person[grep("p3",ROW_labels)]="3"#"M 22.8bmi 30y 185h 78k "
person[grep("p4",ROW_labels)]="4"#"M 23.6bmi 25y 182h 78k "
person[grep("p5",ROW_labels)]="5"#"M 23.0bmi 26y 183h 77k "
person[grep("p6",ROW_labels)]="6"#"F 18.4bmi 23y 165h 50k "
person[grep("p7",ROW_labels)]="7"#"F 20.4bmi 21y 167h 57k "
person[grep("p8",ROW_labels)]="8"#"M 24.5bmi 24y 175h 75k "

activity=rep(0,length(ROW_labels))
group.activity=rep("zz",length(ROW_labels))
group.sex=rep("M",length(ROW_labels))
group.sex[which((person==1)|(person==2)|(person==6)|(person==7))]="F"
for (act in 1:19){
  if (act<10){
    pos=grep(paste0("a0",act),ROW_labels)
    activity[pos]=act
    
  }else{
    pos=grep(paste0("a",act),ROW_labels)
    activity[pos]=act
    
  }
  
}
#groups of activities
# sitting (A1), 
# standing (A2), 
# lying on back and on right side (A3 and A4), 
# ascending and descending stairs (A5 and A6), 
# standing in an elevator still (A7) 
# and moving around in an elevator (A8), 
# walking in a parking lot (A9), 
# walking on a treadmill with a speed of 4 km/h (in flat and 15 deg inclined positions) (A1 
#                                                                                        0 and A11), 
# running on a treadmill with a speed of 8 km/h (A12), 
# exercising on a stepper (A13), 
# exercising on a cross trainer (A14), 
# cycling on an exercise bike in horizontal and vertical positions (A15 and A16), 
# rowing (A17),
# jumping (A18), 
# and playing basketball (A19).

group.activity[which(activity==1)]="01_sitting"
group.activity[which(activity==2)]="02_standing"
group.activity[which(activity==3)]="03_lying on Back"
group.activity[which(activity==4)]="04_lying on right"
group.activity[which(activity==5)]="05_Ascending stairs"
group.activity[which(activity==6)]="06_Descending stairs"
group.activity[which(activity==7)]="07_Stan_in_elev_"
group.activity[which(activity==8)]="08_Moving in elevator"
group.activity[which(activity==9)]="09_walk in park.lot"
group.activity[which(activity==10)]="10_walk 4kmH flat"
group.activity[which(activity==11)]="11_walk 4kmH 15 deg"
group.activity[which(activity==12)]="12_runn 8kmH"
group.activity[which(activity==13)]="13_stepper"
group.activity[which(activity==14)]="14_Cross trainer"
group.activity[which(activity==15)]="15_Cycl hor pos"
group.activity[which(activity==16)]="16_Cycl ver pos"
group.activity[which(activity==17)]="17_Rowing"
group.activity[which(activity==18)]="18_Jumping"
group.activity[which(activity==19)]="19_Play basket"
###################################################################



s_act=c(9)
allact=FALSE
s_person =c(1:8)


selection=SelACT_and_persons(sel_activity = s_act,
                             sel_person = s_person,
                             activity,person)
v_ACT=activity[selection]
v_PER=person[selection]
GR_act_sel=group.activity[selection]




OnlyACC=c(1:3,10:12,19:21,28:30,37:39)
OnlyGYR=c(4:6,13:15,22:24,31:33,40:42)
OnlyMAG=c(7:9,16:18,25:27,34:36,43:45)
list_of_vars=c(1:45)
list_of_vars=sort(c(OnlyACC,OnlyGYR))
# list_of_vars=OnlyMAG

data=A[selection,list_of_vars]

#source('C:/Users/Antonio/Desktop/HISTO_R_PACK/SOM_appl/Kohonen_maps_fun.R')
#Rcpp::sourceCpp('CPP.cpp')
# [1] "TO_xacc" "TO_yacc" "TO_zacc" "TO_xgyr" "TO_ygyr" "TO_zgyr" "TO_xmag" "TO_ymag" "TO_zmag"
# [10] "RA_xacc" "RA_yacc" "RA_zacc" "RA_xgyr" "RA_ygyr" "RA_zgyr" "RA_xmag" "RA_ymag" "RA_zmag"
# [19] "LA_xacc" "LA_yacc" "LA_zacc" "LA_xgyr" "LA_ygyr" "LA_zgyr" "LA_xmag" "LA_ymag" "LA_zmag"
# [28] "RL_xacc" "RL_yacc" "RL_zacc" "RL_xgyr" "RL_ygyr" "RL_zgyr" "RL_xmag" "RL_ymag" "RL_zmag"
# [37] "LL_xacc" "LL_yacc" "LL_zacc" "LL_xgyr" "LL_ygyr" "LL_zgyr" "LL_xmag" "LL_ymag" "LL_zmag"

# #select activity and select persons

# group.activity[which(activity==1)]="01_sitting"
# group.activity[which(activity==2)]="02_standing"
# group.activity[which(activity==3)]="03_lying on Back"
# group.activity[which(activity==4)]="04_lying on right"
# group.activity[which(activity==5)]="05_Ascending stairs"
# group.activity[which(activity==6)]="06_Descending stairs"
# group.activity[which(activity==7)]="07_Stan. in elev."
# group.activity[which(activity==8)]="08_Moving in elevator"
# group.activity[which(activity==9)]="09_walk in park.lot"
# group.activity[which(activity==10)]="10_walk 4km/h flat"
# group.activity[which(activity==11)]="11_walk 4km/h 15 deg"
# group.activity[which(activity==12)]="12_runn 8km/h"
# group.activity[which(activity==13)]="13_stepper"
# group.activity[which(activity==14)]="14_Cross trainer"
# group.activity[which(activity==15)]="15_Cycl hor pos"
# group.activity[which(activity==16)]="16_Cycl ver pos"
# group.activity[which(activity==17)]="17_Rowing"
# group.activity[which(activity==18)]="18_Jumping"
# group.activity[which(activity==19)]="19_Play basket"

## PERSONS
# person[grep("p1",ROW_labels)]="1"#"F 21.8bmi 25y 170h 63k "
# person[grep("p2",ROW_labels)]="2"#"F 20.6bmi 20y 162h 54k "
# person[grep("p3",ROW_labels)]="3"#"M 22.8bmi 30y 185h 78k "
# person[grep("p4",ROW_labels)]="4"#"M 23.6bmi 25y 182h 78k "
# person[grep("p5",ROW_labels)]="5"#"M 23.0bmi 26y 183h 77k "
# person[grep("p6",ROW_labels)]="6"#"F 18.4bmi 23y 165h 50k "
# person[grep("p7",ROW_labels)]="7"#"F 20.4bmi 21y 167h 57k "
# person[grep("p8",ROW_labels)]="8"#"M 24.5bmi 24y 175h 75k "