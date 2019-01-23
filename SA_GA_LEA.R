### tsp.R file ###

#### import package ####
library(TSP) # load TSP package
library(RCurl) # load RCurl package
source("oea.R") # load ordered evolutionary algorithm

file_list = c("qa194.csv","uy734.csv","lu980.csv","rw1621.csv","mu1979.csv")
MD_list = list()
N_list = c()
for (i in 1:length(file_list)) {
  #### import data ####
  Data=read.table(file_list[i],sep=" ")
  Data=Data[,3:2] # longitude and latitude
  names(Data)=c("x","y") # x and y labels
  N_list[i]=nrow(Data) # number of cities
  
  #### calculate data distance ####
  D=dist(Data,upper=TRUE)
  D[1:length(D)]=round(D[1:length(D)])
  MAXIT=250000
  Methods=c("SANN_qa","SANN_uy","SANN_lu","SANN_rw","SANN_mu",
            "LEA_qa","LEA_uy","LEA_lu","LEA_rw","LEA_mu",
            "GA_qa","GA_uy","GA_lu","GA_rw","GA_mu") # comparison of 3 methods
  RES=matrix(nrow=MAXIT,ncol=length(Methods))
  MD_list[[i]] =as.matrix(D)
}
rm(Data)


#### define eval function ####
# overall distance of a tour (evaluation function):
#輸入#主要分析函數
tour=function(s){ 
  # compute tour length:
  EV<<-EV+1 # increase evaluations
  s=c(s,s[1]) # start city is also end city
  res=0
  #計算總旅程的距離
  for(i in 2:length(s)) res=res+MD[s[i],s[i-1]] #去MD矩陣中找距離
  # store memory with best values:
  if(res<BEST) BEST<<-res
  if(EV<=MAXIT) F[EV]<<-BEST
  # only for hybrid method:
  # return tour
  return(res)
}

## LEA function ##
# 計算移動城市的index長度，用在LEA裡面
mindex=function(i,dir,s=NULL,N=length(s)){
  res=i+dir #positive or negative jump
  if(res<1){
    res=N+res
  }else if(res>N){
    res=res-N
  }
  return(res)
}

#先進行一次local imp tour，用在LEA裡面
local_imp_tour=function(s,p=NA){
  # local search
  N=length(s); ALL=1:N
  if(is.na(p)) p=sample(ALL,1) # select random position
  I=setdiff(ALL,p) #排除掉P
  # current distance: p to neighbors(找下一個要移動的城市)
  pprev=mindex(p,-1,N=N); pnext=mindex(p,1,N=N)
  dpcur=MD[s[pprev],s[p]]+MD[s[p],s[pnext]] #上一個到這個城市加上這個到下一個的距離
  # new distance if p is remove to another position:
  dpnew=MD[s[pprev],s[pnext]] #上一個到下一個的距離
  
  # search for best insertion position for p:
  ibest=0;best=-Inf
  
  for(i in I){ # extra cycle that increases computation
    inext=mindex(i,1,N=N);iprev=mindex(i,-1,N=N)
    if(inext==p) inext=pnext
    if(iprev==p) iprev=pprev
    # dinew: new distance p to neighbors if p inserted:
    # current i distance without p:
    if(i<p) {
      dinew=MD[s[iprev],s[p]]+MD[s[p],s[i]]
      dicur=MD[s[iprev],s[i]]
    }else{
      dinew=MD[s[i],s[p]]+MD[s[p],s[inext]]
      dicur=MD[s[i],s[inext]]
    }
    # difference between current tour and new one:
    dif=(dicur+dpcur)-(dinew+dpnew)
    if(dif>0 && dif>best){ # improved solution
      best=dif
      ibest=i
    }
  }
  if(ibest>0) # insert p in i
    s=insertion(s,p=p,i=ibest)
  return(list(eval=tour(s),solution=s))
}

#### Output ####
# SANN: #模擬退火法
for (r in 1:length(file_list)) {
  cat(paste0("SANN run:",r,"\n"))
  set.seed(9530) # for replicability
  N = N_list[r]
  s=sample(1:N,N) # initial solution
  MD = MD_list[[r]]
  EV=0; BEST=Inf; F=rep(NA,MAXIT) # reset these vars.
  C=list(maxit=MAXIT,temp=500,trace=TRUE,REPORT=MAXIT)
  PTM=proc.time() # start clock
  SANN=optim(s,fn=tour,gr=insertion,method="SANN",control=C)
  sec=(proc.time()-PTM)[3] # get seconds elapsed
  cat("time elapsed:",sec,"\n")
  RES[,r]=F
  cat(paste0("Finish SANN run:",r,"\n"))
}
# LEA : 拉馬克演化法
for (r in 1:length(file_list)) {
  # Lamarckian EA (LEA): (使用的eval function不同)
  cat("LEA run:\n")
  set.seed(911) # for replicability
  EV=0; BEST=Inf; F=rep(NA,MAXIT) # reset these vars.
  pSize=10;iters=ceiling((MAXIT-pSize)/(pSize-1))
  PTM=proc.time() # start clock
  N = N_list[r]
  MD = MD_list[[r]]
  LEA=oea(size=N,popSize=pSize,iters=iters,evalFunc=local_imp_tour,
          crossfunc=ox,mutfunc=displacement,REPORT=iters,elitism=1)
  sec=(proc.time()-PTM)[3] # get seconds elapsed
  cat("time elapsed:",sec,"\n")
  RES[,r+5]=F
}
# GA : 基因演化法
for (r in 1:length(file_list)) {
  cat(paste0("GA run:",r,"\n"))
  N = N_list[r]
  MD = MD_list[[r]]
  set.seed(12345) # for replicability
  EV=0; BEST=Inf; F=rep(NA,MAXIT) # reset these vars.
  pSize=8;iters=ceiling((MAXIT-pSize)/(pSize-1))
  PTM=proc.time() # start clock
  OEA=oea(size=N,popSize=pSize,iters=iters,evalFunc=tour,crossfunc
          =ox,mutfunc=insertion,mutationChance=0.6,REPORT=iters,elitism=1)
  sec=(proc.time()-PTM)[3] # get seconds elapsed
  cat("time elapsed:",sec,"\n")
  cat(paste0("Finish GA run:",r,"\n"))
  RES[,r+10]=F
}


# create PDF with comparison:
pdf("lea-all.pdf",paper="special")
par(mar=c(4.0,4.0,0.1,0.1))
X=seq(1,MAXIT,length.out=200)
ylim=c(min(RES)-50,max(RES))
plot(X,RES[X,6],ylim=ylim,type="l",lty=5,lwd=2,xlab="evaluations
     ",ylab="tour distance")
lines(X,RES[X,7],type="l",lty=4,lwd=2)
lines(X,RES[X,8],type="l",lty=3,lwd=2)
lines(X,RES[X,9],type="l",lty=2,lwd=2)
lines(X,RES[X,10],type="l",lty=1,lwd=2)
legend("topright",Methods[c(6:10)],lwd=2,lty=5:1)
dev.off()
