### oea.R file ###
### mutation operators:
#S 都是輸入 parents
#exchange : 隨機選兩個進行swap
exchange=function(s,N=length(s)){ 
  p=sample(1:N,2) # select two positions
  temp=s[p[1]] # swap values
  s[p[1]]=s[p[2]]
  s[p[2]]=temp
  return(s)
}
#insertion : 隨機選一個數字，對另外一個位置進行插入
insertion=function(s,N=length(s),p=NA,i=NA){
  if(is.na(p)) p=sample(1:N,1) # select a position
  I=setdiff(1:N,p) # ALL except p
  if(is.na(i)) i=sample(I,1) # select random place
  if(i>p) i=i+1 # need to produce a change
  I1=which(I<i) # first part
  I2=which(I>=i) # last part
  s=s[c(I[I1],p,I[I2])] # new solution
  return(s)
}
#displacement : 隨機選一串數字，對另外一個位置進行插入
displacement=function(s,N=length(s)){
  p=c(1,N)
  # select random tour different than s
  while(p[1]==1&&p[2]==N) p=sort(sample(1:N,2))
  I=setdiff(1:N,p[1]:p[2]) # ALL except p
  i=sample(I,1) # select random place
  I1=which(I<i) # first part
  I2=which(I>=i) # last part
  s=s[c(I[I1],p[1]:p[2],I[I2])]
  return(s)
}

### crossover operators:
# partially matched crossover (PMX) operator:
# m is a matrix with 2 parent x ordered solutions
# m 都是輸入兩個 parents 的矩陣 (n*2)
# pmx 是交換中間的，並比較中間的連結關係做衝突處理
pmx=function(m){
  N=ncol(m)
  p=sample(1:N,2) # two cutting points
  c=m # children
  for(i in p[1]:p[2]){
    # rearrange:
    c[1,which(c[1,]==m[2,i])]=c[1,i]
    # crossed section:
    c[1,i]=m[2,i]
    # rearrange:
    c[2,which(c[2,]==m[1,i])]=c[2,i]
    # crossed section:
    c[2,i]=m[1,i]
  }
  return(c)
}
# order crossover (OX) operator:
# m is a matrix with 2 parent x ordered solutions
# ox 是保留中間的，並運用另一個parents的其他數字按照順序填入
ox=function(m){
  N=ncol(m)
  p=sort(sample(1:N,2)) # two cutting points
  c=matrix(rep(NA,N*2),ncol=N)
  # keep selected section:
  c[,p[1]:p[2]]=m[,p[1]:p[2]]
  # rotate after cut 2 (p[2]):
  I=((p[2]+1):(p[2]+N))
  I=ifelse(I<=N,I,I-N)
  a=m[,I]
  # fill remaining genes:
  a1=setdiff(a[2,],c[1,p[1]:p[2]])
  a2=setdiff(a[1,],c[2,p[1]:p[2]])
  I2=setdiff(I,p[1]:p[2])
  c[,I2]=rbind(a1,a2)
  return(c)
}



### order (representation) evolutionary algorithm:
# adapted version of rbga.bin that works with ordered vectors,
# accepts used defined mutation and crossover operators and
# accepts a Lamarckian evolution if evalFunc returns a list
# note: assumes solution with values from the range 1,2,...,size

# Order Evolutionary Algorithm 序列基因演化法
oea=function(size=10,suggestions=NULL,popSize=200,iters=100,
             mutationChance=NA,
             elitism=NA,evalFunc=NULL,
             crossfunc=NULL,mutfunc=mutfunc,REPORT=0){
  if(is.na(mutationChance)) { mutationChance=0.5 }
  if(is.na(elitism)) { elitism=floor(popSize/5)}
  # population initialization: #母體建立
  population=matrix(nrow=popSize,ncol=size)
  #print(population)
  if(!is.null(suggestions)){
    suggestionCount=dim(suggestions)[1]
    for(i in 1:suggestionCount){
      population[i, ] = suggestions[i, ]
    }
    I=(suggestionCount+1):popSize ### new code
  }else {
    I=1:popSize ### new code
  }
  #test
  #print(population)
  # print(I)
  # print(size)
  for(child in I){
    ### new code
    #population 每列 從 1:維度D抽D個
    population[child,]=sample(1:size,size) ### new code
  } 
  #print(population)
  # evaluate population:
  evalVals = rep(NA, popSize)
  # main GA cycle:
  for(iter in 1:iters){
    # evaluate population
    for(object in 1:popSize){
      ### new code
      #print(population[object,])
      EF = evalFunc(population[object,]) #把194個城市丟進去eval func求解
      if(is.list(EF)){ # Lamarckian change of solution
        population[object,]=EF$solution
        evalVals[object] = EF$eval #將求出的解提出來
      }else {
        evalVals[object]=EF
      }
      ### end of new code
    }
    sortedEvaluations=sort(evalVals,index=TRUE)
    if(REPORT>0 && (iter%%REPORT==0||iter==1)){
      cat(iter,"best:",sortedEvaluations$x[1],"mean:",mean(sortedEvaluations$x),"\n") 
    }
    sortedPopulation=matrix(population[sortedEvaluations$ix,],ncol=size)
    
    # check elitism: #挑選出elitism的解
    newPopulation=matrix(nrow=popSize,ncol=size)
    if(elitism>0) # applying elitism:
      newPopulation[1:elitism,]=sortedPopulation[1:elitism,]
    ### very new code inserted here : ###
    # roulette wheel selection of remaining individuals
    others=popSize-elitism
    prob=(max(sortedEvaluations$x)-sortedEvaluations$x+1)
    prob=prob/sum(prob) # such that sum(prob)==1
    # crossover with half of the population (if !is.null)
    if(!is.null(crossfunc)) # need to crossover
      half=round(others/2)
    else half=0 # no crossover
    
    #雜交演化 main function
    if(!is.null(crossfunc)){
      for(child in seq(1,half,by=2)){
        selIDs=sample(1:popSize,2,prob=prob)
        parents=sortedPopulation[selIDs, ]
        if(child==half)
          newPopulation[elitism+child,]=crossfunc(parents)[1,] # 1st child
        else
          newPopulation[elitism+child:(child+1),]=crossfunc(parents) # two children
      }
    }
    # 剩下的做突變 
    # mutation with remaining individuals
    for(child in (half+1):others){
      selIDs=sample(1:popSize,1,prob=prob)
      newPopulation[elitism+child,]=mutfunc(sortedPopulation[selIDs,])
    }
    ### end of very new code ###
    #print(population)
    population=newPopulation # store new population
  } # end of GA main cycle
  result=list(type="ordered chromosome",size=size,
              popSize=popSize, iters=iters,population=population,
              elitism=elitism, mutationChance=mutationChance,
              evaluations=evalVals)
  return(result)
}