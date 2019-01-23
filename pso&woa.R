# 載入粒子群聚法的套件

# install.packages("pso")
# install.packages("metaheuristicOpt")
library("pso")
library("metaheuristicOpt")

# 設定環境變數

setwd(dir = "./")
path = c("qa194.csv", "uy734.csv", "lu980.csv", "rw1621.csv", "mu1979.csv")

# 定義目標函數

fitness <- function(X)
{
  total = 0.
  
  # Largest Value Mapping
  
  LVM = rank(X, ties.method= "first")
  
  # 計算解函數值
  
  for(i in 1:length(LVM))
  {
    if(i + 1 < length(LVM))
    {
      p1 = LVM[i + 0]
      p2 = LVM[i + 1]
      total = total + sqrt((data[p1, "x"] - data[p2, "x"])^2 + (data[p1, "y"] - data[p2, "y"])^2)
    }
  }
  return(total)
}

# 結果矩陣宣告

r = 1
result.f = array(dim = c(5, 5, 2))
result.t = array(dim = c(5, 5, 2))

# 正式進行實驗
i = 1
for(i in (1:length(path)))
{
  data = read.table(path[i], header = FALSE, sep = " ")
  data[, 1] = NULL
  colnames(data) <- c("x", "y")
  
  # 用粒子群聚法求解
  
  beg = proc.time()
  solution = psoptim(rep(NA, nrow(data)), 
                     fitness, 
                     lower = -5, 
                     upper = 5, 
                     control = list(maxit = 100, trace = 1, REPORT = 1, s = 200, w = 6, c.p = 4, c.g = 6))
  end = proc.time()
  result.t[i, r, 1] = (end - beg)[3]
  result.f[i, r, 1] = (solution$value)
  
  # 用鯨魚演算法求解
  
  beg = proc.time()
  solution = WOA(FUN = fitness, 
                 optimType = "MIN", 
                 numVar = nrow(data), 
                 rangeVar = matrix(c(-5, 5), nrow = 2),
                 numPopulation = 50,
                 maxIter = 50) 
  end = proc.time()
  result.t[i, r, 2] = (end - beg)[3]
  result.f[i, r, 2] = (fitness(solution))
}

result.f[,,2]
result.t[,,2]
