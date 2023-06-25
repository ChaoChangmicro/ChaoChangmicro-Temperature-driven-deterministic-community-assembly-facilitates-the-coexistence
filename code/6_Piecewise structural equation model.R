#分段 SEM，详情 ?psem
#这里，对于每个独立的响应方程，直接使用简单线性回归，确定响应关系
#其它情况，如果已知变量间的某种非线性关系，即可以使用非线性模型
library(piecewiseSEM)
dat <- read.csv("sem.csv",row.names = 1,header = T,sep=",")
keeley_psem <- psem(
  lm(multidiversity ~ cohesion+stochasticity+pH+latitude, data = dat),
  lm(cohesion ~stochasticity+pH+latitude, data = dat),
  lm(stochasticity ~pH+tem+latitude, data = dat),
  lm(pH ~tem+latitude, data = dat),
  lm(tem ~latitude, data = dat),
  data = dat)
summary(keeley_psem)

keeley_psem1 <- psem(
  lm(multidiversity ~ cohesion+pH+latitude, data = dat),
  lm(cohesion ~stochasticity+pH, data = dat),
  lm(stochasticity ~pH+tem, data = dat),
  lm(pH ~latitude, data = dat),
  lm(tem ~latitude, data = dat),
  data = dat)
summary(keeley_psem1)
A <- lm(cohesion ~stochasticity+pH+latitude, data = dat)
summary(A)
keeley_psem2 <- psem(
  lm(multidiversity ~ cohesion+pH+latitude, data = dat),
  lm(cohesion ~stochasticity+pH+latitude, data = dat),
  lm(stochasticity ~pH+tem, data = dat),
  lm(pH ~latitude, data = dat),
  lm(tem ~latitude, data = dat),
  data = dat)
summary(keeley_psem2)
keeley_psem3 <- psem(
  lm(multidiversity ~ cohesion+pH+latitude, data = dat),
  lm(cohesion ~stochasticity+pH, data = dat),
  lm(stochasticity ~pH+tem, data = dat),
  lm(pH ~latitude, data = dat),
  lm(tem ~latitude, data = dat),
  data = dat)
summary(keeley_psem3)
AIC(keeley_psem2,keeley_psem3)



plot(keeley_psem2)

dSep(keeley_psem, .progressBar = FALSE)
fisherC(keeley_psem2)
coefs(keeley_psem)
#coefs(keeley.sem2,intercepts=T)#截距也可以给出
rsquared(keeley_psem)#给出每个内生变量的R2

