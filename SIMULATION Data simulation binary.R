library(xlsx)
library(MASS)#mvrnorm
library(matrixcalc)#ispositivedefinite


###simulation of binary data###
for (x in 0:99) {
  set.seed(999+x)
  #simulate the parameters of the error structure
  r12_r=rnorm(1,mean=0.20,sd=0.01)
  r13_r=rnorm(1,mean=0.30,sd=0.01)
  r23_r=rnorm(1,mean=0.20,sd=0.01)
  corrmat=matrix(c(1,r12_r,r13_r,r12_r,1,r23_r,r13_r,r23_r,1),nrow=3)
  
  #simulate the regression coefficients
  sex_r=rnorm(3,mean=0,sd=0.1)
  age_r=rnorm(3,mean=0,sd=0.01)
  r2_r=rnorm(3,mean=0,sd=0.01)
  r3_r=rnorm(3,mean=0,sd=0.01)
  r4_r=rnorm(3,mean=0,sd=0.01)
  r5_r=rnorm(3,mean=0,sd=0.01)
  r6_r=rnorm(3,mean=0,sd=0.01)
  
  #simulate the latent variables
  latent = mvrnorm(n = 2000, mu = c(0,0,0), Sigma = corrmat)
  mydata=as.data.frame(latent)
  
  ##add three predictors
  mydata$sex = rbinom(2000,size=1,prob=0.5)
  mydata$age = runif(2000, min = 18, max = 80)
  mydata$region = round(runif(2000, min = 0.5, max = 5.5))
  
  #add effect sex,age and region to the latent variable
  head(mydata)
  mydata$ys1=mydata[,1]+sex_r[1]*mydata$sex+age_r[1]*mydata$age+
    r2_r[1]*(mydata$region==2)+
    r3_r[1]*(mydata$region==3)+
    r4_r[1]*(mydata$region==4)+
    r5_r[1]*(mydata$region==5)+
    r6_r[1]*(mydata$region==6)
  
  mydata$ys2=mydata[,2]+sex_r[2]*mydata$sex+age_r[2]*mydata$age+
    r2_r[2]*(mydata$region==2)+
    r3_r[2]*(mydata$region==3)+
    r4_r[2]*(mydata$region==4)+
    r5_r[2]*(mydata$region==5)+
    r6_r[2]*(mydata$region==6)
  
  mydata$ys3=mydata[,3]+sex_r[3]*mydata$sex+age_r[3]*mydata$age+
    r2_r[3]*(mydata$region==2)+
    r3_r[3]*(mydata$region==3)+
    r4_r[3]*(mydata$region==4)+
    r5_r[3]*(mydata$region==5)+
    r6_r[3]*(mydata$region==6)
  
  mydata$y1=(mydata$ys1>mean(mydata$ys1))
  mydata$y2=(mydata$ys2>mean(mydata$ys2))
  mydata$y3=(mydata$ys3>mean(mydata$ys3))
  
  mydata=mydata[,-c(1:3)]
  
  mydata$region_2=(mydata$region==2)
  mydata$region_3=(mydata$region==3)
  mydata$region_4=(mydata$region==4)
  mydata$region_5=(mydata$region==5)
  mydata$y1b=ifelse(mydata$y1=="TRUE",1,0)
  mydata$y2b=ifelse(mydata$y2=="TRUE",1,0)
  mydata$y3b=ifelse(mydata$y3=="TRUE",1,0)
  
  new_x <- paste("simulated_", x, ".xlsx", sep = "")
  write.xlsx(mydata, new_x, row.names = F)
  print(x)
}
