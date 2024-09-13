library(xlsx)
library(MASS)#mvrnorm
library(matrixcalc)#ispositivedefinite
setwd('C:/Users/u0118563/OneDrive - KU Leuven/Projecten/MV probit/Data/Simulation study')

#categorize function
categorize = function(cut1, cut2,cut3, latent) {
  category = ifelse(latent < cut1, 0, ifelse(latent < cut2, 1, ifelse(latent < cut3,2,3)))
  return(category)
}

###############start simulation#################################################
for (x in 1:100) {
  set.seed(x)
  m = 50
  #sample variance-covariance matrices
  cmat_female = matrix(rWishart(n=1,df=m,Sigma=diag(m)),nrow=m)
  cmat_male = matrix(rWishart(n=1,df=m,Sigma=diag(m)),nrow=m)
  
  #sample the responses
  mean_female = rnorm(n = m, mean = 0.5, sd = 1)
  mean_male = mean_female + rnorm(n = m, mean = 0.1, sd = 0.1)
  
  latent_female = mvrnorm(n = 1000, mu = mean_female, Sigma = cmat_female)
  latent_male = mvrnorm(n = 1000, mu = mean_male, Sigma = cmat_male)
  
  sex = c(rep(0, 1000), rep(1, 1000))
  mydata = cbind(rbind(latent_female, latent_male), sex)
  mydata = as.data.frame(mydata)
  
  ##add two extra predictors
  mydata$age = runif(2000, min = 18, max = 80)
  mydata$region = round(runif(2000, min = 0.5, max = 5.5))
  
  #add effect age and region to the latent variable
  mydata[, 1:m] = mydata[, 1:m] + rnorm(1, mean = 0, sd = 0.01) * mydata$age
  mydata[, 1:m] = mydata[, 1:m] + rnorm(1, mean = 0, sd = 0.1) * (mydata$region ==
                                                                    2)
  mydata[, 1:m] = mydata[, 1:m] + rnorm(1, mean = 0, sd = 0.1) * (mydata$region ==
                                                                    3)
  mydata[, 1:m] = mydata[, 1:m] + rnorm(1, mean = 0, sd = 0.1) * (mydata$region ==
                                                                    4)
  mydata[, 1:m] = mydata[, 1:m] + rnorm(1, mean = 0, sd = 0.1) * (mydata$region ==
                                                                    5)
  
  #select thresholds
  means = unlist(lapply(mydata[, 1:m], mean))
  var = unlist(lapply(mydata[, 1:m], var))
  threshold_1 = means - runif(n = m, min = 0.76, max = 1.5)*sqrt(var)
  threshold_2 = means + runif(n = m, min = -0.75, max = 0.75)*sqrt(var)
  threshold_3 = means + runif(n = m, min = 0.76, max = 1.5)*sqrt(var)
  
  i = 1
  j = 1
  for (i in 1:m) {
    cat = c()
    for (j in 1:2000) {
      cat = c(cat, categorize(threshold_1[i], threshold_2[i], threshold_3[i], mydata[j, i]))
    }
    mydata[, m + 3 + i] = cat
  }
  names(mydata) = c(paste('z_', 1:m, sep = ''),
                    'sex',
                    'age',
                    'region',
                    paste('y_', 1:m, sep = ''))
  mydata$sex = as.factor(mydata$sex)
  mydata$region = as.factor(mydata$region)
  mydata=mydata[,-c(1:m)]
  new_x <- paste("simulated_", x, ".xlsx", sep = "")
  write.xlsx(mydata, new_x, row.names = F)
  print(x)
}

