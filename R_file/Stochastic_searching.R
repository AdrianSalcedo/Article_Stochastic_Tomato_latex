

#################################################################################
 ##### Generating uniform random values for each parameter on its interval #####
#################################################################################
sampling <- function(){
  beta_p <- runif(1,0,1)#runif(1,0,0.1)
  r_1 <- runif(1,0,2)#runif(1,0.0083,0.013)
  r_2 <- runif(1,0,2)#runif(1,0.0083,0.013)
  b <- runif(1,0.05,0.1)
  beta_v <-runif(1,0,1)#runif(1,0,0.2)
  theta <- runif(1,0,1)
  mu <- runif(1,0,1)
  gamma <- runif(1,0.03,1)
  sigma_L <- runif(1,0,2)#0.005)
  sigma_I <- runif(1,0,2)#0.005)
  sigma_v <- runif(1,0,2)#1)
  
  parameter <- c(beta_p,r_1,r_2,b,beta_v,theta,mu,gamma,sigma_L,sigma_I,sigma_v)
  
  return(parameter)
}
Parameter <-sampling()
while (TRUE){
  beta_p <- Parameter[1]
  r_1 <- Parameter[2]
  r_2 <- Parameter[3]
  b <- Parameter[4]
  beta_v <- Parameter[5]
  theta <- Parameter[6]
  mu <- Parameter[7]
  gamma <- Parameter[8]
  sigma_L <- Parameter[9]
  sigma_I <- Parameter[10]
  sigma_v <- Parameter[11]
  R_d_0 <- (beta_p*beta_v*b)/(gamma*(b+r_1)*r_2)
  cond1 <- R_d_0 > 1
  cond2 <- TRUE#beta_p > sigma_L^2
  cond3 <- TRUE#r_2 > sigma_I^2
  cond4 <- (beta_p^2)/(2*sigma_L^2)+2*beta_p-r_1+(r_2^2)/(2*sigma_I^2) < 0
  cond5 <-TRUE #beta_v > sigma_v^2
  cond6 <- (beta_v^2)/(2*sigma_v^2)+beta_v-gamma+theta*mu < 0
  if (isTRUE(cond1 & cond2 & cond3 & cond4 & cond5 & cond6) == TRUE)
    break
  else {
    Parameter <-sampling()
  }    
}
print(R_d_0)
print(Parameter)

