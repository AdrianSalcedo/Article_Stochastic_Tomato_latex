library(yuima)
library(ggplot2)
library(gridExtra)
library(grid)

library(lattice)
T <- 300
n <- 10000
################################################################################
########################    Parameters    #############################
################################################################################

beta_p <- 0.8
r_1 <- 0.6
r_2 <- 0.6
b <- 0.5
beta_v <- 0.8
theta <- 0.818469652 
mu <- 0.058596601 
gamma <- 0.06
sigma_L <- 0.378173863
sigma_I <- 0.030015876
sigma_v <- 0.254146903
N_v <- mu/gamma

################################################################################
########################    Stochastic R0 and R0    ############################
################################################################################

r<-r_1
R_0 <- sqrt(beta_v * b * beta_p / ( r_2 * ( r_1 + b) * gamma))
print(R_0)
R_s_0 <- (beta_v * beta_p) / (r * gamma)
print(R_s_0)

################################################################################
############################# Initial Condition ################################
################################################################################

xinitial_det <- c(97,1,2,3,4)
xinitial_stc <- xinitial_det
#xinitial_stc <- log(xinitial_det)

N_p <- xinitial_det[1]+xinitial_det[2]+xinitial_det[3]

sol <- c("x","y","z","v","w") # variable for numerical solution

################################################################################
############################# System Dynamic ################################
################################################################################

f1 <- c("-beta_p*x*w/N_v+r_1*y+r_2*z",
        "beta_p*x*w/N_v-(b+r_1)*y",
        "b*y-r_2*z",
        "-beta_p*v*z/N_p-gamma*v+(1-theta)*mu",
        "beta_v*v*z/N_p-gamma*w+theta*mu") # drift vector

g1 <- matrix(c("sigma_L*y","-sigma_L*y","0","0","0",
               "sigma_I*z","0","-sigma_I*z","0","0",
               "0","0","0","-sigma_v*v","-sigma_v*w"),5,3)

################################################################################
############################# Positive drift ###################################
################################################################################
f2 <- c("-(beta_p/N_v)*exp(w)+r_1*exp(y-x)+r_2*exp(z-x)-(1/2)*((sigma_L*exp(y)+
        sigma_I*exp(z))/exp(x))^2",
        "(beta_p/N_v)*exp(x+w-y)-(b+r_1)-(1/2)*sigma_L^2",
        "b*exp(y-z)-r_2-(1/2)*sigma_I^2",
        "-beta_v*exp(z)/(exp(x)+exp(y)+exp(z))-gamma+(1-theta)*mu*exp(-v)-(1/2)*sigma_v^2",
        "beta_v*exp(v+z-w)/(exp(x)+exp(y)+exp(z))-gamma+(theta)*mu*exp(-w)-(1/2)*sigma_v^2"
) # drift vector


g2 <- matrix(c("(sigma_L*exp(y)+sigma_I*exp(z))/exp(x)",
               "-sigma_L",
               "-sigma_I",
               "-sigma_v",
               "-sigma_v"),5,1)

################################################################################
######################### Deterministic Deriva ###############################
################################################################################

g3 <- matrix(c("0",
               "0",
               "0",
               "0",
               "0"),5,1)

################################################################################
########################## Deterministic path ##################################
################################################################################

mymod_det <- setModel(drift=f1,diffusion=g3,solve.variable= sol)

mymod_det.samp <- setSampling(Terminal=T, n=n)
mymod_det <- setYuima(model=mymod_det, sampling=mymod_det.samp)

mymod_det.samp
X <- simulate(mymod_det,xinit= xinitial_det,true.parameter=list(beta_p =beta_p, r_1=r_1,
                                                                r_2=r_2, b=b, beta_v=beta_v,
                                                                theta=theta, mu=mu, gamma = gamma,
                                                                N_v=N_v,sigma_L = sigma_L,
                                                                sigma_I = sigma_I,
                                                                sigma_v =sigma_v
))

################################################################################
########################## Stochastic path #################################
################################################################################

mymod_stc <- setModel(drift=f1,diffusion=g1,solve.variable= sol)




mymod_stc.samp <- setSampling(Terminal=T, n=n)
mymod_stc <- setYuima(model=mymod_stc, sampling=mymod_stc.samp)

mymod_stc.samp
Y <- simulate(mymod_stc,xinit= xinitial_stc,true.parameter=list(beta_p =beta_p, r_1=r_1,
                                                                r_2=r_2, b=b, beta_v=beta_v,
                                                                theta=theta, mu=mu, gamma = gamma,
                                                                N_v=N_v,sigma_L = sigma_L,
                                                                sigma_I = sigma_I,
                                                                sigma_v =sigma_v
),y)

################################################################################
########################### Data recopilation ##################################
################################################################################

Time_det <- X@sampling@grid[[1]]

Df_det <- data.frame(t= Time_det,x = X@data@original.data[,1],
                     y = X@data@original.data[,2], z = X@data@original.data[,3],
                     v = X@data@original.data[,4], w = X@data@original.data[,5])

#write.csv(Df_det, file = "Tomato_data_det_R0_more.csv",row.names = TRUE)

Time_stc <- Y@sampling@grid[[1]]

Df_stc <- data.frame(t= Time_stc,x = Y@data@original.data[,1],
                     y = Y@data@original.data[,2], z = Y@data@original.data[,3],
                     v = Y@data@original.data[,4], w = Y@data@original.data[,5])

#Df_stc <- data.frame(t= Time_stc,x = exp(Y@data@original.data[,1]),
#                     y = exp(Y@data@original.data[,2]), z = exp(Y@data@original.data[,3]),
#                     v = exp(Y@data@original.data[,4]), w = exp(Y@data@original.data[,5]))

#write.csv(Df_stc, file = "Tomato_data_stc_R0_more.csv",row.names = TRUE)

################################################################################
########################## Plot of the data ####################################
################################################################################
plot1 <- ggplot()+ geom_line(data=Df_stc, aes(t,x),linetype = "solid", 
                             color = "blue") + ggtitle("Suceptibles plants")

plot2 <- ggplot()+ geom_line(data=Df_stc, aes(t,y),linetype = "solid",
                             color = "blue") + ggtitle("Latents plants")


plot3 <- ggplot()+ geom_line(data=Df_stc, aes(t,z),linetype = "solid",
                             color = "blue")+ ggtitle("Infected plants")

plot4 <- ggplot()+ geom_line(data=Df_stc, aes(t,v),linetype = "solid",
                             color = "blue") + ggtitle("Susceptibles vectors")

plot5 <- ggplot() + geom_line(data=Df_stc, aes(t,w),linetype = "solid", 
                              color = "blue") + ggtitle("Infected vectors")

grid.arrange(plot1, plot2,plot3,plot4,plot5, nrow = 2,ncol=3)
#################################################################################
plot1 <- ggplot()+ geom_line(data=Df_det, aes(t,x),linetype = "solid", 
                             color = "blue") + ggtitle("Suceptibles plants")

plot2 <- ggplot()+ geom_line(data=Df_det, aes(t,y),linetype = "solid",
                             color = "blue") + ggtitle("Latents plants")


plot3 <- ggplot()+ geom_line(data=Df_det, aes(t,z),linetype = "solid",
                             color = "blue")+ ggtitle("Infected plants")

plot4 <- ggplot()+ geom_line(data=Df_det, aes(t,v),linetype = "solid",
                             color = "blue") + ggtitle("Susceptibles vectors")

plot5 <- ggplot() + geom_line(data=Df_det, aes(t,w),linetype = "solid", 
                              color = "blue") + ggtitle("Infected vectors")

grid.arrange(plot1, plot2,plot3,plot4,plot5, nrow = 2,ncol=3)

################################################################################
################################################################################
################################################################################

plot11 <- ggplot()+ geom_line(data=Df_stc, aes(t,x),linetype = "solid",
                              color = "darkgreen") +
  geom_line(data=Df_det, aes(t,x),linetype = "solid", 
            color = "blue") + #coord_cartesian(ylim = c(0, 100)) +
  ggtitle("Suceptibles plants")

plot12 <- ggplot()+ geom_line(data=Df_stc, aes(t,y),linetype = "solid",
                              color = "orange2")+
  geom_line(data=Df_det, aes(t,y),linetype = "solid",
            color = "blue")  + coord_cartesian(ylim = c(0, 100))  + 
  ggtitle("Latents plants")

plot13 <- ggplot()+ geom_line(data=Df_stc, aes(t,z),linetype = "solid",
                              color = "darkred") +
  geom_line(data=Df_det, aes(t,z),linetype = "solid",
            color = "blue")  + coord_cartesian(ylim = c(0, 100))  + 
  ggtitle("Infected plants") 

plot14 <- ggplot()+ geom_line(data=Df_stc, aes(t,v),linetype = "solid",
                              color = "darkorchid4") +
  geom_line(data=Df_det, aes(t,v),linetype = "solid",
            color = "blue")  + coord_cartesian(ylim = c(0, 5))  + 
  ggtitle("Susceptibles vectors")

plot15 <- ggplot() + geom_line(data=Df_stc, aes(t,w),linetype = "solid",
                               color = "orange") +
  geom_line(data=Df_det, aes(t,w),linetype = "solid", color = "blue") +
  ggtitle("Infected vectors")  + coord_cartesian(ylim = c(0, 5))
Law_cons <- data.frame(t=Df_stc["t"],LC = Df_stc["x"]+Df_stc["y"]+Df_stc["z"])


plot16 <- ggplot() + geom_line(data=Law_cons, aes(t,x),linetype = "solid", 
                               color = "orange") + coord_cartesian(ylim = c(0,120))

grid.arrange(plot11, plot12,plot13,plot14,plot15,plot16, nrow = 2,ncol=3)


