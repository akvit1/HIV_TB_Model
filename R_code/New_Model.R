

library(deSolve)

#New model includes: Suceptible (X), Latently Infected (L), Tuberculosis (T), 
# chemoprophylaxed (C), and recovered (R) compartments for HIV negative people. 
#Suceptible (Xhiv), Latently Infected (Lhiv), Tuberculosis (Thiv), 
# chemoprophylaxed (Chiv), and recovered (Rhiv) compartments for HIV positive people

#parameters: 
# HIV negative:
# lambda = per capita force of infection
# pi = rate of people being born and entering the suceptible compartment
# u = rate of people dying in each compartment
# ut = mortaility rate for tuberculosis
# delta = rate of chemoprophylaxis
# theta = rate of being treated of TB
# v = rate of developing TB as a latently infected person
# p = proportion of people who develop TB directly
# (1-p) proportion of people who become latently infected
# h = rate of people getting HIV
#
# HIV positive:
# lambda_hiv = per capita force of infection
# r = the rate at which HIV+ people who recovered from TB become sucpetible to TB again
# u_hiv = rate of people dying in each compartment
# ut_hiv = mortaility rate for tuberculosis
# delta_hiv = rate of chemoprophylaxis
# theta_hiv = rate of being treated of TB
# v_hiv = rate of developing TB as a latently infected person

# beta in the original model is called lambda in this model

new_model<-function(t, y, param) {
    # N is equal to the total number of people in the model
   # N <- y["X"] + y["L"] + y["T"] + y["C"] + y["R"] + y["Xhiv"] + y["Lhiv"] + y["Thiv"] + y["Chiv"] + y["Rhiv"] 
    
    # HIV negative people
    dX <-  param["pi"] - param["lambda"]*y["T"]*y["X"] - param["u"]*y["X"] - param["h"]*y["X"]
    
    # param p is the proportion of people who develop TB directly
    dL <- (1 - param["p"])*param["lambda"]*y["T"]*y["X"] - (param["v"] + param["u"] + param["h"] + param["delta"])*y["L"]
    
    dC <- param["delta"]*y["L"] - param["u"]*y["C"] - param["h"]*y["C"]
    
    dT <- param["v"]*y["L"] + param["p"]*param["lambda"]*y["T"]*y["X"] - (param["u"] + param["ut"] + param["theta"] + param["h"])*y["T"]
    
    #Total Non HIV TB Infection Counter
    dTn <- param["v"]*y["L"] + param["p"]*param["lambda"]*y["T"]*y["X"]# - (param["u"] + param["ut"] + param["theta"] + param["h"])*y["T"]
    
    #count number of TB deaths in non hiv
    dTBdeth <- param["ut"]*y["T"]
    
    dR <-  param["theta"]*y["T"] - param["u"]*y["R"] - param["h"]*y["R"]
    
    #HIV positive people
    
    dXhiv <-  param["h"]*y["X"] + param["r"]*y["Rhiv"] - param["lambda_hiv"]*y["Xhiv"] - param["u_hiv"]*y["Xhiv"]
   
    
    dLhiv <- param["h"]*y["L"] + (1 - param["p"])*param["lambda_hiv"]*y["Thiv"]*y["Xhiv"] - (param["v_hiv"] + param["u_hiv"]  + param["delta_hiv"])*y["Lhiv"]
    
    dChiv <- param["h"]*y["C"] + param["delta_hiv"]*y["Lhiv"] - param["u_hiv"]*y["Chiv"]
    
    dThiv <- param["h"]*y["T"] + param["v_hiv"]*y["Lhiv"] + param["p"]*param["lambda_hiv"]*y["Thiv"]*y["Xhiv"] - (param["u_hiv"] + param["ut_hiv"] + param["theta_hiv"])*y["Thiv"]
    
    #HIV+ TB counter
    dTnhiv <- param["h"]*y["T"] + param["v_hiv"]*y["Lhiv"] + param["p"]*param["lambda_hiv"]*y["Thiv"]*y["Xhiv"]# - (param["u_hiv"] + param["ut_hiv"] + param["theta_hiv"])*y["Thiv"]
    
    #count number of TB deaths in HIV+
    dTBhivdeth <- param["ut_hiv"]*y["Thiv"]
    
    dRhiv <-  param["h"]*y["R"] + param["theta_hiv"]*y["Thiv"] - param["u_hiv"]*y["Rhiv"] - param["r"]*y["Rhiv"]
    
    return(list(c(dX, dL, dC, dT, dTn, dTBdeth, dR, dXhiv, dLhiv, dChiv, dThiv, dTnhiv, dTBhivdeth, dRhiv)))
}



run_new_model<-function(lambda = 0.0000468, pi = 4400, u = 0.0222, ut = 0.139, delta = 0, theta = 0, v = 0.00256,
                        p = 0.05, h = 0.01, lambda_hiv = 0.000067, r = 0.2, u_hiv = 0.044, ut_hiv = 0.222,
                        delta_hiv = 0, theta_hiv = 0, v_hiv = 0.00256,
    
                         initial.state= c(X=200000, L=0, C=0, T = 1, Tn = 1,dTBdeth = 0, R=0, Xhiv = 0, Lhiv = 0, Chiv = 0,
                                          Thiv = 0, dTnhiv = 0, dTBhivdeth=0,  Rhiv = 0),
                         max.time = MAX.TIME,
                         freq.dependent=FALSE) {
    
    #Is our model frequency dependent?
    
                             
    #If the model is frequency dependent we modify beta
    #based on the total populations size
    lambda.divisor <- ifelse(freq.dependent,
                           
                           initial.state["X"]+initial.state["L"]+initial.state["C"]+initial.state["T"]+initial.state["Tn"]  +initial.state["TBdeth"]+ initial.state["R"] + initial.state["Xhiv"]+initial.state["Lhiv"]+initial.state["Chiv"]+initial.state["Thiv"] +initial.state["Tnhiv"]+ initial.state["TBhivdeth"]+ initial.state["Rhiv"], 1)
    
    
    
    #create the parameter vector.
    ###Here is a change (#4)
    # the equal signs are so that you can call the names of points in the vector
    
    # what is beta/beta.divisor? Why do you need the divisor???
    
    #Should we include the proportion p as one of the parameters???
    
    param <- c(lambda=lambda/lambda.divisor, pi=pi, u=u, ut=ut, delta=delta, theta=theta, v=v, p=p, h=h,
               lambda_hiv = lambda_hiv, r=r, u_hiv = u_hiv, delta_hiv=delta_hiv, theta_hiv=theta_hiv,
               v_hiv=v_hiv, ut_hiv=ut_hiv)
    
    #Sequence of times at which we want estimates..
    #   ..here we say daily until max.time
    times <- seq(0,max.time,1)
    
    #Run the ODE-Solver
 
    sir.output <- lsoda(initial.state, times, new_model, param)
    
    
    return(sir.output)
}


once<-run_new_model(max.time=10, theta = 0.1, delta = 0, theta_hiv = 0.1, delta_hiv = 0)
once

##############################################################################################

# HIV negative counter

counter1 <- as.data.frame(matrix(0, ncol = 10, nrow = 11))
names(counter1)<- c("trt_only", "theta", "chemo_only", "untreated", "trt_only_HIV", "theta_HIV", "chemo_only_HIV", "untreated_HIV", "TB death", "HIV TB death")
counter1



for (i in seq(0, 1, 0.1)) {
    defaultrun<-run_new_model(max.time=1000, theta = i, delta = 0.0, theta_hiv = i, delta_hiv = 0.0, lambda = 0.0000219488, lambda_hiv = 0.000119488)
    #cases with the control program/# cases with no control program
    
    counter1[1+i*10,1]<-(defaultrun[1000,6])
    counter1[1+i*10,2]<-i
    counter1[1+i*10,5]<-(defaultrun[1000,13]) # Tn HIV
    counter1[1+i*10,6]<-i # theta HIV
    
    counter1[1+i*10,9]<-(defaultrun[1000,7])  #number of TB deaths
    counter1[1+i*10,10]<-(defaultrun[1000,14])  #number of HIV TB deaths
    
}


counter1

# With no treatment, chemo ONLY
for (i in seq(0, 1, 0.1)) {
    defaultrun<-run_new_model(max.time=1000, theta = 0.0, delta = i, theta_hiv = 0.0, delta_hiv = i, lambda = 0.0000219488, lambda_hiv = 0.000119488)
    #cases with the control program/# cases with no control program
    counter1[1+i*10,3]<-(defaultrun[1000,6])
    counter1[1+i*10,7]<-(defaultrun[1000,13])
    
}

# NO TREATMENTS


for (i in seq(0, 1, 0.1)) {
    defaultrun<-run_new_model(max.time=1000, theta = 0.0, delta = 0.0, theta_hiv = 0.0, delta_hiv = 0.0)
    #cases with the control program/# cases with no control program
    counter1[1+i*10,4]<-(defaultrun[1000,6]) #untreated people
    counter1[1+i*10,8]<-(defaultrun[1000,13]) # untreated people HIV
    
    
}



counter1

# FIGURE 2C PLOT
#red is treatment
#blue is chemo


# HIV NEGATIVE PEOPLE

plot(counter1[,1]/counter1[,4]~counter1[,2], type = "l", col = "red", main = "Treatment (red) Chemo (blue) HIV Negative")

lines(counter1[,3]/counter1[,4]~counter1[,2], type = "l", col = "blue")



#HIV POSITIVE PEOPLE
plot(counter1[,5]/counter1[,8]~counter1[,6], type = "l", col = "red", main = "Treatment (red) Chemo (blue) HIV Positive")

lines(counter1[,7]/counter1[,8]~counter1[,6], type = "l", col = "blue")


# TB Deaths of HIV negative people

plot(counter1[,9]/counter1[,4]~counter1[,2], type = "l", col = "red", main = "HIV- (red) and HIV+ (blue) TB Death Rate")
lines(counter1[,10]/counter1[,8]~counter1[,6], type = "l", col = "blue")



############## default plot


fullrun<-run_new_model(max.time=10)
plot(fullrun)
legend(200, 0.7, fullrun)


