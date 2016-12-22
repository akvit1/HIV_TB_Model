
library(deSolve)



#original model includes: Suceptible (X), Latently Infected (L), Tuberculosis (T), 
# chemoprophylaxed (C), and treated (E) compartments. 

#parameters: 
# pi = rate of people being born and entering the suceptible compartment
# u = rate of people dying in each compartment
# ut = mortaility rate for tuberculosis
# delta = rate of chemoprophylaxis
# theta = rate of being treated of TB
# v = rate of developing TB as a latently infected person
# p = proportion of people who develop TB directly
# (1-p) proportion of people who become latently infected
# beta = per capita force of infection
#



original_model<-function(t, y, param) {
   # N is equal to the total number of people in the model
    N <- y["X"] + y["L"] + y["T"] + y["C"] + y["E"]
    
    
    dX <-  param["pi"] - param["beta"]*y["T"]*y["X"] - param["u"]*y["X"]
    

    
    
    # prop p is the proportion of people who develop TB directly
    dL <- (1 - param["p"])*param["beta"]*y["T"]*y["X"] - (param["v"] + param["u"] + param["delta"])*y["L"]
    
    dC <- param["delta"]*y["L"] - param["u"]*y["C"]
    
    dT <- param["v"]*y["L"] + param["p"]*param["beta"]*y["T"]*y["X"] - (param["u"] + param["ut"] + param["theta"])*y["T"]
    
    # counter of the number of people who have TB, with no treatment and no chemo
    dTn <- param["v"]*y["L"] + param["p"]*param["beta"]*y["T"]*y["X"]# - (param["u"] + param["ut"])*y["T"]
    
    
    
    
    # NOT USING THESE TWO
    dLn <- (1 - param["p"])*param["beta"]*y["X"]*y["T"] - (param["v"] + param["u"])*y["L"]
    dTnod <- param["v"]*y["L"] + param["p"]*param["beta"]*y["X"]*y["T"] - (param["theta"])*y["T"]
    
    
    
    
    dE <-  param["theta"]*y["T"] - param["u"]*y["E"]
    
   
    
   
    return(list(c(dX, dL, dC, dT, dTn, dLn, dE, dTnod)))
}





run_orig_model<-function(beta = 0.0000468, pi = 4400, u = 0.0222, ut = 0.139, delta, theta, v = 0.00256,
                        p = 0.05,
                        
                        initial.state= c(X=200000, L=0, C=0, T = 1, Tn = 1, Ln=0, E=0, Tnod = 1),
                        max.time = MAX.TIME,
                        freq.dependent=FALSE) {
    
   
    beta.divisor <- ifelse(freq.dependent,
                           
                           initial.state["X"]+initial.state["L"]+initial.state["C"]+initial.state["T"]+initial.state["Tn"]+initial.state["Ln"] + initial.state["E"]+ initial.state["Tnod"], 1)
    
    
    
  
    
    #Should we include the proportion p as one of the parameters???
    
    param <- c(beta=beta/beta.divisor, pi=pi, u=u, ut=ut, delta=delta, theta=theta, v=v, p=p)
    
    #Sequence of times at which we want estimates..
    #   ..here we say daily until max.time
    times <- seq(0,max.time,1)
    
    #Run the ODE-Solver
    
    #Change "Original_model" here to run something different
    
    sir.output <- lsoda(initial.state, times, original_model, param)
    
    
    return(sir.output)
    
}


once<-run_orig_model(max.time=10, theta = 0.1, delta = 0)
once

counter1 <- as.data.frame(matrix(0, ncol = 4, nrow = 11))
names(counter1)<- c("trt_only", "prop_treated", "chemo_only", "untreated")
counter1

# NO CHEMO, TREATMENT ONLY!
for (i in seq(0, 1, 0.1)) {
    defaultrun<-run_orig_model(max.time=1000, theta = i, delta = 0.0, beta = 0.0000219488)
    #cases with the control program/# cases with no control program
    
    counter1[1+i*10,1]<-(defaultrun[1000,6])
    counter1[1+i*10,2]<-i
}


# With no treatment, chemo ONLY
for (i in seq(0, 1, 0.1)) {
    defaultrun<-run_orig_model(max.time=1000, theta = 0.0, delta = i, beta = 0.0000219488)
    #cases with the control program/# cases with no control program
    counter1[1+i*10,3]<-(defaultrun[1000,6])
    

}

# NO TREATMENTS


for (i in seq(0, 1, 0.1)) {
    defaultrun<-run_orig_model(max.time=1000, theta = 0.0, delta = 0.0)
    #cases with the control program/# cases with no control program
    counter1[1+i*10,4]<-(defaultrun[1000,6])
    
}



counter1
# FIGURE 2C PLOT
#red is treatment
#blue is chemo
    
plot(counter1[,1]/counter1[,4]~counter1[,2], type = "l", col = "red", main = "Treatment (red) Chemo (blue)")

lines(counter1[,3]/counter1[,4]~counter1[,2], type = "l", col = "blue")







par(mfrow=c(1,1))
plot(defaultrun)

legend(0, 25000, defaultrun)








