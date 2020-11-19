# Installing and Initializing packages
packages <- c('dplyr',
              'ggplot2',
              'deSolve'
)
install.packages(packages)
library(dplyr)
library(ggplot2)
library(deSolve)

# Defining the equations

sir_model <- function(time, state,parametres){
  with(as.list(c(state,parametres)),
       {dS <- -beta*I*S
       dI <- beta*I*S-gamma*I
       dR <- gamma*I
       return(list(c(dS,dI,dR)))
       })
}

# Defining initial values for the variables
state <- c(S =0.99999, I = 1/100000, R = 0)

# Defining the points in time where to calculate variable values

time_values <- seq(0,180)

# Defining the range of the normalized transmission rate(ntrrange)
ntrrange <- c(3.0,2.8,2.5,2.2,2.0,1.8,1.6)

# Defining a function that solves the SIR model for different Normalized Transmission Rates
sirModelSolutions <- function(ntrrange,state, time_values, sir_model){
  
  # Defining parameter values as per the paper
  
  gamma = (1/8)
  implied_beta = 0 #Initalizing the beta 'holder' variable
  parameter_values <- c(beta = 0, gamma = 0)
  sir_values <- list() # to store the different sir models for each ntrrange
  
  
  for(i in 1:length(ntrrange)){ 
    # Defining the parametre values for each Normalized Transmission Rate
    implied_beta = ntrrange[i]*gamma
    parameter_values <- c(beta = implied_beta, gamma = gamma)
    
    # Solving the ODE for each parametre set
    sir_values[[i]] <- ode( y = state,
                            time = time_values,
                            func = sir_model,
                            parms = parameter_values
                            
    )
    i <- i + 1}
  return(sir_values)
}

#####################################################################################################################
# SCENARIO 1: NO 'RESISTANT' AGENTS AT THE START                                                                    #
#####################################################################################################################

# Calling the above function

sir_models_1 <- sirModelSolutions(ntrrange = ntrrange, time_values = time_values, state = state, sir_model = sir_model)

# Calculating the cumulative burden for each model

CB <- c()
for(i in 1:length(sir_models_1)){
  CB <- 1 - sir_models_1[[i]][,2]
  sir_models_1[[i]] <- cbind(sir_models_1[[i]],CB)
  CB <- c()
  i <- i + 1
}

# Plotting the graphs

colors <- c('NTR 3.0'='blue','NTR 2.8'='red','NTR 2.5'='yellow','NTR 2.2'='purple','NTR 2.0'='green','NTR 1.8'='light blue','NTR 1.6'='maroon')
# Plots of the fraction of active infections for each Normalized Transmission Rate
ggplot() +
  geom_line(data = as.data.frame(sir_models_1[[1]]), aes(x = time, y = I, color = 'NTR 3.0')) +
  geom_line(data = as.data.frame(sir_models_1[[2]]), aes(x = time, y = I,color = 'NTR 2.8')) +
  geom_line(data = as.data.frame(sir_models_1[[3]]), aes(x = time, y = I, color = 'NTR 2.5')) +
  geom_line(data = as.data.frame(sir_models_1[[4]]), aes(x = time, y = I, color = 'NTR 2.2')) +
  geom_line(data = as.data.frame(sir_models_1[[5]]), aes(x = time, y = I, color = 'NTR 2.0')) +
  geom_line(data = as.data.frame(sir_models_1[[6]]), aes(x = time, y = I, color = 'NTR 1.8')) +
  geom_line(data = as.data.frame(sir_models_1[[7]]), aes(x = time, y = I, color = 'NTR 1.6')) +
  labs(
    x = "Time (Days)",
    y = "Fraction of Infected",
    color = "Legend") +
  scale_color_manual(values = colors)+
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )


# Plots of the Cumulative burden for each model
ggplot() +
  geom_line(data = as.data.frame(sir_models_1[[1]]), aes(x = time, y = CB, color = 'NTR 3.0')) +
  geom_line(data = as.data.frame(sir_models_1[[2]]), aes(x = time, y = CB,color = 'NTR 2.8')) +
  geom_line(data = as.data.frame(sir_models_1[[3]]), aes(x = time, y = CB, color = 'NTR 2.5')) +
  geom_line(data = as.data.frame(sir_models_1[[4]]), aes(x = time, y = CB, color = 'NTR 2.2')) +
  geom_line(data = as.data.frame(sir_models_1[[5]]), aes(x = time, y = CB, color = 'NTR 2.0')) +
  geom_line(data = as.data.frame(sir_models_1[[6]]), aes(x = time, y = CB, color = 'NTR 1.8')) +
  geom_line(data = as.data.frame(sir_models_1[[7]]), aes(x = time, y = CB, color = 'NTR 1.6')) +
  labs(
    x = "Time (Days)",
    y = "Cumulative burden",
    color = "Legend") +
  scale_color_manual(values = colors)+
  theme(
    legend.position = c(0, 1),
    legend.justification = c("left", "top")
  )
#####################################################################################################################
# SCENARIO 2:  INITIAL 'RECOVERED' ARE 1/3 OF POPULATION                                                            #
#####################################################################################################################

# Defining initial values for the variables
state <- c(S = 0.6666567, I = 1/100000, R = 1/3) # The state var has been altered

# Calling the above function

sir_models_2 <- sirModelSolutions(ntrrange = ntrrange, time_values = time_values, state = state, sir_model = sir_model)

# Calculating the cumulative burden for each model

CB <- c()
for(i in 1:length(sir_models_2)){
  CB <- 1 - sir_models_2[[i]][,2]
  sir_models_2[[i]] <- cbind(sir_models_2[[i]],CB)
  CB <- c()
  i <- i + 1
}

# Plotting the graphs
colors <- c('NTR 3.0'='blue','NTR 2.8'='red','NTR 2.5'='yellow','NTR 2.2'='purple','NTR 2.0'='green','NTR 1.8'='light blue','NTR 1.6'='maroon')
# Plots of the fraction of active infections for each Normalized Transmission Rate
ggplot() +
  geom_line(data = as.data.frame(sir_models_2[[1]]), aes(x = time, y = I, color = 'NTR 3.0')) +
  geom_line(data = as.data.frame(sir_models_2[[2]]), aes(x = time, y = I,color = 'NTR 2.8')) +
  geom_line(data = as.data.frame(sir_models_2[[3]]), aes(x = time, y = I, color = 'NTR 2.5')) +
  geom_line(data = as.data.frame(sir_models_2[[4]]), aes(x = time, y = I, color = 'NTR 2.2')) +
  geom_line(data = as.data.frame(sir_models_2[[5]]), aes(x = time, y = I, color = 'NTR 2.0')) +
  geom_line(data = as.data.frame(sir_models_2[[6]]), aes(x = time, y = I, color = 'NTR 1.8')) +
  geom_line(data = as.data.frame(sir_models_2[[7]]), aes(x = time, y = I, color = 'NTR 1.6')) +
  labs(
    x = "Time (Days)",
    y = "Fraction of Infected",
    color = "Legend") +
  scale_color_manual(values = colors)+
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top")
  )
# Plots of the Cumulative burden for each model

ggplot() +
  geom_line(data = as.data.frame(sir_models_2[[1]]), aes(x = time, y = CB, color = 'NTR 3.0')) +
  geom_line(data = as.data.frame(sir_models_2[[2]]), aes(x = time, y = CB,color = 'NTR 2.8')) +
  geom_line(data = as.data.frame(sir_models_2[[3]]), aes(x = time, y = CB, color = 'NTR 2.5')) +
  geom_line(data = as.data.frame(sir_models_2[[4]]), aes(x = time, y = CB, color = 'NTR 2.2')) +
  geom_line(data = as.data.frame(sir_models_2[[5]]), aes(x = time, y = CB, color = 'NTR 2.0')) +
  geom_line(data = as.data.frame(sir_models_2[[6]]), aes(x = time, y = CB, color = 'NTR 1.8')) +
  geom_line(data = as.data.frame(sir_models_2[[7]]), aes(x = time, y = CB, color = 'NTR 1.6')) +
  labs(
    x = "Time (Days)",
    y = "Cumulative burden",
    color = "Legend") +
  scale_color_manual(values = colors)+
  theme(
    legend.position = c(0, 1),
    legend.justification = c("left", "top")
  )
################################################################################################################
# SCENARIO 3: TEMPORARY MITIGATION                                                                             #
################################################################################################################

# When mitigation is applied, normalized transmission rate changes from 2.5 to 1.5
no_mitigation <- seq(0,180)
state <- c(S =0.99999, I = 1/100000, R = 0)

mitigation_before20 <- seq(0,20)
mitigation_20to50 <- seq(21,50)
mitigation_after50 <- seq(51,180)

mitigation_before50 <- seq(0,50)
mitigation_50to80 <- seq(51,80)
mitigation_after80 <- seq(81,180)

# No mitigation
sir_model_no_mitigation <- sirModelSolutions(ntrrange = 2.5, time_values = time_values, state = state, sir_model = sir_model)
sir_model_no_mitigation <- as.data.frame(sir_model_no_mitigation)
CB <- 1 -sir_model_no_mitigation$S
sir_model_no_mitigation <- cbind(sir_model_no_mitigation,CB)

# Mitigation after 20 days
mitigationImplementor <- function(premitigation,mitigation, postmitigation){
  preMitigation_SIR <- sirModelSolutions(ntrrange = 2.5, time_values = premitigation, state = state, sir_model = sir_model)
  preMitigation_SIR <- as.data.frame(preMitigation_SIR)
  
  state <- c(S = preMitigation_SIR$S[nrow(preMitigation_SIR)],I = preMitigation_SIR$I[nrow(preMitigation_SIR)], R = preMitigation_SIR$R[nrow(preMitigation_SIR)])
  Mitigation_SIR <- sirModelSolutions(ntrrange = 1.5, time_values = mitigation, state = state, sir_model = sir_model)
  Mitigation_SIR <- as.data.frame(Mitigation_SIR)
  
  state <- c(S = Mitigation_SIR$S[nrow(Mitigation_SIR)],I=Mitigation_SIR$I[nrow(Mitigation_SIR)], R = Mitigation_SIR$R[nrow(Mitigation_SIR)])
  postMitigation_SIR <- sirModelSolutions(ntrrange = 2.5, time_values = postmitigation, state = state, sir_model = sir_model)
  postMitigation_SIR <- as.data.frame(postMitigation_SIR)
  
  FullMitigationSIR <- rbind(preMitigation_SIR, Mitigation_SIR, postMitigation_SIR)
  return(FullMitigationSIR)
}

# Implementing temporary mitigation after 20 days

TempMit20to50 <- mitigationImplementor( premitigation = mitigation_before20, mitigation = mitigation_20to50, postmitigation = mitigation_after50)
CB <- 1 -TempMit20to50$S
TempMit20to50 <- cbind(TempMit20to50,CB)
# Implementing temporary mitigation after 50 days

TempMit50to80 <- mitigationImplementor( premitigation = mitigation_before50, mitigation = mitigation_50to80, postmitigation = mitigation_after80)
CB <- 1 -TempMit50to80$S
TempMit50to80 <- cbind(TempMit50to80,CB)

# Plotting the fraction of active infections
colors <- c('No Mitigation' = 'blue', 'Mitigation days 20-50' = 'red', 'Mitigation days 50-80' = 'orange')
ggplot() +
  geom_line(data = sir_model_no_mitigation, aes(x = time, y = I, color = 'No Mitigation')) +
  geom_line(data = TempMit20to50, aes(x = time, y = I, color = 'Mitigation days 20-50')) +
  geom_line(data = TempMit50to80, aes(x = time, y = I, color = 'Mitigation days 50-80')) +
  labs(
    x = "Time (Days)",
    y = "Fraction of Infected",
    color = "Legend") +
  scale_color_manual(values = colors)+
  theme(
    legend.position = c(1,1 ),
    legend.justification = c("right", "top")
  )
# Plot of Cumulative Burden
ggplot() +
  geom_line(data = sir_model_no_mitigation, aes(x = time, y = CB, color = 'No Mitigation')) +
  geom_line(data = TempMit20to50, aes(x = time, y = CB, color = 'Mitigation days 20-50')) +
  geom_line(data = TempMit50to80, aes(x = time, y = CB, color = 'Mitigation days 50-80')) +
  labs(
    x = "Time (Days)",
    y = "Cumulative Burden",
    color = "Legend") +
  scale_color_manual(values = colors)+
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom")
  )
#################################################################################################################
#   FORECASTING A PANDEMIC: STRUCTURAL APPROACH                                                                 #
#################################################################################################################

install.packages('covid19.analytics')
library(covid19.analytics)
US.data <-  covid19.data(case = 'ts-deaths-US')

# Processing the data
# Removing unwanted variables
US.data <- US.data[,-c(1:4)]
# Aggregating all US locations on each day
US.data <- colSums(US.data)
US.data <- as.data.frame(US.data)
dates <- as.Date(rownames(US.data))
US.data <- cbind(dates,US.data)
# Extracting the data from March 8th to April 8th
US.data <- select(filter(US.data, dates >= '2020-04-08' & dates <= '2020-05-08'), names(US.data))

# Estimating the average daily deaths parameter 'a'  using a linear model
t <- seq(1, nrow(US.data))
US.data <- cbind(t,US.data)
lm.model <- lm(US.data$US.data ~ US.data$t)
a <- as.numeric(lm.model$coefficients[2])

# Plotting the Cumulative deaths and the model
predicted <- data.frame(cumdeaths_pred = predict(lm.model), t=US.data$t)

colors <- c('model' = 'orange', 'data' = 'blue')
ggplot()+
  geom_line(data = US.data, aes( x = t, y = US.data, color = 'data'))+
  geom_line(data = predicted, aes(x = t, y = cumdeaths_pred, color = 'model'))+
  labs(
    x = "Time (Days)",
    y = "Cumulative Deaths",
    color = "Legend") +
  scale_color_manual(values = colors)+
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("right", "top"))


# Calculating S,I and R values using the structural approach.(Formulas in pg. 20)
# Normalizing the number of deaths
Dt <- US.data$US.data/328200000
a <- a/328200000
US.data <- cbind(Dt, US.data)

gamma <- 1/8
vee <- c(0.005,0.01,0.002) # Fatality rate

structural_approach <- function(data, gamma,vee,a){
  struc_models <- list()
  S <- I <- R <- 0 
  for(i in 1:length(vee)){
    v <- vee[i]
    I <- (1/(v*gamma))*a
    R <- (1/v)*data$Dt
    S <- 1 - R - I
    NTR <- 1/S
    struc_models[[i]] <- data.frame(t = data$t, S = S, I = I, R = R, NTR = NTR)
    
    i <- i + 1
  }
  return(struc_models)
}
structuralModels <- structural_approach(data = US.data, vee = vee, gamma = gamma, a = a)

# Plotting the fraction of agents still susceptible
colors <- c('Fatality rate 0.005' = 'blue', 'Fatality rate 0.01' = 'red', 'Fatality rate 0.002' = 'orange')
ggplot() +
  geom_line(data = structuralModels[[1]], aes(x = t, y = S, color = 'Fatality rate 0.005'))+
  geom_line(data = structuralModels[[2]], aes(x = t, y = S, color= 'Fatality rate 0.01'))+
  geom_line(data = structuralModels[[3]], aes(x = t, y = S, color = 'Fatality rate 0.002'))+
  labs(
    x = "Time (Days)",
    y = "Fraction of Susceptibles",
    color = "Legend") +
  scale_color_manual(values = colors)+
  theme(
    legend.position = c(0.05, 0.05),
    legend.justification = c("left", "bottom"))
# Plots of Normalized Transmission Rates for each fatality rate
ggplot() +
  geom_line(data = structuralModels[[1]], aes(x = t, y = NTR, color = 'Fatality rate 0.005'))+
  geom_line(data = structuralModels[[2]], aes(x = t, y = NTR, color= 'Fatality rate 0.01'))+
  geom_line(data = structuralModels[[3]], aes(x = t, y = NTR, color = 'Fatality rate 0.002'))+
  labs(
    x = "Time (Days)",
    y = "Normalized Transmission Rates",
    color = "Legend") +
  scale_color_manual(values = colors)+
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"))

# Using the final values in the model, forecasts for the next 30 days are made

# Defining initial values for the variables
state <- c(S =  0.9425, I = 0.0448 ,  R = 0.0477) # The state var has been altered
# Calling the above function

sir_models_3 <- sirModelSolutions(ntrrange = ntrrange, time_values = time_values, state = state, sir_model = sir_model)

# Calculating the cumulative burden for each model

CB <- c()
for(i in 1:length(sir_models_3)){
  CB <- 1 - sir_models_3[[i]][,2]
  sir_models_3[[i]] <- cbind(sir_models_3[[i]],CB)
  CB <- c()
  i <- i + 1
}

# Plotting the graphs
colors <- c('NTR 3.0'='blue','NTR 2.8'='red','NTR 2.5'='yellow','NTR 2.2'='purple','NTR 2.0'='green','NTR 1.8'='light blue','NTR 1.6'='maroon')
# Plots of the fraction of active infections for each Normalized Transmission Rate
ggplot() +
  geom_line(data = as.data.frame(sir_models_3[[1]]), aes(x = time, y = I, color = 'NTR 3.0')) +
  geom_line(data = as.data.frame(sir_models_3[[2]]), aes(x = time, y = I,color = 'NTR 2.8')) +
  geom_line(data = as.data.frame(sir_models_3[[3]]), aes(x = time, y = I, color = 'NTR 2.5')) +
  geom_line(data = as.data.frame(sir_models_3[[4]]), aes(x = time, y = I, color = 'NTR 2.2')) +
  geom_line(data = as.data.frame(sir_models_3[[5]]), aes(x = time, y = I, color = 'NTR 2.0')) +
  geom_line(data = as.data.frame(sir_models_3[[6]]), aes(x = time, y = I, color = 'NTR 1.8')) +
  geom_line(data = as.data.frame(sir_models_3[[7]]), aes(x = time, y = I, color = 'NTR 1.6')) +
  labs(
    x = "Time (Days)",
    y = "Fraction of Infected",
    color = "Legend") +
  scale_color_manual(values = colors)+
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )


# Plots of the Cumulative burden for each model
ggplot() +
  geom_line(data = as.data.frame(sir_models_3[[1]]), aes(x = time, y = CB, color = 'NTR 3.0')) +
  geom_line(data = as.data.frame(sir_models_3[[2]]), aes(x = time, y = CB,color = 'NTR 2.8')) +
  geom_line(data = as.data.frame(sir_models_3[[3]]), aes(x = time, y = CB, color = 'NTR 2.5')) +
  geom_line(data = as.data.frame(sir_models_3[[4]]), aes(x = time, y = CB, color = 'NTR 2.2')) +
  geom_line(data = as.data.frame(sir_models_3[[5]]), aes(x = time, y = CB, color = 'NTR 2.0')) +
  geom_line(data = as.data.frame(sir_models_3[[6]]), aes(x = time, y = CB, color = 'NTR 1.8')) +
  geom_line(data = as.data.frame(sir_models_3[[7]]), aes(x = time, y = CB, color = 'NTR 1.6')) +
  labs(
    x = "Time (Days)",
    y = "Cumulative burden",
    color = "Legend") +
  scale_color_manual(values = colors)+
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom")
  )
