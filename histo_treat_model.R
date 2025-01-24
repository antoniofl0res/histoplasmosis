### Histoplasmosis treatment model ###

L = dim(post_seq)[1] #length of vector parameters

PopN = 1e6
pos1 = rbinom(L, PopN, post_seq$Pa)
pos2 = rbinom(L, pos1 , post_seq$Pa2)



p = list(
  Tr = rbeta(L, 90, 10),  #proprotion treated
  E = rbeta(L, 45, 5),  # Treatment effect
  A = rbeta(L, 6, 44),  # Risk of adverse events
  Ma = rbeta(L, 15, 35), # Risk of death (AHD)
  Mh = rbeta(L, 20, 30) # Risk of death (histoplasmosis)
)

summary(as.bayesboot(p))


dh2 = data.frame(TP = rbinom(L, pos2, post_seq$PPV2_pos),
                FP = rbinom(L, pos2, (1 - post_seq$PPV2_pos)),
                TN = rbinom(L, (PopN - pos1), post_seq$NPV) + rbinom(L, pos1-pos2, post_seq$NPV2_pos),
                FN = rbinom(L, (PopN - pos1), (1-post_seq$NPV)) + rbinom(L, pos1-pos2, (1 - post_seq$NPV2_pos)))

dh1 = data.frame(TP = rbinom(L, pos1, post_seq$PPV),
                 FP = rbinom(L, pos1, (1 - post_seq$PPV)),
                 TN = rbinom(L, (PopN - pos1), post_seq$NPV),
                 FN = rbinom(L, (PopN - pos1), (1-post_seq$NPV)))



histo <- function(dataset, params = p) {
 
  output = data.frame(Disease = numeric(L), Positive = numeric(L), Treated = numeric(L), 
                      Overtreat = numeric(L), Undertreat = numeric(L),
                      Adverse_events = numeric(L), 
                      Alive = numeric(L), Died = numeric(L))
  output$Disease = dataset$TP + dataset$FN
  output$Positive <- (dataset$TP + dataset$FP)
  output$Treated <- (dataset$TP + dataset$FP) * params$Tr
  output$Overtreat <- dataset$FP * params$Tr
  output$Undertreat <- dataset$FN
  output$Undertreat_mortality <- dataset$FN * params$Mh
  output$Adverse_events <- output$Treated * params$A
  output$excess_AE <- output$Overtreat * params$A
  output$Alive <- dataset$TP * (params$Tr * params$E + (1 - params$Tr) * (1 - params$Mh)) +
    dataset$FN * (1 - params$Mh) + (dataset$FP + dataset$TN)  * (1 - params$Ma)
   
  output$Died <- dataset$TP * (params$Tr * (1 - params$E) + (1 - params$Tr) * params$Mh) +
    dataset$FN * params$Mh + (dataset$FP + dataset$TN) * params$Ma
  
  output = round(output)
  output = output |> rowwise() |> mutate(Total = Alive + Died,
                                         Treat_prop = Treated/Total,
                                         overtreat_AE_prop = excess_AE/Treated,
                                         undertreat_mort_rate = Undertreat_mortality/Total*1000,
                                         undertreat_mort_prop = Undertreat_mortality/Died,
                                         mortality = Died/(Alive+Died))
  
  return(output)
}

histo.sum = function(df) {
  
  t(round(rbind(expected = colMeans(df), hdi(df)),3))
}


T1 = histo(dataset = dh1) # Treatment decision based on single test (LFA)

histo.sum(T1)


T2 = histo(dataset = dh2) # Treatment decision based on sequential testing
histo.sum(T2)


# Treatment #

plotPost(T1$Treat_prop-T2$Treat_prop, compVal = 0)

# Unnecessary harm 
### how much more toxicity a single test may cause #

plotPost(T1$Treat_prop/ T2$Treat_prop, compVal = 0)
plotPost(1/(T1$Treat_prop - T2$Treat_prop), compVal = 0)

plotPost(T1$Overtreat/T2$Overtreat, compVal = 0)
plotPost(T1$overtreat_AE_prop/T2$overtreat_AE_prop, compVal = 0) # Risk difference
plotPost(1/(T1$overtreat_AE_prop-T2$overtreat_AE_prop), xlab="NNH") # NNH - number need to harm


### mortality attributable to underdiagnosis

plotPost(T2$undertreat_mort_rate/T1$undertreat_mort_rate, compVal = 1)
plotPost(1/(T2$undertreat_mort_prop - T1$undertreat_mort_prop), xlab="NNH - no treatment")


### mortality

### Cost

cost_first_test = 15
cost_second_test = 30
cost_treatment = 240
cost_mngtm_AE = 60


cost2 = (T1$Positive * cost_second_test + 
           T1$Total * cost_first_test + 
           T2$Treated * cost_treatment +
           T2$Adverse_events * cost_mngtm_AE)
cost1 = (T1$Total * cost_first_test + 
           T1$Treated * cost_treatment +
           T1$Adverse_events * cost_mngtm_AE)

cost_per_death_averted = (cost1-cost2)/(T2$Undertreat_mortality-T1$Undertreat_mortality)
summary(as.bayesboot(cost_per_death_averted))

plotPost(cost_per_death_averted, compVal = 0)


#### Markov model ####

ini = c(Non_disease=9000,
        TN=0,
        FP=0,
        Disease=1000, 
        TP=0,
        FN=0,
        Treated=0,
        Untreated=0,
        AdvE=0,
        Recovered=0,
        Died=0) # Initial vector

# Initialize the result matrix to store results for iterations
result <- matrix(0, nrow = 30, ncol = length(ini))
colnames(result) = names(ini)

# Set the initial state as the first row of the result
result[1, ] <- ini


# Stochastic transition matrix
for (i in 2:nrow(result)) {
  
  # Parameters
  
  Se = rbeta(1, 204, 20) 
  Sp = rbeta(1, 956, 44)
  Treat = rbeta(1, 90, 10) #Treatment uptake
  AE = rbeta(1, .9, 99.1) #adverse events (severe)
  M = rbeta(1, 2, 98) # disease-related mortality
  C = rbeta(1, 7.5, 992.5) #mortality from complications
  
  # Create a new random transition matrix for each iteration
  markov <- matrix(0, nrow=11, ncol = 11, 
                   dimnames = list(names(ini), names(ini)))
  markov["Non_disease", "TN"] = Sp
  markov["Non_disease", "FP"] = 1-Sp
  markov["Disease", "TP"] = Se
  markov["Disease", "FN"] = 1-Se
  markov["TN", "TN"] = 1
  markov["TP", "Treated"] = Treat
  markov["TP", "Untreated"] = 1-Treat
  markov["FN", "Untreated"] = 1
  markov["FP", "Treated"] = Treat
  markov["FP", "Recovered"] = 1-Treat
  markov["Treated", "AdvE"] = AE
  markov["Treated", "Recovered"] = 1-AE
  markov["AdvE", "Died"] = C
  markov["AdvE", "Recovered"] = 1-C
  markov["Untreated", 'Died'] = M
  markov["Untreated", 'Untreated'] = 1-M
  markov["Recovered", "Recovered"] = 1
  markov["Died", "Died"] = 1
  
  
  
  # Normalize rows to ensure the sum is 1
  
  #markov <- sweep(markov, 1, rowSums(markov), FUN = "/")
  
  # Multiply the previous state by the transition matrix
  result[i, ] <- result[i - 1, ] %*% markov
}


res = result |> data.frame() |> mutate(Step = 1:n()) |> 
  pivot_longer(1:length(ini), names_to = "state")

ggplot(res |> filter(state!="TN", state!="Non_disease"), aes(x=Step, y=value, color = state)) + 
  geom_line()


### Multiple Markov simulations #####

# Number of iterations for Markov chain (time points)
n_time = 30
# Number of simulations to perform
n_simulations = 1000

ini = c(Non_disease=9900,
        TN=0,
        FP=0,
        Disease=100, 
        TP=0,
        FN=0,
        Treated=0,
        Untreated=0,
        AdvE=0,
        Recovered=0,
        Died=0) # Initial vector


# Create an array to store the results of all simulations
# Dimensions: [time steps, states, simulations]
all_results <- array(0, dim = c(n_time, length(ini), n_simulations),
                     dimnames = list(1:n_time, names(ini), 1:n_simulations))

# Loop for each simulation
for (sim in 1:n_simulations) {
  
  # Initialize the result matrix for each simulation
  result <- matrix(0, nrow = n_time, ncol = length(ini))
  colnames(result) <- names(ini)
  
  # Set the initial state as the first row of the result
  result[1, ] <- ini
  
  
  
  # Stochastic transition matrix loop
  for (i in 2:nrow(result)) {
    
    # Parameters
    Se = rbeta(1, 204, 20) 
    Sp = rbeta(1, 956, 44)
    Treat = rbeta(1, 995, 5) # Treatment uptake
    AE = rbeta(1, 0.9, 99.1) # Adverse events (severe)
    M = rbeta(1, 2, 98) # Disease-related mortality
    C = rbeta(1, 7.5, 992.5) # Mortality from complications
    
    # Create a new random transition matrix for each iteration
    markov <- matrix(0, nrow=11, ncol = 11, 
                     dimnames = list(names(ini), names(ini)))
    markov["Non_disease", "TN"] = Sp
    markov["Non_disease", "FP"] = 1 - Sp
    markov["Disease", "TP"] = Se
    markov["Disease", "FN"] = 1 - Se
    markov["TN", "TN"] = 1
    markov["TP", "Treated"] = Treat
    markov["TP", "Untreated"] = 1 - Treat
    markov["FN", "Untreated"] = 1
    markov["FP", "Treated"] = Treat
    markov["FP", "Recovered"] = 1 - Treat
    markov["Treated", "AdvE"] = AE
    markov["Treated", "Recovered"] = 1 - AE
    markov["AdvE", "Died"] = C
    markov["AdvE", "Recovered"] = 1 - C
    markov["Untreated", 'Died'] = M
    markov["Untreated", 'Untreated'] = 1 - M
    markov["Recovered", "Recovered"] = 1
    markov["Died", "Died"] = 1
    
    # Multiply the previous state by the transition matrix
    result[i, ] <- result[i - 1, ] %*% markov
  }
  
  # Store the result of the current simulation in the array
  all_results[, , sim] <- result
}

# Now `all_results` contains the results of all simulations


df_long = function(df=all_results, col) {df[, col, ] |> data.frame() |> 
    mutate(Step = 1:n()) |> pivot_longer(1:n_simulations, names_to = "Sim")}



plot_sim = function(col, step=n_time){
  df = df_long(col=col)
  df_mean = df[df[["Step"]] == step, "value"] |> 
    summarise(mean=mean(value)/sum(ini)) |> round(digits=3) |> pull(mean)
  
  # Calculate uncertainty interval at the last time step
  df_ui <- df[df[["Step"]] == step, "value"] |> 
    reframe(ui = c(quantile(value/sum(ini), c(.025, .975)))) |> round(digits = 3) |>
    pull(ui)
  
  
  plot = ggplot(df, aes(x=Step, y=value, color=Sim)) +  
           geom_line() + theme(legend.position = "none") +
    labs(title = paste("Simulations:", col),
         subtitle = paste0("Mean: ", df_mean, " (95% UI: ", df_ui[1],", ", df_ui[2], ")"))
  
  return(plot)
}


