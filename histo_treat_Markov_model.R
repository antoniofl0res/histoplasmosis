### Histoplasmosis treatment model ###

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


plot_sim("Died")
