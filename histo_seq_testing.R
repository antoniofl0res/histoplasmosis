### Infering prevalence from a POC test ###

library(tidyverse)
library(gtools)
library(HDInterval)
library(rstan)
library(bayesboot)
library(gt)


### DATA LIST #####

# Store data in a list for passing to Stan models
# Define sample sizes and positive samples for initial and subsequent tests
data_list <- list(
  N = 100, 
  npos = 15,
  N2 = 15,
  npos2 = 12,
  N3 = 0,
  npos3 = 0
)

### STAN SEQUENTIAL MODEL ####

# Fit the sequential Stan model using the provided data
seq_model <- rstan::stan(file="model_seq_pos_retested_binom.stan", 
                         data=data_list, chains=3, cores=(parallel::detectCores()*.8), 
                         iter=6e3, warmup=1e3)



# Extract posterior samples from the sequential model
post_seq <- rstan::extract(seq_model) |> as.data.frame()

# Summarize key statistics from the posterior samples
list(
  N = c(Total=data_list$N, pos1=data_list$npos, 
        pos_retested = data_list$N2, 
        pos2=data_list$npos2,
        neg_retested = data_list$N3,
        pos3 = data_list$npos3),
  expected_PPV = round(mean(post_seq$PPV)*100, 1), 
  hdi_PPV = round(hdi(post_seq$PPV)*100, 1), 
  prob_50 = mean(post_seq$PPV>.5)
)

# Prepare data for display
post_seq = post_seq[, 1:(ncol(post_seq)-1)]

# Generate a table for output display
round((seq.out = as.data.frame(t(rbind(expected=colMeans(post_seq), 
                                       HDInterval::hdi(post_seq))))),3)

data.frame(cbind(parameters=colnames(post_seq), round(seq.out,3))) |> gt()

### POSTERIOR CHECK #####
# Function to check if observed positive tests align with model output

post.check = function(df=post_seq, total, pos_obs, prop="Pa", label="") {
  s = data.frame(samples = rbinom(5e5, total, df[, prop])) 
  
  s |> 
    ggplot()+
    geom_histogram(aes(x=samples, fill="Expected positive samples"), lwd=0.5, binwidth = 1)+
    geom_vline(aes(xintercept = pos_obs, color="Observed positive samples"), lty=2, lwd=.75)+
    scale_color_manual(values = "black")+
    theme_classic()+theme(element_blank())+
    theme(axis.text.y = element_blank()) +  
    labs(fill="", color="",
         title = "Posterior check - observed and expected positive samples",
         subtitle = label,
         x = "# positive samples", y = "",
         caption = paste("Probability of direction:", round(bayestestR::pd(s$samples-pos_obs)$pd,3)))
}

# Check observed vs. expected positive samples for all samples tested in the 1st test
post.check(df=post_seq, total=data_list$N, pos_obs = data_list$npos, label = "All samples tested - 1st test")

# Check observed vs. expected positive samples for oositive samples retested in the 2nd test
post.check(df=post_seq, total=data_list$N2, pos_obs = data_list$npos2, prop = "Pa2", label = "Positive samples retested - 2nd test")

