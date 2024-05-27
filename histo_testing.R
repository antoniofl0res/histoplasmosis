### Infering prevalence from a POC test ###

library(dplyr)
library(tidyr)
library(HDInterval)
library(rethinking)
library(bayesboot)
library(gt)

library(gtsummary)

cmdstanr::cmdstan_path()
set_cmdstan_path("C:/Users/Antonio/Documents/.cmdstan/cmdstan-2.32.2")
Sys.getenv()
Sys.setenv(HOME="C:/Users/Antonio/Documents")
#Sys.setenv(HOME="C:/Users/Antonio/OneDrive - MSF/Documents")

getwd()


#### DATA ####

# Define sample sizes and positive samples for initial and subsequent tests
N = 200  # Total samples
npos = 15  # Positive samples initially

# Initialize variables for second test
N2 = 14  # Number of initially positive samples tested on second test
npos2 = 11  # Positive samples found on second test

N3 = 180  # Number of initially negative samples tested on second test
npos3 = 2  # Positive samples found on second test

### DATA LIST #####

# Store data in a list for passing to Stan models
d = list(
  N = N,
  npos = npos,
  N2 = N2,
  npos2 = npos2,
  N3 = N3,
  npos3 = npos3
)

##### STAN MODEL - ONE TEST ####
# Fit the Stan model for one test using the provided data
hist.mod = stan("mod_one_test.stan", data=d, 
                chains=4, cores=4, iter=7.5e3, warmup=2.5e3)

##### MODEL OUTPUT ####
# Display summary statistics and plots for the model output
precis(hist.mod, prob = .95, depth=3)
plot(precis(hist.mod, prob = .95))

# Uncomment to display traceplot and trankplot
#traceplot(hist.mod)
#trankplot(hist.mod)

# Extract posterior samples from the model
post_hist <- extract.samples(hist.mod) |> as.data.frame()

# Summarize key statistics from the posterior samples
list(
  N = c(Total=N, pos1=d$npos),
  expected_PPV = round(mean(post_hist$PPV)*100, 1), 
  hdi_PPV = round(hdi(post_hist$PPV)*100, 1), 
  prob_50 = mean(post_hist$PPV>.5)
)

# Prepare data for display
post_hist = post_hist[, 1:ncol(post_hist)-1]

# Generate a table for output display
(output = round(as.data.frame(t(rbind(expected=colMeans(post_hist), 
                                      HDInterval::hdi(post_hist)))),3))

data.frame(cbind(parameters=colnames(post_hist), output)) |> gt()

### STAN SEQ MODEL ####

# Fit the sequential Stan model using the provided data
seq_model <- stan(file="seq_model.stan", data=d, chains=4, cores=4, 
                  iter=1.2e4, warmup=2e3)

# Uncomment to display traceplot
#traceplot(seq_model)

# Display summary statistics for the sequential model
precis(seq_model, prob=.95)

# Extract posterior samples from the sequential model
post_seq <- rstan::extract(seq_model) |> as.data.frame()

# Summarize key statistics from the posterior samples
list(
  N = c(Total=N, pos1=d$npos, restested = N2, pos2=d$npos2),
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
post.check = function(df=post_seq, total=N, pos_obs=npos, prop="Pa", label="") {
  s = data.frame(samples = rbinom(1e6, total, df[, prop])) 
  
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
         caption = paste("Probability of direction:", round(bayestestR::pd(s$samples-pos_obs),3)))
}

# Check observed vs. expected positive samples for all samples tested in the 1st test
post.check(df=post_seq, total=N, pos_obs = npos, label = "All samples tested - 1st test")

# Check observed vs. expected positive samples for samples retested in the 2nd test
post.check(df=post_seq, total=N2, pos_obs = npos2, prop = "Pa2", label = "Positive samples retested - 2nd test")
