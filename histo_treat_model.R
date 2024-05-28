### Histoplasmosis treatment model ###

L = dim(post_seq)[1] #length of vector parameters

PopN = 1e6
pos1 = rbinom(L, PopN, post_seq$Pa)
pos2 = rbinom(L, pos1 , post_seq$Pa2)
neg = pos1-pos2


p = list(
  Tr = rbeta(L, 90, 10),  #proprotion treated
  E = rbeta(L, 45, 5),  # Treatment effect
  A = rbeta(L, 6, 44),  # Risk of adverse events
  Ma = rbeta(L, 12.5, 37.5), # Risk of death (AHD)
  Mh = rbeta(L, 15, 35) # Risk of death (histoplasmosis)
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
 
  output = data.frame(Positive = numeric(L), Treated = numeric(L), 
                      Overtreat = numeric(L), Undertreat = numeric(L),
                      Adverse_events = numeric(L), 
                      Alive = numeric(L), Died = numeric(L))
  
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

plotPost(T1$Treat_prop - T2$Treat_prop, compVal = 0)

plotPost(T1$overtreat_AE_prop-T2$overtreat_AE_prop, compVal = 0) # Risk difference
plotPost(1/(T1$overtreat_AE_prop-T2$overtreat_AE_prop), xlab="NNH") # NNH - number need to harm


### mortality attributable to underdiagnosis

plotPost(T2$undertreat_mort_rate/T1$undertreat_mort_rate, compVal = 1)
plotPost(1/(T2$undertreat_mort_prop - T1$undertreat_mort_prop), xlab="NNH - no treatment")


### mortality

### Cost

cost_first_test = 15
cost_second_test = 30
cost_treatment = 240


cost2 = (T1$Positive * cost_second_test + T1$Total * cost_first_test + T2$Treated * cost_treatment) / T1$Total
cost1 = (T1$Total * cost_first_test + T1$Treated * cost_treatment) / T1$Total

plotPost(cost2 - cost1, compVal = 0)
