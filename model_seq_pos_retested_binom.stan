data {
  int N;
  int npos;
  int N2;
  int npos2;
}

parameters {
  real<lower=0, upper=1> Pt;
  real<lower=0, upper=1> Se;
  real<lower=0, upper=1> Sp;
  real<lower=0, upper=1> Se2;
  real<lower=0, upper=1> Sp2;
  real<lower=0> alphaSe;
  real<lower=0> betaSe;
  real<lower=0> alphaSp;
  real<lower=0> betaSp;
  real<lower=0> alphaSe2;
  real<lower=0> betaSe2;
  real<lower=0> alphaSp2;
  real<lower=0> betaSp2;
}

transformed parameters {
  real Pa;
  Pa = Pt * Se + (1 - Pt) * (1 - Sp);
}

model {
  // Priors for Pt
  Pt ~ beta(1, 1);

  // Priors for Se and Sp for the first test - MiraVista LFA
  alphaSe ~ gamma(198.0 / 5, 1.0 / 5);
  betaSe ~ gamma(17.0 / 5, 1.0 / 5);
  Se ~ beta(alphaSe, betaSe);

  alphaSp ~ gamma(702.0 / 5, 1.0 / 5);
  betaSp ~ gamma(37.0 / 5, 1.0 / 5);
  Sp ~ beta(alphaSp, betaSp);

  // Priors for Se2 and Sp2 for the second test - IMMY EIA
  alphaSe2 ~ gamma(220.0 / 4, 1.0 / 4);
  betaSe2 ~ gamma(21.0 / 4, 1.0 / 4);
  Se2 ~ beta(alphaSe2, betaSe2);

  alphaSp2 ~ gamma(987.0 / 4, 1.0 / 4);
  betaSp2 ~ gamma(44.0 / 4, 1.0 / 4);
  Sp2 ~ beta(alphaSp2, betaSp2);

  // Likelihood for first test
  npos ~ binomial(N, Pa);

  // Calculate PPV for positive samples
  real TP = Pt * Se;
  real FP = (1 - Pt) * (1 - Sp);
  real PPV = TP / (TP + FP);

  // Likelihood for second test given positive first test
  real Pa2 = PPV * Se2 + (1 - PPV) * (1 - Sp2);
  npos2 ~ binomial(N2, Pa2);
}

generated quantities {
  // True positives, false positives, true negatives, false negatives for the first test
  real TP = Pt * Se;
  real FP = (1 - Pt) * (1 - Sp);
  real TN = (1 - Pt) * Sp;
  real FN = Pt * (1 - Se);
  
  // Positive Predictive Value for the first test
  real PPV = TP / (TP + FP);
  
  // Second test positivity probability given positive first test
  real Pa2 = PPV * Se2 + (1 - PPV) * (1 - Sp2);
  
  // True positives and false positives for the second test
  real TPpos = PPV * Se2;
  real FPpos = (1 - PPV) * (1 - Sp2);
  
  // Positive Predictive Value for the second test given positive first test
  real PPV2_pos = TPpos / (TPpos + FPpos);
}

