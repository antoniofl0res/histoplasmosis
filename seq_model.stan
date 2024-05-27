data{
     int N;
     int npos;
     int N2;
     int npos2;
     int N3;
     int npos3;

}
parameters{
     real<lower=0,upper=1> Pt;
     real<lower=0,upper=1> Se;
     real<lower=0,upper=1> Sp;
     real<lower=0,upper=1> Se2;
     real<lower=0,upper=1> Sp2;
     real<lower=0> alphaSe;
     real<lower=0> betaSe;
     real<lower=0> alphaSp;
     real<lower=0> betaSp;
     real<lower=0> alphaSe2;
     real<lower=0> betaSe2;
     real<lower=0> alphaSp2;
     real<lower=0> betaSp2;
     
     
    

}
model{
     real Pa;

     
    Pt ~ beta( 1 , 1 );
    
    // Sp ~ beta( 771.0 , 42.0); // IMMY EIA
    // Se ~ beta( 156 , 10 );
    alphaSe ~ gamma(182.0/4, 1.0/4);
    betaSe ~ gamma(16.0/4, 1.0/4);
    Se ~ beta( alphaSe , betaSe );
    
    alphaSp ~ gamma(671.0/4, 1.0/4);
    betaSp ~ gamma(37.0/4, 1.0/4);
    Sp ~ beta( alphaSp , betaSp ); // MiraVista LFA
    
    alphaSe2 ~ gamma(204.0/3, 1.0/3);
    betaSe2 ~ gamma(20.0/3, 1.0/3);
    Se2 ~ beta( alphaSe2 , betaSe2 );
    
    alphaSp2 ~ gamma(956.0/3, 1.0/3);
    betaSp2 ~ gamma(44.0/3, 1.0/3);
    Sp2 ~ beta( alphaSp2 , betaSp2 ); // IMMY EIA
    
    Pa = Pt * Se + (1 - Pt) * (1 - Sp);
    npos ~ binomial( N , Pa);
    
    
    
    //TP2 = PPV * Se2;
    //FP2 = (1 - PPV) * (1 - Sp2);
    //TN2 = (1 - PPV) * Sp2;
    //FN2 = PPV * (1 - Se2);
    
   
    
    npos2 ~ binomial( N2 , (Pt*Se*Se2) / Pa + ((1-Pt)*(1-Sp)*(1-Sp2)) / Pa);
    
    //TP3 = (1 - NPV) * Se2;
    //FP3 = (NPV) * (1 - Sp2);
    //TN3 = (NPV) * Sp2;
    //FN3 = (1 - NPV) * (1 - Se2);
    
    
    
    npos3 ~ binomial( N3 , (Pt*(1-Se)*Se2) / (1-Pa) + ((1-Pt)*Sp*(1-Sp2)) / (1-Pa));
   
}
generated quantities{
     real Pa;
 
     real PPV;
     real NPV;
     real PPV2_pos;
     real NPV2_pos;
     real PPV2_neg;
     real NPV2_neg;
     real TP;
     real FP;
     real TN;
     real FN;
     real Pa2;
     real FNpos;
     real TNpos;
     real FPpos;
     real TPpos;
     real Pa3;
     real FNneg;
     real TNneg;
     real FPneg;
     real TPneg;

     
     
    FN = Pt * (1 - Se);
    TN = (1 - Pt) * Sp;
    FP = (1 - Pt) * (1 - Sp);
    TP = Pt * Se;
    
    //Pa = inv_logit(log((TP + FP)/(TN + FN)));
    Pa = TP + FP;
    
    
    PPV = inv_logit(log(TP/FP));
    NPV = inv_logit(log(TN/FN));
    
    
    // Negative samples retested
    //TP3 = (1 - NPV) * Se2; // Pt*(1-Se)*Se2 // +-+
    //FP3 = (NPV) * (1 - Sp2); // (1-Pt)*Sp*(1-Sp2) // --+
    //TN3 = (NPV) * Sp2; // (1-Pt)*Sp*Sp2 // ---
    //FN3 = (1 - NPV) * (1 - Se2); // Pt*(1-Se)*(1-Se2) // +--
    
    // Negative samples retested
    TPneg = (Pt*(1-Se)*Se2) / (1-Pa); // +-+
    FPneg = ((1-Pt)*Sp*(1-Sp2)) / (1-Pa); // --+
    TNneg = ((1-Pt)*Sp*Sp2) / (1-Pa); // ---
    FNneg = (Pt*(1-Se)*(1-Se2)) / (1-Pa); // +--
    
    
    
    //Positive samples retested
    //TP2 = PPV * Se2; // Pt*Se*Se2 // +++
    //FP2 = (1 - PPV) * (1 - Sp2); //(1-Pt)*(1-Sp)*(1-Sp2) // -++
    //TN2 = (1 - PPV) * Sp2; // (1-Pt)*(1-Sp)*Sp2 // -+-
    //FN2 = PPV * (1 - Se2); // Pt*Se*(1-Se2) // ++-
    
    //Positive samples retested
    TPpos = (Pt*Se*Se2) / Pa; // +++
    FPpos = ((1-Pt)*(1-Sp)*(1-Sp2)) / Pa; // -++
    TNpos = ((1-Pt)*(1-Sp)*Sp2) / Pa; // -+-
    FNpos = (Pt*Se*(1-Se2))/ Pa; // ++-
    
    

    PPV2_pos = inv_logit(log(TPpos/FPpos));
    NPV2_pos = inv_logit(log(TNpos/FNpos));
    PPV2_neg = inv_logit(log(TPneg/FPneg));
    NPV2_neg = inv_logit(log(TNneg/FNneg));
    
    Pa3 = inv_logit(log((TPneg + FPneg)/(TNneg + FNneg)));
   
    Pa2 = inv_logit(log((TPpos + FPpos)/(TNpos + FNpos)));
    

   
}

