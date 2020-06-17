#
# Simulation for GCT, GCTm, ZWL, ZWLm
#
library(highD2pop)
library(highDmean)
# Innovations
if (distn == 'Normal') {
    INNOV = function(n,...) return(rnorm(n, 0, 1))
    Signal = 1/8
} else if (distn == 'Gamma') {
    INNOV = function(n,...) return(rgammashift(n,4,2)) 
    Signal = 0.5
} else if (distn == 'Cauchy') {
    INNOV = function(n,...) return(rcauchy(n,0,0.1)) 
    Signal = 1
}

# Proportion of contamination
betaList = seq(0, 1, 0.1)
# Initialzie a storage for all beta values
rejProb = NULL
# Create a file to store the temp result
cat('BETA, pvMy,  pvGCT,  pvMyM,  pvGCTM,', '\n',
    file = paste0('Temp',DEP,'_n',n,'_m', m, '_p',p, 
                  '_',distn,'.txt'), append = F) 

# Loop on portion of contamination 0, 0.1... 1
for (BETA in betaList ){
    # Initalzie storage for 
    Pvalue_My = NULL; Pvalue_GCT = NULL
    Pvalue_MyM = NULL;  Pvalue_GCTM = NULL
    # Generate datasets
    set.seed(BETA*100)
    # Inner loop
    for (k in 1:400){ ## control the number of runs
        DATA <-buildData(n = n,m = m,p = p,
                         muX = rep(0,p),
                         muY = c(rep(Signal, round(p*BETA)), 
                                 rep(0, p - round(p*BETA))),
                         dep = DEP, S = 5, innov = INNOV)
        # Keep the pvalues
        pvMy = zwl_sim(DATA, order = 2 )$pvalue
        pvGCT = GCT.sim(DATA,r=10,smoother="parzen",ntoorderminus = 2)$pvalues
        pvMyM = zwl_sim(DATA, order = 0 )$pvalue
        pvGCTM = GCT.sim(DATA,r=10,smoother="parzen",ntoorderminus = 0)$pvalues
        # Temprorily store the p-values.
        for (q in 1:5){
            cat(BETA, pvMy[q],  pvGCT[q],  pvMyM[q],  pvGCTM[q], '\n',
                file = paste0('Temp',DEP,'_n',n,'_m', m, '_p',p, '_',distn,'.txt'),
                append = T)  
        }
        
        Pvalue_GCT = c(Pvalue_GCT, pvGCT)
        Pvalue_My = c(Pvalue_My, pvMy)
        Pvalue_MyM = c(Pvalue_MyM, pvMyM)
        Pvalue_GCTM = c(Pvalue_GCTM, pvGCTM)
    }
    rejProb = rbind(rejProb, 
                    apply(cbind(Pvalue_My, Pvalue_GCT, Pvalue_MyM, Pvalue_GCTM),
                          2, function(x) mean(x<0.05, rm.na = T)))
}
colnames(rejProb) = c('My', 'GCT', 'MyM', 'GCTM')
row.names(rejProb) = paste0('b=', round(betaList, 2))
write.table(rejProb, file = paste0(DEP,'_n',n,'_m', m, '_p',p, '_',distn,'.txt'))
