## Simulation for Table 5


# Part 1: Helper functions
#
#+-----------------Permutation test function--------------+
# Input: X: n by p matrix, Y: m by p matrix
# Output: p_val: vector of length 6 showing p-values
permTest <- function(X, Y, B = 1000, 
                     alternative = 'two.sided')
{
    n = nrow(X)
    m = nrow(Y)
    p = ncol(X)
    # Bind the data, compute rank
    bind_data = rbind(X, Y)  # n+m by p
    # compute test stat
    observed_stat_vec <- computeTestStat(X, Y)
    # Bootstrap for the pvalues
    stat_mat <- matrix(nrow = B, ncol = length(observed_stat_vec))
    
    for (b in 1:B){
        # indices
        index <- sample(1:(m+n), m+n, replace = T)
        index1 <- index[1:n]
        index2 <- index[-(1:n)]
        # regular version
        X_b <- bind_data[index1,] 
        Y_b <- bind_data[index2,]
        stat_mat[b, ] <- computeTestStat(X_b, Y_b, alternative)
        
        # In fact we only consider two-sided cases 
        if (alternative == 'less'){
            # rejection for lower-tailed
            rejection_mat <- t(t(stat_mat) <= observed_stat_vec)
            p_val <-  colSums(rejection_mat)/B
        }
        else if(alternative == 'greater') {
            # rejection for upper-tailed
            rejection_mat <- t(t(stat_mat) >= observed_stat_vec)
            p_val <-  colSums(rejection_mat)/B
        }
        else{
            # rejection for two-sided
            rejection_mat <- t(t(stat_mat) >= abs(observed_stat_vec)) | 
                t(t(stat_mat) <= -abs(observed_stat_vec))
            p_val <-  colSums(rejection_mat)/B
        }
    }
    p_val
}
#+---------- compute test stat ---------+
# X: n by p matrix
# Output: length-3 vector containing test stats
computeTestStat <- function(X, Y, alternative='two.sided'){
    p = ncol(X)
    n = nrow(X)
    m = nrow(Y)
    
    # Optimized computing t vector
    mean1=drop(matrix(rep(1, n), nrow=1)%*% X/n)
    mean2=drop(matrix(rep(1, m), nrow=1)%*% Y/m)
    var1=drop(matrix(rep(1, n), nrow=1) %*% scale(X, mean1, F)^2/(n-1))
    var2=drop(matrix(rep(1, m), nrow=1) %*% scale(Y, mean2, F)^2/(m-1))
    Sp=sqrt(var1/n+var2/m) 
    t_vector = (mean1-mean2)/Sp
    
    # statistics
    t_sum <- sum(t_vector)  
    t_sq_sum <- as.numeric(t_vector%*%t_vector)
    t_max <- t_vector[which.max(abs(t_vector))]
    t_abs <- sum(abs(t_vector)) 
    # output
    c(t_sum, t_sq_sum, t_max, t_abs)
}  

# Part 2: Run simulation------------------------------------
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
# Save pvalues in a table
number_of_tests = 5
beta_pval_table = matrix(nrow = length(betaList), ncol = number_of_tests) # pv1, pv2..
# Create an empty file to store the temp result
cat('BETA, pvZWL,  pvTsum,  pvTsqsum, pvTmax, pvTabs', '\n',
    file = paste0('Temp_mixed_',DEP,'_n',n,'_m', m, '_p',p,
                  '_',distn,'.txt'), append = F)
number_runs = 2000

# Part 2.1: Mixed signal
# Iterate over the portion of contamination 0, 0.1... 1
for (i in 1:length(betaList) ){
    BETA = betaList[i]
    p_val_mat <- matrix(nrow = number_runs, ncol = number_of_tests )
    set.seed(BETA*100)
    # Inner loop
    for (k in 1:(number_runs/5)){  # 5 datasets produced each time
        # cat(k,' ')
        # DATA: list of 5
        DATA <-buildData(n = n,m = m,p = p,
                         muX = rep(0,p),
                         muY = c(rep(Signal, round(p*BETA)/2), 
                                 rep(-Signal, round(p*BETA)/2), 
                                 rep(0, p - round(p*BETA))),
                         dep = DEP, S = 5, innov = INNOV)
        
        # store the p-values for each run
        for (l in 1:5){
            p_val_per_run  <- c(zwl_test(DATA[[l]]$X, DATA[[l]]$Y)$pvalue,
                                permTest(DATA[[l]]$X, DATA[[l]]$Y))
            p_val_mat[(k-1)*5+l,] <- p_val_per_run
            # Output the pvalues to external file
            cat(BETA, p_val_per_run, '\n',
                file = paste0('Temp_mixed_',DEP,'_n',n,'_m', m, '_p',p, '_',distn,'.txt'),
                append = T)
        }
    }
    beta_pval_table[i,] = apply(p_val_mat, 2, function(x) mean(x<0.05, rm.na = T))
}
colnames(beta_pval_table) = c('ZWL', 'tsum',  'tsqsum',  'tmax', 'tabs')
row.names(beta_pval_table) = paste0('b=', round(betaList, 2))
write.table(beta_pval_table, file = paste0('mixed_',DEP,'_n',n,'_m', m, '_p',p, '_',distn,'.txt'))

# Part 2.2 Pure positive signal
beta_pval_table = matrix(nrow = length(betaList), ncol = number_of_tests) # pv1, pv2..
# Create an empty file to store the temp result
cat('BETA, pvZWL,  pvTsum,  pvTsqsum, pvTmax, pvTabs', '\n',
    file = paste0('Temp_positive_',DEP,'_n',n,'_m', m, '_p',p,
                  '_',distn,'.txt'), append = F)
number_runs = 2000

for (i in 1:length(betaList) ){
    BETA = betaList[i]
    p_val_mat <- matrix(nrow = number_runs, ncol = number_of_tests )
    set.seed(BETA*100)
    # Inner loop
    for (k in 1:(number_runs/5)){  # 5 datasets produced each time
        # cat(k,' ')
        # DATA: list of 5
        DATA <-buildData(n = n,m = m,p = p,
                         muX = rep(0,p),
                         muY = c(rep(Signal, round(p*BETA)), 
                                 rep(0, p - round(p*BETA))),
                         dep = DEP, S = 5, innov = INNOV)
        
        # store the p-values for each run
        for (l in 1:5){
            p_val_per_run  <- c(zwl_test(DATA[[l]]$X, DATA[[l]]$Y)$pvalue,
                                permTest(DATA[[l]]$X, DATA[[l]]$Y))
            p_val_mat[(k-1)*5+l,] <- p_val_per_run
            # Output the pvalues to external file
            cat(BETA, p_val_per_run, '\n',
                file = paste0('Temp_positive_',DEP,'_n',n,'_m', m, '_p',p, '_',distn,'.txt'),
                append = T)
        }
    }
    beta_pval_table[i,] = apply(p_val_mat, 2, function(x) mean(x<0.05, rm.na = T))
}
colnames(beta_pval_table) = c('ZWL', 'tsum',  'tsqsum',  'tmax', 'tabs')
row.names(beta_pval_table) = paste0('b=', round(betaList, 2))
write.table(beta_pval_table, file = paste0('positive_', DEP,'_n',n,'_m', m, '_p',p, '_',distn,'.txt'))




