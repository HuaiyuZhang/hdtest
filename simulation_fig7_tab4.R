# ----------------------------------------------------------
# Factorial Design Simulation
# ----------------------------------------------------------
# Four tests (ZWLm, GCTm, SKK, CQ)
# Four innovation distributions (Normal, Gamma, t3, Cauchy)
# Four dependence structures (IND, LR, WD, SD)
#
# ----------------------------------------------------------
# ----------------------------------------------------------
# Part 1: Run the simulation
# ----------------------------------------------------------
#
if (!requireNamespace("highD2pop")) install.packages("highD2pop")
library(highD2pop)
if (!requireNamespace("highDmean")) install.packages("highDmean")
library(highDmean)


run_simulation <- function(
    # Specify the innovation distribution.
    # Possible choices: 'Normal', 'Gamma', 't3', 'Cauchy'.
    distn,
    # Specify the dependence structure.
    # Possible choices: 'IND', 'LR', 'SD', 'WD'
    DEP,
    # Specify the tests included
    test_type_control,
    # Specify the sample sizes for the two populations and the dimensionality.
    n, m, p)
{
    # Innovations
    if (distn == 'Normal') {
        INNOV <- function(n,...) return(rnorm(n, 0, 1))
    } else if (distn == 'Gamma') {
        INNOV <- function(n,...) return(rgammashift(n,4,2))
    } else if (distn == 'Cauchy') {
        INNOV <- function(n,...) return(rcauchy(n,0,0.1))
    } else if (distn == 't3') {
        INNOV <- function(n,...) return(rt(n,3))
    }
    # Signal strenth lookup table
    signal_strenth_table <- matrix(c(1/8, 1/8, 1/8, 3/8,
                                     1/2, 1/2, 1/2, 3/2,
                                     1, 1, 1, 3,
                                     0.2, 0.25, 0.25, 1),
                                   nrow = 4, ncol = 4, byrow = T,
                                   dimnames = list(c('Normal', 'Gamma', 'Cauchy', 't3'),
                                                   c('IND', 'LR', 'WD', 'SD')))
    # Determine the signal strenth for a specific distribution*dependency combination.
    signal <- signal_strenth_table[distn, DEP]
    # Proportion of contamination
    betaList <- seq(0, 1, 0.1)
    # Initialzie a storage of p-values for all beta values
    rejProb = NULL
    # Create a file to store the temp result
    cat('BETA, pvZWLm,  pvGCTm,  pvCQ,  pvSKK,', '\n',
        file = paste0('TEMP_',DEP, '_',distn,'.txt'), append = F)
    # Loop on portion of contamination 0, 0.1... 1
    for (BETA in betaList ){
        # Initalzie storage for pvalues
        Pvalue_ZWLm <- Pvalue_GCTm <- Pvalue_CQ <- Pvalue_SKK <- NULL
        # Generate datasets
        set.seed(BETA*2000)
        # Inner loop: 400*5=2000 runs
        for (k in 1:400){ ## control the number of runs
            DATA <- buildData(n = n,m = m,p = p,
                              muX = rep(0,p),
                              muY = c(rep(signal, round(p*BETA)),
                                      rep(0, p - round(p*BETA))),
                              dep = DEP, S = 5, innov = INNOV)
            # Initialize collectors for pvalue for temporary output.
            pvZWLm <- pvGCTm <- pvCQ <- pvSKK <- rep(NA, 5)
            if ('ZWLm' %in% test_type_control){
                pvZWLm <- zwl_sim(DATA, order = 0)$pvalue
            }
            if ('GCTm' %in% test_type_control){
                pvGCTm <-  GCT.sim(DATA, r=10, smoother="parzen", ntoorderminus = 0)$pvalues
            }
            if ('CQ' %in% test_type_control){
                pvCQ <-  ChenQin.sim(DATA)[,2]
            }
            if ('SKK' %in% test_type_control){
                pvSKK <-  SKK_sim(DATA)$pvalue
            }
            # Temprorily store the p-values.
            for (q in 1:5){
                cat(BETA, pvZWLm[q],  pvGCTm[q],  pvCQ[q],  pvSKK[q], '\n',
                    file = paste0('TEMP_',DEP, '_',distn,'.txt'),
                    append = T)
            }
            Pvalue_GCTm <- c(Pvalue_GCTm, pvGCTm)
            Pvalue_ZWLm <- c(Pvalue_ZWLm, pvZWLm)
            Pvalue_CQ  <- c(Pvalue_CQ, pvCQ)
            Pvalue_SKK <- c(Pvalue_SKK, pvSKK)
        }
        rejProb <- rbind(rejProb,
                         apply(cbind(Pvalue_ZWLm, Pvalue_GCTm, Pvalue_CQ, Pvalue_SKK),
                               2, function(x) mean(x<0.05, rm.na = T)))
    }
    colnames(rejProb) <- c('ZWLm', 'GCTm', 'CQ', 'SKK')
    row.names(rejProb) <- paste0('b=', round(betaList, 2))
    # Output the rejection proportions
    write.table(rejProb, file = paste0(DEP,'_',distn,'.txt'))
}
# # For example, the Cauchy*SD setting.
# run_simulation(distn = 'Cauchy', DEP = 'SD',
#                test_type_control = c('ZWLm', 'GCTm', 'CQ', 'SKK'),
#                n = 45, m = 60, p = 300)
# For the complete result corresponding to Fig7 in Zhang and Wang (2020),
# run the double loop below
for (distribution in c('Normal', 'Gamma', 'Cauchy', 't3')){
    for (dependency in c('IND', 'LR', 'WD', 'SD')){
        run_simulation(distn = distribution,
                       DEP = dependency,
                       test_type_control = c('ZWLm', 'GCTm', 'CQ', 'SKK'),
                       n = 45, m = 60, p = 300)
    }
}
# ----------------------------------------------------------
# ----------------------------------------------------------
# Part 2: FIGURE 7: Plot the rejection proportions
# ----------------------------------------------------------
#
# Caution: Run below code if the results for all distribution*dependency
#          combinations are ready.
if (!requireNamespace("reshape2")) install.packages("reshape2")
if (!requireNamespace("ggplot2")) install.packages("ggplot2")
library(reshape2)
library(ggplot2)
AllRes <- NULL
distlist <- c('Normal', 'Gamma', 't3','Cauchy' )
deplist <- c('IND', 'LR', 'WD', 'SD')
for (DEP in deplist){
    for (dist in distlist){
        fn <- paste0(DEP,'_',dist,'.txt')
        ## Organize the data: wide to long
        res <- read.table(fn, sep = ' ')
        names(res) <- c('ZWLm','GCTm', 'CQ', 'SKK')
        t1e <- paste0('Size: ', res[1,1]*100, '%(ZWLm), ',
                      res[1,2]*100,'%(GCTm), ',
                      res[1,3]*100,'%(CQ), ',
                      res[1,4]*100,'%(SKK)')
        res2 <- melt(res)
        res2$beta <- paste0(seq(from = 0,to =1, by = 0.1))
        res2$dep <- rep(DEP, nrow(res2))
        res2$t1e <- rep(t1e, nrow(res2))
        res2$dist <- dist
        names(res2) <- c('type', 'power', 'beta', 'dep', 't1e', 'dist')
        AllRes <- rbind(AllRes, res2)
    }
}
# rank the data frame to make the desired display order
AllRes$dep <- factor(AllRes$dep, levels=c('IND',  "LR", "WD", "SD"), ordered=T)
depen <- levels(AllRes$dep)
AllRes$depdist <- factor(paste(AllRes$dep, AllRes$dist),
                         levels=c(outer(depen, distlist, paste)), ordered=T)
# Plot the rejection proportions
ggplot(AllRes,
       aes(x = beta, y = power,
           group = type, colour = type,
           linetype = type) ) +
    scale_color_manual(values=c('blue', 'red', 'green3', 'orange3'))+
    geom_line(size = 0.8) +
    geom_hline(yintercept = 0.05, linetype = 1, colour = 'grey') +
    facet_wrap(~depdist) +
    theme_minimal(base_size = 10) +
    labs( x = expression('Proportion of signal'~ beta),
          y = expression("Proportion of Rejections at level"~ alpha~"=5%"),
          colour='Test: ',
          linetype='Test: ') +
    theme(legend.position = "bottom")+
    theme(axis.text.x = element_text(angle = 90))
# Output the plot
ggsave('Fig7.pdf', width=8, height=8)

# ----------------------------------------------------------
# ----------------------------------------------------------
# Part 3: Make the Table 4 for type I error rates
# ----------------------------------------------------------
#
# Caution: Run below code if the results for all distribution*dependency
#          combinations are ready.
typeI <- AllRes[AllRes$beta<0.01, c('type', 'power', 'dep', 'dist')]
typeI <- typeI[order(typeI$dep),]
rowlabel <- NULL
typeIerr <- NULL
for (type in levels(AllRes$type)){
    myrow <- data.frame()
    rowlab <- NULL
    for (dist in distlist){
        thismany  <-  typeI[(type==typeI$type)&(typeI$dist == dist), ]
        rowlab  <-  rbind(rowlab, thismany[, c('dep', 'dist', 'type')] )
        myrow  <-  c(myrow, thismany$power*100)
    }
    rowlabel <- rbind(rowlabel, t(rowlab))
    typeIerr <- rbind(typeIerr, myrow)
}
# Output latex table
if (!requireNamespace("xtable")) install.packages("xtable")
xtable::xtable(typeIerr, label='typeI', caption='typeI')
# xtable::xtable(rowlabel[(1:4)*3,])
# xtable::xtable(rowlabel[(0:3)*3+1,] )
# xtable::xtable(rowlabel[(0:3)*3+2,])
# END OF CODE
