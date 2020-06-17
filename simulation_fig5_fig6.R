# Part 1: run the simulation and save results

# The simulation for figure 5
distn='Normal';n=45;m=60;p=300;DEP='IND'
source('simulation_fig5_fig6_helper.R') 
distn='Normal';n=45;m=60;p=300;DEP='LR'
source('simulation_fig5_fig6_helper.R') 
distn='Normal';n=45;m=60;p=300;DEP='SD'
source('simulation_fig5_fig6_helper.R') 
distn='Normal';n=45;m=60;p=300;DEP='WD'
source('simulation_fig5_fig6_helper.R') 

# The simulation for figure 6
distn='Normal';n=90;m=120;p=300;DEP='IND'
source('simulation_fig5_fig6_helper.R') 
distn='Normal';n=90;m=120;p=300;DEP='LR'
source('simulation_fig5_fig6_helper.R') 
distn='Normal';n=90;m=120;p=300;DEP='SD'
source('simulation_fig5_fig6_helper.R') 
distn='Normal';n=90;m=120;p=300;DEP='WD'
source('simulation_fig5_fig6_helper.R') 

# ---------------------------------------------
# Part 2:  Plot the results for figure 5 and 6
library(reshape2)
library(ggplot2)
sampleSize = 'n45' # Figure5
sampleSize = 'n90' # Figure6

if (sampleSize == 'n45'){
    fnlist <- dir(pattern = 'n45_m60_p300_Normal.txt$')
    pdf('fig5.pdf')
} else {
    fnlist <- dir(pattern = 'n90_m120_p300_Normal.txt$')
    pdf('fig6.pdf')
}

## Organize the data
AllRes <- NULL
for (fn in fnlist)
{
    DEP <- strsplit(fn, '_')[[1]][1]
    res <- read.table(fn, sep = ' ')
    names(res) <- c('ZWL','GCT', 'ZWLm', 'GCTm')
    t1e <- paste0('Size: ', res[1,1], '(ZWL),', 
                  res[1,2],'(GCT),',
                  res[1,3],'(ZWLm),',
                  res[1,4],'(GCTm).')
    res2 = melt(res)
    res2$beta = paste0(seq(from = 0,to =1, by = 0.1))
    res2$dep = rep(DEP, nrow(res2))
    res2$t1e = rep(t1e, nrow(res2))
    names(res2) = c('type', 'power', 'beta', 'dep', 't1e')
    AllRes <- rbind(AllRes, res2)
}
# Plot
p <-  ggplot(AllRes, aes(x = beta, y = power, 
                         group = type, colour = type, 
                         linetype = type) ) +
    geom_line(size = 0.8) +
    geom_hline(yintercept = 0.05, linetype = 1, colour = 'grey') +
    annotate('text', label = paste("alpha == ", .05), 
             x = 10, y = 0.05, size = 3, parse = TRUE) +
    facet_wrap(~dep + t1e, 
               labeller = label_value) +
    theme_minimal(base_size = 10) +
    labs( x = "Proportion of signal", y = "Proportion of Rejections") +
    theme(legend.position = "bottom")
print(p)          
dev.off()
