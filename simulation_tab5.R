# Simulation for Table 5
# Run the sequence below to get the data for table 5.

distn='Normal';n=45;m=60;p=300;DEP='IND'
source('simulation_tab5_helper.R') 

distn='Normal';n=45;m=60;p=300;DEP='SD'
source('simulation_tab5_helper.R') 

distn='Cauchy';n=45;m=60;p=300;DEP='IND'
source('simulation_tab5_helper.R') 

distn='Cauchy';n=45;m=60;p=300;DEP='SD'
source('simulation_tab5_helper.R') 

# Make the table
t1 <- read.table('mixed_IND_n45_m60_p300_Normal.txt')
t2 <- read.table('mixed_IND_n45_m60_p300_Cauchy.txt')
t3 <- read.table('mixed_SD_n45_m60_p300_Normal.txt')
t4 <- read.table('mixed_SD_n45_m60_p300_Cauchy.txt')
pick <- seq(1, 11, 2)
ind_case <- cbind(t1, t2)[pick,]
sd_case <- cbind(t3, t4)[pick,]
row.names(ind_case) <- c('0', '0.2', '0.4', '0.6', '0.8', '1')
xtable::xtable(ind_case, digit = 3)
row.names(sd_case) <- c('0', '0.2', '0.4', '0.6', '0.8', '1')
xtable::xtable(sd_case, digit = 3)

t1 <- read.table('positive_IND_n45_m60_p300_Normal.txt')
t2 <- read.table('positive_IND_n45_m60_p300_Cauchy.txt')
t3 <- read.table('positive_SD_n45_m60_p300_Normal.txt')
t4 <- read.table('positive_SD_n45_m60_p300_Cauchy.txt')
pick <- seq(1, 11, 2)
ind_case <- cbind(t1, t2)[pick,]
sd_case <- cbind(t3, t4)[pick,]
row.names(ind_case) <- c('0', '0.2', '0.4', '0.6', '0.8', '1')
xtable::xtable(ind_case, digit = 3)
row.names(sd_case) <- c('0', '0.2', '0.4', '0.6', '0.8', '1')
xtable::xtable(sd_case, digit = 3)
