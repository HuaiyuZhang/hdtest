# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install()  # Update packages
# BiocManager::install("topGO")
# BiocManager::install("ALL")
# BiocManager::install("hgu95av2.db")
# BiocManager::install('genefilter')

#+------------------Preprocessing--------------------------------+
# load data
library(ALL)
data(ALL)

# Preprocessing
ind1 <- (substr(ALL$BT,0,1) == 'B') & (ALL$mol.biol == 'BCR/ABL' |ALL$mol.biol == 'NEG') 
ALL1 <- ALL[,ind1]
# filtering: 25% of 79 samples is greater than log2(100)
ind2 <- apply(exprs(ALL1), 1, function(x) sum(x > log2(100)) > 79*0.25)
ALL2 <- ALL1[ind2]
# filtering: IQR is greater than 0.5
ind3 <- apply(exprs(ALL2), 1, function(x) IQR(x) > 0.5)
ALL3 <- ALL2[ind3]

all_gene_names <- rownames(exprs(ALL3))
dat <- data.frame(t(exprs(ALL3)), class=ALL3$mol.biol)
names(dat) <- c(all_gene_names, 'class')

# Match the GO terms.
library(topGO)
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)
gene_list <- rep(0, length(all_gene_names))
names(gene_list) = all_gene_names
topDiffGenes <- function (allScore) 
{
    return(allScore < 0.01)
}

# For Table 1 summary statistics
get_summary <- function(genemap){
    genecount = sapply(GO_gene_map, length)
    list(num_GO = length(genecount),
         min = min(genecount),
         max = max(genecount), 
         mean = mean(genecount),
         median = median(genecount),
         std_dev = sd(genecount))
}


#+-------------------------------------------------------+
#+--------------- Part 1: Test NEG vs BCR/ABL------------+
for (GO_type in c('BP', 'MF', 'CC')){
    GOdata <- new("topGOdata",
                  description = "Simple session", 
                  ontology = GO_type, # choose the ontology
                  allGenes = gene_list, 
                  geneSel = topDiffGenes,
                  nodeSize = 41,
                  annot = annFUN.db, 
                  affyLib = affyLib)
    GO_gene_map = genesInTerm(GOdata)
    # To get the summary statistics in Table 1
    tab1_summary <- get_summary(GO_gene_map)
    write.table(tab1_summary, paste0(GO_type,'_summary.txt'))
    # ready for test
    library(highD2pop)
    library(highDmean)
    N <- length(GO_gene_map)
    res <- matrix(0, nrow = N, ncol = 6)
    res <- data.frame(res)
    colnames(res) <- c('GOterm', 'numGene', 'ZWLm', 'GCTm', 'CQ', 'SKK')
    for (i in 1:N){
        cat(i,"~")
        X <- as.matrix(dat[dat$class=='NEG', GO_gene_map[[i]]])
        Y <- as.matrix(dat[dat$class=='BCR/ABL', GO_gene_map[[i]]])
        res$GOterm[i] <- names(GO_gene_map)[i]
        res$numGene[i] <- length(GO_gene_map[[i]])
        res[i, 3:6] <- c(zwl_test(X,Y)$pvalue, 
                         GCT.test(X,Y, r = 10, ntoorderminus = 0)$pvalue,
                         ChenQin.test(X,Y)$pvalue,
                         SKK_test(X, Y)$pvalue)
    }
    filename = paste0('between_40_node_', GO_type,'.txt')
    write.table(res, filename)
}

#+---------- Table 2: Rejection proportion-----------------+
# Set adjust method for the p-value
adjust_method = 'bonferroni'
# adjust_method = 'BY'
res <- read.table('between_40_node_BP.txt')
res[,3:6] <- apply(res[,3:6], 2, p.adjust, method = adjust_method)
reject_BP <- apply(res[,3:6], 2,function(x) mean(x<0.05))

res <- read.table('between_40_node_MF.txt')
res[,3:6] <- apply(res[,3:6], 2, p.adjust, method = adjust_method)
reject_MF <- apply(res[,3:6], 2,function(x) mean(x<0.05))

res <- read.table('between_40_node_CC.txt')
res[,3:6] <- apply(res[,3:6], 2, p.adjust, method = adjust_method)
reject_CC<- apply(res[,3:6], 2,function(x) mean(x<0.05))
xtable::xtable(round(rbind(reject_BP, reject_MF, reject_CC),3))



#+--------------------- Figure 3 -----------------------+
rm(list = ls())
library(ggplot2)
# Set adjust method for the p-value
adjust_method = 'bonferroni'
res1 <- read.table('between_40_node_BP.txt')
res1[,3:6] <- apply(res1[,3:6], 2, p.adjust, method = adjust_method)
res2 <- read.table('between_40_node_MF.txt')
res2[,3:6] <- apply(res2[,3:6], 2, p.adjust, method = adjust_method)
res3 <- read.table('between_40_node_CC.txt')
res3[,3:6] <- apply(res3[,3:6], 2, p.adjust, method = adjust_method)
dat <- rbind(res1, res2, res3)
ontology_type <- c(rep('BP', nrow(res1)), rep('MF', nrow(res2)), rep('CC', nrow(res3)))
dat <- data.frame(dat, ontology_type)
library(tidyr)
data_long <- gather(dat, test, pvalue, ZWLm:SKK, factor_key=TRUE)

p10 <- ggplot(data_long, 
              aes(x = test, y = pvalue, fill = test)) +
    geom_boxplot(alpha = 0.7,  
                 outlier.colour = "#1F3552", 
                 outlier.shape = 20,
                 width = .8, 
                 notch = TRUE, 
                 position = position_dodge(width = 0)) +
    scale_y_continuous(name = "P-values",
                       limits=c(0, 1)) +
    scale_x_discrete(name = "") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 11),
          legend.position = 'none') +
    facet_grid(. ~ ontology_type)
p10
ggsave('fig3.pdf', p10, width = 11, height = 7)


#+ ------------------------Figure 4: test within the group-----------------------
# 
for (GO_type in c('BP', 'MF', 'CC')){
    GOdata <- new("topGOdata",
                  description = "Simple session", 
                  ontology = GO_type,
                  allGenes = gene_list, 
                  geneSel = topDiffGenes,
                  nodeSize = 41,
                  annot = annFUN.db, 
                  affyLib = affyLib)
    GO_gene_map = genesInTerm(GOdata)
    # ready for test
    library(highD2pop)
    library(highDmean)
    
    N <- length(GO_gene_map) # number of GO terms
    res <- matrix(0, nrow = N, ncol = 6)
    res <- data.frame(res)
    colnames(res) <- c('GOterm', 'numGene', 'ZWLm', 'GCTm', 'CQ', 'SKK')
    for (i in 1: N){
        cat(i,"~")
        XY <- as.matrix(dat[dat$class=='NEG', GO_gene_map[[i]]])
        # 42 observations for the NEG class
        # Split the sample to two groups
        ind <- sample(1:42, size = 21, replace = F)
        X <- XY[ind, ]
        Y <- XY[-ind,]
        res$GOterm[i] <- names(GO_gene_map)[i]
        res$numGene[i] <- length(GO_gene_map[[i]])
        res[i, 3:6] <- c(zwl_test(X,Y)$pvalue,
                         GCT.test(X,Y, r = 10, ntoorderminus = 0)$pvalue,
                         ChenQin.test(X,Y)$pvalue,
                         SKK_test(X, Y)$pvalue)
    }
    # Store data
    filename = paste0('back_test_40_node_', GO_type,'.txt')
    write.table(res, filename)
}


#+----------------- Plot -------------------+
library(ggplot2)
adjust_method = 'bonferroni'
res1 <- read.table('back_test_40_node_BP.txt')
res1[,3:6] <- apply(res1[,3:6], 2, p.adjust, method = adjust_method)
res2 <- read.table('back_test_40_node_MF.txt')
res2[,3:6] <- apply(res2[,3:6], 2, p.adjust, method = adjust_method)
res3 <- read.table('back_test_40_node_CC.txt')
res3[,3:6] <- apply(res3[,3:6], 2, p.adjust, method = adjust_method)
dat <- rbind(res1, res2, res3)
ontology_type <- c(rep('BP', nrow(res1)), 
                   rep('MF', nrow(res2)), 
                   rep('CC', nrow(res3)))
dat <- data.frame(dat, ontology_type)
library(tidyr)
data_long <- gather(dat, test, pvalue, ZWLm:SKK, factor_key=TRUE)
p10 <- ggplot(data_long, aes(x = test, y = pvalue, fill = test)) +
    geom_boxplot(alpha = 0.7,  outlier.colour = "#1F3552", outlier.shape = 20,
                 width = .8, 
                 notch = TRUE, 
                 position = position_dodge(width = 0)) +
    scale_y_continuous(name = "P-values",
                       limits=c(0, 1)) +
    scale_x_discrete(name = "") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 11),
          legend.position = 'none') +
    facet_grid(. ~ ontology_type)
p10
# Out put the graph
ggsave('fig4.pdf',p10, width = 11, height = 7)
