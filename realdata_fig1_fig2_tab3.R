# Create the datafile and GO map files for beocat run
#
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
# Save datafile
saveRDS(dat, "data_file")
# Save GO map
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
    saveRDS(GO_gene_map, paste0("gene_map_",GO_type))
}

# Define a function running the tests given a GO_gene_map
# INPUT: GO_gene_map, data file, random seed
# OUTPUT: p-values for each test
run_tests <- function(GO_gene_map, dat, random_seed = 1, fn_prefix){
    # ready for test
    library(highD2pop)
    library(highDmean)
    # Randomly split
    set.seed(random_seed)
    sample_1_NEG <- sample(1:42, 21)
    sample_2_NEG <- setdiff(1:42, sample_1_NEG)
    sample_1_BCR <- sample(1:37, 18)
    sample_2_BCR <- setdiff(1:37, sample_1_BCR)
    # Initialize the output object
    N <- length(GO_gene_map)
    res <- matrix(0, nrow = N, ncol = 6)
    res <- data.frame(res)
    colnames(res) <- c('GOterm', 'numGene', 'ZWLm', 'GCTm', 'CQ', 'SKK')
    for (k in c(1, 2)){  # Loop over the two splits
        if (k == 1){
            sample_NEG <- sample_1_NEG
            sample_BCR <- sample_1_BCR
            filename <- paste0(fn_prefix, 'sample1_pvals')
        }
        else{
            sample_NEG <- sample_2_NEG
            sample_BCR <- sample_2_BCR
            filename <- paste0(fn_prefix, 'sample2_pvals')
        }
        for (i in 1: N){ # Loop over the GO terms
            # cat(i,"~")
            X <- as.matrix(dat[dat$class=='NEG', GO_gene_map[[i]]])
            Y <- as.matrix(dat[dat$class=='BCR/ABL', GO_gene_map[[i]]])
            X <- X[sample_NEG, ]
            Y <- Y[sample_BCR, ]
            res$GOterm[i] <- names(GO_gene_map)[i]
            res$numGene[i] <- length(GO_gene_map[[i]])
            res[i, 3:6] <- c(zwl_test(X,Y)$pvalue,
                             GCT.test(X,Y, r = 10, ntoorderminus = 0)$pvalue,
                             ChenQin.test(X,Y)$pvalue,
                             SKK_test(X, Y)$pvalue)
        }
        saveRDS(res, file = filename)
    }
}

# RUN 
# Time-consuming part.
# We ran this part on computing cluster
dat <- readRDS('data_file')
for (group in c('BP', 'MF', 'CC')){
    gene_map <- readRDS(paste0('gene_map_', group))
    for (k in 1:100){
        fn_prefix =paste0('res_',group,'_',k)
        run_tests(GO_gene_map = gene_map,
                  dat = dat,
                  random_seed = k,
                  fn_prefix = fn_prefix)
    }
}

# --------------------------------------------------------------------------
# Make Table 3
#
# Split the NEG and BCR/ABL samples to subsets, test the subsets, 
# observe the overlap.
#
# The splits were done by 100 times. This file will compute the 
# average IoU.
#
#+-------------------------------------------------------+
## Compute the overlap: Intersection over union

# Helper function
get_intersection_over_union <- function(pval1, pval2, alpha=0.05, 
                                        adjust_method = 'none')
{
    pval1 <- p.adjust(pval1, method = adjust_method)
    pval2 <- p.adjust(pval2, method = adjust_method)
    flag1 <- which(pval1 < alpha)
    flag2 <- which(pval2 < alpha)
    inter_size <- length(intersect(flag1, flag2))
    union_size <- length(union(flag1, flag2))
    IOU <- inter_size/union_size
    # return(c(inter_size, union_size, IOU))
    return(IOU)
}
adjust_method = 'bonferroni'
# loop over the splits
# GO_type = 'MF'
GO_type = 'BP'
# GO_type = 'CC'
res_matrix <- matrix(0, nrow = 1, ncol = 4)
colnames(res_matrix) <- c('ZWLm','GCTm', 'CQ', 'SKK')
for (k in 1:100){
    fn1 <- paste0('res_',GO_type,'_', k, 'sample1_pvals')
    fn2 <- paste0('res_',GO_type,'_', k, 'sample2_pvals')
    res1 <- readRDS(fn1)
    res2 <- readRDS(fn2)
    
    one_row <- NULL
    for ( i in 3:6) {
        one_row <- c(one_row, 
                     get_intersection_over_union(res1[, i], res2[,i],
                                                 alpha = 0.05,
                                                 adjust_method = adjust_method)
        )
    }
    res_matrix <- rbind(res_matrix, one_row)
}
res_matrix <- res_matrix[-1, ]
apply(res_matrix,2, function(x) round(mean(x, na.rm = T),3))
# --------------------------------------------------------------------------
# Make Figure 2
# Plot the Venn Diagram for one of the split.
#
S = 1
library(VennDiagram)
adjust_method = 'bonferroni'
GO_type = 'BP'
fn1 <- paste0('res_',GO_type,'_',S, 'sample1_pvals')
fn2 <- paste0('res_',GO_type,'_', S, 'sample2_pvals')
res1 <- readRDS(fn1)
res2 <- readRDS(fn2)
# Produce the plots.
alpha = 0.05
result_list = list()
union_size_list = rep(0, 4)
par(mfrow = c(1, 4))
for (i in 3:6){
    column_index = i
    pval1 <- p.adjust(res1[,column_index], method = adjust_method)
    pval2 <- p.adjust(res2[,column_index], method = adjust_method)
    flag1 <- which(pval1 < alpha)
    flag2 <- which(pval2 < alpha)
    inter_size <- length(intersect(flag1, flag2))
    area1 = length(flag1)
    area2 = length(flag2)
    union_size <- length(union(flag1, flag2))
    result_list[[i-2]] <- draw.pairwise.venn(area1 = area1, area2 = area2,
                                             cross.area = inter_size,
                                             lty = rep("blank", 2), 
                                             fill = c("light blue", "pink"),
                                             euler.d = TRUE, 
                                             ext.text= F,
                                             sep.dist = 0.03, 
                                             rotation.degree = 90,
                                             alpha = rep(0.65, 2),
                                             ind=FALSE)
    union_size_list[i-2] <- union_size
}
library(gridExtra)
lay1 <- matrix(rep(1, 36), nrow = 6)
lay2 <- rbind(
    rep(5, 4),
    matrix(rep(2, 16), nrow = 4),
    rep(5, 4)
)
lay3 <- matrix(c(5, 5, 3, 5,5, 5),nrow = 6)
lay4 <- rbind(
    rep(5, 3),
    matrix(rep(4, 9), nrow = 3),
    matrix(rep(5, 6), nrow = 2)
)
lay <- cbind(lay1, lay2, lay3, lay4)

grid.arrange(gTree(children=result_list[[1]]),
             gTree(children=result_list[[2]]),
             gTree(children=result_list[[3]]),
             gTree(children=result_list[[4]]),
             layout_matrix  = lay)
g1 <- arrangeGrob(gTree(children=result_list[[1]]),
                  gTree(children=result_list[[2]]),
                  gTree(children=result_list[[3]]),
                  gTree(children=result_list[[4]]),
                  layout_matrix  = lay)
ggplot2::ggsave(file="fig2.pdf", g1, width = 7, height = 3)


# --------------------------------------------------------------------------
# Make Figure 1
get_intersection_over_union <- function(pval1, pval2, alpha=0.05, 
                                        adjust_method = 'none')
{
    pval1 <- p.adjust(pval1, method = adjust_method)
    pval2 <- p.adjust(pval2, method = adjust_method)
    flag1 <- which(pval1 < alpha)
    flag2 <- which(pval2 < alpha)
    inter_size <- length(intersect(flag1, flag2))
    union_size <- length(union(flag1, flag2))
    IOU <- inter_size/union_size
    inter_size
}

GO_type ='BP'
BYres=NULL; Bonfres=NULL
for(split_group in 1:100){
    fn = paste0('res_', GO_type,'_', split_group) 
    res1= readRDS( paste0(fn, 'sample1_pvals'))
    res2= readRDS( paste0(fn, 'sample2_pvals'))
    BYsig =NULL;  Bonfsig=NULL
    for ( i in 3:6) {
        BYsig=c(BYsig, get_intersection_over_union(res1[, i], res2[,i],      
                                                   alpha = 0.05,      
                                                   adjust_method = 'BY'))    
        Bonfsig=c(Bonfsig,  get_intersection_over_union(res1[, i], res2[,i],      
                                                        alpha = 0.05,      
                                                        adjust_method =  'bonferroni'))    
    }
    BYres=rbind(BYres, BYsig)
    Bonfres=rbind(Bonfres, Bonfsig)
}
write.table(BYres, file=paste0(GO_type,'_', 'BYinterCount.txt') )
write.table(Bonfres, file=paste0(GO_type,'_','BonfinterCount.txt'))
BYres=read.table(paste0(GO_type,'_',	'BYinterCount.txt'), row.names = NULL)	
Bonfres=read.table(paste0(GO_type,'_','BonfinterCount.txt'), row.names = NULL) 
colnames(BYres)=colnames(Bonfres)= c("Adjustment", "ZWLm", "GCTm", "CQ", "SKK" ) 

sig_or_not=function(pval1) p.adjust(pval1, method = 'BY') < 0.05
BP1pv=readRDS(paste0('res_', GO_type,'_1sample1_pvals'))
BP2pv=readRDS(paste0('res_', GO_type,'_1sample2_pvals'))
BP1sig=apply(BP1pv[,3:6], 2, sig_or_not)
BP2sig=apply(BP2pv[,3:6], 2, sig_or_not)
set1GCT=as.character(BP1pv[BP1sig[,2], 1] )
set2GCT=as.character(BP2pv[BP2sig[,2], 1] )
set1CQ=as.character(BP1pv[BP1sig[,3], 1] )
set2CQ=as.character(BP2pv[BP2sig[,3], 1] )
set1SKK=as.character(BP1pv[BP1sig[,4], 1] )
set2SKK=as.character(BP2pv[BP2sig[,4], 1] )

# Make the left panel of figure 1
library(VennDiagram)
v=venn.diagram(
    x = list(set1GCT, set2GCT, set1SKK, set2SKK, set2CQ), 
    category.names = c("Sub-sample 1 GCT" , "Sub-sample 2 GCT", "Sub-sample 1 SKK", "Sub-sample 2 SKK", "Sub-sample 2 CQ"),
    filename = NULL, 
    lwd=4,cex=1.4, cat.cex = 1.4, cat.default.pos = "text",
    cat.dist = c(0.20, 0.25, -0.185, -0.2, 0.21),
    cat.pos = c(0,-0.35, 0,0,0),
    output=TRUE
)
pdf('fig_1_left_panel.pdf')
grid.draw(v);
dev.off()

# Make the right panel of figure 1
library(ggplot2)
library(tidyr)
BYresult=gather(BYres[,-1], key = "method", value = "intersectionCount")
Bonfresult=gather(Bonfres[,-1], key = "method", value = "intersectionCount")
ggplot(BYresult[BYresult$method!='ZWLm',], aes(x=method, y=intersectionCount)) + 
    geom_boxplot(fill="gray")+
    labs(title="Number of jointly identified GO terms from BCR/ABL vs NEG subsamples",
         x="Method", y = paste0("Count (out of",nrow(BP1pv), ")"))+ coord_flip()+
    theme(plot.title = element_text(size=10))
ggsave(paste0('fig_1_right_panel', GO_type, '.pdf'), height=2, width=5.5)




