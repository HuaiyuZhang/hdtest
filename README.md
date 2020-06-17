# hdtest
Reproducible code for high dimensional two sample test (Zhang and Wang, 2020).

These codes are used to reproduce all the tables and figures in the paper (Zhang and Wang, 2020).
The implementation of the new test is available through CRAN R package [highDmean](https://cran.r-project.org/web/packages/highDmean/).

# code structure

The prefix "realdata" and "simulation" indicate the source of the data. 
The "realdata" refers to the acute lymphoblastic leukemia (ALL) microarray gene expressions dataset 
(Chiaretti et al., 2004).
The ALL dataset is used to demonstrate the result consistency evaluation for the differentially expressed 
Gene Ontology (GO) terms and the performance of four competing high-dimensional tests. 
The "simulation" refers to the simulated data that
are used to evaluate the empirical type I error and power of four competing high-dimensional tests. 

# Running time

Most of the codes were run on a computing cluster. 
The following running times are estimated by summing up the computing time on the cluster nodes. 
It may not reflect the actual running time on a desktop with high computing capacity.

* realdata_fig1_fig2_tab3.R (150 hours)
* realdata_tab1_tab2_fig3_fig4.R (2 hours)
* simulation_fig5_fig6.R (24 hours)
* simulation_fig5_fig6_helper.R (a helper file)
* simulation_fig7_tab4.R (100 hours)
* simulation_tab5.R (64 hours)
* simulation_tab5_helper.R (a helper file)