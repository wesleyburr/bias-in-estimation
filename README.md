bias-in-estimation
==================

Code and Results for paper "Bias in Estimation of ... "
by Burr, Takahara and Shin, submitted to Environmetrics, http://onlinelibrary.wiley.com/journal/10.1002/(ISSN)1099-095X.

Uses packages:
* gam: http://cran.r-project.org/web/packages/gam/index.html
* multitaper: http://cran.r-project.org/web/packages/multitaper/index.html
* gplots: http://cran.r-project.org/web/packages/gplots/index.html
* xtable: http://cran.r-project.org/web/packages/xtable/index.html
* MASS: http://cran.r-project.org/web/packages/MASS/index.html
* splines: (base)
* slp: http://github.com/wesleyburr/slp/
* NMMAPSdata: http://www.ihapss.jhsph.edu/data/NMMAPS/R/
* tsinterp: http://github.com/wesleyburr/tsinterp/

To re-generate the tables and figures of the paper, using saved objects
(minimal time commitment), run:

> source("genFig1_2_3_gam.R")
> source("genFig4.R")
> source("genFig5.R")
> source("genFig6.R")
> source("analyzeSims.R")
> source("chicagoAnalysis.R")
> source("genFig8_9_10_gam.R")

In addition, code to re-create the basis objects is contained in several
of these routines as comments.

To re-create the simulated realizations, run
> source("sim_fullModels.R")
although you should expect it to take some time to complete.

All testing was done on R 3.1.0 on an x86_64 Linux machine running
Fedora 20. 
