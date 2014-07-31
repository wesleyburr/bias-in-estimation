bias-in-estimation
==================

Code and Results for paper "Bias in Estimation of ... "
by Burr, Takahara and Shin, submitted to Environmetrics, http://onlinelibrary.wiley.com/journal/10.1002/(ISSN)1099-095X.

Uses packages:
* gam: http://cran.r-project.org/web/packages/gam/index.html
* multitaper: http://cran.r-project.org/web/packages/multitaper/index.html
* xtable: http://cran.r-project.org/web/packages/xtable/index.html
* MASS: http://cran.r-project.org/web/packages/MASS/index.html
* splines: (base)
* slp: http://cran.r-project.org/web/packages/slp/index.html
* NMMAPSdata: http://www.ihapss.jhsph.edu/data/NMMAPS/R/
* tsinterp: http://github.com/wesleyburr/tsinterp/

To re-generate the tables and figures of the paper, using saved objects
(minimal time commitment), run:

````
> source("genFig1_2_3_gam.R")
> source("genFig4.R")
> source("genFig5.R")
> source("genFig6.R")
> source("analyzeSims.R")
> source("chicagoAnalysis.R")
> source("genFig8_9_10_gam.R")
````

In addition, code to re-create the basis objects is contained in several
of these routines as comments.

To re-create the simulated realizations, run
````
> source("sim_fullModels.R")
````
although you should expect it to take some time to complete.

If you are comfortable with makefiles, one is included which will
install all required packages and then run all code but the recreation
of the simulations. 

Under Windows:
````
C:\> make -f Makefile setupWin
C:\> make -f Makefile all
````

Under *nix:
````
[user@system ]$ make -f Makefile setupNix
[user@system ]$ make -f Makefile all
````

Testing was done:
*  x86_64 Linux (Fedora 20) installation, R 3.1.0
*  x86_64 Windows (7, SP1) installation, R 3.1.1
