#
#  Simple make-file, computing figures and tables in order
#
setupWin: 
	R CMD BATCH installCRANpackages.R
	R CMD INSTALL ./packages/*.zip
	rm -vf *.Rout .RData

setupNix:
	R CMD BATCH installCRANpackages.R
	R CMD INSTALL ./packages/*.tar.gz
	rm -vf *.Rout .RData

all: 
	rm -vf ./figures/*.eps 
	rm -vf ./tables/*.tex
	R CMD BATCH genFig1_2_S1.R 
	R CMD BATCH genFig3.R 
	R CMD BATCH genFig4.R 
	R CMD BATCH chicagoAnalysis.R     # saves two tables to ./tables/, one figure to ./figures
	R CMD BATCH analyzeSims.R         # saves four tables to ./tables
	R CMD BATCH genFigS2.R
	rm -vf *.Rout .RData

#  If you want to verify the realizations of the simulations 
#  (or generate a new set) ...
sim: 
	R CMD BATCH sim_fullModels.R
	rm -vf *.Rout .RData


