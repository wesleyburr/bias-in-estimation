#
#  Analysis of Chicago data, 1987 - 2000, using the smoothers
#  developed for the Environmetrics paper
#
#  Full-span data only, using the classic Dominici et al. smooth function
#  link on temperature (tmpd). No humidity, DOW included. Collapsed mortality
#  records (cvd, death) from the NMMAPSdata package for the city. Pollutants 
#  used are pm10tmean and o3tmean, both lag-0 -- just for ease of analysis
#
#  Generates Figure 7 for paper, as well as results in Table 4
#
#  
library("gam")
library("slp")
library("gplots")
library("xtable")

load("./data/chic.RData")

tArr <- 1:length(chic[[1]])
dfGAM <- c(6 * 14, 12 * 14)
#  Code for creation of the 8 basis matrices of interest, as used in sim_fullModels
#  Bases are provided with this repository as "chicagoBases.RDa"
#  > basis1 <- ns(x = tArr, df = dfGAM[1])
#  > basis2 <- ns(x = tArr, df = dfGAM[2])
#  > basis3 <- slp(x = tArr, K = dfGAM[1], naive = TRUE)
#  > basis4 <- slp(x = tArr, K = dfGAM[2], naive = TRUE)
#  > basis5 <- slp(x = tArr, K = dfGAM[1], intercept = TRUE)
#  > basis6 <- slp(x = tArr, K = dfGAM[2], intercept = TRUE)
#  > basis7 <- basis5[, -1] # equiv to running slp with intercept = FALSE, faster
#  > basis8 <- basis6[, -1] #   than recomputing
#  > 
#  > bases <- vector("list", 8)
#  > for(j in 1:8) { bases[[j]] <- get(paste("basis", j, sep="")) }
#  > save(file="./data/chicagoBases.RData", bases)
#

load("./data/chicagoBases.RData")
for(j in 1:8) {
  assign(paste("basis", j, sep = ""), bases[[j]])
}


################################################################################
#
#  Ozone and Death
#
dat <- chic
names(dat)[names(dat) %in% c("death", "o3tmean")] <- c("Y", "X")

commonAll1 <- "Y ~ X + as.factor(dow) + ns(tmpd, df = 3)"
commonAll2 <- "Y ~ X + as.factor(dow)" 

models <- c(paste0(commonAll1, " + basis1"),     paste0(commonAll2, " + basis1"), 
            paste0(commonAll1, " + basis2"),     paste0(commonAll2, " + basis2"),
            paste0(commonAll1, " + basis3"),     paste0(commonAll2, " + basis3"),
            paste0(commonAll1, " + basis4"),     paste0(commonAll2, " + basis4"),
            paste0(commonAll1, " + basis5 - 1"), paste0(commonAll2, " + basis5 - 1"),
            paste0(commonAll1, " + basis6 - 1"), paste0(commonAll2, " + basis6 - 1"),
            paste0(commonAll1, " + basis7"),     paste0(commonAll2, " + basis7"),
            paste0(commonAll1, " + basis8"),     paste0(commonAll2, " + basis8"))

fits1 <- lapply(models, FUN = function(x){ gam(as.formula(x), family = poisson, 
                                               data = dat, na.action = na.omit) })

################################################################################
#
#  Ozone and CVD
#
dat <- chic
names(dat)[names(dat) %in% c("cvd", "o3tmean")] <- c("Y", "X")

commonAll1 <- "Y ~ X + as.factor(dow) + ns(tmpd, df = 3)"
commonAll2 <- "Y ~ X + as.factor(dow)" 

models <- c(paste0(commonAll1, " + basis1"),     paste0(commonAll2, " + basis1"), 
            paste0(commonAll1, " + basis2"),     paste0(commonAll2, " + basis2"),
            paste0(commonAll1, " + basis3"),     paste0(commonAll2, " + basis3"),
            paste0(commonAll1, " + basis4"),     paste0(commonAll2, " + basis4"),
            paste0(commonAll1, " + basis5 - 1"), paste0(commonAll2, " + basis5 - 1"),
            paste0(commonAll1, " + basis6 - 1"), paste0(commonAll2, " + basis6 - 1"),
            paste0(commonAll1, " + basis7"),     paste0(commonAll2, " + basis7"),
            paste0(commonAll1, " + basis8"),     paste0(commonAll2, " + basis8"))

fits2 <- lapply(models, FUN = function(x){ gam(as.formula(x), family = poisson, 
                                               data = dat, na.action = na.omit) })

################################################################################
#
#  PM10 and Death
#
dat <- chic
names(dat)[names(dat) %in% c("death", "pm10tmean")] <- c("Y", "X")

commonAll1 <- "Y ~ X + as.factor(dow) + ns(tmpd, df = 3)"
commonAll2 <- "Y ~ X + as.factor(dow)" 

models <- c(paste0(commonAll1, " + basis1"),     paste0(commonAll2, " + basis1"), 
            paste0(commonAll1, " + basis2"),     paste0(commonAll2, " + basis2"),
            paste0(commonAll1, " + basis3"),     paste0(commonAll2, " + basis3"),
            paste0(commonAll1, " + basis4"),     paste0(commonAll2, " + basis4"),
            paste0(commonAll1, " + basis5 - 1"), paste0(commonAll2, " + basis5 - 1"),
            paste0(commonAll1, " + basis6 - 1"), paste0(commonAll2, " + basis6 - 1"),
            paste0(commonAll1, " + basis7"),     paste0(commonAll2, " + basis7"),
            paste0(commonAll1, " + basis8"),     paste0(commonAll2, " + basis8"))

fits3 <- lapply(models, FUN = function(x){ gam(as.formula(x), family = poisson, 
                                               data = dat, na.action = na.omit) })

################################################################################
#
#  PM10 and CVD
#
dat <- chic
names(dat)[names(dat) %in% c("cvd", "pm10tmean")] <- c("Y", "X")

commonAll1 <- "Y ~ X + as.factor(dow) + ns(tmpd, df = 3)"
commonAll2 <- "Y ~ X + as.factor(dow)" 

models <- c(paste0(commonAll1, " + basis1"),     paste0(commonAll2, " + basis1"), 
            paste0(commonAll1, " + basis2"),     paste0(commonAll2, " + basis2"),
            paste0(commonAll1, " + basis3"),     paste0(commonAll2, " + basis3"),
            paste0(commonAll1, " + basis4"),     paste0(commonAll2, " + basis4"),
            paste0(commonAll1, " + basis5 - 1"), paste0(commonAll2, " + basis5 - 1"),
            paste0(commonAll1, " + basis6 - 1"), paste0(commonAll2, " + basis6 - 1"),
            paste0(commonAll1, " + basis7"),     paste0(commonAll2, " + basis7"),
            paste0(commonAll1, " + basis8"),     paste0(commonAll2, " + basis8"))

fits4 <- lapply(models, FUN = function(x){ gam(as.formula(x), family = poisson, 
                                               data = dat, na.action = na.omit) })

################################################################################
#
#  Consolidate the coefficients for comparison, with their standard errors
#  * S1 - S4, S7, S8 are all normal; S5 and S6 have the intercept in coefficient 
#    slot 11 instead, scaled differently
#
#  * 4 sets of fits, each 16 long, each 8 smoothers and tmpd included
#  * 64 coefficients and standard errors
coefs <- matrix(data = 0, nrow = 4 * 16, ncol = 2)
fits1 <- lapply(fits1, FUN = function(x) { y <- summary.glm(x)$coefficients; y[rownames(y) == "X", 1:2] })[seq(1, 16, 2)]
fits2 <- lapply(fits2, FUN = function(x) { y <- summary.glm(x)$coefficients; y[rownames(y) == "X", 1:2] })[seq(1, 16, 2)]
fits3 <- lapply(fits3, FUN = function(x) { y <- summary.glm(x)$coefficients; y[rownames(y) == "X", 1:2] })[seq(1, 16, 2)]
fits4 <- lapply(fits4, FUN = function(x) { y <- summary.glm(x)$coefficients; y[rownames(y) == "X", 1:2] })[seq(1, 16, 2)]

coefs <- cbind(matrix(unlist(fits1), nrow = 8, ncol = 2, byrow = TRUE),
               matrix(unlist(fits2), nrow = 8, ncol = 2, byrow = TRUE),
               matrix(unlist(fits3), nrow = 8, ncol = 2, byrow = TRUE),
               matrix(unlist(fits4), nrow = 8, ncol = 2, byrow = TRUE))
# rearrange into -6 and -12 sets
coefs <- coefs[c(1,3,5,7,2,4,6,8), ]

#
#  Results for Table 4
# 
nameOrder <- c("NS", "SLP", "SLP2", "SLP3")
row.names(coefs) <- c(paste0("$mathbf{S}$", nameOrder, "-6"),
                      paste0("$mathbf{S}$", nameOrder, "-12"))

tab1 <- xtable(coefs[, 1:4], digits = 3, display = rep("e", 5))
tab2 <- xtable(coefs[, 5:8], digits = 3, display = rep("e", 5))

print.xtable(tab1, file = "./tables/chicagoTableA.tex")
print.xtable(tab2, file = "./tables/chicagoTableB.tex")

################################################################################
#
#  The coefficients: sets of 4 (4 combinations of response-predictor,
#  8 models which are similar (S1,3,5,7 or S2,4,6,8), with/without tmpd
#  = 4 * 8 * 2 = 64 coefficients)
#  * organized into an 8x8 data frame 'coefs'

coefSet <- rbind(coefs[c(1, 7), 1:2], 
                 coefs[c(1, 7), 3:4], 
                 coefs[c(1, 7), 5:6], 
                 coefs[c(1, 7), 7:8])

coefCenter <- coefSet[, 1] 
coefLeft   <- cbind(coefSet[, 1] - 2 * coefSet[, 2],
                    coefCenter - 4e-5)
coefRight  <- cbind(coefCenter + 4e-5, 
                    coefSet[, 1] + 2 * coefSet[, 2])

minX <- min(coefLeft); maxX <- max(coefRight) * 1.1
yLabels <- c("Death/Ozone", "Death/PM10", "CVD/Ozone", "CVD/PM10")

smootherLabels <- c("S-NS-6", "S-SLP2-12")
smootherLabelPos <- coefRight[, 2] 

#
#  Figure 7 - updated July 31, 2014
#
# pdf(file = "figures/chicagoDeathOzoneCI.pdf", width = 9, height = 9)
postscript(file = "figures/fig7-chicagoDeathOzoneCI.eps", width = 9, height = 9,
           horizontal = FALSE, paper = 'special')
par(mar = c(4,4,0.5,0.5))
plot(x = NA, y = NA, ylim = c(0.5,8.5), xlim = c(minX, maxX), ylab = "", xlab ="", yaxt = 'n')
abline(v = 0, lty = 2, col = "grey60")
points(coefCenter, 8:1, pch = 19)
axis(side = 2, at = seq(7.5, 1.5, -2), labels = yLabels)
for(j in 1:8) {
  lines(coefLeft[j, ], rep((8-j)+1, 2))
  lines(coefRight[j, ], rep((8-j)+1, 2))
}
text(x = smootherLabelPos, y = seq(8.2,1.2,-1), labels = smootherLabels)
dev.off()

