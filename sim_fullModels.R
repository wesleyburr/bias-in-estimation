#
#  Simulation: Models
#
#  Simulations for Section 5.1
#  * Model: Y ~ X + smooth function + \epsilon
#  * variations to demonstrate the different smoothers
#
#  Generates realizations which are used as simulations in paper; file
#  analyzeSims.R uses these saved realizations and analyzes them.
#

################################################################################
#
# Two-sided running-mean smoother, uses smoother matrix setup; neh = (lFilt - 1)/2
#
twoSidedSmooth <- function(x, neh) {
    
    # neh = Number Each Half

    stopifnot(is.numeric(neh))
    stopifnot(is.numeric(x), (N <- length(x)) > (2*neh+1))

    lFilt <- 2 * neh + 1
    filt <- c(rep(1/lFilt, lFilt), rep(0, (N + 1 - lFilt)))

    # filter gets applied to the (neh+1)th point ... (N-neh)th point
    # = total of N - lFilt + 1 rows
    filtMat <- matrix(data = rep(filt, N - lFilt + 1)[1:(N * (N - 2 * neh))], ncol = N, nrow = (N - 2 * neh), byrow = TRUE)

    y <- x
    y[1:neh] <- rep(mean(x[1:(neh+1)], na.rm = TRUE), neh)
    y[(N - neh + 1):N] <- rep(mean(x[(N - neh):N], na.rm = TRUE), neh)
    y[(neh+1):(N-neh)] <- (filtMat %*% x) 

    y
}

################################################################################
#
# Create a matrix of sinusoids for setting background noise; all are perfectly
# in phase unless otherwise specified ... could easily set arbitrary phases for 
# more randomness
#
createSineArr <- function(n, dT, nMesh, phse = NULL) {
  tArr <- (1:n)*86400
  sineMat <- matrix(data = 0, nrow = n, ncol = nMesh) 
  freq <- seq(1/nMesh/dT/2, 1/2/dT, 1/nMesh/dT/2) 
  if(is.null(phse)) {
    for(j in 1:nMesh) { sineMat[, j] <- sin(2 * pi * tArr * freq[j]) }
  } else if(is.numeric(phse)) {
     for(j in 1:nMesh) { sineMat[, j] <- sin(2 * pi * tArr * freq[j] + phse[j]) }
  }
  sineMat
}

################################################################################
#
#  Compute metrics and statistics for each realization (e.g., M1, M2, etc.)
#
computeStats <- function(S, Xlp, Xhp, Ylp, Yhp) {
 
  X <- Xlp + Xhp
  Y <- Ylp + Yhp
  oneV <- rep(1, length(X))
  xHat <- X - S %*% X
  yHat <- Y - S %*% Y
  oHat <- oneV - S %*% oneV
   
  M1 <- t(oHat) %*% oHat    # || oHat ||^2
  M2 <- t(oHat) %*% xHat    # oHat * x_{HP}
  xHat2 <- t(xHat) %*% xHat

  alp <- (M1 * xHat2) / (M1 * xHat2 - M2^2)
  gam <- (M2 * t(oHat) %*% yHat) / (M1 * t(xHat) %*% xHat - M2^2)

  fits <- matrix(data = 0, nrow = 12, ncol = 1)

  fits[1]  <- lm(Y ~ X)$coefficients[2]
  fits[2]  <- lm(Yhp ~ X)$coefficients[2]
  fits[3]  <- lm(yHat ~ X)$coefficients[2]
  fits[4]  <- lm(yHat - Yhp ~ X)$coefficients[2]
  fits[5]  <- lm(Y ~ Xhp)$coefficients[2]
  fits[6]  <- lm(Yhp ~ Xhp)$coefficients[2]
  fits[7]  <- lm(yHat ~ Xhp)$coefficients[2]
  fits[8]  <- lm(yHat - Yhp ~ Xhp)$coefficients[2]
  fits[9]  <- lm(Y ~ xHat)$coefficients[2]
  fits[10] <- lm(Yhp ~ xHat)$coefficients[2]
  fits[11] <- lm(yHat ~ xHat)$coefficients[2]
  fits[12] <- lm(yHat - Yhp ~ xHat)$coefficients[2]

  return(list(alp, gam, fits))
}

################################################################################
#
#  Setup for simulations
#

library("gam")  #  auto-loads 'splines'
library("slp")
library("MASS")

N <- 365 * 10
tArr <- 1:N
nSim <- 250    # number of realizations
nMesh <- 2^14
dT <- 86400  # 86,400 seconds in 1 day

sineMat <- createSineArr(n = N, dT = dT, nMesh = nMesh, phse = NULL)
phseArr <- rep(0, nMesh)
freqs <- seq(1/nMesh/dT/2, 1/2/dT, 1/nMesh/dT/2)

yMean <- 1.0
dfGAM <- c(60, 120, 60, 120)
oneV  <- rep(1, N)

# set seed for rnorm() below [through all 4 simulations]
set.seed(23)

#  Code to create bases for various smoothers; these have been precomputed and are
#  provided as part of this reposity
#
#  > basis1 <- ns(x = tArr, df = dfGAM[1])
#  > basis2 <- ns(x = tArr, df = dfGAM[2])
#  > basis3 <- slp(x = tArr, K = dfGAM[3], naive = TRUE)
#  > basis4 <- slp(x = tArr, K = dfGAM[4], naive = TRUE)
#  > basis5 <- slp(x = tArr, K = dfGAM[3], intercept = TRUE)
#  > basis6 <- slp(x = tArr, K = dfGAM[4], intercept = TRUE)
#  > basis7 <- basis5[, -1] # equiv to running slp with intercept = FALSE
#  > basis8 <- basis6[, -1]
#  > bases <- vector("list", 8)
#  > for(j in 1:8) { bases[[j]] <- get(paste("basis", j, sep="")) }
#  > save(file="./data/fullModelsBases.RDa", bases)
#

# Load the 8 basis sets and assign them to variable names used below
load("./data/fullModelsBases.RDa")
for(j in 1:8) {
  assign(paste("basis", j, sep = ""), bases[[j]])
}

S1 <- basis1 %*% ginv(t(basis1) %*% basis1) %*% t(basis1)
S2 <- basis2 %*% ginv(t(basis2) %*% basis2) %*% t(basis2)
S3 <- basis3 %*% t(basis3)  # orthonormal, no inverse needed
S4 <- basis4 %*% t(basis4)
S5 <- basis5 %*% t(basis5)
S6 <- basis6 %*% t(basis6)
S7 <- basis7 %*% t(basis7)
S8 <- basis8 %*% t(basis8)


################################################################################
#
#  Simulation 1: 
#  
#
sim1 <- vector("list", nSim)  # [[1]] is the general stuff, then [[2]] through [[5]] is for the four smoothers
sim1 <- lapply(sim1, FUN = function(x){y <- vector("list", 9)
                                       x <- lapply(y, FUN = function(z){ z <- vector("list", 4) })
                                      })

periodsLPx <- c(183, 75)   
ampsLPx <- c(2, 2)
phasesLPx <- c(0, 0)

periodsLPy <- c(183, 75)   
ampsLPy <- c(2, 2)
phasesLPy <- c(180, 180)

periodsHP <- c(5)
ampsHP <- c(2)
phasesHP <- c(0)

Wcut <- min(which(freqs > 1/60/dT))     # dividing line between "long" and "short"

cat("Simulating: ")
for(k in 1:nSim) {
    if((k %% 25) == 0) { cat(".") }

    # simulate noise background - low frequency random, high frequency shared
    # (this ensures that the 'true' correlation coefficient is consistent)
    Xlp <- sineMat[, 1:(Wcut-1)] %*% rnorm(n = Wcut - 1, sd = sqrt(1/nMesh))
    Ylp <- sineMat[, 1:(Wcut-1)] %*% rnorm(n = Wcut - 1, sd = sqrt(1/nMesh))
    Xhp <- sineMat[, Wcut:nMesh] %*% rnorm(n = (nMesh - Wcut + 1), sd = sqrt(1/nMesh))
    Yhp <- Xhp

    # add low pass sinusoids in
    for(j in 1:length(periodsLPx)) {
      Xlp <- Xlp + ampsLPx[j] * sin(2 * pi * tArr / periodsLPx[j] + phasesLPx[j])
      Ylp <- Ylp + ampsLPy[j] * sin(2 * pi * tArr / periodsLPy[j] + phasesLPy[j])
    }  

    # add high pass sinusoids in
    for(j in 1:length(periodsHP)) {
      Xhp <- Xhp + ampsHP[j] * sin(2 * pi * tArr / periodsHP[j] + phasesHP[j])
      Yhp <- Yhp + ampsHP[j] * sin(2 * pi * tArr / periodsHP[j] + phasesHP[j])
    }  

    X <- Xlp + Xhp
    Y <- Ylp + Yhp

    # add a clip to keep the variance from spinning out of control
    Y[which(Y > quantile(Y, 0.99))] <- quantile(Y, 0.99)
    Y[which(Y < quantile(Y, 0.01))] <- quantile(Y, 0.01)

    sim1[[k]][[9]] <- c(cor(Xlp, Ylp), cor(Xhp, Yhp), cor(X, Y), (t(Xhp) %*% Yhp) / (t(Xhp) %*% Xhp))

    # Find M_1 and M_2 for each smoother; requires finding \hat{x}_HP via each as well,
    # and \hat{1}_HP and \hat{y}_HP
    sim1[[k]][[1]][1:3] <- computeStats(S1, Xlp, Xhp, Ylp, Yhp)
    sim1[[k]][[2]][1:3] <- computeStats(S2, Xlp, Xhp, Ylp, Yhp)
    sim1[[k]][[3]][1:3] <- computeStats(S3, Xlp, Xhp, Ylp, Yhp)
    sim1[[k]][[4]][1:3] <- computeStats(S4, Xlp, Xhp, Ylp, Yhp)
    sim1[[k]][[5]][1:3] <- computeStats(S5, Xlp, Xhp, Ylp, Yhp)
    sim1[[k]][[6]][1:3] <- computeStats(S6, Xlp, Xhp, Ylp, Yhp)
    sim1[[k]][[7]][1:3] <- computeStats(S7, Xlp, Xhp, Ylp, Yhp)
    sim1[[k]][[8]][1:3] <- computeStats(S8, Xlp, Xhp, Ylp, Yhp)

    dat <- as.data.frame(cbind(X, Y + yMean, tArr))
    names(dat) <- c("X", "lY", "tArr")

    modelNS1 <- paste0("lY ~ X + basis1")
    modelNS2 <- paste0("lY ~ X + basis2")
    modelSP1 <- paste0("lY ~ X + basis3")
    modelSP2 <- paste0("lY ~ X + basis4")
    modelSP3 <- paste0("lY ~ X + basis5 - 1")  # remove intercept from model
    modelSP4 <- paste0("lY ~ X + basis6 - 1")  # remove intercept from model
    modelSP5 <- paste0("lY ~ X + basis7")
    modelSP6 <- paste0("lY ~ X + basis8")

    fitNS1 <- gam(as.formula(modelNS1), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitNS2 <- gam(as.formula(modelNS2), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP1 <- gam(as.formula(modelSP1), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP2 <- gam(as.formula(modelSP2), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP3 <- gam(as.formula(modelSP3), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP4 <- gam(as.formula(modelSP4), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP5 <- gam(as.formula(modelSP5), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP6 <- gam(as.formula(modelSP6), family = gaussian, data = dat, 
                  na.action = na.omit)

    sim1[[k]][[1]][4] <- list(summary.glm(fitNS1)$coefficients[1:2])
    sim1[[k]][[2]][4] <- list(summary.glm(fitNS2)$coefficients[1:2])
    sim1[[k]][[3]][4] <- list(summary.glm(fitSP1)$coefficients[1:2])
    sim1[[k]][[4]][4] <- list(summary.glm(fitSP2)$coefficients[1:2])
    sim1[[k]][[5]][4] <- list(c(summary.glm(fitSP3)$coefficients[2] * mean(basis5[, 1]), 
                                           summary.glm(fitSP3)$coefficients[1]))
    sim1[[k]][[6]][4] <- list(c(summary.glm(fitSP4)$coefficients[2] * mean(basis6[, 1]), 
                                           summary.glm(fitSP4)$coefficients[1]))
    sim1[[k]][[7]][4] <- list(summary.glm(fitSP5)$coefficients[1:2])
    sim1[[k]][[8]][4] <- list(summary.glm(fitSP6)$coefficients[1:2])

}
cat("\n")

save(file = "./data/simFullModelsResults1.RData", sim1)

################################################################################
#
#  Simulation 2:
#  * same (uniform) variance for low and high pass (white-ish spectrum), but with
#    correlation below (183), inside (75, 105), and above the bad pass-band area (5, 14)
#  * make all correlation longer than 2 months negative, less than, positive
#
sim2 <- vector("list", nSim)  # [[1]] is the general stuff, then [[2]] through [[5]] is for the four smoothers
sim2 <- lapply(sim2, FUN = function(x){y <- vector("list", 9)
                                       x <- lapply(y, FUN = function(z){ z <- vector("list", 4) })
                                      })

periodsLPx <- c(183, 105, 75, 68)
ampsLPx <- c(2, 2, 2, 2)
phasesLPx <- c(0, 0, 0, 0)

periodsLPy <- c(183, 105, 75, 68)
ampsLPy <- c(2, 2, 2, 2)
phasesLPy <- c(180, 180, 180, 180)  # negative correlation on these points

periodsHP <- c(5, 14)
ampsHP <- c(2, 2)
phasesHP <- c(0, 0)

Wcut <- min(which(freqs > 1/60/dT))   # cut-off of 60 days, 2 months

cat("Simulating: ")
for(k in 1:nSim) {
    if((k %% 25) == 0) { cat(".") }

    # simulate noise background - low frequency random, high frequency shared
    # (this ensures that the 'true' correlation coefficient is consistent)
    Xlp <- sineMat[, 1:(Wcut-1)] %*% rnorm(n = Wcut - 1, sd = (1/nMesh))
    Ylp <- sineMat[, 1:(Wcut-1)] %*% rnorm(n = Wcut - 1, sd = (1/nMesh))
    Xhp <- sineMat[, Wcut:nMesh] %*% rnorm(n = (nMesh - Wcut + 1), sd = (1/nMesh))
    Yhp <- Xhp

    # add low pass sinusoids in
    for(j in 1:length(periodsLPx)) {
      Xlp <- Xlp + ampsLPx[j] * sin(2 * pi * tArr / periodsLPx[j] + phasesLPx[j])
      Ylp <- Ylp + ampsLPy[j] * sin(2 * pi * tArr / periodsLPy[j] + phasesLPy[j])
    }  

    # add high pass sinusoids in
    for(j in 1:length(periodsHP)) {
      Xhp <- Xhp + ampsHP[j] * sin(2 * pi * tArr / periodsHP[j] + phasesHP[j])
      Yhp <- Yhp + ampsHP[j] * sin(2 * pi * tArr / periodsHP[j] + phasesHP[j])
    }  

    X <- Xlp + Xhp
    Y <- Ylp + Yhp

    # add a clip to keep the variance from spinning out of control
    Y[which(Y > quantile(Y, 0.99))] <- quantile(Y, 0.99)
    Y[which(Y < quantile(Y, 0.01))] <- quantile(Y, 0.01)

    sim2[[k]][[9]] <- c(cor(Xlp, Ylp), cor(Xhp, Yhp), cor(X, Y), (t(Xhp) %*% Yhp) / (t(Xhp) %*% Xhp))

    # Find M_1 and M_2 for each smoother; requires finding \hat{x}_HP via each as well,
    # and \hat{1}_HP and \hat{y}_HP
    sim2[[k]][[1]][1:3] <- computeStats(S1, Xlp, Xhp, Ylp, Yhp)
    sim2[[k]][[2]][1:3] <- computeStats(S2, Xlp, Xhp, Ylp, Yhp)
    sim2[[k]][[3]][1:3] <- computeStats(S3, Xlp, Xhp, Ylp, Yhp)
    sim2[[k]][[4]][1:3] <- computeStats(S4, Xlp, Xhp, Ylp, Yhp)
    sim2[[k]][[5]][1:3] <- computeStats(S5, Xlp, Xhp, Ylp, Yhp)
    sim2[[k]][[6]][1:3] <- computeStats(S6, Xlp, Xhp, Ylp, Yhp)
    sim2[[k]][[7]][1:3] <- computeStats(S7, Xlp, Xhp, Ylp, Yhp)
    sim2[[k]][[8]][1:3] <- computeStats(S8, Xlp, Xhp, Ylp, Yhp)

    dat <- as.data.frame(cbind(X, Y + yMean, tArr))
    names(dat) <- c("X", "lY", "tArr")

    modelNS1 <- paste0("lY ~ X + basis1")
    modelNS2 <- paste0("lY ~ X + basis2")
    modelSP1 <- paste0("lY ~ X + basis3")
    modelSP2 <- paste0("lY ~ X + basis4")
    modelSP3 <- paste0("lY ~ X + basis5 - 1")  # remove intercept from model
    modelSP4 <- paste0("lY ~ X + basis6 - 1")  # remove intercept from model
    modelSP5 <- paste0("lY ~ X + basis7")
    modelSP6 <- paste0("lY ~ X + basis8")

    fitNS1 <- gam(as.formula(modelNS1), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitNS2 <- gam(as.formula(modelNS2), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP1 <- gam(as.formula(modelSP1), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP2 <- gam(as.formula(modelSP2), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP3 <- gam(as.formula(modelSP3), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP4 <- gam(as.formula(modelSP4), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP5 <- gam(as.formula(modelSP5), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP6 <- gam(as.formula(modelSP6), family = gaussian, data = dat, 
                  na.action = na.omit)

    sim2[[k]][[1]][4] <- list(summary.glm(fitNS1)$coefficients[1:2])
    sim2[[k]][[2]][4] <- list(summary.glm(fitNS2)$coefficients[1:2])
    sim2[[k]][[3]][4] <- list(summary.glm(fitSP1)$coefficients[1:2])
    sim2[[k]][[4]][4] <- list(summary.glm(fitSP2)$coefficients[1:2])
    sim2[[k]][[5]][4] <- list(c(summary.glm(fitSP3)$coefficients[2] * mean(basis5[, 1]), 
                                           summary.glm(fitSP3)$coefficients[1]))
    sim2[[k]][[6]][4] <- list(c(summary.glm(fitSP4)$coefficients[2] * mean(basis6[, 1]), 
                                           summary.glm(fitSP4)$coefficients[1]))
    sim2[[k]][[7]][4] <- list(summary.glm(fitSP5)$coefficients[1:2])
    sim2[[k]][[8]][4] <- list(summary.glm(fitSP6)$coefficients[1:2])

}
cat("\n")

save(file = "./data/simFullModelsResults2.RData", sim2)

################################################################################
#
#  Simulation 3:
#  * have an exponential decay to the variance, Kolmogorov-like
#  * make phase of LP stuff anti-correlated (cor = -1)
#  * make noise of HP stuff correlated
#  * make all sinusoids anti- or correlated depending on Wcut
#
sim3 <- vector("list", nSim)  # [[1]] is the general stuff, then [[2]] through [[5]] is for the four smoothers
sim3 <- lapply(sim3, FUN = function(x){y <- vector("list", 9)
                                       x <- lapply(y, FUN = function(z){ z <- vector("list", 4) })
                                      })
periodsLPx <- c(183, 105, 75, 68)   
ampsLPx <- c(2, 2, 2, 2)
phasesLPx <- c(0, 0, 0, 0)

periodsLPy <- c(183, 105, 75, 68)  
ampsLPy <- c(2, 2, 2, 2)
phasesLPy <- c(180, 180, 180, 180)

periodsHP <- c(2, 5, 14)
ampsHP <- c(2, 2, 2)
phasesHP <- c(0, 0, 0)

Wcut <- min(which(freqs > 1/60/dT))   # cut-off of 60 days, 2 months

cat("Simulating: ")
for(k in 1:nSim) {
    if((k %% 25) == 0) { cat(".") }

    # noise background has variance (2*exp(-1e6*freqs) + 1) * sqrt(1/nMesh)
    noiseVar <- (2 * exp(-1e6 * freqs) + 1) * sqrt(1 / nMesh)

    # simulate noise background - low frequency random, high frequency shared
    # (this ensures that the 'true' correlation coefficient is consistent)
    Xlp <- sineMat[, 1:(Wcut-1)] %*% rnorm(n = Wcut - 1, sd = noiseVar[1:(Wcut-1)])
    Ylp <- sineMat[, 1:(Wcut-1)] %*% rnorm(n = Wcut - 1, sd = noiseVar[1:(Wcut-1)])
    Xhp <- sineMat[, Wcut:nMesh] %*% rnorm(n = (nMesh - Wcut + 1), sd = noiseVar[Wcut:(nMesh)])
    Yhp <- Xhp

    # add low pass sinusoids in
    for(j in 1:length(periodsLPx)) {
      Xlp <- Xlp + ampsLPx[j] * sin(2 * pi * tArr / periodsLPx[j] + phasesLPx[j])
      Ylp <- Ylp + ampsLPy[j] * sin(2 * pi * tArr / periodsLPy[j] + phasesLPy[j])
    }  

    # add high pass sinusoids in
    for(j in 1:length(periodsHP)) {
      Xhp <- Xhp + ampsHP[j] * sin(2 * pi * tArr / periodsHP[j] + phasesHP[j])
      Yhp <- Yhp + ampsHP[j] * sin(2 * pi * tArr / periodsHP[j] + phasesHP[j])
    }  

    X <- Xlp + Xhp
    Y <- Ylp + Yhp
    
    # add a clip to keep the variance from spinning out of control
    Y[which(Y > quantile(Y, 0.99))] <- quantile(Y, 0.99)
    Y[which(Y < quantile(Y, 0.01))] <- quantile(Y, 0.01)

    sim3[[k]][[9]] <- c(cor(Xlp, Ylp), cor(Xhp, Yhp), cor(X, Y), (t(Xhp) %*% Yhp) / (t(Xhp) %*% Xhp))

    # Find M_1 and M_2 for each smoother; requires finding \hat{x}_HP via each as well,
    # and \hat{1}_HP and \hat{y}_HP
    sim3[[k]][[1]][1:3] <- computeStats(S1, Xlp, Xhp, Ylp, Yhp)
    sim3[[k]][[2]][1:3] <- computeStats(S2, Xlp, Xhp, Ylp, Yhp)
    sim3[[k]][[3]][1:3] <- computeStats(S3, Xlp, Xhp, Ylp, Yhp)
    sim3[[k]][[4]][1:3] <- computeStats(S4, Xlp, Xhp, Ylp, Yhp)
    sim3[[k]][[5]][1:3] <- computeStats(S5, Xlp, Xhp, Ylp, Yhp)
    sim3[[k]][[6]][1:3] <- computeStats(S6, Xlp, Xhp, Ylp, Yhp)
    sim3[[k]][[7]][1:3] <- computeStats(S7, Xlp, Xhp, Ylp, Yhp)
    sim3[[k]][[8]][1:3] <- computeStats(S8, Xlp, Xhp, Ylp, Yhp)

    dat <- as.data.frame(cbind(X, Y + yMean, tArr))
    names(dat) <- c("X", "lY", "tArr")

    modelNS1 <- paste0("lY ~ X + basis1")
    modelNS2 <- paste0("lY ~ X + basis2")
    modelSP1 <- paste0("lY ~ X + basis3")
    modelSP2 <- paste0("lY ~ X + basis4")
    modelSP3 <- paste0("lY ~ X + basis5 - 1")  # remove intercept from model
    modelSP4 <- paste0("lY ~ X + basis6 - 1")  # remove intercept from model
    modelSP5 <- paste0("lY ~ X + basis7")
    modelSP6 <- paste0("lY ~ X + basis8")

    fitNS1 <- gam(as.formula(modelNS1), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitNS2 <- gam(as.formula(modelNS2), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP1 <- gam(as.formula(modelSP1), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP2 <- gam(as.formula(modelSP2), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP3 <- gam(as.formula(modelSP3), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP4 <- gam(as.formula(modelSP4), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP5 <- gam(as.formula(modelSP5), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP6 <- gam(as.formula(modelSP6), family = gaussian, data = dat, 
                  na.action = na.omit)

    sim3[[k]][[1]][4] <- list(summary.glm(fitNS1)$coefficients[1:2])
    sim3[[k]][[2]][4] <- list(summary.glm(fitNS2)$coefficients[1:2])
    sim3[[k]][[3]][4] <- list(summary.glm(fitSP1)$coefficients[1:2])
    sim3[[k]][[4]][4] <- list(summary.glm(fitSP2)$coefficients[1:2])
    sim3[[k]][[5]][4] <- list(c(summary.glm(fitSP3)$coefficients[2] * mean(basis5[, 1]), 
                                           summary.glm(fitSP3)$coefficients[1]))
    sim3[[k]][[6]][4] <- list(c(summary.glm(fitSP4)$coefficients[2] * mean(basis6[, 1]), 
                                           summary.glm(fitSP4)$coefficients[1]))
    sim3[[k]][[7]][4] <- list(summary.glm(fitSP5)$coefficients[1:2])
    sim3[[k]][[8]][4] <- list(summary.glm(fitSP6)$coefficients[1:2])
}
cat("\n")

save(file = "./data/simFullModelsResults3.RData", sim3)

################################################################################
#
#  Simulation 4:
#  * exponential decay variance
#  * make LP correlated and HP correlated, but ...
#  * make LP correlation = 2? All about scale factor - correlation 1, but half the
#    magnitude in X (Y = 2X)
#
sim4 <- vector("list", nSim)  # [[1]] is the general stuff, then [[2]] through [[5]] is for the four smoothers
sim4 <- lapply(sim4, FUN = function(x){y <- vector("list", 9)
                                       x <- lapply(y, FUN = function(z){ z <- vector("list", 4) })
                                      })
periodsLPx <- c(183, 105, 75, 68)   
ampsLPx <- c(2, 2, 2, 2)
phasesLPx <- c(0, 0, 0, 0)

periodsLPy <- c(183, 105, 75, 68)
ampsLPy <- c(4, 4, 4, 4)
phasesLPy <- c(180, 180, 180, 180)

periodsHP <- c(2, 5, 14)
ampsHP <- c(2, 2, 2)
phasesHP <- c(0, 0, 0)

Wcut <- min(which(freqs > 1/60/dT))   # cut-off of 60 days, 2 months

cat("Simulating: ")
for(k in 1:nSim) {
    if((k %% 25) == 0) { cat(".") }

    # noise background has variance (2*exp(-1e6*freqs) + 1) * sqrt(1/nMesh)
    noiseVar <- (2 * exp(-1e6 * freqs) + 1) * sqrt(1 / nMesh)

    # simulate noise background - low frequency random, high frequency shared
    # (this ensures that the 'true' correlation coefficient is consistent)
    Xlp <- sineMat[, 1:(Wcut-1)] %*% rnorm(n = Wcut - 1, sd = noiseVar[1:(Wcut-1)])
    Ylp <- sineMat[, 1:(Wcut-1)] %*% rnorm(n = Wcut - 1, sd = sqrt(2) * noiseVar[1:(Wcut-1)])
    Xhp <- sineMat[, Wcut:nMesh] %*% rnorm(n = (nMesh - Wcut + 1), sd = noiseVar[Wcut:(nMesh)])
    Yhp <- Xhp

    # add low pass sinusoids in
    for(j in 1:length(periodsLPx)) {
      Xlp <- Xlp + ampsLPx[j] * sin(2 * pi * tArr / periodsLPx[j] + phasesLPx[j])
      Ylp <- Ylp + ampsLPy[j] * sin(2 * pi * tArr / periodsLPy[j] + phasesLPy[j])
    }  

    # add high pass sinusoids in
    for(j in 1:length(periodsHP)) {
      Xhp <- Xhp + ampsHP[j] * sin(2 * pi * tArr / periodsHP[j] + phasesHP[j])
      Yhp <- Yhp + ampsHP[j] * sin(2 * pi * tArr / periodsHP[j] + phasesHP[j])
    }  

    X <- Xlp + Xhp
    Y <- Ylp + Yhp

    # add a clip to keep the variance from spinning out of control
    Y[which(Y > quantile(Y, 0.99))] <- quantile(Y, 0.99)
    Y[which(Y < quantile(Y, 0.01))] <- quantile(Y, 0.01)

    sim4[[k]][[9]] <- c(cor(Xlp, Ylp), cor(Xhp, Yhp), cor(X, Y), (t(Xhp) %*% Yhp) / (t(Xhp) %*% Xhp))

    # Find M_1 and M_2 for each smoother; requires finding \hat{x}_HP via each as well,
    # and \hat{1}_HP and \hat{y}_HP
    sim4[[k]][[1]][1:3] <- computeStats(S1, Xlp, Xhp, Ylp, Yhp)
    sim4[[k]][[2]][1:3] <- computeStats(S2, Xlp, Xhp, Ylp, Yhp)
    sim4[[k]][[3]][1:3] <- computeStats(S3, Xlp, Xhp, Ylp, Yhp)
    sim4[[k]][[4]][1:3] <- computeStats(S4, Xlp, Xhp, Ylp, Yhp)
    sim4[[k]][[5]][1:3] <- computeStats(S5, Xlp, Xhp, Ylp, Yhp)
    sim4[[k]][[6]][1:3] <- computeStats(S6, Xlp, Xhp, Ylp, Yhp)
    sim4[[k]][[7]][1:3] <- computeStats(S7, Xlp, Xhp, Ylp, Yhp)
    sim4[[k]][[8]][1:3] <- computeStats(S8, Xlp, Xhp, Ylp, Yhp)

    dat <- as.data.frame(cbind(X, Y + yMean, tArr))
    names(dat) <- c("X", "lY", "tArr")

    modelNS1 <- paste0("lY ~ X + basis1")
    modelNS2 <- paste0("lY ~ X + basis2")
    modelSP1 <- paste0("lY ~ X + basis3")
    modelSP2 <- paste0("lY ~ X + basis4")
    modelSP3 <- paste0("lY ~ X + basis5 - 1")  # remove intercept from model
    modelSP4 <- paste0("lY ~ X + basis6 - 1")  # remove intercept from model
    modelSP5 <- paste0("lY ~ X + basis7")
    modelSP6 <- paste0("lY ~ X + basis8")

    fitNS1 <- gam(as.formula(modelNS1), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitNS2 <- gam(as.formula(modelNS2), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP1 <- gam(as.formula(modelSP1), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP2 <- gam(as.formula(modelSP2), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP3 <- gam(as.formula(modelSP3), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP4 <- gam(as.formula(modelSP4), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP5 <- gam(as.formula(modelSP5), family = gaussian, data = dat, 
                  na.action = na.omit)
    fitSP6 <- gam(as.formula(modelSP6), family = gaussian, data = dat, 
                  na.action = na.omit)

    sim4[[k]][[1]][4] <- list(summary.glm(fitNS1)$coefficients[1:2])
    sim4[[k]][[2]][4] <- list(summary.glm(fitNS2)$coefficients[1:2])
    sim4[[k]][[3]][4] <- list(summary.glm(fitSP1)$coefficients[1:2])
    sim4[[k]][[4]][4] <- list(summary.glm(fitSP2)$coefficients[1:2])
    sim4[[k]][[5]][4] <- list(c(summary.glm(fitSP3)$coefficients[2] * mean(basis5[, 1]), 
                                           summary.glm(fitSP3)$coefficients[1]))
    sim4[[k]][[6]][4] <- list(c(summary.glm(fitSP4)$coefficients[2] * mean(basis6[, 1]), 
                                           summary.glm(fitSP4)$coefficients[1]))
    sim4[[k]][[7]][4] <- list(summary.glm(fitSP5)$coefficients[1:2])
    sim4[[k]][[8]][4] <- list(summary.glm(fitSP6)$coefficients[1:2])

}
cat("\n")

save(file = "./data/simFullModelsResults4.RData", sim4)

