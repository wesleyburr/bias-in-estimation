#
#  Generate Figure 4 for paper
#

library("multitaper")
library("splines")
library("MASS")

# 10 years of data, approximately
time <- 1:3650
Nknot <- 60
test <- ns(time, df=Nknot) # 6 knots per year
spec <- matrix(data=0, nrow=16385, ncol=Nknot)
for(j in 1:Nknot) { spec[,j] <- spec.mtm(test[,j],deltat=86400,nFFT=2*16384,plot=FALSE)$spec }

freqs <- seq(0,0.5,1/2/16384)
theory <- rep(1,min(which(freqs > 6/365))+5)
cutoff <- min(which(freqs > 6/365))+5
theory <- c(theory, rep(1e-21, length(freqs)-length(theory)))
freqs <- freqs*365

# ** only using vectors 1-58, because 59 and 60 have _very_ poor spectral properties
#    and make the theoretical transfer function appear worse than the actual performance
# ** this _helps_ the case of the natural cubic regression splines

# projection matrix has eigenvalues 1 and 0; the eigenvalue of H^T*H gives the magnitude correction factor
etf <- eigen(t(test[, 1:58]) %*% test[, 1:58])$value[1]  
# magnitude transfer function: sum of the individual projection vectors, normalized by etf
tf <- rowSums(spec[, 1:58]) / (etf * 86400)

# estimate via averaged simulation on zero-mean data
S <- test %*% ginv(t(test) %*% test) %*% t(test)
M <- 1000; N <- 3650
spec2 <- matrix(data = 0, nrow = 16385, ncol = M)
for(j in 1:M) {
  simDat <- rnorm(n = N, sd = 2)
  fit <- S %*% simDat
  spec2[, j] <- spec.mtm(fit, deltat = 86400, nFFT = 2*16384, plot = FALSE)$spec /
                spec.mtm(simDat, deltat = 86400, nFFT = 2*16384, plot = FALSE)$spec
}
# average white-noise projections together
tf2 <- rowSums(spec2) / M 

#  smooth transfer functions slightly
moveAvg <- function(x, n = 20) { filter(x, rep(1/n, n), sides=2)}

tf <- moveAvg(tf)
tf2 <- moveAvg(tf2)

PerLab <- c(0:7, 14, 21, 30, 60, 120, 365)
atPerLab <- 365 / PerLab   # already scaled the freq array to be in cycles/year, not cycles/day
yAxisP <- 10^(seq(-15, 1, 1))
yAxisL <- paste0("1e-", 15:1)

#
#  Generate Figure 4
#
pdf(file="../figures/transferFuncCubic.pdf",width=6,height=5)
par(mar=c(4,4,4,0.5))
plot(freqs, tf, type="l", col="black", lwd=2, log="y", xlim=c(0,20),
     ylim=c(1e-14,1.1e0), xlab="Frequency in Cycles/Year", 
     ylab="Magnitude Transfer Functions (TFs)",
     xaxs = "i", xaxt='n', yaxt = 'n')
lines(freqs, theory, type="l", col="red", lty=2, lwd=2)
lines(freqs, tf2, type = "l", col = "grey60", lwd = 2)
axis(side = 3, line = 0, at = atPerLab, labels = PerLab)
axis(side = 2, line = 0, at = yAxisP, labels = yAxisL)
axis(side = 1, line = 0, at = PerLab, labels = PerLab)
mtext("Period in Days", side = 3, line = 2)
text(x = c(1, 1, 1, 10, 6.1), y = c(1e-9, 1e-10, 1e-11, 1e-5, 1e-1), 
     labels = c("Mean TF of", "Basis Vectors", "1-58", "Mean Simulated White Noise TF", "Ideal TF"),
     col = c("black", "black", "black", "grey60", "red"), pos = 4)
dev.off()

