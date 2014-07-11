#
#  Generate Figure 5
#
library("multitaper")
library("splines")
library("MASS")

N <- 3650 # 10 years of data, approximately
time <- 1:N
Nknot <- 60
fftn <- 16384
ns60  <- ns(time, df = Nknot) # 6 knots per year
ns120 <- ns(time, df = Nknot * 2)

specNS60  <- matrix(data = 0, nrow = fftn + 1, ncol = Nknot)
specNS120 <- matrix(data = 0, nrow = fftn + 1, ncol = Nknot * 2)

for(j in 1:Nknot) {
  specNS60[, j]  <- spec.mtm( ns60[, j], deltat = 86400, nFFT = 2 * fftn, plot=FALSE)$spec }
for(j in 1:(Nknot * 2)) {
  specNS120[, j] <- spec.mtm(ns120[, j], deltat = 86400, nFFT = 2 * fftn, plot = FALSE)$spec }

freqs <- seq(0, 0.5, 1/2/fftn)
cutoff <- min(which(freqs > 6/365))+5          # about 6 cycles/year
theory <- rep(1, cutoff)
theory <- c(theory, rep(1e-21, length(freqs)-length(theory)))
freqs <- freqs * 365    # changes freqs from cycles/day to cycles/year

# ** only using vectors 1-58, because 59 and 60 have _very_ poor spectral properties
#    and make the result look even worse than actual performance

# projection matrix has eigenvalues 1 and 0; the eigenvalue of H^T*H gives the magnitude correction factor
etfNS60  <- eigen(t(ns60[, 1:58]) %*% ns60[, 1:58])$value[1]  
etfNS120 <- eigen(t(ns120[, 1:118]) %*% ns120[, 1:118])$value[1]  

# magnitude transfer function: sum of the individual projection vectors, normalized by etf
tfNS60  <- rowSums(specNS60[, 1:58]) / (etfNS60 * 86400)
tfNS120 <- rowSums(specNS120[, 1:118]) / (etfNS120 * 86400)

#
#  Add Slepians the same way
#
slp1 <- multitaper::dpss(n = N, k = Nknot, nw = (Nknot+2)/2)$v # 60 basis vectors, same as option 1 previously, W \approx 120 days
slp2 <- multitaper::dpss(n = N, k = Nknot * 2, nw = ((Nknot * 2 + 2) / 2))$v   # 120 basis vectors, W \approx 60 days

specS1 <- matrix(data=0, nrow=fftn + 1, ncol=Nknot)
specS2 <- matrix(data=0, nrow=fftn + 1, ncol=Nknot * 2)

for(j in 1:Nknot) {
  specS1[, j] <- spec.mtm(slp1[, j], deltat = 86400, nFFT = 2 * fftn, plot = FALSE)$spec }
for(j in 1:(Nknot * 2)) {
  specS2[, j] <- spec.mtm(slp2[, j], deltat = 86400, nFFT = 2 * fftn, plot = FALSE)$spec }

# magnitude transfer function: sum of the individual projection vectors, no normalization -- eigenvalues = 1
tfSLP60  <- rowSums(specS1[, 1:60]) / 86400
tfSLP120 <- rowSums(specS2[, 1:120]) / 86400 

################################################################################
#
#  10 point moving average for each, cleans up plots for comparison
#
moveAvg <- function(x, n = 20) { filter(x, rep(1/n, n), sides=2)}

tfNS60   <- moveAvg(tfNS60)
tfNS120  <- moveAvg(tfNS120)
tfSLP60  <- moveAvg(tfSLP60)
tfSLP120 <- moveAvg(tfSLP120) 

PerLab <- c(0,1,2,3,4,5,6,7,14,21,30,60,120, 365)
atPerLab <- 365/(PerLab)   # already scaled the freq array to be in cycles/year, not cycles/day
yAxisP <- 10^(seq(-15,1,1))
yAxisL <- c("1e-15", "1e-14", "1e-13", "1e-12","1e-11","1e-10","1e-9","1e-8","1e-7","1e-6",
            "1e-5","1e-4","1e-3","1e-2", "1e-1","1", "1e1")
#
#  Generate Figure 5 for paper
#
pdf(file="../figures/compareTransferFunctAll.pdf",width=6,height=5)
par(mar=c(4,4,4,0.5))
plot(freqs, theory, type="l", col="red", lwd=2, lty = 3, log="y", xlim=c(0,20),
     ylim=c(1e-14,1.1e0), xlab="Frequency in Cycles/Year", ylab="Magnitude Transfer Functions (TFs)",
     xaxt = 'n', yaxt = 'n', xaxs = "i")
lines(freqs, tfNS60,   type = "l", col = "black",  lwd = 2)
lines(freqs, tfNS120,  type = "l", col = "black",  lwd = 2, lty = 2)
lines(freqs, tfSLP60,  type = "l", col = "grey60", lwd = 2)
lines(freqs, tfSLP120, type = "l", col = "grey60", lwd = 2, lty = 2)

axis(side = 3, line = 0, at = atPerLab, labels = PerLab)
axis(side = 2, line = 0, at = yAxisP, labels = yAxisL)
axis(side = 1, line = 0, at = PerLab, labels = PerLab)
mtext("Period in Days", side = 3, line = 2)

legend(x = "topright", lty = c(1,2,1,2, 3), col = c("grey60", "grey60", "black", "black", "red"),
       legend = c("slp(60)", "slp(120)", "ns(60)", "ns(120)", "Ideal"))
dev.off()

