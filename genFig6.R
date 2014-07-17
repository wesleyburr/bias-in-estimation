#
#  Generate Figure 6
#
library("multitaper")
library("splines")
library("MASS")

#  3 years of data for "zoomed" in effect
N <- 365 * 3
time <- 1:N
Nknot <- 18

basis1 <- ns(time, df=Nknot) # 6 knots per year
basis2 <- dpss(n = N, k = Nknot, nw = (Nknot+1) / 2)$v

S1 <- basis1 %*% ginv(t(basis1) %*% basis1) %*% t(basis1)
S2 <- basis2 %*% t(basis2)

oneV <- rep(1, N)
oneP1 <- S1 %*% oneV
oneP2 <- S2 %*% oneV

pMax <- max(oneP1, oneP2)

#
#  Generate Figure 6
#
# pdf(file = "figures/gibbsRipples.pdf", width = 6, height = 4)
postscript(file = "figures/gibbsRipples.eps", width = 6, height = 4,
           horizontal = FALSE, paper = 'special')
par(mar = c(4,4,0.5,0.5))
plot(1:N, oneV, type = "l", lwd = 2, xlab = "Time", ylab = "Magnitude", ylim = c(0.75, pMax))
lines(1:N, oneP1, lwd = 3, lty = 1, col = "black")
lines(1:N, oneP2, lwd = 1.5, lty = 1, col = "blue")
text(x = c(170, 500), y = c(0.85, 1.15), pos = 4, labels = c("Projection onto ns(18)", "Projection onto slp(18)"),
     col = c("black", "blue"))
dev.off()


