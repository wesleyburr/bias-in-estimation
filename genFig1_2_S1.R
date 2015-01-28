################################################################################
#
#   Figures 1, 2 for paper and S1 for supplement
#
#   Wesley S. Burr, Glen Takahara and Hwashin H. Shin 
#
#   Last modified: January 27, 2015
#
#   * notes for publication
#   * rendered PDF version of journal uses 11 pt CMbright (Computer Modern)
#   * so we use CM Sans font, 11 pt, true-size graphics, OR
#   * a 6.0in figure fits nicely in the page width of the journal, so we will
#     use package "extrafont" and "fontcm" to get access to Computer Modern Sans
#     for use in our figures, making them true-to-size at 6in, 11points
# 
#   * for some reason (unknown), generating a true-to-size PDF at 6 inches with
#     pointsize = 11 results in text that is too big. Replaced with pointsize = 9
#     as the best integer approximation.
#
################################################################################

################################################################################
#
#  Utility function which leverages discrete prolate spheroidal sequences 
#  to create a digital filter which returns the long timescale portion of the
#  input series. Requires a bandwidth 'w' which is measured against the sampling
#  units of the input series xd; e.g., if xd is measured daily, then w should be
#  also measured as cycles/day. 
#
#  Example Call:
#  > xd <- rnorm(100)   # daily data
#  > xd.low <- dpssFiltNA(N=length(xd), w=7/365, xd=xd)
#
"dpssFiltNA" <- function(N, w, xd) {
  require("multitaper")
  kmax = floor(max(2*N*w - 2, 2))
  dwv <- dpss(N, kmax, N*w)$v

  yk <- rep(NA, kmax)
  for(k in 1:kmax) {
    yk[k] <- sum(dwv[, k] * xd, na.rm=TRUE) 
  }
  xlp <- dwv %*% yk
  xlp
}

################################################################################
#
#  Prepares a collapsed NMMAPSdata database (only possible on Windows, running R 2.9 or earlier)
#
#  Using Commands:
#   > install.packages("NMMAPSdata", 
#   +   contriburl = "http://www.ihapss.jhsph.edu/data/NMMAPS/R/", type = "source")
#   > library("NMMAPSdata")
#   > buildDB(procFunc = collapseEndpoints)
#   > registerDB("collapseEndpoints")
#   > loadCity("chic")
#   > chic <- chic[, c("city", "date", "dow", "death", "cvd", 
#   +                  "tmpd", "o3tmean", "pm10tmean")]
#   > save(file = "chic.RData", chic)
#

#   NMMAPSdata used to have a source version available which could be recompiled
#   for a modern version of R. Unfortunately, the link on the website was removed
#   over a year ago (June 2013). Apparently the database was removed from public
#   access at the request of an external agency. This makes reproducibility more
#   difficult, but it may be possible to get a copy of the package by asking Roger
#   Peng directly. 
#

load("./data/chic.RData")   

load("./data/chicagoBases.RData")
for(j in 1:8) {
  assign(paste("basis", j, sep = ""), bases[[j]])
}

#############################################################################################
#
#  Figures chicagoResidAutocor.pdf, chicagoResidSpectrum.pdf and chicagoResidTimeDomain.pdf
#
#  * DID NOT REMOVE OUTLIER FROM 1995, which influences GAM regression
#

library("gam")
library("multitaper")      
library("extrafont")
loadfonts(device = "postscript")
loadfonts()

#  Code to interpolate PM10 missing values to allow for spectrum estimation
#
#  Note: file chic.RData included with this repository has had this interpolation done already
#
load("./data/chic.RData")
if(!("pm10tmeanGF" %in% names(chic))) {
    library("tsinterp")                       # package for interpolation; available from http://github.com/wesleyburr/tsinterp
    chicPM10 <- chic[, "pm10tmean"]
    chicPM10[1] <- chicPM10[2] <- chicPM10[3] # manually extrapolate first two points; interpolator only works on internal values
    mask <- chicPM10
    mask[!is.na(mask)] <- 1
    nMiss <- length(which(is.na(mask)))
    blocks <- findBlocks(mask)   # findBlocks and linInt are from the \code{tsinterp} package
    pm10tmeanGF <- linInt(chicPM10, blocks)
    time <- seq(1, length(chicPM10), 1)
    chic <- cbind(chic, pm10tmeanGF, time)
    save(file = "chic.RData", chic)
}  # end interpolation

timeAxis <- ISOdate(substr(chic$date, 1, 4), substr(chic$date, 5, 6), substr(chic$date, 7, 8))

# setup gam model: using INTERPOLATED DATA (no dramatic change from non-interpolated; can be 
# confirmed by replacing following pollutant assignment

pollutant <- "pm10tmeanGF"

cause <- "death"
nYr <- length(unique(substr(chic$date, 1, 4)))
df.Temp <- 3 
df.Time <- 14 * 6

# basis6 is a 12-df/year basis, used in Section 5.2
modelFormula <- paste0(cause, " ~ ", pollutant, " + dow + ns(tmpd, df = ", df.Temp, ") + basis6")
modelFormulaOrig <- paste0(cause, " ~ ", pollutant, " + dow + ns(time, df = ", df.Time,") + ns(tmpd, df = ", df.Temp, ")")

fit <- gam(as.formula(modelFormula), family = poisson, data = chic,
               na.action = na.omit)
fitOld <- gam(as.formula(modelFormulaOrig), family = poisson, data = chic, na.action = na.omit)

# autocorrelation of residuals
resid <- fit$residuals
resid.acf <- acf(resid, lag.max=50, plot=FALSE)

resid.old <- fitOld$residuals
resid.acf.old <- acf(fitOld$residuals, lag.max = 50, plot = FALSE)

xaxis <- resid.acf$lag[-1]
n <- length(xaxis)
N <- length(resid)

acfT <- resid.acf$acf[-1]
# Classical confidence interval
cI <- qnorm(0.975) / sqrt(N)

#
# pdf(file="figures/chicagoResidAutocor.pdf", width=6, height=4)
#
#  Figure 1 - modified January 26 to be the old figure 8, modified; compression of figures 1 and 8 together.
#
postscript(file="figures/fig1-chicagoResidAutocor.eps", width=6, height=4,
           horizontal = FALSE, paper = 'special', family = "CM Sans", pointsize = 9)
par(mar=c(4,4,1,1))
plot(xaxis, rep(NA, n), xlab = "Time Lag in Days", ylab = "Autocorrelation Estimate",
     xaxs = "i", xaxt = 'n', ylim = c(-0.10, 0.15))
axis(side = 1, at = c(1,5,10,20,30,40,50), labels = c(1,5,10,20,30,40,50))
rect(1.10, -cI, n-0.05, cI, col = "grey90", border = NA)
abline(h=0)
lines(xaxis, resid.acf$acf[-1], type = "l", col = "black", lwd = 2)
points(xaxis, resid.acf$acf[-1], pch = 18) 
lines(xaxis, resid.acf.old$acf[-1], type = "l", col = "grey30", lty = 2, lwd = 2)
points(xaxis, resid.acf.old$acf[-1], pch = 19, col = "grey30")
legend(x = "topright", legend = c("Model 2 - S-NS-6 Model", "Model 9 - S-SLP2-12 Model", 
       "95% Confidence Interval \n for White Noise"),
       col = c("grey30", "black", "grey90"),
       lty = c(1, 2, 0), 
       lwd = c(2, 2, 0), 
       pch = c(18, 19, 15),
       pt.bg = c("black", "grey30", "grey90"),
       pt.cex = c(1, 1, 1),
       cex = 1)
dev.off()

pdf(file="figures/fig1-chicagoResidAutocor.pdf", width=6, height=4,
    paper = 'special', family = "CM Sans", pointsize = 9)
par(mar=c(4,4,1,1))
plot(xaxis, rep(NA, n), xlab = "Time Lag in Days", ylab = "Autocorrelation Estimate",
     xaxs = "i", xaxt = 'n', ylim = c(-0.10, 0.15))
axis(side = 1, at = c(1,5,10,20,30,40,50), labels = c(1,5,10,20,30,40,50))
rect(1.10, -cI, n-0.05, cI, col = "grey90", border = NA)
abline(h=0)
lines(xaxis, resid.acf$acf[-1], type = "l", col = "black", lwd = 2)
points(xaxis, resid.acf$acf[-1], pch = 18) 
lines(xaxis, resid.acf.old$acf[-1], type = "l", col = "grey30", lty = 2, lwd = 2)
points(xaxis, resid.acf.old$acf[-1], pch = 19, col = "grey30")
legend(x = "topright", legend = c("Model 2 - S-NS-6 Model", "Model 9 - S-SLP2-12 Model", 
       "95% Confidence Interval \n for White Noise"),
       col = c("grey30", "black", "grey90"),
       lty = c(1, 2, 0), 
       lwd = c(2, 2, 0), 
       pch = c(18, 19, 15),
       pt.bg = c("black", "grey30", "grey90"),
       pt.cex = c(1, 1, 3), cex = 1)
dev.off()
embed_fonts("figures/fig1-chicagoResidAutocor.pdf", outfile = "figures/fig1-chicagoResidAutocor_embed.pdf")


# multitaper spectrum estimate of residuals
resid <- fit$residuals
resid.sp <- spec.mtm(resid, deltat=60*60*24, nw=4, k=7, plot=FALSE)  # default units: Hz, cycles/second
resid.Old <- fitOld$residuals
resid.Old.sp <- spec.mtm(resid.Old, deltat=60*60*24, nw=4, k=7, plot=FALSE)
axisPeriod <- c(7,14,21,30,60,120,365.2425)  # periods in days that we are interested in
labelPeriod <- c(7,14,21,30,60,120,365)
axisFreq <- 1/(axisPeriod * 86400) # frequencies in Hz from those periods
axisFreq2 <- c(1:10,15,20,30,40,50)

# pdf(file="figures/fig2-chicagoResidSpectrum.pdf", width=6, height=6)
postscript(file="figures/fig2-chicagoResidSpectrum.eps", width=9, height=7,
           horizontal = FALSE, paper = 'special', family = "CM Sans")
par(mar=c(4,4,4,1))
yLabExp <- expression(paste("Power Spectrum in ppb"^2, "/(cycle/year)", sep = ""))
plot(resid.Old.sp$freq, resid.Old.sp$spec, type="l", log="y", xlab="Frequency in cycles/year", ylab="", 
     xaxs="i", xlim=c(0, 1/(14*86400)), yaxt = 'n', xaxt = 'n', lwd = 2, col = "grey40", lty = 2)
lines(resid.sp$freq, resid.sp$spec, type="l", lwd=2, col = "black")
mtext(side=2, line=2.5, yLabExp)
axis(side=1, at=1 / ((365.2425 * 86400) / axisFreq2), labels=axisFreq2)
axis(side=3, at=axisFreq, labels=labelPeriod)
mtext(side = 3, line = 2.5, "Period in Days")
axis(side=2, at=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000),
     labels = c(1,2,5,"1e1","2e1","5e1","1e2","2e2","5e2","1e3","2e3","5e3","1e4"))
abline(v=1/(60.9*86400), col="red", lty=3, lwd = 2.5)
legend(x = "bottomright", lty = c(2, 1, 3), lwd = c(2,2,2.5), col = c("grey40", "black", "red"), 
       legend = c("Model 2 - S-NS-6", "Model 9 - S-SLP2-12", "Ideal Band-Edge at 6 cycles/year"))
dev.off()

pdf(file="figures/fig2-chicagoResidSpectrum.pdf", width=6, height=4,
    paper = 'special', family = "CM Sans", pointsize = 9)
par(mar=c(4,4,4,1))
yLabExp <- expression(paste("Power Spectrum in ppb"^2, "/(cycle/year)", sep = ""))
plot(resid.Old.sp$freq, resid.Old.sp$spec, type="l", log="y", xlab="Frequency in cycles/year", ylab="", 
     xaxs="i", xlim=c(0, 1/(14*86400)), yaxt = 'n', xaxt = 'n', lwd = 2, col = "grey40", lty = 2)
lines(resid.sp$freq, resid.sp$spec, type="l", lwd=2, col = "black")
mtext(side=2, line=2.5, yLabExp)
axis(side=1, at=1 / ((365.2425 * 86400) / axisFreq2), labels=axisFreq2)
axis(side=3, at=axisFreq, labels=labelPeriod)
mtext(side = 3, line = 2.5, "Period in Days")
axis(side=2, at=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000),
     labels = c(1,2,5,"1e1","2e1","5e1","1e2","2e2","5e2","1e3","2e3","5e3","1e4"))
abline(v=1/(60.9*86400), col="red", lty=3, lwd = 2.5)
legend(x = "bottomright", lty = c(2, 1, 3), lwd = c(2,2,2.5), col = c("grey40", "black", "red"), 
       legend = c("Model 2 - S-NS-6", "Model 9 - S-SLP2-12", "Ideal Band-Edge at 6 cycles/year"))
dev.off()
embed_fonts("figures/fig2-chicagoResidSpectrum.pdf", outfile = "figures/fig2-chicagoResidSpectrum_embed.pdf")

# time domain representation
resid.low <- dpssFiltNA(N=length(resid), w=6/365.2425, xd=resid)
resid.old.low <- dpssFiltNA(N=length(resid.old), w=6/365.2425, xd=resid.old)
yrs <- ISOdate(seq(1987, 2000, 1), rep(1, 8), rep(1, 8))

# pdf(file="figures/figSupp1-chicagoResidTimeDomain.pdf", width=6, height=4)
postscript(file="figures/figSupp1-chicagoResidTimeDomain.eps", width=6, height=4,
           horizontal = FALSE, paper = 'special', family = "CM Sans", pointsize = 9)
par(mar=c(4,4,1,1))
plot(timeAxis, resid, type="l", col="grey80", ylim=c(-0.3,0.5), 
     xlab="Time in Years", ylab="Residuals in log(ppb)", xaxs="i", xaxt = 'n')
axis(side = 1, at = yrs, labels = c(1987, "", "", 1990, "", 1992, "", 1994, "", 1996, "", 1998, "", 2000))
lines(timeAxis, resid.old.low, lwd = 3, lty = 1, type = "l", col = "grey50")
lines(timeAxis, resid.low, lwd=3)
legend(x = "topright", lty = c(1, 1, 1), lwd = c(2, 3, 3), col = c("grey80", "grey50", "black"),
       legend = c("Residuals", "Model 2 - Long Time-scale Residuals", "Model 9 - Long Time-scale Residuals"))
dev.off()

pdf(file="figures/figSupp1-chicagoResidTimeDomain.pdf", width=6, height=4,
    paper = 'special', family = "CM Sans", pointsize = 9)
par(mar=c(4,4,1,1))
plot(timeAxis, resid, type="l", col="grey80", ylim=c(-0.3,0.5), 
     xlab="Time in Years", ylab="Residuals in log(ppb)", xaxs="i", xaxt = 'n')
axis(side = 1, at = yrs, labels = c(1987, "", "", 1990, "", 1992, "", 1994, "", 1996, "", 1998, "", 2000))
lines(timeAxis, resid.old.low, lwd = 3, lty = 1, type = "l", col = "grey50")
lines(timeAxis, resid.low, lwd=3)
legend(x = "topright", lty = c(1, 1, 1), lwd = c(2, 3, 3), col = c("grey80", "grey50", "black"),
       legend = c("Residuals", "Model 2 - Long Time-scale Residuals", "Model 9 - Long Time-scale Residuals"))
dev.off()
embed_fonts("figures/figSupp1-chicagoResidTimeDomain.pdf", outfile = "figures/figSupp1-chicagoResidTimeDomain_embed.pdf")


