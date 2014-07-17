################################################################################
#
#   Figures 8, 9 and 10 for paper  (repeat of 1, 2 and 3, but with SLP2-12)
#
#   Wesley S. Burr, Glen Takahara and Hwashin H. Shin 
#
#   Last modified: June 15, 2014
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
#  Prepares a collapsed NMMAPSdata database
#
#  Using Commands:
#  > install.packages("NMMAPSdata", contriburl = "http://www.ihapss.jhsph.edu/data/NMMAPS/R/", type = "source")
#  > buildDB(procFunc = collapseEndpoints)
#  > registerDB("collapseEndpoints")
#  > loadCity("chic")   # ``chic'' dataframe now available for use

load("./data/chic.RData")   
   # NMMAPSdata has been broken for a year now.
   # It has been reported to Roger Peng several times, but nothing has been done
   # chic.RData was extracted from an old copy of the package and used

load("./data/chicagoBases.RData")
for(j in 1:8) {
  assign(paste("basis", j, sep = ""), bases[[j]])
}

#############################################################################################
#
#  Figures chicagoResidAutocorB.pdf, chicagoResidSpectrumB.pdf and chicagoResidTimeDomainB.pdf
#  * shows residual Barlett autocorrelations from slepian-time GAM applied to 
#    data from Chicago, Il., using timespan 1987--2000. 
#  * shows corresponding power spectrum (multitaper estimated)
#  * shows time-domain representation of power from 0 to 1.92e-7 Hz (60 to infinity days period)
#
#  * DID NOT REMOVE OUTLIER FROM 1995
#

library("gam")
library("multitaper")      
library("tsinterp")       # for linearly interpolating PM10 data for demonstration purposes

# interpolate the PM10 data
chicPM10 <- chic[, "pm10tmean"]
chicPM10[1] <- chicPM10[2] <- chicPM10[3] # manually extrapolate first two points; interpolator only works on internal values
mask <- chicPM10
mask[!is.na(mask)] <- 1
nMiss <- length(which(is.na(mask)))
blocks <- findBlocks(mask)   # findBlocks and linInt are from the \code{tsinterp} package
pm10tmeanGF <- linInt(chicPM10, blocks)
time <- seq(1, length(chicPM10), 1)
chic <- cbind(chic, pm10tmeanGF, time)

timeAxis <- ISOdate(substr(chic$date, 1, 4), substr(chic$date, 5, 6), substr(chic$date, 7, 8))

# setup gam model: using INTERPOLATED DATA (no dramatic change from non-interpolated; can be 
#   confirmed by swapping the following lines re: pollutant

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

resid.acf.old <- acf(fitOld$residuals, lag.max = 50, plot = FALSE)

xaxis <- resid.acf$lag[-1]
n <- length(xaxis)

# pdf(file="figures/chicagoResidAutocorB.pdf", width=6, height=4)
postscript(file="figures/fig8-chicagoResidAutocorB.eps", width=6, height=4,
           horizontal = FALSE, paper = 'special')
par(mar=c(4,4,1,1))
plot(xaxis, rep(NA, n), xlab = "Time Lag in Days", ylab = "Autocorrelation Estimate",
     xaxs = "i", xaxt = 'n', ylim = c(-0.11, 0.18))
axis(side = 1, at = c(1,5,10,20,30,40,50), labels = c(1,5,10,20,30,40,50))
abline(h=0)
lines(xaxis, resid.acf$acf[-1], type = "b", col = "black", lwd = 2)
lines(xaxis, resid.acf.old$acf[-1], type = "b", col = "blue", lty = 2)
legend(x = "topright", col = c("black", "blue"), lty = c(1, 2), 
       legend = c("Updated S-SLP2-12 Model", "Original S-NS-6 Model"))
dev.off()

# multitaper spectrum estimate of residuals

resid <- fit$residuals
resid.sp <- spec.mtm(resid, deltat=60*60*24, nw=4, k=7, plot=FALSE)  # default units: Hz, cycles/second
axisPeriod <- c(7,14,21,30,60,120,365.2425)  # periods in days that we are interested in
labelPeriod <- c(7,14,21,30,60,120,365)
axisFreq <- 1/(axisPeriod * 86400) # frequencies in Hz from those periods
axisFreq2 <- c(1:10,15,20,30,40,50)

# pdf(file="figures/chicagoResidSpectrumB.pdf", width=6, height=6)
postscript(file="figures/fig9-chicagoResidSpectrumB.eps", width=6, height=6,
           horizontal = FALSE, paper = 'special')
par(mar=c(4,4,4,1))
yLabExp <- expression(paste("Power Spectrum in ppb"^2, "/(cycle/year)", sep = ""))
plot(resid.sp$freq, resid.sp$spec, type="l", log="y", xlab="Frequency in cycles/year", ylab="", 
     xaxs="i", xlim=c(0, 1/(7*86400)), yaxt = 'n', xaxt = 'n')
mtext(side=2, line=2.5, yLabExp)
axis(side=1, at=1 / ((365.2425 * 86400) / axisFreq2), labels=axisFreq2)
axis(side=3, at=axisFreq, labels=labelPeriod)
mtext(side = 3, line = 2.5, "Period in Days")
axis(side=2, at=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000),
     labels = c(1,2,5,"1e1","2e1","5e1","1e2","2e2","5e2","1e3","2e3","5e3","1e4"))
abline(v=1/(60*86400), col="blue", lty=2)
dev.off()

# time domain representation
resid.low <- dpssFiltNA(N=length(resid), w=6/365.2425, xd=resid)
yrs <- ISOdate(seq(1987, 2000, 1), rep(1, 8), rep(1, 8))

# pdf(file="figures/chicagoResidTimeDomainB.pdf", width=6, height=4)
postscript(file="figures/fig10-chicagoResidTimeDomainB.eps", width=6, height=4,
           horizontal = FALSE, paper = 'special')
par(mar=c(4,4,1,1))
plot(timeAxis, resid, type="l", col="grey80", ylim=c(-0.3,0.5), 
     xlab="Time in Years", ylab="Residuals in log(ppb)", xaxs="i", xaxt = 'n')
axis(side = 1, at = yrs, labels = c(1987, "", "", 1990, "", 1992, "", 1994, "", 1996, "", 1998, "", 2000))
lines(timeAxis, resid.low, lwd=2)
dev.off()

