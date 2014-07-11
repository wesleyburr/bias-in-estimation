#
#  Uses realizations created by sim_fullModels.R to analyze performance of
#  various smoothers; results are Table 3 in the paper.
#

load("simFullModelsResults1.RDa")
load("simFullModelsResults2.RDa")
load("simFullModelsResults3.RDa")
load("simFullModelsResults4.RDa")

simulations <- vector("list", 4)

# each simX is a list of length nSim, each of the entries is a list of 9 elements,
# each of which is a list of 4. 
simulations[[1]] <- lapply(sim1, FUN = function(x) {
  matrix(data = c(unlist(sapply(x, FUN = function(x) { unlist(x)  })), rep(0, 12)), 
         ncol = 16, byrow = TRUE)
})
simulations[[2]] <- lapply(sim2, FUN = function(x) {
  matrix(data = c(unlist(sapply(x, FUN = function(x) { unlist(x)  })), rep(0, 12)), 
         ncol = 16, byrow = TRUE)
})
simulations[[3]] <- lapply(sim3, FUN = function(x) {
  matrix(data = c(unlist(sapply(x, FUN = function(x) { unlist(x)  })), rep(0, 12)), 
         ncol = 16, byrow = TRUE)
})
simulations[[4]] <- lapply(sim4, FUN = function(x) {
  matrix(data = c(unlist(sapply(x, FUN = function(x) { unlist(x)  })), rep(0, 12)), 
         ncol = 16, byrow = TRUE)
})
save(file="paperSimulations.RDa", simulations)

summaried <- lapply(simulations, FUN = function(x) { Reduce("+", x) / length(x) })
summaried <- lapply(summaried, FUN = function(x) { 
                    row.names(x) <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "X")
                    colnames(x) <- c("alpha", "gamma", paste("lm", 1:12, sep=""), "gam$mean", "gam$beta")
                    x })
summaried <- lapply(summaried, FUN = function(x) {
                    x[, 1] <- x[, 1] - 1.000000 
                    x
                   })
save(file="summariedPaper.RDa", summaried)

library("xtable")  # want rows 1,3,5,7,2,4,6,8 and columns 1,2,13,15,16
# Sim 1
xtable(summaried[[1]][c(1,3,5,7,2,4,6,8), c(1:2, 13, 15:16)], digits = 3, display = rep("e", 6))
# Sim 2
xtable(summaried[[2]][c(1,3,5,7,2,4,6,8), c(1:2, 13, 15:16)], digits = 3, display = rep("e", 6))
# Sim 3
xtable(summaried[[3]][c(1,3,5,7,2,4,6,8), c(1:2, 13, 15:16)], digits = 3, display = rep("e", 6))
# Sim 4
xtable(summaried[[4]][c(1,3,5,7,2,4,6,8), c(1:2, 13, 15:16)], digits = 3, display = rep("e", 6))

