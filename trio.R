library(genio) # Read plink format
library(OpenMx) # Fit model
library(gaston) # Compute GRM

# Read data
dat = read.bed.matrix("simtrio")

# Compute GRM
A = GRM(dat)

# Subset GRM blocks
mid = 1:1000
pid = 1001:2000
oid = 2001:3000
Amm = A[mid, mid]
App = A[pid, pid]
Aoo = A[oid, oid]
Dpm = A[pid, mid] + A[mid, pid]
Dom = A[oid, mid] + A[mid, oid]
Dop = A[oid, pid] + A[pid, oid]

# Set up model
yX = cbind(dat@ped$pheno[oid], 1)
colnames(yX) = c("y", "x")
K = 1000 # Number of trios
modmx = mxModel("trio",
              mxMatrix("Lo", 3, 3, T, sqrt(diag(var(yX[, 1]) / 4, 3)),
                       c("l11", "l21", "l31", "l22", "l32", "l33"), name = "L"), # Cholesky factor of genetic covariance matrix
              mxAlgebra(L %*% t(L), name = "Sg"),
              mxMatrix("Fu", 1, 1, T, var(yX[, 1]) / 4, "se", name = "Se"), # Residual variance
              mxMatrix("Sy", K, K, F, Amm, name = "Amm"),
              mxMatrix("Sy", K, K, F, App, name = "App"),
              mxMatrix("Sy", K, K, F, Aoo, name = "Aoo"),
              mxMatrix("Sy", K, K, F, Dpm, name = "Dpm"),
              mxMatrix("Sy", K, K, F, Dom, name = "Dom"),
              mxMatrix("Sy", K, K, F, Dop, name = "Dop"),
              mxMatrix("Id", K, K, name = "I"),
              mxData(yX, "raw", sort = F),
              # Model implied covariance
              mxAlgebra(Sg[1, 1] * Amm + Sg[2, 2] * App + Sg[3, 3] * Aoo +
                          Sg[2, 1] * Dpm + Sg[3, 1] * Dom + Sg[3, 2] * Dop +
                          se * I, name = "V"),
              mxExpectationGREML("V", dataset.is.yX = T),
              mxFitFunctionGREML())

# Fit the model
modmx = mxRun(modmx)

# Look at results
summary(modmx)
mxEval(Sg, modmx) # Genetic (co)variances