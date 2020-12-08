rm(list = ls())

library(genio)
library(MASS)

set.seed(080318)

setwd("C:/Users/espenmei/OneDrive - Universitetet i Oslo/Prosjekter/GCTA/MFGCTA/BehaviorGenetics/V3/github")

K = 1000
M = 1000
P = runif(M, 0.1, 0.2)
var_g = matrix(c(3, 0.7, 1, 
                 0.7, 2, 0.5,
                 1, 0.5, 1), 3, 3, byrow = T)
var_e = 1

# Sample genotypes
Xm_t = matrix(rbinom(K * M, 1, P), K, M, byrow = T)
Xm_n = matrix(rbinom(K * M, 1, P), K, M, byrow = T)
Xp_t = matrix(rbinom(K * M, 1, P), K, M, byrow = T)
Xp_n = matrix(rbinom(K * M, 1, P), K, M, byrow = T)
Xm = Xm_t + Xm_n
Xp = Xp_t + Xp_n
Xo = Xm_t + Xp_t
X = rbind(Xm, Xp, Xo)

# Simulate phenotypes
Zm = scale(Xm, 2 * P, sqrt(2 * P * (1 - P)))
Zp = scale(Xp, 2 * P, sqrt(2 * P * (1 - P)))
Zo = scale(Xo, 2 * P, sqrt(2 * P * (1 - P)))
L = matrix(1, K, 1)

u = mvrnorm(M, c(0, 0, 0), var_g / M)
e = rnorm(K, 0, sqrt(var_e))
y = L %*% 2 + cbind(Zm, Zp, Zo) %*% c(u) + e

# Fix genetic data
# -----------------------------------------------------
# BIM # One line per SNP
bim = make_bim(n = M)
#bim$chr = paste0("chr", bim$chr)
bim$chr = sample(1:22, M, replace = T)
# Make SNP IDs look like "rs" IDs
bim$id = paste0('rs', bim$id)
bim$posg = bim$posg
bim$pos = bim$pos
# Select randomly between Cs and Gs for the reference alleles
bim$ref = sample(c("C", "G"), M, replace = T)
# Select randomly between As and Ts for the alternative alleles
bim$alt = sample(c("A", "T"), M, replace = T)

# FAM
# Specify the number of individuals
fam = make_fam(n = 3 * K)
fam$fam = rep(1:K, times = 3)
fam$id = fam$id
fam$pat = c(rep(0, 2 * K), fam$id[(K + 1):(2 * K)])
fam$mat = c(rep(0, 2 * K), fam$id[1:K])
# Sex 1 = male, 2 = female
fam$sex = c(rep(2:1, each = K), sample(1:2, K, T))
# Phenotype missing for parents
fam$pheno = c(rep(-9, 2 * K), y)

# BED
Xt = t(X)
dimnames(Xt) = list(bim$id, fam$id)

# Write to file
sim = "simdat"
write_bim(sim, bim)
write_fam(sim, fam)
write_bed(sim, Xt)

# pheno
#pheno = data.frame(fam$fam, fam$id, fam$pheno)
#write.table(pheno, paste0(sim, ".pheno"), quote = F, sep = "\t", row.names = F, col.names = F)
