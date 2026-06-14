## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----lib----------------------------------------------------------------------
library(nmfkc)

## ----data, fig.width = 6, fig.height = 4--------------------------------------
set.seed(123)
r_true <- 3
P <- 24                       # features
n_each <- 20                  # samples per cluster
N <- r_true * n_each          # 60 samples

# basis: each latent factor is a distinct block of features
X <- matrix(0, P, r_true)
X[1:8, 1] <- 1; X[9:16, 2] <- 1; X[17:24, 3] <- 1
X <- X + matrix(runif(P * r_true, 0, 0.1), P, r_true)   # small background

# coefficients: each sample loads mainly on its own cluster's factor
cl_true <- rep(1:r_true, each = n_each)
B <- matrix(runif(r_true * N, 0, 0.2), r_true, N)
for (j in seq_len(N)) B[cl_true[j], j] <- runif(1, 1, 2)

Y <- X %*% B
Y <- Y + matrix(rnorm(P * N, 0, 0.15 * mean(X %*% B)), P, N)  # noise
Y[Y < 0] <- 0

dim(Y)
image(1:N, 1:P, t(Y), xlab = "sample", ylab = "feature",
      main = "Synthetic data (true rank = 3)")

## ----rank, fig.width = 7, fig.height = 6--------------------------------------
rk <- nmfkc.rank(Y, rank = 1:6)
rk$rank.best                       # recommended rank (ECV minimum)
round(rk$criteria, 3)

## ----cv-----------------------------------------------------------------------
ev <- nmfkc.ecv (Y, rank = 1:6)
bv <- nmfkc.bicv(Y, rank = 1:6)   # nfolds = 2 (Owen & Perry) by default
data.frame(rank = 1:6,
           sigma.ecv  = round(ev$sigma, 3),
           sigma.bicv = round(bv$sigma, 3))
cat("ECV  picks rank", which.min(ev$sigma),
    "| bi-CV picks rank", bv$rank[which.min(bv$sigma)], "\n")

## ----consensus, fig.width = 7, fig.height = 5---------------------------------
cs <- nmfkc.consensus(Y, rank = 2:6, nrun = 20, keep.consensus = TRUE)
cs
plot(cs)                               # stability curves

## ----consensus-heatmap, fig.width = 8, fig.height = 6-------------------------
plot(cs, type = "heatmap")   # all ranks

## ----ard, fig.width = 6.5, fig.height = 4.5-----------------------------------
ar <- nmfkc.ard(Y, rank = 8, nrun = 20)   # >=10-20 restarts for a stable mode
ar                                        # see "rank over runs" for confidence
plot(ar)

## ----summary------------------------------------------------------------------
data.frame(
  method = c("nmfkc.rank (ECV)", "nmfkc.ecv", "nmfkc.bicv",
             "nmfkc.consensus (dispersion)", "nmfkc.consensus (PAC)",
             "nmfkc.ard"),
  estimate = c(rk$rank.best,
               which.min(ev$sigma),
               bv$rank[which.min(bv$sigma)],
               cs$rank[which.max(cs$dispersion)],
               cs$rank[which.min(cs$pac)],
               ar$rank),
  true = r_true
)

## ----flow, fig.width = 8, fig.height = 6--------------------------------------
fits <- lapply(1:6, function(q) nmfkc(Y, Q = q, print.dims = FALSE))
fl <- nmf.cluster.flow(fits, reference = 3)
head(fl$clusters)      # rows = individuals, columns = rank, entries = cluster id

