## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----lib----------------------------------------------------------------------
library(nmfkc)

## ----data, fig.width = 6, fig.height = 5.2------------------------------------
set.seed(1)
K <- 3; n_each <- 20; N <- K * n_each       # 60 nodes, 3 communities
block <- rep(1:K, each = n_each)            # true community labels

p_in <- 0.6; p_out <- 0.05
Prob <- matrix(p_out, N, N)
for (k in 1:K) Prob[block == k, block == k] <- p_in
Y <- matrix(rbinom(N * N, 1, Prob), N, N)   # 0/1 adjacency
Y[lower.tri(Y)] <- t(Y)[lower.tri(Y)]        # make symmetric
diag(Y) <- 0

# add three bridge nodes, each also linked to a second community
bridge <- c(20, 40, 60); into <- c(2, 3, 1)
for (i in seq_along(bridge)) {
  b <- bridge[i]; tgt <- which(block == into[i])
  e <- rbinom(length(tgt), 1, 0.45)
  Y[b, tgt] <- pmax(Y[b, tgt], e); Y[tgt, b] <- Y[b, tgt]
}

isSymmetric(Y); sum(Y) / 2                   # symmetric, number of edges

# adjacency matrix in the original (random) node order -- structure is hidden
image(1:N, 1:N, Y[, N:1], col = c("white", "steelblue"),
      xlab = "node", ylab = "node", main = "Adjacency (original order)")

## ----fit----------------------------------------------------------------------
res <- nmfkc.net(Y, rank = 3, type = "bi", nstart = 10, maxit = 500)
res$r.squared                 # fit
head(round(res$X.prob, 2))    # soft membership: rows = nodes, columns = communities
head(res$X.cluster)           # hard label = argmax membership

## ----compare------------------------------------------------------------------
table(true = block, estimated = res$X.cluster)

# the three bridge nodes have genuinely split memberships
round(res$X.prob[bridge, ] * 100, 1)

## ----viz, fig.width = 8, fig.height = 4---------------------------------------
ord <- order(res$X.cluster)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
image(1:N, 1:N, Y[ord, rev(ord)], col = c("white", "steelblue"),
      xlab = "node (reordered)", ylab = "", main = "Adjacency by community")
image(1:N, 1:K, res$X.prob[ord, , drop = FALSE],
      col = hcl.colors(20, "Blues", rev = TRUE),
      xlab = "node (reordered)", ylab = "community", axes = FALSE,
      main = "Soft membership")
axis(1); axis(2, at = 1:K)
par(op)

## ----dot-check, include = FALSE-----------------------------------------------
has_dg <- requireNamespace("DiagrammeR", quietly = TRUE)

## ----dot-soft, eval = has_dg, fig.width = 7, fig.height = 7-------------------
dot_soft <- nmfkc.net.DOT(res,
  Y.label = as.character(1:N), X.label = paste("Community", 1:K),
  Y.title = "Network nodes",   X.title = "Communities",
  layout  = "neato", threshold = 0.25, Y.cluster = "soft")
DiagrammeR::grViz(as.character(dot_soft))

## ----dot-hard, eval = has_dg, fig.width = 7, fig.height = 7-------------------
dot_hard <- nmfkc.net.DOT(res,
  Y.label = as.character(1:N), X.label = paste("Community", 1:K),
  Y.title = "Network nodes",   X.title = "Communities",
  layout  = "neato", threshold = 0.25, Y.cluster = "hard")
DiagrammeR::grViz(as.character(dot_hard))

## ----rank, fig.width = 7, fig.height = 6--------------------------------------
rk <- nmfkc.net.rank(Y, rank = 1:6, type = "bi", nstart = 5)
rk$rank.best
round(rk$criteria, 3)

