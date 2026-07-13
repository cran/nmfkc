## ----setup, include = FALSE---------------------------------------------------
# Data-dependent chunks run only when ade4 (source of the Doubs data) is
# installed, so the vignette still builds on machines without it.
has_ade4 <- requireNamespace("ade4", quietly = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = has_ade4
)

## ----lib----------------------------------------------------------------------
library(nmfkc)

## ----no-ade4, eval = !has_ade4, echo = FALSE, results = "asis"----------------
# cat("> **Note:** this vignette needs the `ade4` package for the Doubs data.",
#     "Install it with `install.packages(\"ade4\")` to run the code below.\n")

## ----data---------------------------------------------------------------------
data(doubs, package = "ade4")

# per-variable min-max to [0,1], then transpose to (variables x sites)
nz <- function(M) t(nmfkc.normalize(as.matrix(M)))
Y1 <- nz(doubs$fish)   # responses: 27 fish species x 30 sites
Y2 <- nz(doubs$env)    # covariates: 11 environment x 30 sites
dim(Y1)
dim(Y2)

## ----fit----------------------------------------------------------------------
fit <- nmf.rrr(Y1, Y2, rank1 = 2, rank2 = 2,
               epsilon = 1e-8, nstart = 20, seed = 1)

# in-sample, column-centered R^2
Y1hat <- fit$X1 %*% fit$C %*% fit$X2 %*% Y2
R2 <- 1 - sum((Y1 - Y1hat)^2) / sum((Y1 - rowMeans(Y1))^2)
round(R2, 3)

## ----resp-groups--------------------------------------------------------------
for (q in 1:ncol(fit$X1))
  cat(sprintf("Resp%d: %s\n", q,
      paste(rownames(Y1)[order(-fit$X1[, q])[1:6]], collapse = ", ")))

## ----cov-groups---------------------------------------------------------------
for (r in 1:nrow(fit$X2))
  cat(sprintf("Cov%d: %s\n", r,
      paste(rownames(Y2)[order(-fit$X2[r, ])[1:5]], collapse = ", ")))

## ----ecv----------------------------------------------------------------------
ecv <- nmf.rrr.ecv(Y1, Y2, rank1 = 1:2, rank2 = 1:2,
                   nfolds = 5, seed = 123)
round(ecv$sigma, 4)

## ----inference----------------------------------------------------------------
inf <- nmf.rrr.inference(fit, Y1, Y2)
co  <- inf$coefficients
print(format(co[order(co$p_value), c("Basis","Covariate","Estimate","SE","z_value","p_value")],
             digits = 3))

## ----theta--------------------------------------------------------------------
round(fit$C, 3)

## ----heatmap, fig.width = 7.5, fig.height = 4---------------------------------
nmf.rrr.heatmap(fit)

