## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## -----------------------------------------------------------------------------
library(lavaan)
data(HolzingerSwineford1939)

d <- HolzingerSwineford1939
index <- complete.cases(d)
d <- d[index, ]

## -----------------------------------------------------------------------------
d$age.rev <- -(d$ageyr + d$agemo / 12)
d$sex.2 <- ifelse(d$sex == 2, 1, 0)
d$school.GW <- ifelse(d$school == "Grant-White", 1, 0)
d[, c("id", "sex", "ageyr", "agemo", "school", "grade")] <- NULL

## -----------------------------------------------------------------------------
library(nmfkc)
d <- nmfkc::nmfkc.normalize(d)

## -----------------------------------------------------------------------------
exogenous_vars <- c("age.rev", "sex.2", "school.GW")
endogenous_vars <- setdiff(colnames(d), exogenous_vars)

Y1 <- t(d[, endogenous_vars])
Y2 <- t(d[, exogenous_vars])

## -----------------------------------------------------------------------------
myepsilon <- 1e-6
Q0 <- 3

res0 <- nmfkc(
  Y = Y1,
  A = Y2,
  rank = Q0,
  epsilon = myepsilon,
  X.L2.ortho = 100
)

M.simple <- res0$X %*% res0$C

## ----eval=FALSE---------------------------------------------------------------
# grid_params <- expand.grid(
#     C1.L1     = c(0,1:9/10,1:10),
#     C2.L1     = c(0,1:9/10,1:10)
# )
# n_iter <- nrow(grid_params)
# mae.cv <- 0*1:n_iter
# 
# for(i in 1:n_iter){
#   if (i %% round(n_iter / 10) == 0) {
#     message(sprintf("Processing... %d%% (%d/%d)", round(i/n_iter*100), i, n_iter))
#   }
#   p <- grid_params[i, ]
#   res.cv <- nmf.ffb.cv(Y1, Y2, rank = Q0,
#                        X.init = res0$X,
#                        X.L2.ortho = 100,
#                        C1.L1 = p$C1.L1,
#                        C2.L1 = p$C2.L1,
#                        seed = 1, epsilon = myepsilon)
#   mae.cv[i] <- res.cv
# }
# 
# f <- data.frame(grid_params,mae.cv)
# f <- f[order(f$mae.cv),]
# head(f,5)
# #    C2.L2.off C1.L1 C2.L1    mae.cv
# #140         0    10   0.6 0.1820841
# #160         0    10   0.7 0.1820843
# #180         0    10   0.8 0.1820877
# #200         0    10   0.9 0.1820907
# #120         0    10   0.5 0.1820908
# print(p <- f[1,])

## -----------------------------------------------------------------------------
p <- list(C1.L1 = 10, C2.L1 = 0.6)

res <- nmf.ffb(
  Y1, Y2,
  rank = Q0,
  X.init = res0$X,
  X.L2.ortho = 100,
  C1.L1 = p$C1.L1,
  C2.L1 = p$C2.L1,
  epsilon = myepsilon
)

## -----------------------------------------------------------------------------
plot(res$objfunc.full, type = "l",
     main = "Objective Function",
     ylab = "Loss", xlab = "Iteration")

## -----------------------------------------------------------------------------
SC.map <- cor(as.vector(res$M.model),
              as.vector(M.simple))

cat("Q=     ", Q0, "\n")
cat("RHO=   ", round(res$XC1.radius, 3), "\n")
cat("AR=    ", round(res$amplification, 3), "\n")
cat("SCmap= ", round(SC.map, 3), "\n")
cat("SCcov= ", round(res$SC.cov, 3), "\n")
cat("MAE=   ", round(res$MAE, 3), "\n")

## -----------------------------------------------------------------------------
res.dot <- nmf.ffb.DOT(
  res,
  weight_scale = 5,
  rankdir = "TB",
  threshold = 0.01,
  fill = FALSE,
  cluster.box = "none"
)

# plot(res.dot)  # requires DiagrammeR

