## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load-data----------------------------------------------------------------
library(nmfkc)
library(nlme)

data(Orthodont)
head(Orthodont)

## ----prepare-Y----------------------------------------------------------------
Y <- matrix(Orthodont$distance, nrow = 4, ncol = 27)
colnames(Y) <- paste0("S", 1:27)
rownames(Y) <- paste("Age", c(8, 10, 12, 14))
Y[, 1:6]

## ----prepare-A----------------------------------------------------------------
male <- ifelse(Orthodont$Sex[seq(1, 108, 4)] == "Male", 1, 0)
A <- rbind(intercept = 1, male = male)
A[, 1:6]

## ----fit-nmfre----------------------------------------------------------------
res <- nmfre(Y, A, rank = 1, prefix = "Trend")

## ----summary------------------------------------------------------------------
summary(res)

## ----diagnostics--------------------------------------------------------------
cat("sigma^2 =", round(res$sigma2, 4), "\n")
cat("tau^2   =", round(res$tau2, 4), "\n")
cat("lambda  =", round(res$lambda, 5), "\n")
cat("dfU     =", round(res$dfU, 2),
    "  dfU/(NQ) =", round(res$dfU.frac, 4), "\n")

## ----inference----------------------------------------------------------------
res <- nmfre.inference(res, Y, A, wild.B = 1000)
summary(res)

## ----plot-growth, fig.height=6------------------------------------------------
age <- c(8, 10, 12, 14)

plot(age, res$XB[, 1], type = "n", ylim = range(Y),
     xlab = "Age (years)", ylab = "Distance (mm)",
     main = "Orthodont: NMF-RE Growth Curves")

# Plot observed data points
for (j in 1:27) {
  pch_j <- ifelse(male[j] == 1, 4, 1)
  points(age, Y[, j], pch = pch_j, col = "gray60")
}

# Plot individual predictions (fixed + random effects)
for (j in 1:27) {
  lines(age, res$XB.blup[, j], col = "steelblue", lty = 3, lwd = 0.8)
}

# Plot population-level fixed effects (two lines: male and female)
for (j in 1:27) {
  lines(age, res$XB[, j], col = "red", lwd = 2)
}

legend("topleft",
  legend = c("Fixed effect (male/female)", "Fixed + Random (BLUP)",
             "Male (observed)", "Female (observed)"),
  lwd = c(2, 1, NA, NA), lty = c(1, 3, NA, NA),
  pch = c(NA, NA, 4, 1),
  col = c("red", "steelblue", "gray60", "gray60"),
  cex = 0.85)

## ----basis-X------------------------------------------------------------------
cat("Basis X (temporal pattern):\n")
print(round(res$X, 4))

## ----coef-C-------------------------------------------------------------------
cat("Coefficient matrix (Theta):\n")
print(round(res$C, 4))

## ----random-effects, fig.height=4---------------------------------------------
barplot(res$U[1, ], names.arg = colnames(Y),
        las = 2, cex.names = 0.7,
        col = ifelse(male == 1, "steelblue", "salmon"),
        main = "Random Effects (U) by Subject",
        ylab = "Random effect value")
legend("topright", fill = c("steelblue", "salmon"),
       legend = c("Male", "Female"), cex = 0.85)

## ----convergence--------------------------------------------------------------
cat("Converged:", res$converged, "\n")
cat("Iterations:", res$iter, "\n")
cat("Stop reason:", res$stop.reason, "\n")

## ----plot-convergence, fig.height=4-------------------------------------------
plot(res, main = "Convergence (marginal NLL)")

## ----residual-analysis, fig.height=4------------------------------------------
residuals <- Y - res$XB.blup
cat("Mean absolute residual (BLUP):", round(mean(abs(residuals)), 4), "\n")
cat("Mean absolute residual (fixed):", round(mean(abs(Y - res$XB)), 4), "\n")

# Fitted vs Observed
plot(as.vector(Y), as.vector(res$XB.blup),
     xlab = "Observed", ylab = "Fitted (BLUP)",
     main = "Observed vs Fitted", pch = 16, col = "steelblue")
abline(0, 1, col = "red", lwd = 2)

## ----compare-models-----------------------------------------------------------
# Standard NMF with covariates (no random effects)
res_fixed <- nmfkc(Y, A = A, rank = 1)

cat("=== Standard NMF (fixed effects only) ===\n")
cat("R-squared:", round(1 - sum((Y - res_fixed$XB)^2) / sum((Y - mean(Y))^2), 4), "\n\n")

cat("=== NMF-RE (fixed + random effects) ===\n")
cat("R-squared (XB):      ", round(res$r.squared.fixed, 4), "\n")
cat("R-squared (XB+blup): ", round(res$r.squared, 4), "\n")
cat("ICC:                 ", round(res$ICC, 4), "\n")

## ----reinference--------------------------------------------------------------
# coefficient table from the inference object built in Section 2.2
res$coefficients[, c("Basis", "Covariate", "Estimate", "SE", "p_value")]

## ----dot-visualization, eval=requireNamespace("DiagrammeR", quietly=TRUE)-----
dot <- nmfkc.DOT(res, type = "YXA", sig.level = 0.05)
plot(dot)

## ----nonneg-------------------------------------------------------------------
res_nn <- nmfre(Y, A, rank = 1, prefix = "Trend", C.signed = FALSE)
res_nn <- nmfre.inference(res_nn, Y, A, wild.B = 500)
cat("All Theta >= 0 :", all(res_nn$C >= 0), "\n")
cat("P-value side   :", res_nn$C.p.side, "\n")
res_nn$coefficients[, c("Basis", "Covariate", "Estimate", "SE", "p_value")]

