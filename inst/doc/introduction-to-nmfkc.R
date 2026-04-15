## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load-package-------------------------------------------------------------
library(nmfkc)

## ----create-data--------------------------------------------------------------
# Rows: Users (U1-U5), Cols: Movies (M1-M4)
# U1, U2, U3 prefer Action movies.
# U4, U5 prefer Romance movies.
Y <- matrix(
  c(5, 4, 1, 1,
    4, 5, 1, 2,
    5, 5, 2, 2,
    1, 2, 5, 4,
    1, 1, 4, 5),
  nrow = 5, byrow = TRUE
)

# Assign names for better interpretation
rownames(Y) <- paste0("User", 1:5)
colnames(Y) <- c("Action1", "Action2", "Romance1", "Romance2")

# Check the data
print(Y)

## ----run-nmfkc----------------------------------------------------------------
# Run NMF with rank = 2
res <- nmfkc(Y, rank = 2, seed = 123)

## ----interpret-X--------------------------------------------------------------
# Each column represents a latent factor (Basis)
res$X

## ----interpret-B--------------------------------------------------------------
# Each row represents a latent factor
res$B

## ----plot-convergence---------------------------------------------------------
plot(res, main = "Convergence Plot")

## ----plot-residual, fig.width=9, fig.height=4---------------------------------
# Visualize Original vs Fitted vs Residuals
nmfkc.residual.plot(Y, res)

## ----create-na----------------------------------------------------------------
Y_missing <- Y
Y_missing["User1", "Action1"] <- NA # Introduce missing value
print(Y_missing)

## ----run-na-------------------------------------------------------------------
res_na <- nmfkc(Y_missing, rank = 2, seed = 123)

## ----impute-na----------------------------------------------------------------
# Extract the predicted value from the fitted matrix XB
predicted_rating <- res_na$XB["User1", "Action1"]
actual_rating <- Y["User1", "Action1"] # The original hidden value (5)

cat(paste0("Actual Rating:    ", actual_rating, "\n"))
cat(paste0("Predicted Rating: ", round(predicted_rating, 2), "\n"))

