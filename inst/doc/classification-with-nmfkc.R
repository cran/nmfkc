## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load-packages------------------------------------------------------------
library(nmfkc)
library(palmerpenguins)

## ----iris-data-prep-----------------------------------------------------------
# 1. Prepare Labels (Y)
label_iris <- iris$Species
Y_iris <- nmfkc.class(label_iris) # One-hot encoding (3 Classes x 150 Samples)
rank_iris <- length(unique(label_iris)) # Number of classes

# 2. Prepare Features (U)
# Normalize features to [0, 1] range and transpose to (Features x Samples)
U_iris <- t(nmfkc.normalize(iris[, -5]))

## ----iris-linear--------------------------------------------------------------
# Use Normalized Features directly as A
A_linear <- U_iris

# Fit Linear NMF-LAB
res_linear <- nmfkc(Y = Y_iris, A = A_linear, rank = rank_iris, seed = 123, prefix = "Class")

# --- Interpretability Check ---
# The matrix C (Q x R) shows the weight of each Feature (columns) for each Class (rows).
# Let's look at the estimated weights:
round(res_linear$C, 2)

## ----iris-linear-acc----------------------------------------------------------
pred_linear <- predict(res_linear, type = "class")
(f_linear <- table(fitted.label = pred_linear, label = label_iris))

# Calculate Accuracy (assuming diagonal correspondence)
acc_linear <- sum(diag(f_linear)) / sum(f_linear)
cat(paste0("Linear Model Accuracy: ", round(acc_linear * 100, 2), "%\n"))

## ----iris-kernel--------------------------------------------------------------
# 1. Optimize Kernel Width (beta)
# Heuristic estimation of beta
res_beta <- nmfkc.kernel.beta.nearest.med(U_iris)

# Cross-validation for fine-tuning (using generated candidates)
cv_res <- nmfkc.kernel.beta.cv(Y_iris, rank = rank_iris, U = U_iris, 
                               beta = res_beta$beta_candidates, plot = FALSE)
best_beta <- cv_res$beta

# 2. Fit Kernel NMF-LAB
A_kernel <- nmfkc.kernel(U_iris, beta = best_beta)
res_kernel <- nmfkc(Y = Y_iris, A = A_kernel, rank = rank_iris, seed = 123, prefix = "Class")

# 3. Prediction and Evaluation
fitted_label <- predict(res_kernel, type = "class")
(f_kernel <- table(fitted.label = fitted_label, label = label_iris))

# Calculate Accuracy
acc_kernel <- sum(diag(f_kernel)) / sum(f_kernel)
cat(paste0("Kernel Model Accuracy: ", round(acc_kernel * 100, 2), "%\n"))

## ----iris-plot-X, fig.height=5------------------------------------------------
image(t(res_kernel$X)[, nrow(res_kernel$X):1], 
      main = "Basis Matrix X (Kernel Model)\nMapping Factors to Species",
      axes = FALSE, col = hcl.colors(12, "YlOrRd", rev = TRUE))

axis(1, at = seq(0, 1, length.out = rank_iris), labels = colnames(res_kernel$X))
axis(2, at = seq(0, 1, length.out = rank_iris), labels = rev(rownames(res_kernel$X)), las = 2)
box()

## ----penguins-data-prep-------------------------------------------------------
# Load and clean data (remove rows with NAs)
d_penguins <- na.omit(palmerpenguins::penguins)

# Prepare Y (Labels)
label_penguins <- d_penguins$species
Y_penguins <- nmfkc.class(label_penguins)

# Prepare U (Features)
U_penguins <- t(nmfkc.normalize(d_penguins[, 3:6]))

## ----penguins-model-fit-------------------------------------------------------
rank_penguins <- length(unique(label_penguins))

# 1. Heuristic beta estimation
best_beta_penguins <- nmfkc.kernel.beta.nearest.med(U_penguins)$beta

# 2. Optimization
A_penguins <- nmfkc.kernel(U_penguins, beta = best_beta_penguins)
res_penguins <- nmfkc(Y = Y_penguins, A = A_penguins, rank = rank_penguins, 
                      seed = 123, prefix = "Class")

## ----penguins-soft-vis, fig.width=8, fig.height=4-----------------------------
# Get probabilistic predictions
probs <- predict(res_penguins, type = "prob")

# Visualize
barplot(probs, col = c("#FF8C00", "#9932CC", "#008B8B"), border = NA,
        main = "Soft Classification Probabilities (Penguins)",
        xlab = "Sample Index", ylab = "Probability")

legend("topright", legend = levels(label_penguins), 
       fill = c("#FF8C00", "#9932CC", "#008B8B"), bg = "white", cex = 0.8)

## ----penguins-acc-------------------------------------------------------------
fitted_label_p <- predict(res_penguins, type = "class")
(f_penguins <- table(Predicted = fitted_label_p, Actual = label_penguins))

acc_p <- sum(diag(f_penguins)) / sum(f_penguins)
cat(paste0("Penguins Accuracy: ", round(acc_p * 100, 2), "%\n"))

