## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7, 
  fig.height = 5
)

## ----load-packages------------------------------------------------------------
library(nmfkc)
library(quanteda)

## ----data-prep----------------------------------------------------------------
# Load the corpus from quanteda
corp <- corpus(data_corpus_inaugural)

# Preprocessing: tokenize, remove stopwords, and punctuation
tok <- tokens(corp, remove_punct = TRUE)
tok <- tokens_remove(tok, pattern = stopwords("en", source = "snowball"))

# Create DFM and filter
df <- dfm(tok)
df <- dfm_select(df, min_nchar = 3) # Remove short words (<= 2 chars)
df <- dfm_trim(df, min_termfreq = 100) # Remove rare words (appearing < 100 times)

# --- CRITICAL STEP ---
# quanteda's DFM is (Documents x Words).
# nmfkc expects (Words x Documents).
# We must transpose the matrix.
d <- as.matrix(df)
# Sort by frequency
index <- order(colSums(d), decreasing=TRUE) 
d <- d[,index]
Y <- t(d)

dim(Y) # Features (Words) x Samples (Documents)

## ----rank-selection, fig.width=7, fig.height=6--------------------------------
# Evaluate ranks from 2 to 6
nmfkc.rank(Y, rank = 2:6)

## ----model-fit----------------------------------------------------------------
rank <- 3
# Set seed for reproducibility
res_std <- nmfkc(Y, rank = rank, seed = 123, prefix = "Topic")

# Check Goodness of Fit (R-squared)
res_std$r.squared

## ----interpret-topics---------------------------------------------------------
# Extract top 10 words for each topic from X.prob (normalized X)
Xp <- res_std$X.prob
for(q in 1:rank){
  message(paste0("----- Featured words on Topic [", q, "] -----"))
  print(paste0(rownames(Xp), "(", rowSums(Y), ") ", round(100*Xp[,q], 1), "%")[Xp[,q]>=0.5])
}

## ----kernel-prep--------------------------------------------------------------
# Covariate: Year of the address
years <- as.numeric(substring(names(data_corpus_inaugural), 1, 4))
U <- t(as.matrix(years))

# Optimize beta (Gaussian Kernel width)
# We test a specific range of betas to locate the minimum CV error.
beta_candidates <- c(0.2, 0.5, 1, 2, 5) / 10000

# Run CV to find the best beta
# Note: We use the same rank (rank=3) as selected above.
cv_res <- nmfkc.kernel.beta.cv(Y, rank = rank, U = U, beta = beta_candidates, plot = FALSE)
best_beta <- cv_res$beta
print(best_beta)

## ----kernel-fit---------------------------------------------------------------
# Create Kernel Matrix
A <- nmfkc.kernel(U, beta = best_beta)

# Fit NMF with Kernel Covariates
res_ker <- nmfkc(Y, A = A, rank = rank, seed = 123, prefix = "Topic")

## ----plot-comparison, fig.width=8, fig.height=8-------------------------------
oldpar <- par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# Prepare Axis Labels (Rounded to integers)
at_points <- seq(1, ncol(Y), length.out = 10)
labels_years <- round(seq(min(years), max(years), length.out = 10))

# 1. Standard NMF (Noisy)
barplot(res_std$B.prob, col = 2:(rank+1), border = NA, xaxt='n',
        main = "Standard NMF: Topic Proportions (Noisy)", ylab = "Probability")
axis(1, at = at_points, labels = labels_years)

# 2. Kernel NMF (Smooth trend)
barplot(res_ker$B.prob, col = 2:(rank+1), border = NA, xaxt='n',
        main = "Kernel NMF: Temporal Topic Evolution (Smooth)", ylab = "Probability")
axis(1, at = at_points, labels = labels_years)

# Legend
legend("topright", legend = paste("Topic", 1:rank), fill = 2:(rank+1), bg="white", cex=0.8)
par(oldpar)

