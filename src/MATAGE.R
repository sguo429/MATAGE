library(data.table)
library(parallel)
library(Matrix)
library(GMMAT)
library(CompQuadForm)

# --- Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)
geno_infile <- args[1]
start_row <- as.numeric(args[2])
num_rows <- as.numeric(args[3])
outfile <- args[4]
null_model_RData <- args[5]
geno_start_col <- as.numeric(args[6])

# Load the null model
load(null_model_RData)

# --- Main Association Function ---
MATAGE <- function(obj, infile, outfile, start_row, num_rows, geno_start_col){
  
  # Part 1: IO
  if (grepl("\\.txt$", infile)) {
    cmd <- sprintf("head -n 1 %s && tail -n +%d %s | head -n %d", infile, start_row + 1, infile, num_rows)
  } else if (grepl("\\.txt\\.gz$", infile)) {
    cmd <- sprintf("gunzip -c %s | { head -n 1; tail -n +%d | head -n %d; }", infile, start_row + 1, num_rows)
  } else {
    stop("The genotype file format is not supported.")
  }

  geno_data_full <- fread(cmd = cmd, header = TRUE, data.table = FALSE)
  
  # Part 2: Data Preparation
  # Match sample IDs between genotype file and null model
  geno_sample_ids <- colnames(geno_data_full)[geno_start_col:ncol(geno_data_full)]
  obj_ids <- unique(obj$id_include)
  common_ids <- geno_sample_ids[geno_sample_ids %in% obj_ids]
  if (any(!obj_ids %in% common_ids)) {
    warning("Some id_include in obj are missing in the sample IDs of the infile!")
  }
  
  # Subset the null model objects to the common samples
  idx <- match(common_ids, obj_ids)
  scaled.residuals <- obj$scaled.residuals[idx]
  Sigma_i <- obj$Sigma_i[idx, idx]
  Sigma_iX <- obj$Sigma_iX[idx, , drop = FALSE]
  
  # Separate variant info from the genotype data
  info_data <- geno_data_full[, 1:(geno_start_col - 1), drop = FALSE]
  
  # Convert all genotype columns to a numeric matrix
  #geno_matrix <- as.matrix(geno_data_full[, c(common_ids)])
  geno_matrix <- do.call(cbind, lapply(geno_data_full[, c(common_ids)], as.numeric))

  # Clean up large intermediate objects
  rm(geno_data_full, geno_sample_ids, obj_ids, common_ids, idx)

  # Check if the model is Linear or Logistic
  is_binary_model <- family$family == "binomial"
  # Check if it is a linear regression model
  is_linear_regression <- !is_binary_model & (length(obj$theta) == 1)

  # Function to compute p-values for a single genotype
  compute_p_values <- function(G, info_row) {
    
    # --- Part 3: Linear Test ---
    # Handle missing values
    if(any(is.na(G))){
        idx <- !is.na(G)
        G <- G[idx]
        scaled.residuals <- scaled.residuals[idx]
        Sigma_i <- Sigma_i[idx, idx]
        Sigma_iX <- Sigma_iX[idx, , drop = FALSE]
    }

    # Summary data of genotype
    mean_G <- mean(G)
    var_G <- var(G)
    n <- length(G)
    G0 <- sort(unique(G))
    r <- length(G0)
    
    if (r == 1){
        return(list(ID = info_row$ID, n = n, knots = r, mean_G = mean_G, var_G = var_G, p_value_l = NA, p_value_nl = NA))
    }

    H_i <-  Sigma_i * obj$theta[1]
    S_G <- as.numeric(crossprod(G, scaled.residuals))
    Gt_res <- S_G * obj$theta[1]
    p <- length(obj$coefficients)
    Gt_Hi_G <- crossprod(G, crossprod(H_i, G))
    Gt_Hi_X <- crossprod(G, Sigma_iX * obj$theta[1])
    Xt_Hi_X_i <- obj$cov / obj$theta[1]
    Gt_Px_G_i <- 1 / as.numeric((Gt_Hi_G - tcrossprod(Gt_Hi_X, tcrossprod(Gt_Hi_X, Xt_Hi_X_i))))
    beta_hat <- Gt_res * Gt_Px_G_i
    dispersion <- ifelse(is_binary_model, 1, ((n - p) * obj$theta[1] - Gt_res^2 * Gt_Px_G_i) / (n - p - 1))
    Var_beta_hat <- dispersion * Gt_Px_G_i
    W <- (beta_hat)^2 / Var_beta_hat
    p_value_l <- ifelse(is_linear_regression, pf(W, 1, n - p - 1, lower.tail = FALSE), pchisq(W, df = 1, lower.tail = FALSE))

    # --- Part 4: Nonlinear Test ---
    if (r == 2) {
        return(list(ID = info_row$ID, n = n, knots = r, mean_G = mean_G, var_G = var_G, p_value_l = p_value_l, p_value_nl = NA))
    }

    if(!is_binary_model){
        # Updating using the new residual variance
        Sigma_i <-  H_i / dispersion
        Sigma_iX <- Sigma_iX * obj$theta[1] / dispersion
        S_G <- Gt_res / dispersion
        scaled.residuals <- scaled.residuals * obj$theta[1] / dispersion
        cov <- Xt_Hi_X_i * dispersion
    } else {
        cov <- obj$cov
    }

    # Generate matrix J
    J <- sparseMatrix(i=1:n, j=match(G, G0), x=1)

    # Prepare to construct matrix Q and R
    h <- diff(G0)
    h1 <- h[c(1:(r-2))]
    h2 <- h[c(2:(r-1))]
    h1_inv <- 1 / h1
    h2_inv <- 1 / h2
    cols <- r - 2

    # Prepare to construct matrix Q
    q1 <- sparseMatrix(
        i = rep(1:cols, each = 2),
        j = as.vector(sapply(1:cols, function(k) k:(k+1))),
        x = rep(c(1, -1), times = cols),
        dims = c(cols, r)
    )

    q2 <- sparseMatrix(
        i = rep(1:cols, each = 2),
        j = as.vector(sapply(1:cols, function(k) c(k+1, k+2))),
        x = rep(c(-1, 1), times = cols),
        dims = c(cols, r)
    )

    # Construct matrix Q
    Q_trans <- q1*h1_inv + q2*h2_inv
    Q <- t(Q_trans)

    # Prepare to construct matrix R
    diag_vals <- rep(1/3, cols)
    off_diag_vals <- rep(1/6, cols - 1)
     
    if (cols > 1){
        r1 <- sparseMatrix(
        i = c(1:cols, 2:cols),
        j = c(1:cols, 1:(cols-1)),
        x = c(diag_vals, off_diag_vals),
        dims = c(cols, cols)
        )

        r2 <- sparseMatrix(
            i = c(1:cols, 1:(cols-1)),
            j = c(1:cols, 2:cols),
            x = c(diag_vals, off_diag_vals),
            dims = c(cols, cols)
        )
    } else {
        r1 <- 1/3
        r2 <- 1/3
    }

    # Construct matrix R
    R <- r1*h1 + r2*h2

    # Cholesky decomposition of R
    C <- chol(R)

    # Generate matrix B
    B0 <- crossprod(Q_trans, tcrossprod(solve(crossprod(Q)), C))
    S_B <- crossprod(B0, crossprod(J, scaled.residuals))

    # Score test for the nonlinear effect
    Jt_Sigma_iX <- crossprod(J, Sigma_iX)
    Jt_Sigma_iJ <- crossprod(J, crossprod(Sigma_i, J))
    JtPJ <- Jt_Sigma_iJ - tcrossprod(Jt_Sigma_iX, tcrossprod(Jt_Sigma_iX, cov))
    JtPJ_B0 <- crossprod(JtPJ, B0)
    GtPB <- crossprod(G0, JtPJ_B0)
    GtPG_i <- solve(crossprod(G0, crossprod(JtPJ, G0)))
    M <- crossprod(GtPG_i, GtPB)
    S_B_hat <- S_B - crossprod(M, S_G)
    V <- crossprod(B0, JtPJ_B0) - crossprod(M, GtPB)
    p_value_nl <- .quad_pval(U = S_B_hat, V = V, method = "davies")
    
    return(list(ID = info_row$ID, n = n, knots = r, mean_G = mean_G, var_G = var_G, p_value_l = p_value_l, p_value_nl = p_value_nl))
  }

  # Loop over the rows of the genotype matrix
  results_list <- mclapply(1:nrow(geno_matrix), function(i) {
    compute_p_values(G = geno_matrix[i, ], info_row = info_data[i, ])
  }, mc.cores = 1)
  
  # Combine and save results
  results <- rbindlist(results_list, use.names=TRUE, fill=TRUE)
  results$p_value_joint <- apply(results[, c("p_value_l", "p_value_nl")], 1, fisher_pval)
  fwrite(results, outfile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="NA")
}

# --- Helper Functions ---
.quad_pval <- function(U, V, method = "davies") {
    Q <- sum(U^2)
    lambda <- eigen(V, only.values = TRUE, symmetric=TRUE)$values
    lambda <- lambda[lambda > 0]
    pval <- .Q_pval(Q, lambda, method = method)
    return(pval)
}


.Q_pval <- function(Q, lambda, method = "davies") {
  if(method == "davies") {
    tmp <- try(suppressWarnings(CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-6)))
    if(inherits(tmp, "try-error") || tmp$ifault > 0 || tmp$Qq <= 1e-5 || tmp$Qq >= 1) method <- "kuonen"
    else return(tmp$Qq)
  }
  if(method == "kuonen") {
    pval <- try(.pKuonen(x = Q, lambda = lambda))
    if(inherits(pval, "try-error") || is.na(pval)) method <- "liu"
    else return(pval)
  }
  if(method == "liu") {
    pval <- try(CompQuadForm::liu(q = Q, lambda = lambda))
    if(inherits(pval, "try-error")) cat("Warning: method \"liu\" failed...\nQ:", Q, "\nlambda:", lambda, "\n")
    else return(pval)
  }
  return(NA)
}


.pKuonen<-function (x, lambda, delta = rep(0, length(lambda)), df = rep(1, length(lambda)))
{
    delta <- delta[lambda != 0]
    df <- df[lambda != 0]
    lambda <- lambda[lambda != 0]
    if(length(lambda) != length(delta)) stop("Error: inconsistent length in lambda and delta!")
    if(length(lambda) != length(df)) stop("Error: inconsistent length in lambda and df!")
    if (length(lambda) == 1) {
        return(pchisq(x/lambda, df = df, ncp = delta, lower.tail = FALSE))
    }
    d <- max(lambda)
    lambda <- lambda/d
    x <- x/d
    k0 <- function(zeta) {
        -sum(df * log(1 - 2 * zeta * lambda))/2 + sum((delta * lambda *
            zeta)/(1 - 2 * zeta * lambda))
    }
    kprime0 <- function(zeta) {
        sapply(zeta, function(zz) {
            sum(((delta + df) * lambda)/(1 - 2 * zz * lambda) + 2 * (delta *
                zz * lambda^2)/(1 - 2 * zz * lambda)^2)
        })
    }
    kpprime0 <- function(zeta) {
        sum((2 * (2 * delta + df) * lambda^2)/(1 - 2 * zeta * lambda)^2 + 8 *
            delta * zeta * lambda^3/(1 - 2 * zeta * lambda)^3)
    }
    if (any(lambda < 0)) {
        lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
    }
    else if (x > sum((df+delta)*lambda)) {
        lmin <- -0.01
    }
    else {
        lmin <- -length(lambda)*max(df+delta)/(2 * x)
    }
    lmax <- min(1/(2 * lambda[lambda > 0])) * 0.99999
    hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin,
        upper = lmax, tol = 1e-08)$root
    w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
    v <- hatzeta * sqrt(kpprime0(hatzeta))
    if (abs(hatzeta) < 1e-04)
        N
    else pnorm(w + log(v/w)/w, lower.tail = FALSE)
}

fisher_pval <- function(p) {
  is.valid.p <- !is.na(p) & p > 0 & p <= 1
  if(sum(is.valid.p) == 0) return(NA)
  p <- p[is.valid.p]
  pchisq(-2*sum(log(p)), df = 2*length(p), lower.tail = FALSE)
}

# --- Run the Analysis ---
MATAGE(null_model, geno_infile, outfile, start_row, num_rows, geno_start_col)
