library(data.table)
library(parallel)
library(Matrix)
library(GMMAT)
library(CompQuadForm)

args <- commandArgs(trailingOnly = TRUE)
geno_infile <- args[1]
start_row <- as.numeric(args[2])
num_rows <- as.numeric(args[3])
outfile <- args[4]
null_model_RData <- args[5]
geno_start_col <- as.numeric(args[6])

# Load the null model
load(null_model_RData)

# Function for conducting the association tests
MATAGE <- function(obj, infile, outfile, start_row, num_rows, geno_start_col){
    cmd <- sprintf("gunzip -c %s | { head -n 1; tail -n +%d | head -n %d; }", infile, start_row + 1, num_rows)
    geno_data <- fread(cmd = cmd, header = TRUE, data.table = FALSE)
    
    geno_sample_ids <- colnames(geno_data)[geno_start_col:ncol(geno_data)]
    obj_ids <- unique(obj$id_include)
    common_ids <- intersect(geno_sample_ids, obj_ids)
    if (any(!obj_ids %in% common_ids)) {
        warning("Some id_include in obj are missing in the sample IDs of the infile!")
    }

    info_cols <- colnames(geno_data)[1:(geno_start_col - 1)]
    geno_data <- geno_data[, c(info_cols, common_ids)]

    idx <- match(common_ids, obj_ids)
    scaled.residuals <- obj$scaled.residuals[idx]
    X <- obj$X[idx, , drop = FALSE]
    Y <- obj$Y[idx]
    Sigma_i <- obj$Sigma_i[idx, idx]
    Sigma_iX <- obj$Sigma_iX[idx, , drop = FALSE]

    rm(common_ids, obj_ids, geno_sample_ids, idx, info_cols)

    # Function to compute p-values for a genotype
    compute_p_values <- function(row) {
        G <- as.vector(as.numeric(row[, geno_start_col:ncol(row)]))

        # Missing values
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
        
        # Return NA p_values if there is only 1 knot
        if (r == 1){
            return(list(ID = row$ID, n = n, knots = r, mean_G = mean_G, var_G = var_G, p_value_l = NA, p_value_nl = NA))
        }

        # Wald test for the linear effect
        H_i <-  Sigma_i * obj$theta[1]
        S_G <- as.numeric(crossprod(G, scaled.residuals))
        Gt_res <- S_G * obj$theta[1]
        p <- length(obj$coefficients)
        Gt_Hi_G <- crossprod(G, crossprod(H_i, G))
        Gt_Hi_X <- crossprod(G, Sigma_iX * obj$theta[1])
        Xt_Hi_X_i <- obj$cov / obj$theta[1]
        Gt_Px_G_i <- 1 / as.numeric((Gt_Hi_G - tcrossprod(Gt_Hi_X, tcrossprod(Gt_Hi_X, Xt_Hi_X_i))))
        beta_hat <- Gt_res * Gt_Px_G_i
        res_var <- ((n - p) * obj$theta[1] - Gt_res^2 * Gt_Px_G_i) / (n - p - 1)
        Var_beta_hat <- res_var * Gt_Px_G_i
        W <- (beta_hat)^2 / Var_beta_hat
        p_value_l <- pf(W, 1, n - p - 1, lower.tail = FALSE)


        # Return NA non-linear p_values if there are only 2 knots
        if (r == 2) {
            return(list(ID = row$ID, n = n, knots = r, mean_G = mean_G, var_G = var_G, p_value_l = p_value_l, p_value_nl = NA))
        }

        # Update using the new residual variance
        Sigma_i <-  H_i / res_var
        Sigma_iX <- Sigma_iX * obj$theta[1] / res_var
        S_G <- Gt_res / res_var
        scaled.residuals <- scaled.residuals * obj$theta[1] / res_var
        cov <- Xt_Hi_X_i * res_var

        # Generate J matrix
        J <- sparseMatrix(i=1:n, j=match(G, G0), x=1)

        # Prepare to generate matrix Q and R
        h <- diff(G0)
        h1 <- h[c(1:(r-2))]
        h2 <- h[c(2:(r-1))]
        h1_inv <- 1 / h1
        h2_inv <- 1 / h2
        cols <- r - 2

        # Prepare to generate matrix Q
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

        # Generate matrix Q
        Q_trans <- q1*h1_inv + q2*h2_inv
        Q <- t(Q_trans)

        # Prepare to generate matrix R
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

        # Generate matrix R
        R <- r1*h1 + r2*h2

        # Cholesky decomposition of R
        C <- chol(R)

        # Generate matrix B
        B0 <- crossprod(Q_trans, tcrossprod(solve(crossprod(Q)), C))
        B <- J %*% B0
        S_B <- crossprod(B, scaled.residuals)

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

        return(list(ID = row$ID, n = n, knots = r, mean_G = mean_G, var_G = var_G, p_value_l = p_value_l, p_value_nl = p_value_nl))
    }

    # Compute test results for all the genotypes
    results_list <- mclapply(1:nrow(geno_data), function(i) compute_p_values(geno_data[i, ]), mc.cores = 1)

    # Combine results into a data frame 
    results <- rbindlist(results_list, use.names=TRUE, fill=TRUE)

    # Calculate the joint p-value
    results$p_value_joint <- apply(results[, c("p_value_l", "p_value_nl")], 1, fisher_pval)

    # Save the results
    fwrite(results, outfile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="NA")
}


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

# Conduct the association tests for selected genotypes
MATAGE(null_model, geno_infile, outfile, start_row, num_rows, geno_start_col)
