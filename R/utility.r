############################################################
# Utility functions for simulation studies
#
# This file provides auxiliary functions used to reproduce the
# simulation studies in the paper. These functions support data
# generation, censoring mechanism construction, and summary of empirical size and power.
############################################################

#' Set the covariate dimension
#'
#' Computes the covariate dimension used in the simulation studies as a function
#' of the sample size.
#'
#' @param n An integer specifying the sample size.
#'
#' @return An integer giving the covariate dimension.
#'
#' @details
#' The dimension is calculated as
#' \deqn{p = \mathrm{round}\{\exp((2n/3)^{0.4}) + 100\}.}
#'
#' @examples
#' set_p(90)
#' set_p(120)
#' set_p(150)
#'
#' @export
set_p <- function(n) {
    p <- round(exp((2/3*n)^0.4) + 100, 0)
    return(p)
}

#' Reformat simulation settings for summary tables and plots
#'
#' Reformats a simulation setting list for displaying summarized simulation
#' results in tables and figures.
#'
#' @param setting A list containing simulation settings, typically including
#'   `n`, `error`, and other simulation factors.
#'
#' @return A reformatted setting list. The sample size `n` is replaced by
#'   `np`, which records the pair `(n, p)`, and selected error distribution
#'   names are relabeled for display.
#'
#' @details
#' The covariate dimension `p` is computed from `n` by `set_p()`. The element
#' `n` is then replaced by `np = "(n,p)"`. If the setting list contains an
#' `error` element, `"lnorm"` is relabeled as `"log-norm"` and `"mixnorm"` is
#' relabeled as `"mix-norm"`.
#'
#' @examples
#' setting <- list(
#'   SNR = seq(0, 1, 0.5),
#'   n = c(90, 120),
#'   error = c("norm", "lnorm", "mixnorm")
#' )
#' reformat_setting(setting)
#'
#' @export
reformat_setting <- function(setting) {
    # Calculate p according to the given relationship
    n <- setting$n
    p <- set_p(n)
    np <- paste0("(", n, ",", p, ")")
    
    # Replace n by np
    name_set <- names(setting)
    name_set[which(name_set == "n")] <- "np"

    setting$np <- np
    setting$n <- NULL

    setting <- setting[name_set]
    
    # Format error names for display
    if ("error" %in% names(setting)) {
        setting$error <- gsub("lnorm", "log-norm", setting$error)
        setting$error <- gsub("mixnorm", "mix-norm", setting$error)
    }
    
    return(setting)
}

#' Count simulation settings
#'
#' Counts the total number of simulation setting combinations defined by a
#' setting list.
#'
#' @param setting A list containing simulation settings. Each element should be
#'   a vector of candidate values for one simulation factor.
#' @param exclude_method A logical value indicating whether the element named
#'   `"method"` should be excluded from the count. Default is `TRUE`.
#'
#' @return An integer giving the total number of simulation setting combinations.
#'
#' @details
#' This function is used to determine the number of simulation jobs. The
#' `"method"` element is usually excluded because it records the testing methods
#' evaluated under each setting, rather than defining separate data-generating
#' settings.
#'
#' @examples
#' setting <- list(
#'   SNR = seq(0, 1, 0.5),
#'   n = c(90, 120),
#'   error = c("norm", "t3"),
#'   method = c("CR")
#' )
#' count_settings(setting)
#'
#' @export
count_settings <- function(setting, exclude_method = TRUE) {
    if (exclude_method && "method" %in% names(setting)) {
        setting <- setting[names(setting) != "method"]
    }
    prod(sapply(setting, length))
}

#' Assign simulation setting indices
#'
#' Generates all combinations of indices for a collection of simulation setting
#' factors.
#'
#' @param n_args An integer vector giving the number of candidate values for each
#'   simulation setting factor.
#'
#' @return A matrix whose columns correspond to simulation jobs and whose rows
#'   correspond to simulation setting factors. Each column gives one combination
#'   of setting indices.
#'
#' @details
#' This function is used to map a job number to a specific combination of
#' simulation settings. For example, if `n_args = c(2, 3)`, then the function
#' returns all index combinations from `{1, 2} x {1, 2, 3}`.
#'
#' @examples
#' assig(c(2, 3))
#'
#' @export
assig <- function(n_args) 
{
    cargs <- vector("list", length(n_args))
    for(i in 1:length(n_args)) cargs[[i]] = 1:n_args[i]
    t(expand.grid(cargs))
}

#' Build a moving-average covariance matrix
#'
#' Constructs a banded covariance matrix based on a moving-average structure for
#' simulation studies.
#'
#' @param rho A numeric vector of moving-average coefficients.
#' @param p An integer specifying the dimension of the covariance matrix.
#'
#' @return A `p` by `p` covariance matrix generated from the moving-average
#'   structure.
#'
#' @details
#' The vector `rho` is first normalized so that the resulting covariance matrix
#' has diagonal entries equal to one. The off-diagonal entries are determined by
#' the inner products of shifted moving-average coefficients, leading to a banded
#' covariance structure.
#'
#' @examples
#' rho <- runif(10)
#' Sigma <- build_ma_sigma(rho, p = 50)
#'
#' @importFrom stats runif
#' @export
build_ma_sigma <- function(rho, p)
{
    rho = rho/sqrt(sum(rho^2))
    Sigma.x = matrix(0, p, p)
    diag(Sigma.x) = 1
    n.rho = length(rho)
    for (i in 1:(n.rho-1)) {
        tmp = rho[(1+i):n.rho]%*% rho[1:(n.rho-i)]
        diag(Sigma.x[-c(1:i),]) = tmp
        diag(Sigma.x[,-c(1:i)]) = tmp
    }
    return(Sigma.x)
}

#' Build an autoregressive covariance matrix
#'
#' Constructs an autoregressive covariance matrix for simulation studies.
#'
#' @param rho A numeric value specifying the autoregressive correlation
#'   parameter.
#' @param p An integer specifying the dimension of the covariance matrix.
#'
#' @return A `p` by `p` covariance matrix with entries
#'   \eqn{\Sigma_{ij} = \rho^{|i-j|}}.
#'
#' @details
#' This function is used to generate covariance matrices under the AR structure
#' in the simulation studies.
#'
#' @examples
#' Sigma <- build_ar_sigma(rho = 0.5, p = 10)
#'
#' @export
build_ar_sigma <- function(rho, p)
{
    Sigma.x = outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y)))
    return(Sigma.x)
}

#' Generate regression coefficients
#'
#' Generates a regression coefficient vector for simulation studies under
#' sparse, dense, random-location, or span-based alternatives.
#'
#' @param p An integer specifying the dimension of the coefficient vector.
#' @param s An integer specifying the number of nonzero coefficients. If missing
#'   or larger than `p`, it is reset to 5.
#' @param snr A numeric value specifying the signal-to-noise ratio. Required
#'   when `normalize_by_snr = TRUE`.
#' @param Sigma A covariance matrix used to normalize the coefficient vector
#'   according to the specified signal-to-noise ratio. Required when
#'   `normalize_by_snr = TRUE`.
#' @param Omat An optional matrix used to generate coefficients in the span of
#'   its columns when `loc = "span"`.
#' @param loc A character string specifying the location pattern of the
#'   coefficients. Options are `"fix"`, `"random"`, and `"span"`.
#' @param value A character string specifying the values of the nonzero
#'   coefficients. Options are `"rnorm"`, `"runif"`, and `"one"`.
#' @param normalize_by_snr A logical value indicating whether the coefficient
#'   vector should be rescaled to match the specified signal-to-noise ratio.
#'   Default is `TRUE`.
#' @param seed An integer specifying the random seed. Default is 1.
#'
#' @return A numeric vector of regression coefficients with length `p`.
#'
#' @details
#' When `loc = "fix"`, the first `s` entries are nonzero. When `loc = "random"`,
#' `s` nonzero entries are randomly selected. When `loc = "span"`, the
#' coefficient vector is generated from the column space of `Omat`.
#'
#' If `normalize_by_snr = TRUE`, the coefficient vector is rescaled so that
#' \eqn{\beta^\top \Sigma \beta} equals the specified `snr`. When `snr = 0`,
#' the function returns the zero vector.
#'
#' @examples
#' Sigma <- diag(10)
#' beta <- generate_beta(
#'   p = 10, s = 3, snr = 1, Sigma = Sigma,
#'   loc = "fix", value = "one"
#' )
#'
#' @importFrom stats rnorm runif
#' @export
generate_beta <- function(
    p, s, snr = NULL, Sigma = NULL, Omat = NULL,
    loc = c("fix", "random", "span"), value = c("rnorm", "runif", "one"),
    normalize_by_snr = TRUE, seed = 1) {
    set.seed(seed)
    if (normalize_by_snr) {
        flag <- is.null(snr) || is.null(Sigma)
        if (flag) stop("please provide snr and Sigma!")
        if (snr == 0) {
            beta = rep(0,p)
            return(beta)
        }
    }

    if (missing(s)) {
        warning("s is missing, set s as 5!")
        s <- 5
    }
    if (s > p) {
        warning("s > p, set s as 5!")
        s <- 5
    }
    loc <- match.arg(loc)
    value <- match.arg(value)

    if (loc == "fix") {
        if (value == "rnorm") {
            beta <- c(rnorm(s), rep(0, p - s))
        } else if (value == "runif") {
            beta <- c(runif(s), rep(0, p - s))
        } else {
            beta <- c(rep(1, s), rep(0, p - s))
        }
    } else if (loc == "random") {
        beta <- rep(0, p)
        if (value == "rnorm") {
            beta[sample(1:p, s)] <- rnorm(s)
        } else if (value == "runif") {
            beta[sample(1:p, s)] <- runif(s)
        } else {
            beta[sample(1:p, s)] <- 1
        }
    } else if (loc == "span") {
        if (is.null(Omat)) stop("Omat is missing!")
        if (nrow(Omat) != p) stop("row number of Omat should be p!")
        beta <- Omat %*% rnorm(ncol(Omat))
    }

    if (normalize_by_snr) {
        bSb <- sum(beta * (Sigma %*% beta))
        beta <- sqrt(snr / bSb) * beta
    }

    return(beta)
}

#' Generate random errors
#'
#' Generates random error terms from a specified distribution for simulation
#' studies.
#'
#' @param n An integer specifying the sample size.
#' @param error A character string specifying the error distribution. Available
#'   options are `"norm"`, `"t3"`, `"lnorm"`, and `"mixnorm"`.
#' @param normalize A logical value indicating whether the generated errors
#'   should be standardized to have mean zero and variance one when applicable.
#'   Default is `TRUE`.
#'
#' @return A numeric vector of generated random errors with length `n`.
#'
#' @details
#' This function is used to generate error terms in the simulation studies.
#' The supported distributions include the standard normal distribution,
#' Student's t distribution, log-normal distribution, and a contaminated
#' mixture normal distribution. When `normalize = TRUE`, the generated errors
#' are standardized whenever the corresponding mean and variance are finite.
#'
#' @examples
#' e1 <- generate_error(n = 100, error = "norm")
#' e2 <- generate_error(n = 100, error = "t3")
#' e3 <- generate_error(n = 100, error = "lnorm")
#' e4 <- generate_error(n = 100, error = "mixnorm")
#'
#' @importFrom stats rt rnorm runif rlnorm
#' @export
generate_error <- function(n, error, normalize = TRUE) {

    # Check whether a package is installed. If not, install it from CRAN.
    require_or_install <- function(pkg) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            install.packages(pkg)
        }
    }

    if (substr(error, 1, 1) == "t") {
        # Student's t distribution.
        # The variance exists when v > 2.
        v <- as.numeric(substr(error, 2, nchar(error)))
        e <- stats::rt(n, v)

        if (normalize && v != 2) {
            e <- e / sqrt(v / (v - 2))
        }

        return(e)
    }

    if (error == "norm") {
        # Standard normal distribution.
        e <- stats::rnorm(n)
        return(e)
    }

    if (error == "mixnorm") {
        # Mixture normal distribution.
        ratio <- 0.9
        sigma2_1 <- 1
        sigma2_2 <- 100

        tmp <- stats::runif(n)
        n1 <- sum(tmp > ratio)

        e <- rep(0, n)
        e[tmp > ratio]  <- stats::rnorm(n1, sd = sqrt(sigma2_2))
        e[tmp <= ratio] <- stats::rnorm(n - n1, sd = sqrt(sigma2_1))

        if (normalize) {
            e <- e / sqrt(sigma2_2 * (1 - ratio) + sigma2_1 * ratio)
        }

        return(e)
    }

    if (error == "lnorm") {
        # Log-normal distribution.
        mu <- 0
        sigma2 <- 1

        e <- stats::rlnorm(n, meanlog = mu, sdlog = sqrt(sigma2))

        if (normalize) {
            e <- (e - exp(mu + sigma2 / 2)) /
                sqrt((exp(sigma2) - 1) * exp(2 * mu + sigma2))
        }

        return(e)
    }

    stop("Unsupported error distribution: ", error)
}

#' Choose the censoring constant for a target censoring rate
#'
#' Estimates the censoring constant used to achieve a desired censoring rate in
#' simulation studies.
#'
#' @param n An integer specifying the Monte Carlo sample size used to estimate
#'   the censoring constant.
#' @param beta A numeric vector of regression coefficients.
#' @param error A character string specifying the error distribution.
#' @param mu.x A numeric vector specifying the mean of the covariates. If missing,
#'   it is set to a zero vector.
#' @param Sigma.x A covariance matrix for the covariates. If missing, it is set
#'   to the identity matrix.
#' @param design_distribution A character string specifying the distribution used
#'   to generate the covariates. Options are `"norm"`, `"chi"`, and `"unif"`.
#' @param normalize_error A logical value indicating whether the generated errors
#'   should be standardized. Default is `TRUE`.
#' @param cenRate A numeric value specifying the target censoring rate. Default
#'   is `0.2`.
#' @param Gamma.x An optional matrix satisfying
#'   \eqn{\Gamma_x \Gamma_x^\top = \Sigma_x}, used to generate covariates from
#'   a factor model.
#'
#' @return A numeric value giving the censoring constant.
#'
#' @details
#' This function generates a large auxiliary sample and estimates the constant
#' added to the censoring time so that the simulated data approximately achieve
#' the target censoring rate. The censoring time is generated from a chi-square
#' distribution with one degree of freedom.
#'
#' @examples
#' Sigma <- diag(5)
#' beta <- rep(1, 5)
#' c0 <- control_c_cersor(
#'   n = 1000, beta = beta, error = "norm",
#'   mu.x = rep(0, 5), Sigma.x = Sigma,
#'   cenRate = 0.2
#' )
#'
#' @importFrom stats rnorm runif rchisq quantile
#' @importFrom mnormt rmnorm
#' @export
control_c_cersor <- function(
    n, beta, error, mu.x, Sigma.x,
    design_distribution = c("norm", "chi", "unif"),
    normalize_error = TRUE, cenRate = 0.2, Gamma.x = NULL)
{
    # error distirbuion ========================================================
    e <- generate_error(n, error, normalize_error)

    p <- length(beta)
    if (missing(mu.x)) {
        mu.x = rep(0, p)
    }
    if (missing(Sigma.x)) {
        Sigma.x = diag(p)
    }

    # x from a factor model, x = Gamma.x %*% z
    design_distribution <- match.arg(design_distribution)
    if (design_distribution == "norm") {
        if (is.null(Gamma.x)) {
            x <- rmnorm(n, mean = mu.x, varcov = Sigma.x)
        } else {
            x <- matrix(rnorm(n * p), n, p) %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "chi") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            x <- (matrix(rnorm(n * p), n, p)^2 - 1) * sqrt(0.5)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "unif") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            x <- (matrix(runif(n * p), n, p) - 0.5) * sqrt(12)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    }

    xb <- x %*% beta + e

    C <- rchisq(n, 1)
    c_cersor <- quantile(xb - C, 1 - cenRate)

    return(c_cersor)
}

#' Generate simulated transformation-model data
#'
#' Generates simulated data under nonparametric transformation models with
#' optional censoring.
#'
#' @param n An integer specifying the sample size.
#' @param beta A numeric vector of regression coefficients.
#' @param censor A character string specifying the censoring mechanism. Options
#'   are `"non"` for uncensored data, `"random"` for random right censoring, and
#'   `"right"` for fixed right censoring.
#' @param tran A character string specifying the transformation function.
#'   Options are `"lm"`, `"log"`, and `"boxcox"`.
#' @param error A character string specifying the error distribution.
#' @param mu.x A numeric vector specifying the mean of the covariates. If missing,
#'   it is set to a zero vector.
#' @param Sigma.x A covariance matrix for the covariates. If missing, it is set
#'   to the identity matrix.
#' @param design_distribution A character string specifying the distribution used
#'   to generate the covariates. Options are `"norm"`, `"chi"`, and `"unif"`.
#' @param normalize_error A logical value indicating whether the generated errors
#'   should be standardized. Default is `TRUE`.
#' @param c_cersor A numeric value used as the censoring constant for random
#'   right censoring. If `NULL`, it is set to zero with a warning.
#' @param cenRate A numeric value specifying the target censoring rate for fixed
#'   right censoring. Default is `0.2`.
#' @param Gamma.x An optional matrix satisfying
#'   \eqn{\Gamma_x \Gamma_x^\top = \Sigma_x}, used to generate covariates from
#'   a factor model.
#' @param lambda_boxcox A numeric value specifying the Box-Cox transformation
#'   parameter. Default is `1`.
#'
#' @return A list containing:
#' \item{x}{The generated covariate matrix.}
#' \item{y}{The observed response vector.}
#' \item{beta}{The regression coefficient vector used to generate the data.}
#' \item{status}{The censoring indicator, where `TRUE` indicates an observed
#' event and `FALSE` indicates right censoring.}
#'
#' @details
#' The latent linear predictor is generated as
#' \deqn{x_i^\top \beta + \epsilon_i.}
#' The observed response is then obtained by applying the specified
#' transformation function. For random right censoring, the censoring time is
#' generated from a chi-square distribution with one degree of freedom plus the
#' censoring constant `c_cersor`. For fixed right censoring, the censoring
#' threshold is chosen from the empirical quantile of the latent response.
#'
#' @examples
#' Sigma <- diag(5)
#' beta <- rep(0, 5)
#' dat <- generator(
#'   n = 50, beta = beta, censor = "random", tran = "lm",
#'   error = "norm", mu.x = rep(0, 5), Sigma.x = Sigma,
#'   design_distribution = "norm", c_cersor = 1
#' )
#'
#' @importFrom stats rnorm runif rchisq quantile
#' @importFrom mnormt rmnorm
#' @export
generator <- function(
    n, beta, censor = c("non", "random", "right"), tran = c("lm", "log", "boxcox"),
    error, mu.x, Sigma.x, design_distribution = c("norm", "chi", "unif"),
    normalize_error = TRUE, c_cersor = NULL, cenRate = 0.2, Gamma.x = NULL,
    lambda_boxcox = 1)
{
    # error distirbuion ========================================================
    e <- generate_error(n, error, normalize_error)

    p <- length(beta)
    if (missing(mu.x)) {
        mu.x = rep(0, p)
    }
    if (missing(Sigma.x)) {
        Sigma.x = diag(p)
    }

    # x from a factor model, x = Gamma.x %*% z
    design_distribution <- match.arg(design_distribution)
    if (design_distribution == "norm") {
        if (is.null(Gamma.x)) {
            x <- rmnorm(n, mean = mu.x, varcov = Sigma.x)
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            x <- matrix(rnorm(n * p), n, p) %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "chi") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled chi-square distribution,
            # which has zero mean and unit variance
            x <- (matrix(rnorm(n * p), n, p)^2 - 1) * sqrt(0.5)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "unif") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled uniform distribution,
            # which has zero mean and unit variance
            x <- (matrix(runif(n * p), n, p) - 0.5) * sqrt(12)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    }

    xb <- x %*% beta + e

    # transformation ===========================================================
    tran <- match.arg(tran)
    if (tran == "lm") g <- function(x) x
    if (tran == "log") g <- function(x) exp(x)
    if (tran == "boxcox") {
        if (lambda_boxcox == 0) {
            g <- function(x) exp(x)
        } else {
            g <- function(x) (lambda_boxcox*x + 1.0)^(1.0/lambda_boxcox)
        }
    }

    if (censor == "non") {
        y <- g(xb)
        status <- rep(1, n)
    } else if (censor == "random") {
        # random censor
        if (is.null(c_cersor)) {
            warning("c_cersor is null! We recommend using the 'control_c_cersor' function to calculate it to provides desired censor rate.")
            c_cersor = 0
        }
        C <- rchisq(n, 1) + c_cersor
        status <- xb < C
        y <- g(pmin(xb, C))
    } else if (censor == "right") {
        # right censor
        C <- quantile(xb, 1 - cenRate)
        status <- xb < C
        y <- g(pmin(xb, C))
    }

    # In the transformation model, y may be infinity.
    y[which(y == "NaN" | y == "Inf")] <- max(y[which(y != "NaN" & y != "Inf")])

    list(x      = x,
         y      = y,
         beta   = beta,
         status = status)
}

#' Summarize simulation results
#'
#' Reads saved p-value matrices from simulation outputs and computes empirical
#' rejection rates for all simulation settings.
#'
#' @param setting A list containing simulation settings. The last element should
#'   contain the names of the methods.
#' @param na.rm A logical value indicating whether missing values should be
#'   removed when computing empirical rejection rates. Default is `FALSE`.
#' @param folder A character string specifying the folder where simulation
#'   result files are saved. Default is `"result"`.
#' @param missing_scenario A logical value indicating whether missing result
#'   files should be recorded as zero. Default is `FALSE`.
#'
#' @return An array containing empirical rejection rates across all simulation
#'   settings and methods.
#'
#' @details
#' For each simulation setting, this function reads the corresponding `.rds`
#' file containing p-values from Monte Carlo replications. The empirical
#' rejection rate is computed as the proportion of p-values less than or equal
#' to 0.05. Under the null hypothesis, this quantity corresponds to the
#' empirical Type-I error rate; under alternatives, it corresponds to empirical
#' power.
#'
#' The file names are assumed to be generated by concatenating the setting
#' values with hyphens and using the `.rds` extension.
#'
#' @examples
#' setting <- list(
#'   SNR = c(0, 1),
#'   n = c(90, 120),
#'   error = c("norm", "t3"),
#'   method = c("CR")
#' )
#' # res <- summary_results(setting, folder = "Results_simulation")
#'
#' @export
summary_results <- function(setting, na.rm = FALSE, folder = "result", missing_scenario = FALSE) {
    n_args   = sapply(setting, length)
    nmethod  = n_args[length(n_args)]
    n_args   = n_args[-length(n_args)]
    nset     = length(setting)-1
    jobs     = assig(n_args)
    ntc      = prod(n_args)
    neach    = 1000
    methods <- setting[[length(setting)]]

    if (missing_scenario) {
        message("Results of missing scenario with be recorded as 0.")
    }

    paramList = vector("list", nset)
    size_all = c()
    for (number in 1:ntc) {
        id <- jobs[, number]

        filename <- paste0(folder, "/")
        for (j in 1:nset) {
            paramList[[j]] <- setting[[j]][id[j]]
            filename <- paste0(filename, paramList[[j]], "-")
        }
        filename <- paste0(filename, ".rds")

        if (missing_scenario && (!file.exists(filename))) {
            size_k <- rep(0, nmethod)
            size_all <- c(size_all, size_k)
            next
        }

        output = readRDS(filename)
        size_k <- colMeans(output <= 0.05, na.rm = na.rm)
        if (length(size_k) < nmethod) {
            size_k0 <- rep(0, nmethod)
            names(size_k0) <- methods
            size_k0[colnames(output)] <- size_k
            size_k <- size_k0
        }

        size_all = c(size_all, size_k)

        if (number %% 50 == 0) cat(paste0(number,", "))
    }

    setting0 = setting
    setting0 = setting0[c(length(setting), 1:(length(setting)-1))]
    size_Res = array(size_all, dim = c(nmethod, n_args), dimnames = setting0)
    size_Res = aperm(size_Res, c(c(1:nset)+1,1))
    return(size_Res)
}

#' Compute reverse cumulative products
#'
#' Computes cumulative products of a numeric vector in reverse order.
#'
#' @param vec A numeric vector.
#'
#' @return A numeric vector of the same length as `vec`. The first element is
#'   the last element of `vec`, and each subsequent element is the cumulative
#'   product of elements from the end of `vec` moving backward.
#'
#' @details
#' For a vector `vec = c(v1, v2, ..., vn)`, this function returns
#' `c(vn, vn * v_{n-1}, ..., vn * v_{n-1} * ... * v1)`.
#'
#' @examples
#' revcumprod(c(2, 3, 4))
#'
#' @export
revcumprod <- function(vec) {
    n = length(vec) 
    if (n == 1) return(vec)
    res = rep(0, n) 
    res[1] = vec[n]
    for (i in 2:n) {
        res[i] = res[i-1] * vec[n+1-i]
    }
    return(res)
}

#' Convert a result array to a data frame
#'
#' Converts summarized simulation results stored in an array into a data frame
#' for plotting and table construction.
#'
#' @param ResArray An array containing summarized simulation results, such as
#'   empirical rejection rates.
#' @param valuesList A list containing the simulation setting values
#'   corresponding to the dimensions of `ResArray`.
#' @param nameVec An optional character vector specifying the column names used
#'   for the setting variables in the returned data frame. If `NULL`, the names
#'   of `valuesList` are used.
#'
#' @return A data frame containing the array values and the corresponding
#'   simulation setting labels.
#'
#' @details
#' This function first reformats the simulation settings using
#' `reformat_setting()`, replacing `n` by `(n,p)` when applicable. It then
#' expands the array into a long-format data frame, where each row corresponds
#' to one combination of simulation settings and one summarized result.
#'
#' The variable `np`, if present, is converted to a factor to preserve the order
#' of sample-size and dimension pairs in tables and plots.
#'
#' @examples
#' setting <- list(
#'   SNR = c(0, 1),
#'   n = c(90, 120),
#'   error = c("norm", "t3"),
#'   method = c("CR")
#' )
#' ResArray <- array(runif(8), dim = c(2, 2, 2))
#' df <- array2df(ResArray, setting)
#'
#' @importFrom stats setNames
#' @export
array2df <- function(ResArray, valuesList, nameVec = NULL) {
    valuesList <- reformat_setting(valuesList)
    if (is.null(nameVec)) nameVec = names(valuesList)
    df <- data.frame(values = c(ResArray))
    dims = dim(ResArray)
    eachVec = c(1, cumprod(dims))
    timesVec = c(rev(revcumprod(dims))[-1], 1)
    for (i in 1:length(valuesList)) {
        df[nameVec[i]] = rep(valuesList[[i]], times=timesVec[i], each = eachVec[i])
    }

    df$np <- factor(df$np, levels = as.character(valuesList$np))
    return(df)
}

#' Plot empirical power curves for one simulation setting
#'
#' Generates a faceted empirical power plot from summarized simulation results.
#'
#' @param df A data frame containing summarized simulation results. It should
#'   include columns such as `SNR` or `beta`, `values`, `method`, `rowfeature`,
#'   and `colfeature`.
#' @param setting A list containing simulation settings. The element `method`
#'   is used to determine the number of methods and plotting colors.
#' @param col A character vector specifying colors for different methods.
#' @param sizeline A logical value indicating whether to add a horizontal line
#'   at 0.05. Default is `TRUE`.
#' @param myylab A character string specifying the y-axis label. Default is
#'   `"Power"`.
#'
#' @return A `ggplot` object showing empirical power curves.
#'
#' @details
#' If the first element of `setting` is named `"SNR"`, the x-axis is the
#' signal-to-noise ratio. Otherwise, the x-axis is `beta`. The plot is faceted
#' by `rowfeature` and `colfeature`, and different testing methods are shown
#' using different line types and colors.
#'
#' @examples
#' # p <- power_plot(df, setting)
#'
#' @importFrom ggplot2 ggplot aes xlab ylab geom_line geom_hline facet_grid
#' @importFrom ggplot2 labeller label_parsed scale_fill_manual theme_bw theme
#' @importFrom ggplot2 element_text element_blank guide_legend guides
#' @importFrom stringr str_replace
#' @importFrom grid unit
#' @export
power_plot <- function(
    df, setting, col = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628"), 
    sizeline = TRUE, myylab = "Power"
) {
  library(ggplot2)
  col = col[1:length(setting$method)]
  
  if (names(setting)[1] == "SNR") {
    bp2 = bp2 <- ggplot(df, aes(x=SNR, y=values, linetype=method, colour=method)) + xlab("SNR")
  } else {
    bp2 = bp2 <- ggplot(df, aes(x=beta, y=values, linetype=method, colour=method)) + xlab(expression(beta))
  }
  
  custom_labeller <- function(labels) {
    # Replacing hyphens with a space or other character
    labels <- str_replace(labels, "-", " ")
    labels
  }
  
  bp2 <- bp2 +
    
    geom_line(size = 1, aes(color=method)) + 
    
    # facet_grid(rowfeature ~ colfeature, scales = "fixed", labeller=labeller(colfeature = custom_labeller)) +  
    facet_grid(rowfeature ~ colfeature, scales = "fixed", labeller=labeller(colfeature = label_parsed)) +
    
    scale_fill_manual(values=col)+
    
    theme_bw() +
    
    ylim(0,1) + 
    
    theme(plot.title = element_text(size=25,face = "bold",hjust=0.5),
          axis.title.x = element_text(size=20,face = "bold"),
          axis.text.x = element_text(size=15,face = "bold"),
          axis.title.y = element_text(size=25,face = "bold"),
          axis.text.y = element_text(size=25,face = "bold"),
          legend.position="bottom",
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          #strip.background = element_blank(),
          # strip.background = element_rect(fill = "grey"),
          strip.text = element_text(size=25,face = "bold"),
          # legend.title=element_text(size=25,face = "bold"),
          legend.title=element_blank(),
          legend.key.width= unit(1.5, 'cm'),
          legend.text=element_text(size=25,face = "bold")) + 
    
    ylab(myylab) + guides(color = guide_legend(nrow = 1))
  
  if (sizeline) {
    bp2 = bp2 + geom_hline(aes(yintercept = 0.05), color = "black", linetype="solid",size = 1)
  }
  
  bp2
}

#' Plot summarized simulation results
#'
#' Generates empirical power plots from a summarized simulation result data frame.
#'
#' @param df A data frame containing summarized simulation results. It should
#'   include the empirical rejection rate column `values` and columns
#'   corresponding to the simulation settings.
#' @param setting A list containing simulation settings. The list is reformatted
#'   internally using `reformat_setting()`.
#' @param plotfeature A character vector specifying the simulation factors used
#'   to split the results into different plot files.
#' @param rowfeature A character vector specifying the simulation factors used
#'   to define facet rows. If `NULL`, the function uses an existing
#'   `rowfeature` column in `df`.
#' @param colfeature A character vector specifying the simulation factors used
#'   to define facet columns. If `NULL`, the function uses an existing
#'   `colfeature` column in `df`.
#' @param removespace A logical value indicating whether spaces should be
#'   removed from combined facet labels. Default is `FALSE`.
#' @param na.rm A logical value indicating whether rows with missing `values`
#'   should be removed before plotting. Default is `TRUE`.
#' @param plot_method An optional character vector specifying the methods to be
#'   included in the plots. If `NULL`, all methods in `setting$method` are used.
#' @param plotfile A character string specifying the folder where plots are
#'   saved. Default is `"Plot/"`.
#' @param plotfile_prefix A character string specifying the prefix of saved plot
#'   files. Default is `"Power_"`.
#' @param returnplot A logical value indicating whether to return the generated
#'   plot objects and file information. Default is `FALSE`.
#' @param ... Additional arguments passed to `power_plot()`.
#'
#' @return If `returnplot = TRUE`, a list containing:
#' \item{plotList}{A list of generated `ggplot` objects.}
#' \item{filenames}{A list of file name prefixes for the generated plots.}
#' \item{width}{The plot width used in `ggsave()`.}
#' \item{height}{The plot height used in `ggsave()`.}
#' If `returnplot = FALSE`, the function saves the plots and returns `NULL`
#' invisibly.
#'
#' @details
#' This function is used by the simulation summary scripts to generate empirical
#' power curves. It first prepares display labels, constructs facet variables,
#' optionally filters methods and removes missing values, and then creates one
#' plot for each combination of `plotfeature`. Each plot is saved as a PDF file.
#'
#' @examples
#' # plot_data(
#' #   df = df,
#' #   setting = setting,
#' #   plotfeature = c("cenRate", "tran"),
#' #   rowfeature = "np",
#' #   colfeature = "error"
#' # )
#'
#' @importFrom ggplot2 ggsave
#' @importFrom latex2exp TeX
#' @importFrom stats na.omit
#' @export
plot_data <- function(
    df, setting, plotfeature, rowfeature = NULL, colfeature = NULL, removespace = FALSE, 
    na.rm = TRUE, plot_method = NULL, plotfile = "Plot/", plotfile_prefix = "Power_", 
    returnplot = FALSE, ...
) {
  library(ggplot2)
  library(latex2exp)
  ifelse(!dir.exists(plotfile), dir.create(plotfile), FALSE)
  setting <- reformat_setting(setting)
  df$error <- factor(
        df$error,
        levels = c("norm", "t3", "log-norm", "mix-norm"),
        labels = c(
            parse(text = TeX("norm", bold = TRUE)),
            parse(text = TeX("$t_3$", bold = TRUE)),
            parse(text = TeX("log$\\;$norm", bold = TRUE)),
            parse(text = TeX("mixture$\\;$normals", bold = TRUE))
        )
    )
  
  if (is.null(plot_method)) {
    plot_method = setting$method
  } else {
    df <- df[df$method %in% plot_method, ]
  }
  df$method <- factor(df$method, levels=plot_method)
  
  if (is.null(rowfeature)) {
    if (!("rowfeature" %in% colnames(df))) stop("rowfeature is not in df, and not provided!")
    message("rowfeature is NULL, use provided rowfeature!")
  } else {
    if (length(rowfeature) == 1) {
      df$rowfeature <- df[,rowfeature]
    } else {
      n_rowfeature = length(rowfeature)
      tmp = paste0( rowfeature[1], " = ", df[,rowfeature[1]] )
      for (i in 2:n_rowfeature) {
        tmp = paste0(tmp, ", ", rowfeature[i], " = ", df[,rowfeature[i]])
      }
      if (removespace) tmp = gsub(" ", "", tmp)
      df$rowfeature <- paste0("'", tmp, "'")
    }
  }
  
  if (is.null(colfeature)) {
    if (!("colfeature" %in% colnames(df))) stop("colfeature is not in df, and not provided!")
    message("colfeature is NULL, use provided colfeature!")
  } else {
    if (length(colfeature) == 1) {
      df$colfeature <- df[,colfeature]
    } else {
      n_colfeature = length(colfeature)
      tmp = paste0( colfeature[1], " = ", df[,colfeature[1]] )
      for (i in 2:n_colfeature) {
        tmp = paste0(tmp, ", ", colfeature[i], " = ", df[,colfeature[i]])
      }
      if (removespace) tmp = gsub(" ", "", tmp)
      df$colfeature <- paste0("'", tmp, "'")
    }
  }
  
  # df$situation <- paste0("'rho = ", df$rho, ", CR = ", df$cenRate, "'")
  # # df$situation <- factor(df$situation, levels=unique(df$situation))
  # df$situation <- factor(df$situation, labels=c(
  #     'rho = 0, CR = 0.1'  =parse(text=TeX('$\\rho$ = 0, CR = 0.1')), 
  #     'rho = 0.3, CR = 0.1'=parse(text=TeX('$\\rho$ = 0.3, CR = 0.1')),
  #     'rho = 0, CR = 0.2'  =parse(text=TeX('$\\rho$ = 0, CR = 0.2')),
  #     'rho = 0.3, CR = 0.2'=parse(text=TeX('$\\rho$ = 0.3, CR = 0.2'))))
  if (na.rm && any(is.na(df$values))) {
    message("NA values are removed!")
    df <- df[!is.na(df$values), ]
  }
  
  prefix = paste0(plotfile, plotfile_prefix)
  
  mywidth = length(unique(df$colfeature))*4
  myheight = length(unique(df$rowfeature))*4
  
  
  
  n_args   = sapply(setting, length)[plotfeature]
  ntc      = prod(n_args) 
  jobs     = assig(n_args)
  plotList = vector("list", ntc)
  plot_rgs = vector("list", ntc)
  for (number in 1:ntc) {
    id = jobs[, number]
    df0 = df
    filename = prefix
    for (i in 1:length(plotfeature)) {
      pfvalue = setting[[plotfeature[i]]][id[i]]
      df0 = df0[df0[,plotfeature[i]] == pfvalue, ]
      # filename = paste0(filename, plotfeature, "_", pfvalue)
      filename = paste0(filename, pfvalue, "_")
    }
    plot_rgs[[number]] <- filename
    # filename = paste0(filename, ".eps")
    plotList[[number]] <- power_plot( df0, setting, ... )
    # ggsave(file=filename, plot = p1, width = 14.9, height = 12, units = "in", dpi = 1000)
    # filename = paste0(filename, ".png")
    # ggsave(file=filename, plot = plotList[[number]], width = mywidth, height = myheight, units = "in")
    # filename = paste0(filename, ".eps")
    # ggsave(file=filename, plot = plotList[[number]], width = mywidth, height = myheight, device = "eps")
    filename = paste0(filename, ".pdf")
    ggsave(file=filename, plot = plotList[[number]], width = mywidth, height = myheight, units = "in")
  }
  
  if (returnplot) {
    return(
      list(
        plotList = plotList,
        filenames = plot_rgs,
        width = mywidth, 
        height = myheight
      )
    )
  }
}

#' Summarize empirical Type-I error rates
#'
#' Extracts and formats empirical Type-I error rates from summarized simulation
#' results.
#'
#' @param df A data frame containing summarized simulation results, typically
#'   generated by `array2df()`. It should contain a column named `values` and
#'   columns corresponding to the simulation settings.
#' @param setting A list containing simulation settings. The first element is
#'   assumed to be the signal strength variable, such as `SNR`.
#' @param colfeature A character vector specifying the simulation factors to be
#'   displayed across columns.
#' @param rowfeature An optional character vector specifying the simulation
#'   factors to be displayed across rows. If `NULL`, all remaining factors except
#'   `colfeature` and `method` are used.
#' @param table_method An optional character vector specifying the methods to be
#'   included in the table. If `NULL`, all methods in `setting$method` are used.
#' @param include_feature A logical value indicating whether row feature labels
#'   should be included in the returned table. Default is `FALSE`.
#'
#' @return A matrix or data frame containing formatted empirical Type-I error
#'   rates. Missing values are displayed as blank entries.
#'
#' @details
#' This function extracts simulation results under the null hypothesis by
#' selecting rows for which the first setting variable, usually `SNR`, equals
#' zero. It then arranges empirical Type-I error rates into a table according to
#' the specified row and column features. The output is formatted to three
#' decimal places for direct use in simulation summary tables.
#'
#' @examples
#' # size_table <- summary_size(
#' #   df = df,
#' #   setting = setting,
#' #   colfeature = "tran",
#' #   include_feature = TRUE
#' # )
#'
#' @export
summary_size <- function(
    df, setting, colfeature, rowfeature = NULL, table_method = NULL,
    include_feature = FALSE
) {
    setting <- reformat_setting(setting)
    df <- df[df[, names(setting)[1]] == 0, ]
    features <- names(setting)
    features <- features[-1]
    features <- features[-length(features)]

    n_args <- sapply(setting, length)[features]
    jobs   <- assig(n_args)

    if (is.null(table_method)) {
        table_method <- setting$method
    } else {
        df <- df[df$method %in% table_method, ]
    }

    if (is.null(rowfeature)) {
        rowfeature <- rev(setdiff(features, colfeature))
    }
    nmethod  <- length(table_method)
    rfeature <- match(rowfeature, features)
    cfeature <- match(colfeature, features)
    if (length(rfeature) == 1) { rweight = 1 } else { rweight = c(rev(cumprod(rev(n_args[rfeature])))[-1], 1) }
    if (length(cfeature) == 1) { cweight = 1 } else { cweight = c(rev(cumprod(rev(n_args[cfeature])))[-1], 1) }
    get_rid <- function(idrf, rweight) return( (idrf-1)%*%rweight+1 )
    get_cid <- function(idcf, cweight, nmethod) { tmp = (idcf-1)%*%cweight; return(c((tmp*nmethod+1):(tmp*nmethod+nmethod))) }

    ResTable = matrix(NA, prod(n_args[rfeature]), prod(n_args[cfeature])*nmethod)
    for (number in 1:prod(n_args)) {
        id = jobs[, number]
        df0 = df
        for (i in 1:length(features)) {
            pfvalue = setting[[features[i]]][id[i]]
            df0 = df0[df0[,features[i]] == pfvalue, ]
        }

        rid = get_rid( id[rfeature], rweight )
        cid = get_cid( id[cfeature], cweight, nmethod )

        ResTable[rid, cid] = df0$values
    }

    ResTable <- apply(ResTable, c(1,2), function(x) ifelse(is.na(x), "     ", sprintf("%.3f", x)))

    if (include_feature) {
        rfeature <- rev(rowfeature)
        n_rfeature <- n_args[rfeature]
        n_row <- prod(n_rfeature)
        df_rfeature <- matrix(nrow = n_row, ncol = length(rfeature))
        each_i <- 1
        for (i in 1:length(rfeature)) {
            n_r_i <- length(setting[[rfeature[i]]])
            df_rfeature[, i] = rep(
                rep(setting[[rfeature[i]]], each = each_i), n_row / (each_i * n_r_i))
            each_i <- each_i * n_r_i
        }
        colnames(df_rfeature) <- rfeature
        df_rfeature <- as.data.frame(df_rfeature)[,rev(rfeature)]
        ResTable <- cbind(df_rfeature, ResTable)
    }

    return(ResTable)
}

#' Choose the censoring constant for covariate-dependent censoring
#'
#' Estimates the censoring constant used to achieve a desired censoring rate
#' under a covariate-dependent censoring mechanism.
#'
#' @param n An integer specifying the Monte Carlo sample size used to estimate
#'   the censoring constant.
#' @param beta A numeric vector of regression coefficients.
#' @param error A character string specifying the error distribution.
#' @param mu.x A numeric vector specifying the mean of the covariates. If missing,
#'   it is set to a zero vector.
#' @param Sigma.x A covariance matrix for the covariates. If missing, it is set
#'   to the identity matrix.
#' @param design_distribution A character string specifying the distribution used
#'   to generate the covariates. Options are `"norm"`, `"chi"`, and `"unif"`.
#' @param normalize_error A logical value indicating whether the generated errors
#'   should be standardized. Default is `TRUE`.
#' @param cenRate A numeric value specifying the target censoring rate. Default
#'   is `0.2`.
#' @param Gamma.x An optional matrix satisfying
#'   \eqn{\Gamma_x \Gamma_x^\top = \Sigma_x}, used to generate covariates from
#'   a factor model.
#'
#' @return A numeric value giving the censoring constant.
#'
#' @details
#' This function is similar to `control_c_cersor()`, but the censoring time
#' depends on the first covariate through \eqn{C_i = X_{i1}^2 + c_0}. The
#' constant \eqn{c_0} is estimated from a large auxiliary sample so that the
#' simulated data approximately achieve the target censoring rate.
#'
#' @examples
#' Sigma <- diag(5)
#' beta <- rep(1, 5)
#' c0 <- control_c_dependent_censoring(
#'   n = 1000, beta = beta, error = "norm",
#'   mu.x = rep(0, 5), Sigma.x = Sigma,
#'   cenRate = 0.2
#' )
#'
#' @importFrom stats rnorm runif quantile
#' @importFrom mnormt rmnorm
#' @export
control_c_dependent_censoring <- function(
    n, beta, error, mu.x, Sigma.x,
    design_distribution = c("norm", "chi", "unif"),
    normalize_error = TRUE, cenRate = 0.2, Gamma.x = NULL)
{
    # error distirbuion ========================================================
    e <- generate_error(n, error, normalize_error)

    p <- length(beta)
    if (missing(mu.x)) {
        mu.x = rep(0, p)
    }
    if (missing(Sigma.x)) {
        Sigma.x = diag(p)
    }

    # x from a factor model, x = Gamma.x %*% z
    design_distribution <- match.arg(design_distribution)
    if (design_distribution == "norm") {
        if (is.null(Gamma.x)) {
            x <- rmnorm(n, mean = mu.x, varcov = Sigma.x)
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            x <- matrix(rnorm(n * p), n, p) %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "chi") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled chi-square distribution,
            # which has zero mean and unit variance
            x <- (matrix(rnorm(n * p), n, p)^2 - 1) * sqrt(0.5)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "unif") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled uniform distribution,
            # which has zero mean and unit variance
            x <- (matrix(runif(n * p), n, p) - 0.5) * sqrt(12)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    }

    xb <- x %*% beta + e

    C <- x[,1]^2
    c_cersor <- quantile(xb - C, 1 - cenRate)

    return(c_cersor)
}

#' Generate covariate-dependent right-censored data
#'
#' Generates simulated data under a nonparametric transformation model with a
#' covariate-dependent right-censoring mechanism.
#'
#' @param n An integer specifying the sample size.
#' @param beta A numeric vector of regression coefficients.
#' @param censor A character string specifying the censoring mechanism. Options
#'   are `"non"` for uncensored data, `"random"` for covariate-dependent random
#'   right censoring, and `"right"` for fixed right censoring.
#' @param tran A character string specifying the transformation function.
#'   Options are `"lm"`, `"log"`, and `"boxcox"`.
#' @param error A character string specifying the error distribution.
#' @param mu.x A numeric vector specifying the mean of the covariates. If missing,
#'   it is set to a zero vector.
#' @param Sigma.x A covariance matrix for the covariates. If missing, it is set
#'   to the identity matrix.
#' @param design_distribution A character string specifying the distribution used
#'   to generate the covariates. Options are `"norm"`, `"chi"`, and `"unif"`.
#' @param normalize_error A logical value indicating whether the generated errors
#'   should be standardized. Default is `TRUE`.
#' @param c_cersor A numeric value used as the censoring constant for
#'   covariate-dependent random right censoring. If `NULL`, it is set to zero
#'   with a warning.
#' @param cenRate A numeric value specifying the target censoring rate for fixed
#'   right censoring. Default is `0.2`.
#' @param Gamma.x An optional matrix satisfying
#'   \eqn{\Gamma_x \Gamma_x^\top = \Sigma_x}, used to generate covariates from
#'   a factor model.
#' @param lambda_boxcox A numeric value specifying the Box-Cox transformation
#'   parameter. Default is `1`.
#'
#' @return A list containing:
#' \item{x}{The generated covariate matrix.}
#' \item{y}{The observed response vector.}
#' \item{beta}{The regression coefficient vector used to generate the data.}
#' \item{status}{The censoring indicator, where `TRUE` indicates an observed
#' event and `FALSE` indicates right censoring.}
#'
#' @details
#' The latent response is generated as
#' \deqn{x_i^\top \beta + \epsilon_i.}
#' For covariate-dependent random right censoring, the censoring time is
#' generated as
#' \deqn{C_i = X_{i1}^2 + c_0,}
#' where \eqn{c_0} is the censoring constant supplied by `c_cersor`. The observed
#' response is obtained by applying the specified transformation to the minimum
#' of the latent response and censoring time.
#'
#' @examples
#' Sigma <- diag(5)
#' beta <- rep(0, 5)
#' dat <- generator_dependent_censoring(
#'   n = 50, beta = beta, censor = "random", tran = "lm",
#'   error = "norm", mu.x = rep(0, 5), Sigma.x = Sigma,
#'   design_distribution = "norm", c_cersor = 1
#' )
#'
#' @importFrom stats rnorm runif quantile
#' @importFrom mnormt rmnorm
#' @export
generator_dependent_censoring <- function(
    n, beta, censor = c("non", "random", "right"), tran = c("lm", "log", "boxcox"),
    error, mu.x, Sigma.x, design_distribution = c("norm", "chi", "unif"),
    normalize_error = TRUE, c_cersor = NULL, cenRate = 0.2, Gamma.x = NULL,
    lambda_boxcox = 1)
{
    # error distirbuion ========================================================
    e <- generate_error(n, error, normalize_error)

    p <- length(beta)
    if (missing(mu.x)) {
        mu.x = rep(0, p)
    }
    if (missing(Sigma.x)) {
        Sigma.x = diag(p)
    }

    # x from a factor model, x = Gamma.x %*% z
    design_distribution <- match.arg(design_distribution)
    if (design_distribution == "norm") {
        if (is.null(Gamma.x)) {
            x <- rmnorm(n, mean = mu.x, varcov = Sigma.x)
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            x <- matrix(rnorm(n * p), n, p) %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "chi") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled chi-square distribution,
            # which has zero mean and unit variance
            x <- (matrix(rnorm(n * p), n, p)^2 - 1) * sqrt(0.5)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "unif") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled uniform distribution,
            # which has zero mean and unit variance
            x <- (matrix(runif(n * p), n, p) - 0.5) * sqrt(12)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    }

    xb <- x %*% beta + e

    # transformation ===========================================================
    tran <- match.arg(tran)
    if (tran == "lm") g <- function(x) x
    if (tran == "log") g <- function(x) exp(x)
    if (tran == "boxcox") {
        if (lambda_boxcox == 0) {
            g <- function(x) exp(x)
        } else {
            g <- function(x) (lambda_boxcox*x + 1.0)^(1.0/lambda_boxcox)
        }
    }

    if (censor == "non") {
        y <- g(xb)
        status <- rep(1, n)
    } else if (censor == "random") {
        # random censor
        if (is.null(c_cersor)) {
            warning("c_cersor is null! We recommend using the 'control_c_cersor' function to calculate it to provides desired censor rate.")
            c_cersor = 0
        }
        C <- x[,1]^2 + c_cersor
        status <- xb < C
        y <- g(pmin(xb, C))
    } else if (censor == "right") {
        # right censor
        C <- quantile(xb, 1 - cenRate)
        status <- xb < C
        y <- g(pmin(xb, C))
    }

    # In the transformation model, y may be infinity.
    y[which(y == "NaN" | y == "Inf")] <- max(y[which(y != "NaN" & y != "Inf")])

    list(x      = x,
         y      = y,
         beta   = beta,
         status = status)
}


#' Zhong and Chen global test
#'
#' Implements the high-dimensional global test of Zhong and Chen (2011) for
#' uncensored data.
#'
#' @param X A numeric covariate matrix with rows corresponding to observations
#'   and columns corresponding to covariates.
#' @param y A numeric response vector.
#' @param maxabsy A numeric threshold used to truncate extremely large absolute
#'   response values. Default is `exp(10)`.
#'
#' @return A list containing:
#' \item{Tn}{The standardized test statistic.}
#' \item{pvals}{The one-sided p-value of the test.}
#'
#' @details
#' This function calls the C routine `_ZC_Test` to compute the test statistic and
#' variance components. If the response contains extremely large values, it is
#' truncated at `maxabsy` to avoid numerical instability.
#'
#' @references
#' Zhong, P.-S. and Chen, S. X. (2011). Tests for high-dimensional regression
#' coefficients with factorial designs. \emph{Journal of the American
#' Statistical Association}, 106(493), 260--274.
#'
#' @examples
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' res <- zctest(X, y)
#'
#' @importFrom stats pnorm
#' @useDynLib hdcrtRepro, .registration = TRUE
#' @export
zctest = function(X, y, maxabsy = exp(10))
{
    param     = dim(X)
    if (max(abs(y)) > maxabsy) {
        message(paste0("The maximum absolute value of y is larger than ", maxabsy, "! Please check whether data fits our assumptions!"))
        y[y>maxabsy]   = maxabsy
        y[y< -maxabsy] = -maxabsy
    }
    
    fit <- .Call("_ZC_Test", 
            as.numeric(t(X)), 
            as.numeric(y), 
            as.integer(param) )

    tmp = sqrt(2*fit$tr_sigma)/nrow(X)
    Tn  = fit$test / (fit$sigma2*tmp)
    
    result = list(
        Tn = Tn, 
        pvals = 1-pnorm(Tn)
    )
    return(result)
}


#' Zhong and Chen global test
#'
#' Implements the high-dimensional global test of Zhong and Chen (2011) for
#' uncensored data.
#'
#' @param X A numeric covariate matrix with rows corresponding to observations
#'   and columns corresponding to covariates.
#' @param y A numeric response vector.
#' @param maxabsy A numeric threshold used to truncate extremely large absolute
#'   response values. Default is `exp(10)`.
#'
#' @return A list containing:
#' \item{Tn}{The standardized test statistic.}
#' \item{pvals}{The one-sided p-value of the test.}
#'
#' @details
#' This function calls the C routine `_ZC_Test` to compute the test statistic and
#' variance components. If the response contains extremely large values, it is
#' truncated at `maxabsy` to avoid numerical instability.
#'
#' @references
#' Zhong, P.-S. and Chen, S. X. (2011). Tests for high-dimensional regression
#' coefficients with factorial designs. \emph{Journal of the American
#' Statistical Association}, 106(493), 260--274.
#'
#' @examples
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' res <- zctest(X, y)
#'
#' @importFrom stats pnorm
#' @useDynLib hdcrtRepro, .registration = TRUE
#' @export
cgztest = function(X, y, maxabsy = exp(10))
{
    param     = dim(X)
    if (max(abs(y)) > maxabsy) {
        message(paste0("The maximum absolute value of y is larger than ", maxabsy, "! Please check whether data fits our assumptions!"))
        y[y>maxabsy]   = maxabsy
        y[y< -maxabsy] = -maxabsy
    }
    
    fit <- .Call("_CGZ_Test", 
            as.numeric(t(X)), 
            as.numeric(y), 
            as.integer(param) )

    tmp = sqrt(2*fit$tr_sigma)/nrow(X)
    y0 = y-mean(y)
    X0 = t(t(X)-colMeans(X))
    sigma2 <- VAR_RCV(y0, X0)$sigma2
    
    Tn  = fit$test / (sigma2*tmp)
    
    result = list(
        Tn = Tn, 
        pvals = 1-pnorm(Tn)
    )
    return(result)
}


#' Feng, Zou, and Wang rank-based global test
#'
#' Implements the rank-based high-dimensional global test of Feng et al. (2013)
#' for uncensored data.
#'
#' @param X A numeric covariate matrix with rows corresponding to observations
#'   and columns corresponding to covariates.
#' @param y A numeric response vector used to rank the observations.
#'
#' @return A list containing:
#' \item{Tn}{The standardized test statistic.}
#' \item{pvals}{The one-sided p-value of the test.}
#'
#' @details
#' The observations are first ordered according to the response `y`. The ordered
#' covariate matrix is then passed to the C routine `_FZW_Test` to compute the
#' rank-based test statistic and variance component.
#'
#' @references
#' Feng, L., Zou, C., and Wang, Z. (2013). Rank-based score tests for
#' high-dimensional regression coefficients. \emph{Electronic Journal of
#' Statistics}, 7, 2131--2149.
#'
#' @examples
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' res <- fzwtest(X, y)
#'
#' @importFrom stats pnorm rnorm
#' @useDynLib hdcrtRepro, .registration = TRUE
#' @export
fzwtest = function(X, y)
{
    param   = dim(X)
    y.order = order(y)
    # y       = y[y.order]
    X       = X[y.order, ]
    
    fit <- .Call("_FZW_Test", 
            as.numeric(t(X)), 
            as.integer(param) )

    tmp = sqrt(2*fit$tr_sigma)/nrow(X)
    Tn  = fit$test / tmp
    
    result = list(
        Tn = Tn, 
        pvals = 1-pnorm(Tn)
    )
    return(result)
}

#' Choose censoring constants for double censoring
#'
#' Estimates the censoring constants used to achieve desired left and right
#' censoring rates in double-censored simulation studies.
#'
#' @param n An integer specifying the Monte Carlo sample size used to estimate
#'   the censoring constants.
#' @param beta A numeric vector of regression coefficients.
#' @param error A character string specifying the error distribution.
#' @param mu.x A numeric vector specifying the mean of the covariates. If missing,
#'   it is set to a zero vector.
#' @param Sigma.x A covariance matrix for the covariates. If missing, it is set
#'   to the identity matrix.
#' @param design_distribution A character string specifying the distribution used
#'   to generate the covariates. Options are `"norm"`, `"chi"`, and `"unif"`.
#' @param normalize_error A logical value indicating whether the generated errors
#'   should be standardized. Default is `TRUE`.
#' @param cenRate A numeric value specifying the total target censoring rate.
#'   The left and right censoring rates are each set to `cenRate / 2`.
#'   Default is `0.2`.
#' @param Gamma.x An optional matrix satisfying
#'   \eqn{\Gamma_x \Gamma_x^\top = \Sigma_x}, used to generate covariates from
#'   a factor model.
#'
#' @return A numeric vector of length two. The first element is the censoring
#'   constant for right censoring, and the second element is the censoring
#'   constant for left censoring.
#'
#' @details
#' This function generates a large auxiliary sample and estimates two constants
#' for double censoring. The censoring variable is generated from a chi-square
#' distribution with one degree of freedom. The constants are chosen from
#' empirical quantiles so that the simulated data approximately achieve the
#' desired left and right censoring rates.
#'
#' @examples
#' Sigma <- diag(5)
#' beta <- rep(1, 5)
#' c0 <- control_c_double_censoring(
#'   n = 1000, beta = beta, error = "norm",
#'   mu.x = rep(0, 5), Sigma.x = Sigma,
#'   cenRate = 0.2
#' )
#'
#' @importFrom stats rnorm runif rchisq quantile
#' @importFrom mnormt rmnorm
#' @export
control_c_double_censoring <- function(
    n, beta, error, mu.x, Sigma.x,
    design_distribution = c("norm", "chi", "unif"),
    normalize_error = TRUE, cenRate = 0.2, Gamma.x = NULL)
{
    # error distirbuion ========================================================
    e <- generate_error(n, error, normalize_error)

    p <- length(beta)
    if (missing(mu.x)) {
        mu.x = rep(0, p)
    }
    if (missing(Sigma.x)) {
        Sigma.x = diag(p)
    }

    # x from a factor model, x = Gamma.x %*% z
    design_distribution <- match.arg(design_distribution)
    if (design_distribution == "norm") {
        if (is.null(Gamma.x)) {
            x <- rmnorm(n, mean = mu.x, varcov = Sigma.x)
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            x <- matrix(rnorm(n * p), n, p) %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "chi") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled chi-square distribution,
            # which has zero mean and unit variance
            x <- (matrix(rnorm(n * p), n, p)^2 - 1) * sqrt(0.5)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "unif") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled uniform distribution,
            # which has zero mean and unit variance
            x <- (matrix(runif(n * p), n, p) - 0.5) * sqrt(12)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    }

    xb <- x %*% beta + e

    C <- rchisq(n, 1)
    c_cersor1 <- quantile(xb - C, 1 - cenRate/2)
    c_cersor2 <- quantile(xb - C, cenRate/2)

    return(c(c_cersor1, c_cersor2))
}

#' Generate double-censored transformation-model data
#'
#' Generates simulated data under a nonparametric transformation model with
#' double censoring.
#'
#' @param n An integer specifying the sample size.
#' @param beta A numeric vector of regression coefficients.
#' @param censor A character string specifying the censoring mechanism. Options
#'   are `"non"` for uncensored data, `"random"` for random double censoring,
#'   and `"right"` for fixed right censoring.
#' @param tran A character string specifying the transformation function.
#'   Options are `"lm"`, `"log"`, and `"boxcox"`.
#' @param error A character string specifying the error distribution.
#' @param mu.x A numeric vector specifying the mean of the covariates. If missing,
#'   it is set to a zero vector.
#' @param Sigma.x A covariance matrix for the covariates. If missing, it is set
#'   to the identity matrix.
#' @param design_distribution A character string specifying the distribution used
#'   to generate the covariates. Options are `"norm"`, `"chi"`, and `"unif"`.
#' @param normalize_error A logical value indicating whether the generated errors
#'   should be standardized. Default is `TRUE`.
#' @param c_cersor A numeric vector of length two containing the censoring
#'   constants for random double censoring. The first element is used for right
#'   censoring and the second for left censoring.
#' @param cenRate A numeric value specifying the target censoring rate for fixed
#'   right censoring. Default is `0.2`.
#' @param Gamma.x An optional matrix satisfying
#'   \eqn{\Gamma_x \Gamma_x^\top = \Sigma_x}, used to generate covariates from
#'   a factor model.
#' @param lambda_boxcox A numeric value specifying the Box-Cox transformation
#'   parameter. Default is `1`.
#'
#' @return A list containing:
#' \item{x}{The generated covariate matrix.}
#' \item{y}{The observed response vector after double censoring and
#' transformation.}
#' \item{beta}{The regression coefficient vector used to generate the data.}
#' \item{status_right}{The right-censoring indicator, where `TRUE` indicates
#' that the latent response is not right-censored.}
#' \item{status_left}{The left-censoring indicator, where `TRUE` indicates that
#' the latent response is not left-censored.}
#'
#' @details
#' The latent response is generated as
#' \deqn{x_i^\top \beta + \epsilon_i.}
#' For random double censoring, a chi-square random variable is generated and
#' shifted by two censoring constants to form the right- and left-censoring
#' thresholds. The observed latent response is
#' \deqn{\max\{\min(x_i^\top \beta + \epsilon_i, C_{ri}), C_{li}\},}
#' and the specified transformation is then applied to obtain the observed
#' response.
#'
#' @examples
#' Sigma <- diag(5)
#' beta <- rep(0, 5)
#' dat <- generator_double_censoring(
#'   n = 50, beta = beta, censor = "random", tran = "lm",
#'   error = "norm", mu.x = rep(0, 5), Sigma.x = Sigma,
#'   design_distribution = "norm", c_cersor = c(1, -1)
#' )
#'
#' @importFrom stats rnorm runif rchisq quantile
#' @importFrom mnormt rmnorm
#' @export
generator_double_censoring <- function(
    n, beta, censor = c("non", "random", "right"), tran = c("lm", "log", "boxcox"),
    error, mu.x, Sigma.x, design_distribution = c("norm", "chi", "unif"),
    normalize_error = TRUE, c_cersor = NULL, cenRate = 0.2, Gamma.x = NULL,
    lambda_boxcox = 1)
{
    # error distirbuion ========================================================
    e <- generate_error(n, error, normalize_error)

    p <- length(beta)
    if (missing(mu.x)) {
        mu.x = rep(0, p)
    }
    if (missing(Sigma.x)) {
        Sigma.x = diag(p)
    }

    # x from a factor model, x = Gamma.x %*% z
    design_distribution <- match.arg(design_distribution)
    if (design_distribution == "norm") {
        if (is.null(Gamma.x)) {
            x <- rmnorm(n, mean = mu.x, varcov = Sigma.x)
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            x <- matrix(rnorm(n * p), n, p) %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "chi") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled chi-square distribution,
            # which has zero mean and unit variance
            x <- (matrix(rnorm(n * p), n, p)^2 - 1) * sqrt(0.5)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "unif") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled uniform distribution,
            # which has zero mean and unit variance
            x <- (matrix(runif(n * p), n, p) - 0.5) * sqrt(12)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    }

    xb <- x %*% beta + e

    # transformation ===========================================================
    tran <- match.arg(tran)
    if (tran == "lm") g <- function(x) x
    if (tran == "log") g <- function(x) exp(x)
    if (tran == "boxcox") {
        if (lambda_boxcox == 0) {
            g <- function(x) exp(x)
        } else {
            g <- function(x) (lambda_boxcox*x + 1.0)^(1.0/lambda_boxcox)
        }
    }

    if (censor == "non") {
        y <- g(xb)
        status <- rep(1, n)
    } else if (censor == "random") {
        # random censor
        if (is.null(c_cersor)) {
            warning("c_cersor is null! We recommend using the 'control_c_cersor' function to calculate it to provides desired censor rate.")
            c_cersor = 0
        }
        C <- rchisq(n, 1)
        C1 <- C + c_cersor[1]
        C2 <- C + c_cersor[2]
        status_right <- xb < C1
        status_left <- xb > C2
        y <- g(pmax(pmin(xb, C1), C2))
    } else if (censor == "right") {
        # right censor
        C <- quantile(xb, 1 - cenRate)
        status <- xb < C
        y <- g(pmin(xb, C))
    }

    # In the transformation model, y may be infinity.
    y[which(y == "NaN" | y == "Inf")] <- max(y[which(y != "NaN" & y != "Inf")])

    list(x      = x,
         y      = y,
         beta   = beta,
         status_right = status_right,
         status_left = status_left)
}

#' Choose the censoring constant for partial-test simulations
#'
#' Estimates the censoring constant used to achieve a desired right-censoring
#' rate in simulation studies for the partial test with control factors.
#'
#' @param n An integer specifying the Monte Carlo sample size used to estimate
#'   the censoring constant.
#' @param beta A numeric vector of regression coefficients for all covariates.
#' @param q An integer specifying the number of control variables. This argument
#'   is included for consistency with the partial-test simulation setting.
#' @param error A character string specifying the error distribution.
#' @param mu.x A numeric vector specifying the mean of the covariates. If missing,
#'   it is set to a zero vector.
#' @param Sigma.x A covariance matrix for the covariates. If missing, it is set
#'   to the identity matrix.
#' @param design_distribution A character string specifying the distribution used
#'   to generate the covariates. Options are `"norm"`, `"chi"`, and `"unif"`.
#' @param normalize_error A logical value indicating whether the generated errors
#'   should be standardized. Default is `TRUE`.
#' @param cenRate A numeric value specifying the target censoring rate. Default
#'   is `0.2`.
#' @param Gamma.x An optional matrix satisfying
#'   \eqn{\Gamma_x \Gamma_x^\top = \Sigma_x}, used to generate covariates from
#'   a factor model.
#'
#' @return A numeric value giving the censoring constant.
#'
#' @details
#' This function generates a large auxiliary sample and estimates the constant
#' added to the right-censoring time so that the simulated data approximately
#' achieve the target censoring rate. The censoring time is generated from a
#' chi-square distribution with one degree of freedom.
#'
#' @examples
#' Sigma <- diag(10)
#' beta <- rep(1, 10)
#' c0 <- control_c_right_cersor_partial(
#'   n = 1000, beta = beta, q = 5, error = "norm",
#'   mu.x = rep(0, 10), Sigma.x = Sigma,
#'   cenRate = 0.2
#' )
#'
#' @importFrom stats rnorm runif rchisq quantile
#' @importFrom mnormt rmnorm
#' @export
control_c_right_cersor_partial <- function(
    n, beta, q, error, mu.x, Sigma.x,
    design_distribution = c("norm", "chi", "unif"),
    normalize_error = TRUE, cenRate = 0.2, Gamma.x = NULL)
{
    # error distirbuion ========================================================
    e <- generate_error(n, error, normalize_error)

    p <- length(beta)
    if (missing(mu.x)) {
        mu.x = rep(0, p)
    }
    if (missing(Sigma.x)) {
        Sigma.x = diag(p)
    }

    # x from a factor model, x = Gamma.x %*% z
    design_distribution <- match.arg(design_distribution)
    if (design_distribution == "norm") {
        if (is.null(Gamma.x)) {
            x <- rmnorm(n, mean = mu.x, varcov = Sigma.x)
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            x <- matrix(rnorm(n * p), n, p) %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "chi") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled chi-square distribution,
            # which has zero mean and unit variance
            x <- (matrix(rnorm(n * p), n, p)^2 - 1) * sqrt(0.5)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    } else if (design_distribution == "unif") {
        if (is.null(Gamma.x)) {
            stop("not support for now!")
        } else {
            # Gamma.x Gamma.x^T = Sigma.x
            # z from a scaled uniform distribution,
            # which has zero mean and unit variance
            x <- (matrix(runif(n * p), n, p) - 0.5) * sqrt(12)
            x <- x %*% t(Gamma.x) + rep(mu.x, each = n)
        }
    }

    xb <- x %*% beta + e

    C <- rchisq(n, 1)
    c_cersor <- quantile(xb - C, 1 - cenRate)

    return(c_cersor)
}

#' Estimate coefficients by smoothed partial rank regression
#'
#' Estimates regression coefficients using the smoothed partial rank procedure
#' and selects the model by a modified BIC criterion.
#'
#' @param y A numeric vector of observed survival times.
#' @param x A numeric covariate matrix with rows corresponding to observations
#'   and columns corresponding to covariates.
#' @param status A numeric or logical vector of censoring indicators, where
#'   nonzero or `TRUE` values indicate observed events.
#' @param gamma A numeric value controlling the penalty strength in the modified
#'   BIC criterion. Default is `0.2`.
#'
#' @return A numeric vector of estimated regression coefficients.
#'
#' @details
#' This function calls `hdcrt::sprfabs()` to compute a solution path for the
#' smoothed partial rank regression. The final model is selected by maximizing
#' a modified BIC criterion involving the event count, model degrees of freedom,
#' and the combinatorial penalty `lchoose(p, df)`.
#'
#' @examples
#' # betahat <- sprest(y, x, status)
#'
#' @importFrom hdcrt sprfabs
#' @export
sprest <- function(y, x, status, gamma = 0.2) {
    ns <- sum(status)
    n  <- length(y)
    p <- ncol(x)
    nmax <- as.integer(0.05*p) + 10

    res <- hdcrt::sprfabs(y, x, status, nmax = nmax, maxIter = p*10, message = FALSE)

    mbic <- 2*n*log(-res$loss) - res$df*log(ns) + 2*gamma*lchoose(p, res$df)
    opt <- which.max(mbic)
    betahat <- as.numeric(res$theta[,opt,drop=FALSE])

    return(betahat)
}

#' Yang, Guo, and Zhu score-based test
#'
#' Applies the score-function-based test of Yang, Guo, and Zhu (2024) by first
#' adjusting for control variables and then applying a high-dimensional global
#' test to the residuals.
#'
#' @param y A numeric response vector.
#' @param u A numeric matrix of control covariates.
#' @param x A numeric matrix of covariates of interest to be tested.
#'
#' @return A list returned by `zctest()`, containing the standardized test
#'   statistic and p-value.
#'
#' @details
#' This function is used as a benchmark method in the partial-test simulation
#' study. It fits a Gaussian lasso regression of `y` on the control covariates
#' `u` using `glmnet::cv.glmnet()`. The residuals are then used as the response
#' for testing the overall effect of the covariates of interest `x`.
#'
#' @references
#' Yang, W., Guo, X., and Zhu, L. (2024). Score function-based tests for
#' ultrahigh-dimensional linear models. \emph{Electronic Journal of Statistics},
#' 18(2), 4407--4458.
#'
#' @examples
#' # res <- ygztest(y, u, x)
#'
#' @importFrom glmnet cv.glmnet
#' @export
ygztest <- function(y, u, x) {
    fitglm  <- cv.glmnet(u, y, family = "gaussian", type.measure = "mse")
    betahat <- coef(fitglm)
    resids  <- y - (betahat[1] + u %*% betahat[-1])
    ygz     <- zctest(x, resids)

    return(ygz)
}

#' Estimate coefficients for the real data analysis
#'
#' Estimates regression coefficients using smoothed partial rank regression for
#' the real data analysis.
#'
#' @param y A numeric vector of observed survival times.
#' @param x A numeric covariate matrix with rows corresponding to observations
#'   and columns corresponding to covariates.
#' @param status A numeric or logical vector of censoring indicators, where
#'   nonzero or `TRUE` values indicate observed events.
#' @param gamma A numeric value controlling the penalty strength in the modified
#'   BIC criterion. If `NULL`, the last model on the solution path is used.
#'   Default is `0.4`.
#' @param nmax An integer specifying the maximum number of variables selected
#'   along the solution path. Default is `200`.
#'
#' @return A list containing:
#' \item{fit}{The fitted object returned by `hdcrt::sprfabs()`.}
#' \item{betahat}{The normalized estimated coefficient vector.}
#'
#' @details
#' This function is used in the TCGA SKCM real data analysis. It computes a
#' smoothed partial rank regression solution path using `hdcrt::sprfabs()`.
#' When `gamma` is not `NULL`, the final model is selected by maximizing a
#' modified BIC criterion involving the number of observed events, model degrees
#' of freedom, and the combinatorial penalty `lchoose(p, df)`. The selected
#' coefficient vector is normalized to have Euclidean norm one.
#'
#' @examples
#' # res <- sprest_realdata(y, x, status)
#'
#' @importFrom hdcrt sprfabs
#' @export
sprest_realdata <- function(y, x, status, gamma = 0.4, nmax = 200) {
    ns <- sum(status)
    n  <- length(y)
    p <- ncol(x)

    fit <- hdcrt::sprfabs(
        y, x, status, eps = 0.001, nmax = nmax,
        maxIter = 10000, message = FALSE
    )

    if (is.null(gamma)) {
        opt <- ncol(fit$theta)
    } else {
        mbic <- 2*n*log(-fit$loss) - fit$df*log(ns) + 2*gamma*lchoose(p, fit$df)
        opt <- which.max(mbic)
    }
    betahat <- as.numeric(fit$theta[,opt,drop=FALSE])
    betahat <- betahat / norm(betahat, "2")

    return(list(fit = fit, betahat = betahat))
}

#' Pathway-level partial test
#'
#' Applies the proposed partial test to evaluate the overall effect of genes in
#' a specified pathway while adjusting for clinical variables and genes outside
#' the pathway.
#'
#' @param pathway A character vector containing the gene names in the pathway of
#'   interest. These names should match column names of `x`.
#' @param y A numeric vector of observed survival times.
#' @param x A numeric matrix containing clinical variables and gene expression
#'   variables. Rows correspond to observations and columns correspond to
#'   variables.
#' @param status A numeric or logical vector of censoring indicators, where
#'   nonzero or `TRUE` values indicate observed events.
#' @param eps A numeric value specifying the convergence tolerance used in
#'   `hdcrt::sprfabs()`. Default is `0.001`.
#' @param gamma A numeric value controlling the penalty strength in the modified
#'   BIC criterion. Default is `0.4`.
#'
#' @return A list returned by `hdcrt::hdcrpt()`, containing the pathway-level
#'   partial test statistic and p-value.
#'
#' @details
#' Genes in `pathway` are treated as variables of interest. All remaining
#' variables in `x` are treated as control variables. The control effect is first
#' estimated by smoothed partial rank regression using `hdcrt::sprfabs()`, with
#' the model selected by a modified BIC criterion. The fitted linear predictor is
#' then used as the adjustment term in the proposed partial test with the sigmoid
#' kernel.
#'
#' This function is used in the TCGA SKCM pathway-level analysis.
#'
#' @examples
#' # res <- hdcrpt_pathway(
#' #   pathway = pathway_genes,
#' #   y = y,
#' #   x = cbind(Es, Xs),
#' #   status = status
#' # )
#'
#' @importFrom hdcrt sprfabs hdcrpt
#' @export
hdcrpt_pathway <- function(pathway, y, x, status, eps = 0.001, gamma = 0.4) {
    ns <- sum(status)
    n  <- length(y)
    pu <- ncol(x) - length(pathway)
    
    u <- x[, setdiff(colnames(x), pathway)]
    z <- x[, pathway, drop = FALSE]

    fit <- hdcrt::sprfabs(
        y, u, status, eps = eps, maxIter = 10000, message = FALSE
    )

    mbic <- 2*n*log(-fit$loss) - fit$df*log(ns) + 2*gamma*lchoose(pu, fit$df)
    opt <- which.max(mbic)
    betahat <- as.numeric(fit$theta[,opt,drop=FALSE])
    betahat <- betahat / norm(betahat, "2")

    uhat <- u %*% as.vector(betahat)

    cr_sigmoid <- hdcrt::hdcrpt(
        z, y, status, uhat, "sigmoid", covariate_dependence = TRUE
    )

    return(cr_sigmoid)
}





























































#' Random projection global test
#'
#' Implements a random projection-based global test for high-dimensional linear
#' regression models.
#'
#' @param x A numeric covariate matrix with rows corresponding to observations
#'   and columns corresponding to covariates.
#' @param y A numeric response vector.
#' @param rho A numeric value specifying the projection ratio. Default is `0.4`.
#' @param Pk An optional projection matrix. If `NULL`, a random Gaussian
#'   projection matrix is generated internally when projection is needed.
#'
#' @return A list containing:
#' \item{Tn}{The standardized test statistic.}
#' \item{pvals}{The p-value of the test.}
#' \item{rho}{The projection ratio used in the test.}
#'
#' @details
#' This function implements a random projection-based global test for testing
#' whether the regression coefficients of `x` are zero. If the number of
#' covariates is smaller than `rho * n`, the test is computed without random
#' projection. Otherwise, the covariates are projected to a lower-dimensional
#' space using either the user-provided projection matrix `Pk` or a randomly
#' generated Gaussian projection matrix.
#'
#' The core test statistic is computed by C routines registered in the package
#' and then standardized using the normal approximation.
#'
#' @examples
#' x <- matrix(rnorm(100 * 50), 100, 50)
#' y <- rnorm(100)
#' res <- pvalrp(x, y, rho = 0.4)
#' res$pvals
#'
#' @importFrom stats rnorm pnorm
#' @useDynLib hdcrtRepro, .registration = TRUE
#' @export
pvalrp <- function(x, y, rho = 0.4, Pk = NULL){
	n = nrow(x)
	p = ncol(x)
	if(is.null(n)) n = length(x)
	if(is.null(p)) p = 1


	if(p < rho*n){
		dims = c(n,p)
		Tn 	<- .Call("RPtest0",
					as.numeric(x),
					as.numeric(y),
					as.integer(dims),
                    PACKAGE = "hdcrtRepro"
				)
		rho = p/n
	}
	else{
		if(is.null(Pk)){
			k 	= ceiling(rho*n)
			Pk 	= matrix(rnorm((k*p),0,1), nrow = p, ncol = k)
		}
		else{
			k 	= ncol(Pk)
			rho = k/n
		}
		dims = c(n,p,k)
		Tn 	<- .Call("RPtest",
					as.numeric(x),
					as.numeric(y),
					as.numeric(Pk),
					as.integer(dims),
                    PACKAGE = "hdcrtRepro"
				)
	}
	Tk = (Tn-1) / sqrt(2/(n*rho*(1-rho)))

	pvalue = 1 - pnorm(Tk)
	return(list(Tn = Tk, pvals = pvalue, rho = rho))

}









#' Guo and Chen high-dimensional generalized linear model test
#'
#' Implements the high-dimensional global test of Guo and Chen (2016) for
#' generalized linear regression models.
#'
#' @param x A numeric covariate matrix with rows corresponding to observations
#'   and columns corresponding to covariates of interest.
#' @param y A numeric response vector.
#' @param z An optional numeric matrix of control covariates. If `NULL`, no
#'   control covariates are adjusted for. Default is `NULL`.
#' @param family A character string specifying the generalized linear model
#'   family. Available options are `"gaussian"`, `"binomial"`, and `"poisson"`.
#'   Default is `"gaussian"`.
#' @param resids An optional numeric vector of residuals. If provided, these
#'   residuals are used directly in the test. Default is `NULL`.
#' @param psi An optional numeric vector of weights used in the test statistic.
#'   If `NULL`, it is set to a vector of ones. Default is `NULL`.
#'
#' @return A list containing:
#' \item{Tn}{The standardized test statistic.}
#' \item{pvals}{The p-value of the test.}
#'
#' @details
#' This function tests the global null hypothesis that the coefficients of the
#' covariates of interest `x` are zero in a generalized linear regression model.
#' If `resids` is not provided, residuals are computed internally. When control
#' covariates `z` are provided, a generalized linear model of `y` on `z` is first
#' fitted using `glm()`, and response residuals are extracted. When `z = NULL`,
#' simple null residuals are used according to the specified family.
#'
#' The core test statistic is computed by a C routine registered in the package.
#'
#' @references
#' Guo, B. and Chen, S. X. (2016). Tests for high dimensional generalized linear
#' models. \emph{Journal of the Royal Statistical Society: Series B}, 78,
#' 1079--1102.
#'
#' Chen, J., Li, Q., and Chen, H. Y. (2023). Testing generalized linear models
#' with high-dimensional nuisance parameters. \emph{Biometrika}, 110, 83--99.
#'
#' @examples
#' x <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' res <- pvalgc(x, y, family = "gaussian")
#' res$pvals
#'
#' @importFrom stats glm residuals pnorm rnorm
#' @useDynLib hdcrtRepro, .registration = TRUE
#' @export
pvalgc <- function(x, y, z = NULL, family = "gaussian", resids = NULL, psi = NULL){
	# High-dimensional testing of coefficient in linear regressions.
	if(!(family %in% c('gaussian', 'binomial','poisson'))){
		stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
	}
	n 	= length(y)
	p 	= ifelse(is.null(ncol(x)), 1, ncol(x))

	if(is.null(resids)){
		if(is.null(z)){
			if(family=='gaussian'){
				resids  = y
			}
			else if(family == 'binomial'){
				resids  = y - 0.5
			}
			else if(family == 'poisson'){
				resids  = y - 1
			}
			else{
				stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
			}
		}
		else{
			z 		= z
			fitglm 	= glm(y~z, family = family)
			resids 	= residuals(fitglm, type = "response")
		}
	}
	if(is.null(psi)){
		psi     = rep(1, n)
		ispsi   = 0
	}
	else{
		ispsi = 1
	}



	dims 	= c(n, p, ispsi)
	Tn		= .Call("GCtest_",
					as.numeric(x),
					as.numeric(resids),
					as.numeric(psi),
					as.integer(dims)
			)
	pvals 	= pnorm(Tn, lower.tail = F)

	return(list(Tn = Tn, pvals = pvals))
}












#' Chen, Li, and Chen high-dimensional nuisance test
#'
#' Implements a high-dimensional global test for generalized linear models in
#' the presence of high-dimensional control variables.
#'
#' @param x A numeric covariate matrix with rows corresponding to observations
#'   and columns corresponding to covariates of interest.
#' @param y A numeric response vector.
#' @param z An optional numeric matrix of high-dimensional control covariates.
#'   If `NULL`, no control covariates are adjusted for. Default is `NULL`.
#' @param family A character string specifying the generalized linear model
#'   family. Available options are `"gaussian"`, `"binomial"`, and `"poisson"`.
#'   Default is `"gaussian"`.
#' @param resids An optional numeric vector of residuals. If provided, these
#'   residuals are used directly in the test. Default is `NULL`.
#' @param psi An optional numeric vector of weights used in the test statistic.
#'   If `NULL`, it is set to a vector of ones. Default is `NULL`.
#'
#' @return A list containing:
#' \item{Tn}{The absolute value of the standardized test statistic.}
#' \item{pvals}{The two-sided p-value of the test.}
#'
#' @details
#' This function tests the global null hypothesis that the coefficients of the
#' covariates of interest `x` are zero while adjusting for possibly
#' high-dimensional control covariates `z`.
#'
#' If `resids` is not provided, residuals are computed internally. 
#' The nuisance effect of `z` is estimated by cross-validated lasso 
#' using `glmnet::cv.glmnet()`. If `z = NULL`, null residuals
#' are constructed according to the specified family.
#'
#' The core test statistic is computed by a C routine registered in the package.
#'
#' @references
#' Guo, B. and Chen, S. X. (2016). Tests for high dimensional generalized linear
#' models. \emph{Journal of the Royal Statistical Society: Series B}, 78,
#' 1079--1102.
#'
#' Chen, J., Li, Q., and Chen, H. Y. (2023). Testing generalized linear models
#' with high-dimensional nuisance parameters. \emph{Biometrika}, 110, 83--99.
#'
#' @examples
#' x <- matrix(rnorm(100 * 20), 100, 20)
#' z <- matrix(rnorm(100 * 30), 100, 30)
#' y <- rnorm(100)
#' res <- pvalclc(x, y, z = z, family = "gaussian")
#' res$pvals
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats pnorm rnorm
#' @useDynLib hdcrtRepro, .registration = TRUE
#' @export
pvalclc <- function(x, y, z=NULL, family = "gaussian", resids = NULL, psi = NULL){
	# High-dimensional testing of coefficient in linear regressions in presence of high-dimensional control factors.
	if(!(family %in% c('gaussian', 'binomial','poisson'))){
		stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
	}
	n 	= length(y)
	p 	= ifelse(is.null(ncol(x)), 1, ncol(x))

	if(is.null(resids)){
		if(is.null(z)){
			if(family=='gaussian'){
				resids  = y
			}
			else if(family == 'binomial'){
				resids  = y - 0.5
			}
			else if(family == 'poisson'){
				resids  = y - 1
			}
			else{
				stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
			}
		}
		else {
			z 		= z
			fitglm 	= cv.glmnet(z, y, family = family, type.measure="mse")
			betahat = coef(fitglm)
			mu 		= betahat[1] + z %*% betahat[-1]

			if(family=='gaussian'){
				resids	<- y - mu
			}
			else if(family == 'binomial'){
				resids	<- y - 1/(1+exp(-mu))
			}
			else if(family == 'poisson'){
				resids	<- y - exp(mu)
			}
			else{
				stop("family must be one of {'gaussian', 'binomial', 'poisson'} !")
			}
		}
	}
	if(is.null(psi)){
		psi     = rep(1, n)
		ispsi   = 0
	}
	else{
		ispsi = 1
	}


	dims 	= c(n, p, ispsi)
	Tn		= .Call("GCtest_",
				as.numeric(x),
				as.numeric(resids),
				as.numeric(psi),
				as.integer(dims)
			)
	pvals 	= pnorm(abs(Tn), lower.tail = F)

	return(list(Tn = abs(Tn), pvals = pvals))
}



























#' Refitted cross-validation variance estimation
#'
#' Estimates the residual variance in a high-dimensional linear model using
#' refitted cross-validation.
#'
#' @param y A numeric response vector of length `n`.
#' @param x A numeric design matrix with rows corresponding to observations and
#'   columns corresponding to covariates.
#'
#' @return A list containing:
#' \item{sigma2}{The estimated residual variance.}
#'
#' @details
#' This function implements the refitted cross-validation variance estimator of
#' Fan et al. (2012). The data are split into two parts. A lasso model is fitted
#' on one part, and the selected variables are refitted on the other part to
#' estimate the residual variance. The procedure is then repeated in the reverse
#' direction, and the two variance estimates are averaged.
#'
#' @references
#' Fan, J., Guo, S., and Hao, N. (2012). Variance estimation using refitted
#' cross-validation in ultrahigh dimensional regression. \emph{Journal of the
#' Royal Statistical Society: Series B}, 74, 37--65.
#'
#' @examples
#' n <- 80
#' p <- 100
#' beta <- c(sqrt(0.1 / p) * rep(1, p / 2), rep(0, p / 2))
#' eps <- rnorm(n)
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x %*% beta + eps
#' fit <- VAR_RCV(y, x)
#' fit$sigma2
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats rnorm
#' @export
VAR_RCV <- function(y,x){
  if(is.null(x)) stop("x must not be NA")
  if(is.null(y)) stop("y must not be NA")
  n = nrow(x)
  p = ncol(x)
  if(is.null(p)) p = 1
  half = ceiling(n/2)
  x1 = x[1:half,]
  y1 = y[1:half]
  x2 = x[-c(1:half),]
  y2 = y[-c(1:half)]

  fit.cv = cv.glmnet(x1,y1,family="gaussian")
  ind = which.min(fit.cv$cvm)
  fits = fit.cv$glmnet.fit
  if(ind==1){
    sigmahat1 = sum(y2^2)/half
  } else{
    ind1 = which(abs(fits$beta[,ind])>0)
    if(length(ind1)>half/2){
      betasort = sort(abs(fits$beta[ind1,ind]),decreasing =TRUE,index.return=T)
      ind1 = ind1[betasort$ix[1:ceiling(half/2)]]
    }
    xm12 = x2[,ind1]
    pm12 = xm12%*%solve(t(xm12)%*%xm12)%*%t(xm12)
    sigmahat1 = (sum(y2^2)-t(y2)%*%pm12%*%y2)/(half-length(ind1))
  }

  fit.cv = cv.glmnet(x2,y2,family="gaussian")
  ind = which.min(fit.cv$cvm)
  fits = fit.cv$glmnet.fit
  if(ind==1){
    sigmahat2 = sum(y1^2)/(n-half)
  } else{
    ind1 = which(abs(fits$beta[,ind])>0)
    if(length(ind1)>half/2){
      betasort = sort(abs(fits$beta[ind1,ind]),decreasing=TRUE,index.return=T)
      ind1 = ind1[betasort$ix[1:ceiling(half/2)]]
    }
    xm21 = x1[,ind1]
    pm21 = xm21%*%solve(t(xm21)%*%xm21)%*%t(xm21)
    sigmahat2 = (sum(y1^2)-t(y1)%*%pm21%*%y1)/(n-half-length(ind1))
  }
  sigmahat = (sigmahat1+sigmahat2)/2
  return(list(sigma2=sigmahat))
}






