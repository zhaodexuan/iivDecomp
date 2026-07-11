#' Decompose IIV Using EM Mixture Model (M4)
#'
#' Core function implementing the proposed variance decomposition framework.
#' Separates observed intraindividual variability into true construct-related
#' variability, measurement error, IER contamination, and context effects
#' using an expectation-maximization (EM) mixture model on multilevel
#' model residuals.
#'
#' @param data A data frame from \code{\link{iiv_generate}} with columns
#'   id, time, item, y, ys (continuous pre-discretization values).
#' @param R Number of Likert scale points (default 5). Used if ys column
#'   is unavailable.
#' @param max_iter Maximum EM iterations (default 100).
#' @param tol Convergence tolerance for EM (default 1e-5).
#' @param use_ys Use continuous-scale ys column for sigma2_eps estimation
#'   (default TRUE). If FALSE, uses Likert correction formula.
#'
#' @return A list of class "iivDecomp" with components:
#'   \item{s2_tau}{Estimated true IIV variance}
#'   \item{s2_eps}{Estimated measurement error variance}
#'   \item{s2_ier}{Estimated IER-induced variance}
#'   \item{s2_gamma}{Estimated context effect variance}
#'   \item{s2_between}{Estimated between-person variance}
#'   \item{s2_within}{Total within-person residual variance}
#'   \item{p_ier}{Estimated proportion of IER-contaminated observations}
#'   \item{em_iter}{Number of EM iterations used}
#'   \item{em_converged}{Logical: did EM converge?}
#'   \item{em_loglik}{Final EM log-likelihood}
#'   \item{components}{Data frame of component estimates with proportions}
#'
#' @examples
#' dat <- iiv_generate(N = 200, T = 30, K = 8, sigma_tau = 0.3,
#'   phi = 0.5, p_IER = 0.15, seed = 123)
#' res <- iiv_decompose(dat)
#' print(res)
#' summary(res)
#'
#' @export
iiv_decompose <- function(data, R = 5, max_iter = 100, tol = 1e-5,
                          use_ys = TRUE) {

  K <- length(unique(data$item))
  if (K < 3) warning("K >= 3 is recommended for identifiability")

  # Stage 1: Multilevel model
  fit_res <- iiv_fit_lmer(data)
  s2_w <- fit_res$s2_w
  s2_b <- fit_res$s2_b
  da <- fit_res$d

  # Stage 2: EM mixture on residuals
  if (!is.null(fit_res$fit)) {
    resid_pt <- residuals(fit_res$fit, type = "response")
  } else {
    resid_pt <- da$y - rep(tapply(da$y, da$id, mean, na.rm = TRUE),
                            table(da$id))
  }

  n <- length(resid_pt)
  madv <- mad(resid_pt, constant = 1, na.rm = TRUE)
  narrow <- abs(resid_pt) < 2 * madv

  pi1 <- .clamp(mean(narrow, na.rm = TRUE), 0.1, 0.9)
  pi2 <- 1 - pi1
  s1 <- var(resid_pt[narrow], na.rm = TRUE)
  if (is.na(s1) || s1 <= 0) s1 <- var(resid_pt, na.rm = TRUE) * 0.5
  s2 <- var(resid_pt[!narrow], na.rm = TRUE)
  if (is.na(s2) || s2 <= s1) s2 <- s1 * 3

  em_iter <- 0
  em_converged <- FALSE
  em_loglik <- -Inf

  for (it in 1:max_iter) {
    d1 <- dnorm(resid_pt, 0, sqrt(s1)) * pi1
    d2 <- dnorm(resid_pt, 0, sqrt(s2)) * pi2
    denom <- d1 + d2
    denom[denom < 1e-15] <- 1e-15
    g1 <- d1 / denom
    g2 <- 1 - g1

    pi1_new <- .clamp(mean(g1, na.rm = TRUE), 0.05, 0.95)
    pi2_new <- 1 - pi1_new
    s1_new <- max(0.001, sum(g1 * resid_pt^2, na.rm = TRUE) /
                    max(sum(g1, na.rm = TRUE), 1))
    s2_new <- max(s1_new * 1.5, sum(g2 * resid_pt^2, na.rm = TRUE) /
                    max(sum(g2, na.rm = TRUE), 1))

    ll <- sum(log(d1 + d2), na.rm = TRUE)

    if (abs(pi1_new - pi1) < tol &&
        abs(s1_new - s1) / max(s1, 0.01) < tol) {
      em_converged <- TRUE
      em_iter <- it
      em_loglik <- ll
      break
    }
    if (!is.na(ll) && !is.na(em_loglik) && abs(ll - em_loglik) < 1e-6) {
      em_converged <- TRUE
      em_iter <- it
      em_loglik <- ll
      break
    }

    pi1 <- pi1_new; pi2 <- pi2_new
    s1 <- s1_new; s2 <- s2_new
    em_loglik <- ll
    em_iter <- it
  }
  if (em_iter >= max_iter - 1) em_converged <- FALSE

  p_ier <- pi2

  # Stage 3: sigma2_eps estimation
  if (use_ys && "ys" %in% names(data)) {
    data$pt <- paste(data$id, data$time, sep = "_")
    s2_eps <- as.numeric(median(tapply(data$ys, data$pt, var, na.rm = TRUE),
                                na.rm = TRUE))
  } else {
    # Likert correction formula
    data$pt <- paste(data$id, data$time, sep = "_")
    likert_var <- median(tapply(data$y, data$pt, var, na.rm = TRUE),
                         na.rm = TRUE)
    s2_eps <- max(0.01, likert_var * (R^2 - 1) / 12)
  }
  if (is.na(s2_eps) || s2_eps <= 0) s2_eps <- 0.01

  # Stage 4: Decompose
  s2_tau <- max(0.01, s1 - s2_eps / K)
  s2_ier <- max(0, s2 - s1 - s2_eps / K)
  s2_gamma <- max(0, var(tapply(da$y, da$time, mean, na.rm = TRUE),
                         na.rm = TRUE))

  # Component summary
  s2_total <- s2_tau + s2_eps + s2_ier + s2_gamma
  comp <- data.frame(
    Component = c("True IIV", "Measurement Error", "IER", "Context"),
    Symbol = c("sigma2_tau", "sigma2_eps", "sigma2_ier", "sigma2_gamma"),
    Variance = c(s2_tau, s2_eps, s2_ier, s2_gamma),
    Proportion = c(s2_tau, s2_eps, s2_ier, s2_gamma) / s2_total,
    stringsAsFactors = FALSE
  )

  result <- list(
    s2_tau = s2_tau, s2_eps = s2_eps, s2_ier = s2_ier,
    s2_gamma = s2_gamma, s2_between = s2_b, s2_within = s2_w,
    p_ier = p_ier,
    em_iter = em_iter, em_converged = em_converged, em_loglik = em_loglik,
    components = comp,
    fit = fit_res$fit
  )
  class(result) <- "iivDecomp"
  result
}

#' @export
print.iivDecomp <- function(x, ...) {
  cat("IIV Variance Decomposition Results\n")
  cat("==================================\n")
  cat(sprintf("True IIV (sigma2_tau):        %.4f\n", x$s2_tau))
  cat(sprintf("Measurement Error (sigma2_eps): %.4f\n", x$s2_eps))
  cat(sprintf("IER Contamination (sigma2_ier): %.4f\n", x$s2_ier))
  cat(sprintf("Context Effects (sigma2_gamma): %.4f\n", x$s2_gamma))
  cat(sprintf("Between-Person (sigma2_b):     %.4f\n", x$s2_between))
  cat(sprintf("\nEstimated IER proportion: %.1f%%\n", x$p_ier * 100))
  cat(sprintf("EM iterations: %d (%s)\n", x$em_iter,
      ifelse(x$em_converged, "converged", "not converged")))
  cat("\nComponent proportions:\n")
  print(x$components, row.names = FALSE)
}

#' @export
summary.iivDecomp <- function(object, ...) {
  object$components
}

#' Compare Multiple IIV Estimation Methods
#'
#' Runs M1 (aggregation), M2 (multilevel), M3a (longstring screening),
#' M3b (Mahalanobis screening), and M4 (EM decomposition) on the same
#' dataset and returns a comparison table.
#'
#' @param data A data frame from \code{\link{iiv_generate}}.
#' @param R Number of Likert scale points (default 5).
#'
#' @return A data frame comparing IIV estimates across methods.
#'
#' @examples
#' dat <- iiv_generate(N = 200, T = 30, K = 8, sigma_tau = 0.3,
#'   phi = 0.5, p_IER = 0.15, seed = 123)
#' iiv_compare(dat)
#'
#' @export
iiv_compare <- function(data, R = 5) {
  K <- length(unique(data$item))

  # M1: Aggregation
  m1 <- iiv_aggregate(data)

  # M2: Multilevel
  m2_res <- iiv_fit_lmer(data)
  m2 <- m2_res$s2_w

  # M3a: Longstring + refit
  ls <- iiv_detect_longstring(data, R)
  if (any(ls$flag)) {
    bad <- ls[ls$flag, c("id", "time")]
    keep <- rep(TRUE, nrow(data))
    for (r in seq_len(nrow(bad))) {
      keep[data$id == bad$id[r] & data$time == bad$time[r]] <- FALSE
    }
    m3a_res <- iiv_fit_lmer(data[keep, ])
  } else {
    m3a_res <- m2_res
  }

  # M3b: Mahalanobis + refit
  md <- iiv_detect_mahalanobis(data)
  if (any(md$flag)) {
    bad <- md[md$flag, c("id", "time")]
    keep <- rep(TRUE, nrow(data))
    for (r in seq_len(nrow(bad))) {
      keep[data$id == bad$id[r] & data$time == bad$time[r]] <- FALSE
    }
    m3b_res <- iiv_fit_lmer(data[keep, ])
  } else {
    m3b_res <- m2_res
  }

  # M4: Decomposition
  m4_res <- iiv_decompose(data, R)

  data.frame(
    Method = c("M1 Aggregation", "M2 Multilevel", "M3a Longstring",
               "M3b Mahalanobis", "M4 Decomposition"),
    sigma2_IIV = c(m1["s2_within"], m2, m3a_res$s2_w, m3b_res$s2_w,
                   m4_res$s2_tau),
    Notes = c("All IIV = noise", "All IIV = signal",
              "Longstring pre-screen", "Mahalanobis pre-screen",
              "EM decomposed"),
    stringsAsFactors = FALSE
  )
}
