#' Fit Multilevel Model for IIV Estimation
#'
#' Fits a two-level random-intercept model to person-time aggregated data
#' as the first stage of IIV decomposition.
#'
#' @param data A data frame from \code{\link{iiv_generate}} or similar structure.
#' @param ycol Column name for the response variable (default "y").
#'
#' @details Requires the \pkg{lme4} package (declared in \code{Imports}). The
#'   model is fitted with \code{lme4::lmer()} so that no explicit
#'   \code{library(lme4)} call is needed by the user. For single-cluster
#'   input (e.g., decompositions run separately per participant), the
#'   random-intercept model is not separately identified; the function
#'   deliberately uses a variance-decomposition fallback
#'   (\code{var(y) - var(person means)}) in that case. For multi-cluster
#'   input, if \code{lmer()} fails to fit, a warning is issued and the
#'   fallback is used; the returned \code{method} element records which
#'   estimator was applied ("lmer" or "fallback").
#'
#' @return A list with components: s2_b (between-person variance),
#'   s2_w (within-person residual variance), fit (lmer model object or NULL),
#'   d (aggregated data frame), method ("lmer" or "fallback").
#'
#' @examples
#' dat <- iiv_generate(N = 50, T = 30, K = 4, sigma_tau = 0.3, phi = 0.5, seed = 123)
#' fit <- iiv_fit_lmer(dat)
#' fit$s2_w
#'
#' @export
iiv_fit_lmer <- function(data, ycol = "y") {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required by iivDecomp for multilevel model fitting. ",
         "Install it with install.packages('lme4').")
  }

  fml <- as.formula(paste(ycol, "~ id + time"))
  da_all <- data[!is.na(data[[ycol]]), c("id", "time", ycol)]
  if (nrow(da_all) < 10) stop("Not enough non-missing observations for model fitting")

  da <- aggregate(fml, da_all, mean)
  da <- da[order(da$id, da$time), ]

  fallback <- function(warn = FALSE) {
    if (warn) {
      warning("lme4::lmer() failed to fit; using variance-decomposition fallback estimator.")
    }
    bv <- tryCatch(var(tapply(da[[ycol]], da$id, mean, na.rm = TRUE), na.rm = TRUE),
                   error = function(e) NA_real_)
    tv <- var(da[[ycol]], na.rm = TRUE)
    wv <- if (is.na(bv)) max(0.01, tv) else max(0.01, tv - bv)
    list(s2_b = bv, s2_w = wv, fit = NULL, d = da, method = "fallback")
  }

  # Single cluster: random-intercept model is not separately identified.
  # Deliberate fallback (operative path for person-by-person decompositions).
  if (length(unique(da$id)) < 2) return(fallback(warn = FALSE))

  fit <- tryCatch(
    lme4::lmer(as.formula(paste(ycol, "~ 1 + (1|id)")), data = da, REML = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(fallback(warn = TRUE))

  vc <- as.data.frame(lme4::VarCorr(fit))
  s2_b <- vc$vcov[1]
  s2_w <- attr(lme4::VarCorr(fit), "sc")^2
  if (is.na(s2_w) || s2_w <= 0) s2_w <- var(residuals(fit))

  list(s2_b = s2_b, s2_w = s2_w, fit = fit, d = da, method = "lmer")
}
