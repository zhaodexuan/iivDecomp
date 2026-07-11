#' Fit Multilevel Model for IIV Estimation
#'
#' Fits a two-level random-intercept model to person-time aggregated data
#' as the first stage of IIV decomposition.
#'
#' @param data A data frame from \code{\link{iiv_generate}} or similar structure.
#' @param ycol Column name for the response variable (default "y").
#'
#' @return A list with components: s2_b (between-person variance),
#'   s2_w (within-person residual variance), fit (lmer model object),
#'   d (aggregated data frame).
#'
#' @examples
#' dat <- iiv_generate(N = 50, T = 30, K = 4, sigma_tau = 0.3, phi = 0.5, seed = 123)
#' fit <- iiv_fit_lmer(dat)
#' fit$s2_w
#'
#' @export
iiv_fit_lmer <- function(data, ycol = "y") {
  fml <- as.formula(paste(ycol, "~ id + time"))
  da_all <- data[!is.na(data[[ycol]]), c("id", "time", ycol)]
  if (nrow(da_all) < 10) stop("Not enough non-missing observations for model fitting")

  da <- aggregate(fml, da_all, mean)
  da <- da[order(da$id, da$time), ]

  fit <- tryCatch(
    lmer(as.formula(paste(ycol, "~ 1 + (1|id)")), data = da, REML = TRUE),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    tv <- var(da[[ycol]], na.rm = TRUE)
    bv <- var(tapply(da[[ycol]], da$id, mean, na.rm = TRUE), na.rm = TRUE)
    wv <- max(0.01, tv - bv)
    return(list(s2_b = bv, s2_w = wv, fit = NULL, d = da))
  }

  vc <- as.data.frame(VarCorr(fit))
  s2_b <- vc$vcov[1]
  s2_w <- attr(VarCorr(fit), "sc")^2
  if (is.na(s2_w) || s2_w <= 0) s2_w <- var(residuals(fit))

  list(s2_b = s2_b, s2_w = s2_w, fit = fit, d = da)
}
