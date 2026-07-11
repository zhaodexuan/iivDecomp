#' Detect IER via Longstring Analysis
#'
#' Flags person-time combinations where responses show excessive repetition,
#' indicating insufficient effort responding (IER).
#'
#' @param data A data frame from \code{\link{iiv_generate}}.
#' @param R Number of Likert scale points (default 5).
#'
#' @return A data frame with columns id, time, y (max count),
#'   and flag (logical: suspected IER).
#'
#' @examples
#' dat <- iiv_generate(N = 50, T = 30, K = 4, sigma_tau = 0.3,
#'   phi = 0.5, p_IER = 0.15, seed = 123)
#' ls <- iiv_detect_longstring(dat)
#' table(ls$flag)
#'
#' @export
iiv_detect_longstring <- function(data, R = 5) {
  ls <- aggregate(y ~ id + time, data, function(x) {
    if (all(is.na(x))) return(0)
    max(table(x, useNA = "no"))
  })
  ls$flag <- ls$y >= max(2, R - 1)
  ls
}

#' Detect IER via Mahalanobis Distance
#'
#' Flags person-time combinations that are multivariate outliers in item
#' response space, indicating random careless responding.
#'
#' @param data A data frame from \code{\link{iiv_generate}}.
#' @param alpha Significance threshold for outlier detection (default 0.001).
#'
#' @return A data frame with columns id, time, md (Mahalanobis distance),
#'   and flag (logical: suspected IER).
#'
#' @examples
#' dat <- iiv_generate(N = 50, T = 30, K = 4, sigma_tau = 0.3,
#'   phi = 0.5, p_IER = 0.15, IER_type = "random", seed = 123)
#' md <- iiv_detect_mahalanobis(dat)
#' table(md$flag)
#'
#' @export
iiv_detect_mahalanobis <- function(data, alpha = 0.001) {
  d_complete <- data[!is.na(data$y), ]
  pw <- reshape(d_complete[, c("id", "time", "item", "y")],
                v.names = "y", timevar = "item",
                idvar = c("id", "time"), direction = "wide")
  icols <- grep("^y\\.", names(pw))

  if (length(icols) < 2) {
    return(data.frame(id = integer(), time = integer(),
                      md = numeric(), flag = logical()))
  }

  X <- as.matrix(pw[, icols, drop = FALSE])
  X <- X[complete.cases(X), , drop = FALSE]
  if (nrow(X) < 3) {
    return(data.frame(id = integer(), time = integer(),
                      md = numeric(), flag = logical()))
  }

  mu <- colMeans(X, na.rm = TRUE)
  S <- cov(X, use = "pairwise.complete.obs")
  md <- tryCatch(
    mahalanobis(X, mu, S),
    error = function(e) mahalanobis(X, mu, diag(diag(S)))
  )

  K <- length(icols)
  threshold <- qchisq(1 - alpha, df = K)
  flagged <- pw[md > threshold, c("id", "time")]
  flagged$md <- md[md > threshold]
  flagged$flag <- TRUE

  # Add non-flagged
  clean <- pw[md <= threshold, c("id", "time")]
  clean$md <- md[md <= threshold]
  clean$flag <- FALSE

  rbind(flagged, clean)
}

#' Classical Aggregation Method (M1)
#'
#' Estimates IIV by treating all within-person variance as noise
#' (classical test theory approach).
#'
#' @param data A data frame from \code{\link{iiv_generate}}.
#'
#' @return A named numeric vector with s2_total, s2_within, and ICC.
#'
#' @examples
#' dat <- iiv_generate(N = 50, T = 30, K = 4, sigma_tau = 0.3, phi = 0.5, seed = 123)
#' iiv_aggregate(dat)
#'
#' @export
iiv_aggregate <- function(data) {
  d <- data[!is.na(data$y), ]
  pm <- aggregate(y ~ id, d, mean)
  tv <- var(d$y)
  bv <- var(pm$y)
  c(s2_total = tv, s2_within = tv - bv, ICC = bv / tv)
}
