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
#' @return A data frame with one row per complete person-time (after reshaping
#'   to wide format): columns id, time, md (Mahalanobis distance), and flag
#'   (logical: suspected IER). Returns a zero-row data frame when fewer than
#'   three complete person-times are available, and returns all-FALSE flags
#'   when no observation exceeds the threshold (i.e., the empty-flag case is
#'   handled without error).
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

  empty_result <- data.frame(id = integer(), time = integer(),
                             md = numeric(), flag = logical())

  if (length(icols) < 2) {
    return(empty_result)
  }

  X <- as.matrix(pw[, icols, drop = FALSE])
  keep <- complete.cases(X)
  if (sum(keep) < 3) {
    return(empty_result)
  }

  # Keep distances aligned with rows: compute on complete cases only,
  # then restrict the person-time table to the same rows.
  X <- X[keep, , drop = FALSE]
  pw <- pw[keep, , drop = FALSE]

  mu <- colMeans(X, na.rm = TRUE)
  S <- cov(X, use = "pairwise.complete.obs")
  md <- tryCatch(
    mahalanobis(X, mu, S),
    error = function(e) mahalanobis(X, mu, diag(diag(S)))
  )

  threshold <- qchisq(1 - alpha, df = length(icols))
  data.frame(id = pw$id, time = pw$time, md = md, flag = md > threshold)
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
