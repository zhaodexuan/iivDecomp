#' Generate Simulated Intensive Longitudinal Data
#'
#' Creates simulated data with known true parameters for testing the
#' IIV decomposition framework. Supports AR(1) true IIV, measurement error,
#' IER contamination, context effects, and missing data.
#'
#' @param N Number of participants.
#' @param T Number of time points per participant.
#' @param K Number of items per time point.
#' @param sigma_tau True IIV variance (marginal).
#' @param phi Autocorrelation parameter (0 < phi < 1).
#' @param sigma_eps Measurement error standard deviation (default 0.3).
#' @param sigma_u Between-person standard deviation (default 1.0).
#' @param sigma_gamma Context effect standard deviation (default 0).
#' @param p_context Proportion of time points with context shifts (default 0.20).
#' @param p_IER Proportion of person-times contaminated by IER (default 0).
#' @param IER_type Type of IER: "random", "fixed", or "invariant" (default "random").
#' @param R Number of Likert scale points (default 5).
#' @param burn_in Number of burn-in time points for AR(1) (default 20).
#' @param seed Random seed for reproducibility (optional).
#' @param missing_type Missing data mechanism: "none", "MCAR", or "MAR" (default "none").
#' @param missing_rate Proportion of missing observations (default 0).
#'
#' @return A data frame with columns: id, time, item, y (observed Likert response),
#'   ys (continuous pre-discretization response), tau_true (true IIV value),
#'   ier (logical: IER-contaminated), missing (logical: missing).
#'
#' @examples
#' dat <- iiv_generate(N = 50, T = 30, K = 4, sigma_tau = 0.3, phi = 0.5, seed = 123)
#' head(dat)
#'
#' @export
iiv_generate <- function(N, T, K, sigma_tau, phi,
                         sigma_eps = 0.3, sigma_u = 1.0,
                         sigma_gamma = 0, p_context = 0.20,
                         p_IER = 0, IER_type = c("random", "fixed", "invariant"),
                         R = 5, burn_in = 20, seed = NULL,
                         missing_type = c("none", "MCAR", "MAR"),
                         missing_rate = 0) {

  IER_type <- match.arg(IER_type)
  missing_type <- match.arg(missing_type)
  if (!is.null(seed)) set.seed(seed)

  mu_i <- rnorm(N, 0, sigma_u)
  nttl <- T + burn_in
  innov_sd <- sqrt(sigma_tau * (1 - phi^2))

  # Vectorized AR(1): loop over time only
  innov <- matrix(rnorm(N * nttl, 0, innov_sd), nrow = N, ncol = nttl)
  tau <- matrix(0, nrow = N, ncol = nttl)
  tau[, 1] <- rnorm(N, 0, sqrt(sigma_tau))
  for (j in 2:nttl) tau[, j] <- phi * tau[, j - 1] + innov[, j]
  tau <- tau[, (burn_in + 1):nttl, drop = FALSE]

  nr <- N * T * K
  id <- rep(1:N, each = T * K)
  tm <- rep(rep(1:T, each = K), N)
  it <- rep(1:K, N * T)

  tl <- as.vector(t(tau))[rep(1:(N * T), each = K)]
  ml <- rep(mu_i, each = T * K)
  ys <- ml + tl + rnorm(nr, 0, sigma_eps)

  # Context effects
  if (sigma_gamma > 0) {
    nb <- max(1, round(p_context * T / 3))
    for (b in 1:nb) {
      ts <- sample(1:max(1, T - 2), 1)
      te <- min(ts + sample(2:4, 1), T)
      ys[tm >= ts & tm <= te] <- ys[tm >= ts & tm <= te] + rnorm(1, 0, sigma_gamma)
    }
  }

  y <- pmax(1, pmin(R, round(ys + (R + 1) / 2)))  # center latent on the response scale before discretization
  ier_mask <- rep(FALSE, nr)

  # IER injection
  if (p_IER > 0) {
    pt <- expand.grid(id = 1:N, time = 1:T)
    ni <- round(p_IER * nrow(pt))
    for (idx in sample(1:nrow(pt), ni)) {
      rows <- which(id == pt$id[idx] & tm == pt$time[idx])
      if (IER_type == "random") {
        y[rows] <- sample(1:R, length(rows), replace = TRUE)
      } else if (IER_type == "fixed") {
        y[rows] <- sample(c(1, R), 1)
      } else {
        y[rows] <- ceiling(R / 2)
      }
      ier_mask[rows] <- TRUE
    }
  }

  # Missing data
  missing_mask <- rep(FALSE, nr)
  if (missing_type != "none" && missing_rate > 0) {
    if (missing_type == "MCAR") {
      missing_mask <- runif(nr) < missing_rate
    } else if (missing_type == "MAR") {
      prob_miss <- 1 - pnorm(tl, mean = 0, sd = sd(tl))
      prob_miss <- missing_rate * prob_miss / mean(prob_miss)
      prob_miss <- pmax(0, pmin(0.95, prob_miss))
      missing_mask <- runif(nr) < prob_miss
    }
    y[missing_mask] <- NA
  }

  data.frame(id = id, time = tm, item = it, y = y, ys = ys,
             tau_true = tl, ier = ier_mask, missing = missing_mask,
             stringsAsFactors = FALSE)
}

#' @keywords internal
.clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))
