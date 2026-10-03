# Test that iivDecomp functions run without errors

test_that("iiv_generate produces correct dimensions", {
  dat <- iiv_generate(N = 20, T = 10, K = 4, sigma_tau = 0.3, phi = 0.5, seed = 123)
  expect_equal(nrow(dat), 20 * 10 * 4)
  expect_true(all(c("id", "time", "item", "y", "ys", "tau_true", "ier") %in% names(dat)))
})

test_that("iiv_decompose returns correct structure", {
  dat <- iiv_generate(N = 50, T = 30, K = 4, sigma_tau = 0.3, phi = 0.5, seed = 123)
  res <- iiv_decompose(dat)
  expect_s3_class(res, "iivDecomp")
  expect_true(res$s2_tau > 0)
  expect_true(res$s2_eps > 0)
  expect_true(res$p_ier >= 0 && res$p_ier <= 1)
})

test_that("iiv_compare returns all methods", {
  dat <- iiv_generate(N = 50, T = 30, K = 4, sigma_tau = 0.3, phi = 0.5, seed = 123)
  comp <- iiv_compare(dat)
  expect_equal(nrow(comp), 5)
  expect_equal(comp$Method[1], "M1 Aggregation")
})

test_that("iiv_detect_longstring works", {
  dat <- iiv_generate(N = 20, T = 10, K = 4, sigma_tau = 0.3, phi = 0.5,
                       p_IER = 0.15, IER_type = "fixed", seed = 123)
  ls <- iiv_detect_longstring(dat, R = 5)
  expect_true("flag" %in% names(ls))
})

test_that("iiv_detect_mahalanobis works", {
  dat <- iiv_generate(N = 20, T = 10, K = 4, sigma_tau = 0.3, phi = 0.5,
                       p_IER = 0.15, IER_type = "random", seed = 123)
  md <- iiv_detect_mahalanobis(dat)
  expect_true("flag" %in% names(md))
})

test_that("iiv_detect_mahalanobis is empty-safe and aligned on clean data", {
  dat <- iiv_generate(N = 50, T = 30, K = 4, sigma_tau = 0.3, phi = 0.5, seed = 123)
  md <- iiv_detect_mahalanobis(dat)
  expect_equal(nrow(md), 50 * 30)
  expect_false(any(md$flag))            # clean data: nothing should be flagged
  expect_equal(nrow(md[md$flag, ]), 0)  # flagged-subset ops must not error
})
