# Find a nice example
library(bbbq)
library(testthat)

for (i in seq(1, 1000)) {
  message(i, ": ", Sys.time())
  f <- bbbq::create_random_tmh(n_aas = 30)
  expect_true(pureseqtmr::is_tmh(f))
  h <- sample(bbbq::get_mhc1_haplotypes(), size = 1)
  p <- calc_p_det_tmh_mut(f, h, n_adjancent_sequences = 5, percentile = 0.05)
  if (p > 0.0 && p < 1.0) {
    cat(f)
    cat(h)
    cat(p)
    break
  }
}
# 1: 2020-06-17 09:29:39

# Six minutes
# Percentile 0.1
# CMIGFVIYLLFGFILSLMCVFVLFVILILI HLA-B*27:05 0.3216847
