variance_twostage <- function(c_h, k_h, W_h, sigma2_bh, sigma2_wh) {
  sum(W_h^2 * ((sigma2_bh / c_h) + (sigma2_wh / (c_h * k_h))))
}