noquantile <- function(x, probs) {
  # Estimates quantiles based on NO estimator introduced by 
  # (Navruz and Özdemir, 2020) at probabilites in probs.
  # (Navruz, G., & Özdemir, A. F. (2020). A new quantile estimator 
  # with weights based on a subsampling approach. British Journal of 
  # Mathematical and Statistical Psychology, 73(3), 506-521.)
  # implementation: https://aakinshin.net/posts/navruz-ozdemir-quantile-estimator/
  x <- sort(x)
  n <- length(x)
  if (n <= 2) return(quantile(x, probs))
  sapply(probs, function(p) {
    B <- function(x) dbinom(x, n, p)
    (B(0) * 2 * p + B(1) * p) * x[1] + B(0) * (2 - 3 * p) * x[2] - B(0) * (1 - p) * x[3] +
      sum(sapply(1:(n-2), function(i) (B(i) * (1 - p) + B(i + 1) * p) * x[i + 1])) -
      B(n) * p * x[n - 2] + B(n) * (3 * p - 1) * x[n - 1] + (B(n - 1) * (1 - p) + B(n) * (2 - 2 * p)) * x[n]
  })
}