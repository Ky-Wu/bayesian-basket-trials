UniformBasketEfficacy <- function(n_b, y, p0, a = 1, b = 1) {
  stopifnot(length(n_b) == length(y))
  vapply(seq_along(y), function(i) {
    a_new <- a + y[i]
    b_new <- b + n_b[i] - y[i]
    pbeta(p0, a_new, b_new, lower.tail = FALSE)
  }, numeric(1))
}
