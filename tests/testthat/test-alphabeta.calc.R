test_that("alpha is correct!", {
  mu <- 0.8
  var <- 0.1
  alpha <- (mu * (((mu*(1-mu))/var) - 1))
  expect_equal(alphabeta.calc(mu, var)[1], alpha)
})

test_that("beta is correct!", {
  mu <- 0.8
  var <- 0.1
  beta <- ((1-mu) * (((mu*(1-mu))/var) - 1))
  expect_equal(alphabeta.calc(mu, var)[2], beta)
})
