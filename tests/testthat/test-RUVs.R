context("Test RUVs")

test_that(
  "Case A.1: missing parameters",
  {
    mat = matrix(rnorm(1200), nrow = 12)
    M = mat
    ctl = letters[1:12]
    expect_error(RUVphospho(M = M, ctl = ctl))
    expect_error(RUVphospho(mat = mat, ctl = ctl))
    expect_error(RUVphospho(mat = mat, M = M))
  }
)





