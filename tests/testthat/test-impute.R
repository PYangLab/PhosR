context("Test impute")

############################################
# scImpute
############################################

test_that(
  "Case A.1: missing parameters",
  {
    mat = matrix(rnorm(1200), nrow = 100)
    grps = 1:12
    percent = 0.3
    expect_error(scImpute(grps = grps, percent))
    expect_error(scImpute(mat = mat, grps, percent))
    expect_error(scImpute(mat, grps))
  }
)

# Error handling
test_that(
  "Case A.2: Error handling",
  {
    mat = matrix(rnorm(1200), nrow = 100)
    grps = 1:12
    percent = 0.3
    expect_error(scImpute(mat = mat, grps = 1:14, percent))
    expect_error(scImpute(mat = mat, grps, -1))
    expect_error(scImpute(mat = mat, grps, 50))
  }
)

# Parsed from scImpute
test_that(
  "Case A.3: expected number of rows imputed",
  {
    mat_s = matrix(1:1000, ncol = 10)
    mat_s3 = mat_s2 = mat_s

    tot_row = 100

    mat_s2[1:10, 1:3] = NA # imputation here
    mat_s2[11:20,1:5] = NA # no imputation here

    mat_s3[1:10, 1:2] = NA # imputation here

    percent_s = 0.7

    result1 = stImp(mat_s, percent_s)
    result2 = stImp(mat_s2, percent_s)
    result3 = stImp(mat_s3, percent_s)

    expect_identical(which(rowSums(!is.na(result1)) != ncol(mat_s)), integer(0))
    expect_identical(which(rowSums(!is.na(result2)) != ncol(mat_s2)), 11:20) # expect 11:20 to be not imputed
    expect_identical(which(rowSums(!is.na(result3)) != ncol(mat_s3)), integer(0))
  }
)

############################################
# tImpute
############################################

# Error handling
test_that(
  "Case B.1: Error handling",
  {
    expect_error(tImpute())
  }
)

# Expected output
test_that(
  "Case B.2: Expected output",
  {
    m = s = 1
    mat_34 = matrix(1:1200, ncol = 12)

    mat_34[1:10, 1:3] = NA
    mat_34_1 = tImpute(mat_34, m, s)
    expect_identical(all(!apply(mat_34, 2, is.na)), FALSE)
    expect_identical(all(!apply(mat_34_1, 2, is.na)), TRUE)
    expect_identical(dim(mat_34_1), dim(mat_34), TRUE)

  }
)

############################################
# ptImpute
############################################
# Error handling
test_that(
  "Case C.1: Error handling",
  {
    set.seed(123)
    mat = matrix(rnorm(1200), nrow = 100)
    percent = 0.5
    expect_error(ptImpute(mat1 = mat, percent1 = percent, percent2 = percent))
    expect_error(ptImpute(mat2 = mat, percent1 = percent, percent2 = percent))
    expect_error(ptImpute(mat1 = mat, mat2 = mat, percent2 = percent))
    expect_error(ptImpute(mat1 = mat, mat2 = mat, percent1 = percent))

  }
)

test_that(
  "Case C.2: Expected output",
  {
    mat = matrix(1:1200, ncol = 12)

    mat1 = mat[,1:6]
    mat2 = mat[,7:12]
    percent1 = 0.5
    percent2 = 0.5

    tmp = ptImpute(mat1, mat2, percent1, percent2)
    expect_identical(any(apply(tmp, 2, is.na)), FALSE)
    expect_identical(dim(tmp), dim(cbind(mat1, mat2)))
  }
)










