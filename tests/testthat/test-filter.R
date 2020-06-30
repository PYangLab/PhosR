context("Test filter")



###################################################
# selectGrps
###################################################

# Missing parameters
test_that(
  "Case A.1: missing parameters",
  {
    set.seed(123)
    mat = matrix(rnorm(1200), nrow = 100)
    n = 1
    expect_error(selectGrps(grps = grps, percent = percent))
    expect_error(selectGrps(mat = mat, percent = percent))
    expect_error(selectGrps(mat = mat, grps = grps))
  }
)

# Error handling
test_that(
  "Case A.2: Error handling",
  {
    set.seed(123)
    mat = matrix(rnorm(1200), nrow = 100) # 12 columns
    n = 1
    # Different size of grps and mat columns
    expect_error(selectGrps(mat = mat, grps = 1:10, n = n))
    expect_error(selectGrps(mat = mat, grps = 1:15, n = n))
    
    # Invalid percent values
    expect_error(selectGrps(mat = mat, grps = 1:12, percent = -1))
    expect_error(selectGrps(mat = mat, grps = 1:12, percent = 4))
  }
)

# Expected output
test_that(
  "Case A.3: expect filtered output",
  {
    mat = matrix(1:1200, ncol = 12)
    grps = c(rep("A", 6), rep("B", 6))
    
    mat2 = mat
    mat2[51:100, c(1:4,7:11)] = NA
    mat_expected = mat[1:50,]
    expect_identical(selectGrps(mat2, grps, 0.5, 1), mat_expected)
  }
)

test_that(
  "Case A.4: expect unfiltered output",
  {
    mat = matrix(1:1200, ncol = 12)
    grps = c(rep("A", 6), rep("B", 6))
    
    mat2 = mat
    mat2[51:100, c(1:3,7:9)] = NA
    expect_identical(selectGrps(mat2, grps, 0.5, 1), mat2)
  }
)



###################################################
# selectOverallPercent
###################################################

# Missing parameters
test_that(
  "Case B.1: missing parameters",
  {
    set.seed(123)
    mat = matrix(rnorm(1200), nrow = 100)
    n = 1
    expect_error(selectOverallPercent(percent = percent))
    expect_error(selectOverallPercent(mat = mat))
  }
)

# Error handling
test_that(
  "Case B.2: Error handling",
  {
    set.seed(123)
    mat = matrix(rnorm(1200), nrow = 100) # 12 columns
    
    # Invalid percent values
    expect_error(selectOverallPercent(mat = mat, percent = -1))
    expect_error(selectOverallPercent(mat = mat, percent = 4))
  }
)

test_that(
  "Case B.3: expect filtered output",
  {
    percent = 0.5
    percent_n = NULL
    n = 6
    n_n = NULL
    
    mat = matrix(1:1200, ncol = 12)
    
    mat1 = mat
    mat1[51:100, c(1:4,7:11)] = NA
    
    mat_expected = mat[1:50,]
    expect_identical(selectOverallPercent(mat1, percent, n_n), mat_expected)
    expect_identical(selectOverallPercent(mat1, percent_n, n), mat_expected)
    
  }
)


test_that(
  "Case B.4: expect unfiltered output",
  {
    percent = 0.1
    percent_n = NULL
    n = 1
    n_n = NULL
    
    mat = matrix(1:1200, ncol = 12)
    
    mat1 = mat
    mat1[51:100, c(1:4,7:11)] = NA
    
    
    mat_expected = mat1
    expect_identical(selectOverallPercent(mat1, percent, n_n), mat_expected)
    expect_identical(selectOverallPercent(mat1, percent_n, n), mat_expected)
    
  }
)

###################################################
# selectOverallPercent
###################################################
# Error handling
test_that(
  "Case C.1: missing parameters",
  {
    set.seed(123)
    mat = matrix(rnorm(1200), nrow = 100)
    timepoint = 1:12
    order = 1:12
    percent = 0.5
    expect_error(selectTimes(timepoint = timepoint, order = order, percent = percent))
    expect_error(selectTimes(mat = mat, order = order, percent = percent))
    expect_error(selectTimes(mat = mat, timepoint = timepoint, percent = percent))
    expect_error(selectTimes(mat = mat, timepoint = timepoint, order = order))
  }
)

test_that(
  "Incorrect input",
  {
    set.seed(123)
    mat = matrix(rnorm(1200), nrow = 100)
    timepoint = 1:12
    order = 1:12
    percent = 0.5
    
    # w > length(order)
    expect_error(selectTimes(mat, timepoint, order, percent, w = 1:15))
    # incorrect percent
    expect_error(selectTimes(mat, timepoint, order, percent = -1))
    expect_error(selectTimes(mat, timepoint, order, percent = 5))
  }
)






