context("Test pathwayAnalysis")

# loading data
alter = "greater"
set.seed(123)
geneStats = rnorm(11)
names(geneStats) = letters[c(5:10, 18:22)]
geneSet = names(geneStats)
annotation = list(
  pathway1 = letters[1:10],
  pathway2 = tail(letters, 10),
  pathway3 = letters[15:20],
  pathway4 = LETTERS[c(1:3,4,20)]
)
universe = c(letters, LETTERS)



###################################################
# pathwayOverrepresent
###################################################

# Error handling
test_that(
  "Case B.1: Missing parameters",
  {
    expect_error(patwhayOverrepresent(annotation = annotation, universe = universe))
    expect_error(patwhayOverrepresent(geneSet = geneSet, universe = universe))
    expect_error(patwhayOverrepresent(geneSet = geneSet, annotation = annotation))
  }
)

# Expected output
test_that(
  "Case B.2: Expected output",
  {
    alter1 = "greater"
    alter2 = "less"
    alter3 = "two.sided"
    rNames = paste("pathway", 1:4, sep = "")
    cNames = c("pvalue", "# of substrates", "substrates")
    result1 = matrix(
      c(
        c("0.00318849739409564", "0.0250735382249623", "0.100641107506833", "1"),
        c("6", "5", "3", "0"),
        c(paste(letters[5:10], collapse = "|"), paste(letters[18:22], collapse = "|"), paste(letters[18:20], collapse = "|"), "")
      ),
      nrow = 4
    )
    colnames(result1) = cNames
    rownames(result1) = rNames
    result2 = matrix(
      c(
        c("0.288345338135254", "0.98575515312508", "0.996811502605904", "0.999768940556022"),
        rev(c("6", "5", "3", "0")),
        rev(c(paste(letters[5:10], collapse = "|"), paste(letters[18:22], collapse = "|"), paste(letters[18:20], collapse = "|"), ""))
      ),
      nrow = 4
    )
    colnames(result2) = cNames
    rownames(result2) = rev(rNames)
    result3 = matrix(
      c(
        c("0.00318849739409564", "0.0250735382249623", "0.100641107506833", "0.571378551420568"),
        c("6", "5", "3", "0"),
        c(paste(letters[5:10], collapse = "|"), paste(letters[18:22], collapse = "|"), paste(letters[18:20], collapse = "|"), "")
      ),
      nrow = 4
    )
    colnames(result3) = cNames
    rownames(result3) = rNames


    expect_identical(result1, pathwayOverrepresent(geneSet, annotation, universe, alter1))
    expect_identical(result2, pathwayOverrepresent(geneSet, annotation, universe, alter2))
    expect_identical(result3, pathwayOverrepresent(geneSet, annotation, universe, alter3))
  }
)


###################################################
# pathwayRankBasedEnrichment
###################################################

# Error handling
test_that(
  "Case C.1: Missing parameters",
  {
    expect_error(pathwayRankBasedEnrichment(annotation = annotation))
    expect_error(pathwayRankBasedEnrichment(geneStats = geneStats))
  }
)

# Expected output
test_that(
  "Case 3.4: Expected output",
  {
    alter1 = "greater"
    alter2 = "less"
    alter3 = "two.sided"
    rNames = paste("pathway", 1:4, sep = "")
    cNames = c("pvalue", "# of substrates", "substrates")
    result1 = matrix(
      c(
        c("0.164502164502164", "0.876623376623377", NA, NA),
        c("6", "5", NA, NA),
        c(paste(letters[5:10], collapse = ";"), paste(letters[18:22], collapse = ";"), NA, NA)
      ),
      nrow = 4
    )
    colnames(result1) = cNames
    rownames(result1) = rNames
    result2 = matrix(
      c(
        c("0.164502164502164", "0.876623376623377", NA, NA),
        c("5", "6", NA, NA),
        c(paste(letters[18:22], collapse = ";"), paste(letters[5:10], collapse = ";"), NA, NA)
      ),
      nrow = 4
    )
    colnames(result2) = cNames
    rownames(result2) = c(paste("pathway", c("2", "1", "3", "4"), sep = ""))
    result3 = matrix(
      c(
        c("0.329004329004329", "0.329004329004329", NA, NA),
        c("6", "5", NA, NA),
        c(paste(letters[5:10], collapse = ";"), paste(letters[18:22], collapse = ";"), NA, NA)
      ),
      nrow = 4
    )
    colnames(result3) = cNames
    rownames(result3) = rNames


    expect_identical(result1, pathwayRankBasedEnrichment(geneStats, annotation, alter1))
    expect_identical(result2, pathwayRankBasedEnrichment(geneStats, annotation, alter2))
    expect_identical(result3, pathwayRankBasedEnrichment(geneStats, annotation, alter3))
  }
)


















