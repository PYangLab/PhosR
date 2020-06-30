#' Site- and condition-specific (sc) impute
#'
#' Impute the missing values for a phosphosite across replicates within a single
#' condition (or treatment)
#' if there are n or more quantified values of that phosphosite in that
#' condition.
#'
#' @usage scImpute(mat, percent, grps)
#'
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to replicates within a condition.
#' @param percent a percent from 0 to 1, specifying the percentage of quantified
#' values in any treatment group.
#' @param grps a string specifying the grouping (replciates).
#'
#' @return An imputed matrix
#'
#' @examples
#'
#' data('phospho.cells.Ins.sample')
#' grps = gsub('_[0-9]{1}', '', colnames(phospho.cells.Ins))
#' phospho.cells.Ins.filtered <- selectGrps(phospho.cells.Ins, grps, 0.5, n=1)
#'
#' set.seed(123)
#' phospho.cells.Ins.impute <-
#'     scImpute(phospho.cells.Ins.filtered,
#'     0.5,
#'     grps)[,colnames(phospho.cells.Ins.filtered)]
#'
#' @export
scImpute <- function(mat, percent, grps) {
    if (missing(mat)) {
        stop("Parameter mat is missing!")
    }
    if (missing(percent)) {
        stop("Parameter percent is percent!")
    }
    if (missing(grps)) {
        stop("Parameter grps is missing!")
    }
    if (length(grps) != ncol(mat)) {
        stop("Length of vector grps must be equal to number of columns in mat")
    }
    if ((percent < 0) || (percent > 1)) {
        stop("Parameter percent must be a numeric value between 0 and 1")
    }

    tmp <- lapply(split(seq_len(ncol(mat)), grps), function(i) mat[,
        i])
    mat.imputed <- do.call(cbind, lapply(tmp, stImp, percent = percent))[,
        colnames(mat)]
    return(mat.imputed)
}


stImp <- function(mat, percent) {
    for (i in seq_len(nrow(mat))) {
        idx <- which(!is.na(mat[i, ]))
        if (length(idx)/ncol(mat) >= percent) {
            ms <- mean(mat[i, ], na.rm = TRUE)
            sds <- stats::sd(mat[i, ], na.rm = TRUE)
            nid <- which(is.na(mat[i, ]))
            mat[i, nid] <- stats::rnorm(length(nid), mean = ms,
                sd = sds)
        }
    }
    return(mat)
}

#' Tail-based impute
#'
#' Tail-based imputation approach as implemented in Perseus.
#'
#' @usage tImpute(mat, m, s)
#'
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to samples.
#' @param m a numeric number for controlling mean downshifting.
#' @param s a numeric number for controlling standard deviation of downshifted
#' sampling values.
#'
#' @return An imputed matrix
#'
#' @examples
#'
#' data('phospho.cells.Ins.sample')
#' grps = gsub('_[0-9]{1}', '', colnames(phospho.cells.Ins))
#' phospho.cells.Ins.filtered <- selectGrps(phospho.cells.Ins, grps, 0.5, n=1)
#'
#' set.seed(123)
#' phospho.cells.Ins.impute <- tImpute(phospho.cells.Ins.filtered)
#'
#' @export
tImpute <- function(mat, m = 1.6, s = 0.6) {
    if (missing(mat)) {
        stop("Paramter mat is missing!")
    }

    ms <- colMeans(mat, na.rm = TRUE)
    sds <- apply(mat, 2, stats::sd, na.rm = TRUE)
    mat.impute <- mat
    for (i in seq_len(ncol(mat.impute))) {
        r <- stats::rnorm(n = sum(is.na(mat.impute[, i])), mean = (ms[i] -
            sds[i] * m), sd = (sds[i] * s))
        mat.impute[which(is.na(mat.impute[, i])), i] <- r
    }
    return(mat.impute)
}



#' Paired-tail (pt) based impute
#'
#' Impute the missing values for mat2 using tail imputation approach if mat1 has
#' more than percent1 (percentage) of quantified values
#' and mat2 has less than percent2 (percentage) quantified values,
#' and vice versa if paired is set to be true. That is if mat2 has percentage of
#' quantified values more than percent1 and mat1 has percentage quantified
#' values less than percent2.
#'
#' @usage ptImpute(mat1, mat2, percent1, percent2, m, s, paired)
#'
#' @param mat1 a matrix with rows correspond to phosphosites and columns
#' correspond to replicates within treatment1.
#' @param mat2 a matrix with rows correspond to phosphosites and columns
#' correspond to replicates within treatment2.
#' @param percent1 a percent indicating minimum quantified percentages required
#' for considering for imputation.
#' @param percent2 a percent indicating minimum quantified percentages required
#' for considering for imputation.
#' @param m a numeric number of for controlling mean downshifting.
#' @param s a numeric number of for controlling standard deviation of
#' downshifted sampling values.
#' @param paired a flag indicating whether to impute for both treatment1 and
#' treatment2 (default) or treatment2 only (if paired=FALSE).
#'
#' @return An imputed matrix
#'
#' @examples
#'
#' data('phospho.cells.Ins.sample')
#' grps = gsub('_[0-9]{1}', '', colnames(phospho.cells.Ins))
#' phospho.cells.Ins.filtered <- selectGrps(phospho.cells.Ins, grps, 0.5, n=1)
#'
#' set.seed(123)
#' phospho.cells.Ins.impute <-
#'     scImpute(
#'     phospho.cells.Ins.filtered,
#'     0.5,
#'     grps)[,colnames(phospho.cells.Ins.filtered)]
#'
#' set.seed(123)
#' phospho.cells.Ins.impute[,1:5] <- ptImpute(phospho.cells.Ins.impute[,6:10],
#' phospho.cells.Ins.impute[,1:5], percent1 = 0.6, percent2 = 0, paired = FALSE)
#'
#' @export
#'
ptImpute <- function(mat1, mat2, percent1, percent2, m = 1.6,
    s = 0.6, paired = TRUE) {
    if (missing(mat1))
        stop("Paramter mat1 is missing!")
    if (missing(mat2))
        stop("Paramter mat2 is missing!")
    if (missing(percent1))
        stop("Paramter percent1 is missing!")
    if (missing(percent2))
        stop("Paramter percent2 is missing!")

    # impute for mat2
    idx1 <- which((rowSums(!is.na(mat1))/nrow(mat1)) >= percent1 &
        (rowSums(!is.na(mat2))/nrow(mat2)) <= percent2)
    print(paste("idx1:", length(idx1)))

    ms <- colMeans(mat2, na.rm = TRUE)
    sds <- apply(mat2, 2, stats::sd, na.rm = TRUE)
    if (length(idx1) > 0) {
        for (i in seq_len(ncol(mat2))) {
            mat2[idx1, i] <- stats::rnorm(length(idx1), mean = (ms[i] -
                sds[i] * m), sd = (sds[i] * s))
        }
    }

    if (paired == TRUE) {
        # impute for mat1
        idx2 <- which((rowSums(!is.na(mat2))/nrow(mat2)) >= percent1 &
            (rowSums(!is.na(mat1))/nrow(mat1)) <= percent2)
        print(paste("idx2:", length(idx2)))

        ms <- colMeans(mat1, na.rm = TRUE)
        sds <- apply(mat1, 2, stats::sd, na.rm = TRUE)
        if (length(idx2) > 0) {
            for (i in seq_len(ncol(mat1))) {
                mat1[idx2, i] <- stats::rnorm(length(idx2), mean = (ms[i] -
                    sds[i] * m), sd = (sds[i] * s))
            }
        }
        return(cbind(mat1, mat2))
    } else {
        return(mat2)
    }
}
