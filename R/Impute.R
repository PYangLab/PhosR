#' Site- and condition-specific (sc) impute
#'
#' Impute the missing values for a phosphosite across replicates within a single
#' condition (or treatment)
#' if there are n or more quantified values of that phosphosite in that
#' condition.
#'
#' @usage scImpute(mat, percent, grps, assay)
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to replicates within a condition.
#' @param percent a percent from 0 to 1, specifying the percentage of quantified
#' values in any treatment group.
#' @param grps a string specifying the grouping (replciates).
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return An imputed matrix. If param \code{mat} is a PhosphoExperiment 
#' object, a PhosphoExperiment object will be returned.
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
#' # for PhosphoExperiment Object
#' data('phospho.cells.Ins.pe')
#' grps = gsub('_[0-9]{1}', '', colnames(phospho.cells.Ins.pe))
#' phospho.cells.Ins.filtered <- selectGrps(phospho.cells.Ins.pe, grps, 
#'     0.5, n=1)
#'
#' set.seed(123)
#' phospho.cells.Ins.impute <-
#'     scImpute(phospho.cells.Ins.filtered,
#'     0.5,
#'     grps)[,colnames(phospho.cells.Ins.filtered)]
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' 
#' @export
scImpute <- function(mat, percent, grps, assay = NULL) {
    if (missing(mat)) {
        stop("Parameter mat is missing!")
    }
    if (missing(percent)) {
        stop("Parameter percent is missing!")
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
    
    pe = FALSE
    if (methods::is(mat, "PhosphoExperiment")) {
        pe = TRUE
        mat.orig = mat
        if (is.null(assay)) {
            mat = SummarizedExperiment::assay(mat)
        } else {
            mat = SummarizedExperiment::assay(mat, assay)
        }
    }
    tmp <- lapply(split(seq_len(ncol(mat)), grps), function(i) mat[,
        i])
    mat.imputed <- do.call(cbind, lapply(tmp, stImp, percent = percent))[,
        colnames(mat)]
    
    if (pe) {
        SummarizedExperiment::assay(mat.orig, "imputed") = mat.imputed
        mat.imputed = mat.orig
    }
    
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
#' @usage tImpute(mat, m, s, assay)
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples.
#' @param m a numeric number for controlling mean downshifting.
#' @param s a numeric number for controlling standard deviation of downshifted
#' sampling values.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return An imputed matrix. If param \code{mat} is a SummarizedExperiment 
#' object, a SummarizedExperiment object will be returned.
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
#' # For PhosphoExperiment Object
#' data('phospho.cells.Ins.pe')
#' grps = gsub('_[0-9]{1}', '', colnames(phospho.cells.Ins.pe))
#' phospho.cells.Ins.filtered <- selectGrps(phospho.cells.Ins.pe, grps, 
#'     0.5, n=1)
#'
#' set.seed(123)
#' phospho.cells.Ins.impute <- tImpute(phospho.cells.Ins.filtered)
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom methods is
#'
#' @export
tImpute <- function(mat, m = 1.6, s = 0.6, assay = NULL) {
    if (missing(mat)) {
        stop("Paramter mat is missing!")
    }
    
    mat.orig = mat
    pe = FALSE
    if (methods::is(mat, "PhosphoExperiment")) {
        if (is.null(assay)) {
            mat = SummarizedExperiment::assay(mat)
        } else {
            mat = SummarizedExperiment::assay(mat, assay)
        }
        pe = TRUE
    }

    ms <- colMeans(mat, na.rm = TRUE)
    sds <- apply(mat, 2, stats::sd, na.rm = TRUE)
    mat.impute <- mat
    for (i in seq_len(ncol(mat.impute))) {
        r <- stats::rnorm(n = sum(is.na(mat.impute[, i])), mean = (ms[i] -
            sds[i] * m), sd = (sds[i] * s))
        mat.impute[which(is.na(mat.impute[, i])), i] <- r
    }
    
    if (pe) {
        SummarizedExperiment::assay(mat.orig, "imputed") = mat.impute
        mat.impute = mat.orig
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
#' @usage ptImpute(
#'     mat1, 
#'     mat2, 
#'     percent1, 
#'     percent2, 
#'     m = 1.6, 
#'     s = 0.6, 
#'     paired = TRUE, 
#'     verbose = TRUE,
#'     assay
#' )
#'
#' @param mat1 a matrix (or PhosphoExperiment object) with rows correspond to
#'  phosphosites and columns correspond to replicates within treatment1.
#' @param mat2 a matrix (or PhosphoExperiment object) with rows correspond to
#'  phosphosites and columns correspond to replicates within treatment2.
#' @param percent1 a percent indicating minimum quantified percentages required
#' for considering for imputation.
#' @param percent2 a percent indicating minimum quantified percentages required
#' for considering for imputation.
#' @param m a numeric number of for controlling mean downshifting.
#' @param s a numeric number of for controlling standard deviation of
#' downshifted sampling values.
#' @param paired a flag indicating whether to impute for both treatment1 and
#' treatment2 (default) or treatment2 only (if paired=FALSE).
#' @param verbose Default to \code{TRUE} to show messages during the progress.
#' All messages will be suppressed if set to \code{FALSE}
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#' 
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
#' phospho.cells.Ins.impute[,seq(6)] <- 
#'     ptImpute(phospho.cells.Ins.impute[,seq(7,12)],
#' phospho.cells.Ins.impute[,seq(6)], percent1 = 0.6, percent2 = 0, 
#'     paired = FALSE)
#' 
#' 
#' # For PhosphoExperiment objects
#' # mat = PhosphoExperiment(
#' #     assay = phospho.cells.Ins.impute,
#' #     colData = S4Vectors::DataFrame(
#' #         groups = grps
#' #     )
#' # )
#' # SummarizedExperiment::assay(mat)[,seq(6)] <- 
#' #     ptImpute(SummarizedExperiment::assay(mat)[,seq(7,12)],
#' #         SummarizedExperiment::assay(mat)[,seq(6)], percent1 = 0.6, 
#' #         percent2 = 0, paired = FALSE)
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowData colData
#' @importFrom methods is
#' 
#' @export
#'
ptImpute <- function(mat1, mat2, percent1, percent2, m = 1.6,
    s = 0.6, paired = TRUE, verbose = TRUE, assay = NULL) {
    if (missing(mat1))
        stop("Paramter mat1 is missing!")
    if (missing(mat2))
        stop("Paramter mat2 is missing!")
    if (missing(percent1))
        stop("Paramter percent1 is missing!")
    if (missing(percent2))
        stop("Paramter percent2 is missing!")
    
    pe = FALSE
    mat1.orig = mat1
    mat2.orig = mat2
    
    if (methods::is(mat1, "PhosphoExperiment") && 
            methods::is(mat2, "PhosphoExperiment")) {
        pe = TRUE
        if (is.null(assay)) {
            mat1 = SummarizedExperiment::assay(mat1)
            mat2 = SummarizedExperiment::assay(mat2)
        } else {
            mat1 = SummarizedExperiment::assay(mat1, assay)
            mat2 = SummarizedExperiment::assay(mat2, assay)
        }
    }
    # matrix before imputation
    mat1.raw = mat1
    mat2.raw = mat2
    
    # impute for mat2
    idx1 <- which((rowSums(!is.na(mat1))/ncol(mat1)) >= percent1 &
        (rowSums(!is.na(mat2))/ncol(mat2)) <= percent2)
    
    if (verbose)
        message(paste("idx1:", length(idx1)))
    
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
        # If paired = TRUE, estimate idx2 before imputing mat2
        idx2 <- which((rowSums(!is.na(mat2.raw))/ncol(mat2.raw)) >= percent1 &
                (rowSums(!is.na(mat1.raw))/ncol(mat1.raw)) <= percent2)
        if (verbose)
            message(paste("idx2:", length(idx2)))

        ms <- colMeans(mat1, na.rm = TRUE)
        sds <- apply(mat1, 2, stats::sd, na.rm = TRUE)
        if (length(idx2) > 0) {
            for (i in seq_len(ncol(mat1))) {
                mat1[idx2, i] <- stats::rnorm(length(idx2), mean = (ms[i] -
                    sds[i] * m), sd = (sds[i] * s))
            }
        }
        if (pe) {
            mat = cbind(mat1.orig, mat2.orig)
            if (is.null(assay)) {
                SummarizedExperiment::assay(mat, withDimnames = FALSE) = 
                    cbind(mat1, mat2)
            } else {
                SummarizedExperiment::assay(mat, assay, withDimnames = FALSE) = 
                    cbind(mat1, mat2)
            }
            
        } else {
            mat = cbind(mat1, mat2)
        }
        
        return(mat)
    } else {
        if (pe) {
            if (is.null(assay)) {
                SummarizedExperiment::assay(mat2.orig, withDimnames = FALSE) = 
                    mat2
            } else {
                SummarizedExperiment::assay(mat2.orig, assay, 
                    withDimnames = FALSE) = mat2            
            }
            mat2 = mat2.orig
        }
        return(mat2)
    }
}
