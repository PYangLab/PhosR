#' @title Select by treatment groups (replicate block)
#'
#' @description Select phosphosites that have been quantified in a given
#' percentage of treatment groups (e.g. 0.75 as 3 out of 4 replicates)
#' in n groups.
#'
#' @author Pengyi Yang, Taiyun Kim
#'
#' @usage selectGrps(mat, grps, percent, n)
#'
#' @param mat a matrix (SummarizedExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples in replicates for different 
#' treatments. 
#' @param grps a string specifying the grouping (replicates).
#' @param percent a percent from 0 to 1, specifying the percentage of quantified
#'  values in any treatment group.
#' @param n an integer indicating n or more replicates pass the percentage
#' filtering for a phosphosite to be included.
#'
#' @return a filtered matrix with at least 'percent' quantification
#' in one or more conditions. If an input \code{mat} is a SummarizedExperiment 
#' object, filtered SummarizedExperiment object will be returned.
#'
#' @examples
#' data('phospho.cells.Ins.sample')
#' grps = gsub('_[0-9]{1}', '', colnames(phospho.cells.Ins))
#'
#' phospho.cells.Ins.filtered <- selectGrps(phospho.cells.Ins, grps, 0.5, n=1)
#' 
#' # For SummarizedExperiment objects
#' # mat = SummarizedExperiment::SummarizedExperiment(
#' #     assay = phospho.cells.Ins,
#' #     colData = S4Vectors::DataFrame(
#' #         groups = grps
#' #     )
#' # )
#' # phospho.cells.Ins.filtered = selectGrps(mat, mat$groups, 0.5, n=1)
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#'
#' @export
#'
selectGrps <- function(mat, grps, percent, n = 1) {
    if (missing(mat))
        stop("Parameter mat is missing!")
    if (missing(grps))
        stop("Parameter grps is missing!")
    if (missing(percent))
        stop("Parameter percent is missing!")

    if (length(grps) != ncol(mat))
        stop("Length of vector grps must be equal to number of columns in mat")
    if ((percent < 0) || (percent > 1))
        stop("Parameter percent must be a numeric value between 0 and 1")
    
    se = FALSE
    if (methods::is(mat, "SummarizedExperiment")) {
        mat.orig = mat
        mat = SummarizedExperiment::assay(mat)
        se = TRUE
    }
    

    # split the matrix by groups and organise them to a list
    tmp <- lapply(split(seq_len(ncol(mat)), grps), function(i) mat[,
        i])

    test <- do.call(cbind, lapply(tmp, function(x) {
        rowSums(!is.na(x))/ncol(x) >= percent
    }))
    
    sel = rowSums(test) >= n
    if (se) {
        mat.filtered <- SummarizedExperiment::SummarizedExperiment(
            mat[sel,],
            colData = SummarizedExperiment::colData(mat.orig),
            rowData = SummarizedExperiment::rowData(mat.orig)[sel,,drop = FALSE]
        )
    } else {
        mat.filtered <- mat[sel, ]
    }
    
    return(mat.filtered)
}


#' @title selectTimes
#'
#' @usage selectTimes(mat, timepoint, order, percent, w)
#'
#' @param mat a matrix (or SummarizedExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples in replicates for different 
#' treatments.
#' @param timepoint a timepoint as factor with a length equal to the
#' number of columns of mat.
#' @param order a vector specifying the order of timepoints.
#' @param percent a percent (decimal) from 0 to 1, to filter phosphosites with
#' with missing value larger than percent per timepoint.
#' @param w a timepoint window for selection of phosphosites to remove.
#'
#' @return a filtered matrix. If param \code{mat} is a SummarizedExperiment 
#' object, a SummarizedExperiment object will be returned.
#'
#' @examples
#' data("phospho_liverInsTC_RUV_sample")
#' timepoint = gsub("(.*)(\\d+[ms])(.*)", "\\2",
#'                 colnames(phospho.liver.Ins.TC.ratio.RUV))
#' timepoint[which(timepoint == "0m")] = "0s"
#' timepoint = factor(timepoint)
#' timepointOrder = c("0s", "5s", "1m", "2m", "3m", "4m", "6m")
#'
#' # For demonstration purpose, we introduce missing value at 0s
#' table(timepoint)
#'
#' phospho.liver.Ins.TC.sim = phospho.liver.Ins.TC.ratio.RUV
#' rmId = which(timepoint == "0s")
#'
#' # We replace the values to NA for the first 26 (~60%) of the '0s' samples
#' # for the first 100 phosphosite as NA
#' phospho.liver.Ins.TC.sim[seq(100),rmId[seq(26)]] = NA
#'
#' phospho.liver.Ins.TC.sim = selectTimes(phospho.liver.Ins.TC.sim,
#'                                     timepoint, timepointOrder, 0.5,
#'                                     w = length(table(timepoint)))
#' 
#' # For SummarizedExperiment objects
#' # mat = SummarizedExperiment::SummarizedExperiment(
#' #     assay = phospho.liver.Ins.TC.sim,
#' #     colData = S4Vectors::DataFrame(
#' #         timepoint = timepoint
#' #     )
#' # )
#' # phospho.liver.Ins.TC.sim = selectTimes(mat, mat$timepoint, timepointOrder, 
#' #       0.5, w = length(table(mat$timepoint)))
#' 
#' # Before filtering
#' dim(phospho.liver.Ins.TC.ratio.RUV)
#' # After filtering
#' dim(phospho.liver.Ins.TC.sim)
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' 
#' @export
#'
selectTimes <- function(mat, timepoint, order, percent, w = 1) {
    if (missing(mat))
        stop("Parameter mat is missing!")
    if (missing(timepoint))
        stop("Parameter timepoint is missing!")
    if (missing(order))
        stop("Parameter order is missing!")
    if (missing(percent))
        stop("Parameter percent is missing!")
    if ((percent < 0) || (percent > 1))
        stop("Parameter percent must be a numeric value between 0 and 1")
    
    
    se = FALSE
    if (methods::is(mat, "SummarizedExperiment")) {
        mat.orig = mat
        mat = SummarizedExperiment::assay(mat)
        se = TRUE
    }
    
    # split the matrix by groups and organise them to a list
    tmp <- lapply(split(seq_len(ncol(mat)), timepoint), function(i) mat[,
        i,drop = FALSE])

    test <- do.call(cbind, lapply(tmp, function(x) {
        rowSums(!is.na(x))/ncol(x) >= percent
    }))[, order]

    if ((w == 1) | (w == length(order))) {
        sel = rowSums(test) >= w
        mat.filtered <- mat[sel, ]
    } else if (w > length(order)) {
        stop("w is greater than the number of time points.
Please try a smaller w.")
    } else {
        idx <- length(order) - w + 1
        sel = Reduce(union, lapply(seq_len(idx), function(i) {
            which(rowSums(test[, seq(i, i+w-1),drop = FALSE]) == w)
        }))
        mat.filtered <- mat[sel, ]
    }
    if (se) {
        mat.filtered = SummarizedExperiment::SummarizedExperiment(
            mat.filtered,
            colData = SummarizedExperiment::colData(mat.orig),
            rowData = SummarizedExperiment::rowData(mat.orig)[sel,,drop = FALSE]
        )
    }
    
    return(mat.filtered)
}



#' @title Select phosphosite by percentage of quantification
#'
#' @description Select phosphosites that have been quantified in more than a
#' given percentage of samples
#'
#' @usage selectOverallPercent(mat, percent=NULL, n=NULL)
#'
#' @param mat a matrix (or SummarizedExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples in replicates for different 
#' treatments.
#' @param percent a percent from 0 to 1, specifying the percentage of quantified
#' values in across all samples for retaining a phosphosite for subsequent
#' analysis.
#' @param n an integer indicating n or more quantified values required for
#' retaining a phosphosite for subsequent analysis.
#'
#' @return a filtered matrix
#'
#' @examples
#'
#' data('phospho.cells.Ins.sample')
#'
#' phospho.cells.Ins.filtered <- selectOverallPercent(phospho.cells.Ins, 0.5)
#'
#' # Before filtering
#' dim(phospho.cells.Ins)
#' # After filtering
#' dim(phospho.cells.Ins.filtered)
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' 
#' @export
#'
selectOverallPercent <- function(mat, percent = NULL, n = NULL) {
    mat.filtered <- mat
    se = FALSE
    if (methods::is(mat, "SummarizedExperiment")) {
        mat.orig = mat
        se = TRUE
        mat.filtered = SummarizedExperiment::assay(mat)
    }

    if (missing(mat))
        stop("Parameter mat is missing!")
    if ((!is.null(percent)) && ((percent < 0) || (percent > 1)))
        stop("Parameter percent must be a numeric value between 0 and 1")
    if (is.null(percent) & is.null(n))
        stop("specify either percentage of number of quantified values for a
given phosphosite to be retained.")

    if (!is.null(percent) & is.null(n)) {
        sel = rowSums(!is.na(mat))/ncol(mat) >= percent
        mat.filtered <- mat[sel, ]
    } else if (is.null(percent) & !is.null(n)) {
        sel = rowSums(!is.na(mat)) >= n
        mat.filtered <- mat[sel, ]
    } else {
        sel = (rowSums(!is.na(mat))/ncol(mat) >=
                percent) & (rowSums(!is.na(mat)) >= n)
        mat.filtered <- mat[sel, ]
    }
    if (se) {
        mat.filtered = SummarizedExperiment::SummarizedExperiment(
            mat.filtered,
            colData = SummarizedExperiment::colData(mat.orig),
            rowData = SummarizedExperiment::rowData(mat.orig)[sel,,drop = FALSE]
        )
    }
    
    return(mat.filtered)
}
