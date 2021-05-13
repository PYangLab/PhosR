#' @title Select by treatment groups (replicate block)
#'
#' @description Select phosphosites that have been quantified in a given
#' percentage of treatment groups (e.g. 0.75 as 3 out of 4 replicates)
#' in n groups.
#'
#' @author Pengyi Yang, Taiyun Kim
#'
#' @usage selectGrps(mat, grps, percent, n, assay)
#'
#' @param mat a matrix (PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples in replicates for different 
#' treatments. 
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#' @param grps a string specifying the grouping (replicates).
#' @param percent a percent from 0 to 1, specifying the percentage of quantified
#'  values in any treatment group.
#' @param n an integer indicating n or more replicates pass the percentage
#' filtering for a phosphosite to be included.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return a filtered matrix (or a PhosphoExperiment Oject) with at least 
#' 'percent' quantification in one or more conditions. If an input \code{mat} is 
#' a SummarizedExperiment object, filtered SummarizedExperiment object will be 
#' returned.
#'
#' @examples
#' data('phospho.cells.Ins.sample')
#' grps = gsub('_[0-9]{1}', '', colnames(phospho.cells.Ins))
#' phospho.cells.Ins.filtered <- selectGrps(phospho.cells.Ins, grps, 0.5, n=1)
#' 
#' # For PhosphoExperiment object
#' data('phospho.cells.Ins.pe')
#' grps = gsub('_[0-9]{1}', '', colnames(phospho.cells.Ins.pe))
#' phospho.cells.Ins.filtered <- selectGrps(phospho.cells.Ins.pe, grps, 0.5, n=1)
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#'
#' @export
#'
selectGrps <- function(mat, grps, percent, n = 1, assay = NULL) {
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
    
    mat.filtered = mat
    
    if (methods::is(mat, "PhosphoExperiment")) {
        if (is.null(assay)) {
            mat = SummarizedExperiment::assay(mat)
        } else {
            mat = SummarizedExperiment::assay(mat, assay)
        }
    }
    

    # split the matrix by groups and organise them to a list
    tmp <- lapply(split(seq_len(ncol(mat)), grps), function(i) mat[,
        i])

    test <- do.call(cbind, lapply(tmp, function(x) {
        rowSums(!is.na(x))/ncol(x) >= percent
    }))
    
    sel = rowSums(test) >= n
    
    return(mat.filtered[sel,])
}


#' @title selectTimes
#'
#' @usage selectTimes(mat, timepoint, order, percent, w, assay)
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples in replicates for different 
#' treatments.
#' @param timepoint a timepoint as factor with a length equal to the
#' number of columns of mat.
#' @param order a vector specifying the order of timepoints.
#' @param percent a percent (decimal) from 0 to 1, to filter phosphosites with
#' with missing value larger than percent per timepoint.
#' @param w a timepoint window for selection of phosphosites to remove.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
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
#' # For PhosphoExperiment objects
#' # mat = PhosR::PhosphoExperiment(
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
selectTimes <- function(mat, timepoint, order, percent, w = 1, assay = NULL) {
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
    
    
    mat.orig = mat
    if (methods::is(mat, "PhosphoExperiment")) {
        if (is.null(assay)) {
            mat = SummarizedExperiment::assay(mat)
        } else {
            mat = SummarizedExperiment::assay(mat, assay)
        }
    }
    
    # split the matrix by groups and organise them to a list
    tmp <- lapply(split(seq_len(ncol(mat)), timepoint), function(i) mat[,
        i,drop = FALSE])

    test <- do.call(cbind, lapply(tmp, function(x) {
        rowSums(!is.na(x))/ncol(x) >= percent
    }))[, order]

    if ((w == 1) | (w == length(order))) {
        sel = rowSums(test) >= w
        mat.filtered <- mat.orig[sel, ]
    } else if (w > length(order)) {
        stop("w is greater than the number of time points.
Please try a smaller w.")
    } else {
        idx <- length(order) - w + 1
        sel = Reduce(union, lapply(seq_len(idx), function(i) {
            which(rowSums(test[, seq(i, i+w-1),drop = FALSE]) == w)
        }))
        mat.filtered <- mat.orig[sel, ]
    }
    
    
    return(mat.filtered)
}



#' @title Select phosphosite by percentage of quantification
#'
#' @description Select phosphosites that have been quantified in more than a
#' given percentage of samples
#'
#' @usage selectOverallPercent(mat, percent, n, assay)
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples in replicates for different 
#' treatments.
#' @param percent a percent from 0 to 1, specifying the percentage of quantified
#' values in across all samples for retaining a phosphosite for subsequent
#' analysis.
#' @param n an integer indicating n or more quantified values required for
#' retaining a phosphosite for subsequent analysis.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
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
selectOverallPercent <- function(mat, percent = NULL, n = NULL, assay = NULL) { 


    if (missing(mat))
        stop("Parameter mat is missing!")
    if ((!is.null(percent)) && ((percent < 0) || (percent > 1)))
        stop("Parameter percent must be a numeric value between 0 and 1")
    if (is.null(percent) & is.null(n))
        stop("specify either percentage of number of quantified values for a
given phosphosite to be retained.")
    
    mat.orig = mat
    if (methods::is(mat, "PhosphoExperiment")) {
        if (is.null(assay)) {
            mat = SummarizedExperiment::assay(mat)
        } else {
            mat = SummarizedExperiment::assay(mat, assay)
        }
    }
    
    if (!is.null(percent) & is.null(n)) {
        sel = rowSums(!is.na(mat))/ncol(mat) >= percent
        mat.filtered <- mat.orig[sel, ]
    } else if (is.null(percent) & !is.null(n)) {
        sel = rowSums(!is.na(mat)) >= n
        mat.filtered <- mat.orig[sel, ]
    } else {
        sel = (rowSums(!is.na(mat))/ncol(mat) >=
                percent) & (rowSums(!is.na(mat)) >= n)
        mat.filtered <- mat.orig[sel, ]
    }
    
    
    return(mat.filtered)
}


#' @title Select phosphosites by localisation score
#'
#' @description Select phosphosites with a localisation score higher than the 
#' pre-defined probability score (default score = 0.75)
#'
#' @usage selectLocalisedSites(mat, loc=NULL, prob = 0.75)
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows corresponding to 
#' phosphosites and columns corresponding to samples in replicates for different 
#' treatments.
#' @param loc a vector of localisation scores
#' @param prob a percent from 0 to 1, specifying the localisation probability 
#' of quantified values in across all samples for retaining a phosphosite for 
#' subsequent analysis.
#'
#' @return a filtered matrix
#'
#' @examples
#'
#' data('phospho.cells.Ins.pe')
#' ppe <- phospho.cells.Ins.pe
#' ppe_mat <- as.data.frame(SummarizedExperiment::assay(ppe))
#' # Before filtering
#' dim(ppe)
#' dim(ppe_mat)
#' 
#' # Generate arbitrary localisation probabilities for each phosphosite
#' set.seed(2020)
#' localisation_scores <- round(rnorm(nrow(ppe), 0.8, 0.05), 2)
#' table(localisation_scores >= 0.75)
#' 
#' # Filter
#' Localisation(ppe) <- localisation_scores
#' ppe_filtered <- selectLocalisedSites(ppe, prob=0.75)
#' ppe_mat_filtered <- selectLocalisedSites(ppe_mat, loc=localisation_scores, 
#'      prob=0.75)
#' 
#' # After filtering
#' dim(ppe_filtered)
#' dim(ppe_mat_filtered)
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' 
#' @export
#'
#'
 
 selectLocalisedSites <- function(mat, loc=NULL, prob = 0.75) { 
    
    if (missing(mat))
        stop("Parameter mat is missing!")
    if ((!is.null(prob)) && ((prob < 0) || (prob > 1)))
        stop("Parameter prob must be a numeric value between 0 and 1")
    if ((is.null(loc)) && is.null(Localisation(mat)) && 
        is.null(Localisation(mat)) && sum(is.na(Localisation(mat))) > 0)
         stop("Some or all localisation scores are missing")
    if ((!is.null(loc)) && (length(loc) != nrow(mat)))
         stop("Length of loc should equal to the number of rows in mat")
     
    mat.orig = mat
    if (methods::is(mat, "PhosphoExperiment")) {
            loc = Localisation(mat)
    } else {
        if (is.null(loc))
            stop("Localisation probabilty scores missing")
    }
    
    sel = loc >= prob
    mat.filtered <- mat.orig[sel,]
    
    
    return(mat.filtered)
}
