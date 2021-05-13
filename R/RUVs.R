#' @title RUV for phosphoproteomics data normalisation
#'
#' @description This is a wrapper implementation of RUVIII for phosphoproteomics
#' data normalisation. This function will call
#' tailImpute function to impute all the missing values (if there is any) in the
#' phosphoproteomics data for
#' applying RUVIII. It will then return the normalised values for quantified
#' phosphosites and remove imputed values.
#'
#' @param mat a matrix (or PhosphoExperiment object) with rows correspond to 
#' phosphosites and columns correspond to samples.
#' @param M is the design matrix as defined in RUVIII.
#' @param ctl is the stable phosphosites (or negative controls as defined in
#' RUVIII).
#' @param k is the number of unwanted factors as defined in RUVIII.
#' @param m a numeric number for controlling mean downshifting.
#' @param s a numeric number for controlling standard deviation of downshifted
#' sampling values.
#' @param keepImpute a boolean to keep the missing value in the returned matrix.
#' @param ... additional parameters that may be passed to RUVIII.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return A normalised matrix.
#'
#' @examples
#'
#' data('phospho_L6_ratio_pe')
#' data('SPSs')
#'
#' grps = gsub('_.+', '', colnames(phospho.L6.ratio.pe))
#'
#' L6.sites = paste(sapply(GeneSymbol(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  ";",
#'                  sapply(Residue(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  sapply(Site(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  ";", sep = "")
#'
#' # Construct a design matrix by condition
#' design = model.matrix(~ grps - 1)
#'
#' # phosphoproteomics data normalisation using RUV
#' ctl = which(L6.sites %in% SPSs)
#' phospho.L6.ratio.RUV = RUVphospho(
#'     SummarizedExperiment::assay(phospho.L6.ratio.pe, "Quantification"), 
#'     M = design, k = 3, ctl = ctl)
#' 
#' @importFrom ruv RUVIII
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay
#'
#' @aliases RUVproteome
#'
#' @export
RUVphospho <- function(mat, M, ctl, k = NULL, m = 1.6, s = 0.6,
    keepImpute = FALSE, assay = NULL, ...) {
    if (missing(mat))
        stop("Parameter mat is missing!")
    if (missing(M))
        stop("Parameter M is missing!")
    if (missing(ctl))
        stop("Parameter ctl is missing!")
    
    pe = FALSE
    if (methods::is(mat, "PhosphoExperiment")) {
        pe = TRUE
        mat.orig = mat
        if (is.null(assay)) {
            mat = as.matrix(SummarizedExperiment::assay(mat))
        } else {
            mat = as.matrix(SummarizedExperiment::assay(mat, assay))
        }
    }
    mat = RUV(mat = mat, M = M, ctl = ctl, k = k, m = m, s = s,
        keepImpute = keepImpute, ...)
    if (pe) {
        SummarizedExperiment::assay(mat.orig,"normalised",withDimnames=FALSE) = 
            mat
        mat = mat.orig
    }
    mat
}


#' @importFrom ruv RUVIII
#' @export RUVproteome
RUVproteome <- function(mat, M, ctl, k = NULL, m = 1.8, s = 0.3,
    keepImpute = FALSE, assay = NULL, ...) {
    pe = FALSE
    if (methods::is(mat, "PhosphoExperiment")) {
        mat.orig = mat
        if (is.null(assay)) {
            mat = as.matrix(SummarizedExperiment::assay(mat))
        } else {
            mat = as.matrix(SummarizedExperiment::assay(mat, assay))
        }
        pe = TRUE
    }
    
    mat = RUV(mat = mat, M = M, ctl = ctl, k = k, m = m, s = s,
        keepImpute = keepImpute, ...)
    
    if (pe) {
        SummarizedExperiment::assay(mat.orig,"normalised",withDimnames=FALSE) = 
            mat
        mat = mat.orig
    }
    mat
}

# An internal wrapper function for calling RUVIII
RUV <- function(mat, M, ctl, k = NULL, m = m, s = s, keepImpute = keepImpute,
    ...) {
    mat.complete <- c()
    if (sum(is.na(mat)) > 0) {
        mat.complete <- t(tImpute(mat, m = m, s = s))
    } else {
        mat.complete <- t(mat)
    }

    mat.ruv <- t(RUVIII(Y = mat.complete, M = M, k = k, ctl = ctl,
        return.info = FALSE, ...))

    if ((sum(is.na(mat)) > 0) & (keepImpute == FALSE)) {
        mat.ruv[which(is.na(mat))] <- NA
        return(mat.ruv)
    } else {
        return(mat.ruv)
    }
}
