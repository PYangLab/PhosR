#' @title A set of function for data QC plot
#'
#' @usage plotQC(mat, cols, labels, panel, ...)
#'
#' @param mat A p by n matrix, where p is the number of phosphosites and n is
#' the number of samples.
#' @param cols A vector of colours to be used in the plot. The length should be
#' equal to the columns of the mat.
#' @param labels A vector of sample names. Used the label points in PCA plot
#' (panel=4)
#' @param panel A numeric value (0-4) to choose the plot type. See description
#' for details.
#' @param ... Plotting parameters for base plots
#'
#' @description
#' The `panel` parameter allows different type of visualisation for output
#' object from PhosR.
#' `panel = 0` is used to create a 2*2 panel of plots including the following.
#' `panel = 1` is used to visualise percentage of quantification after
#' imputataion.
#' `panel = 2` is used to visualise dendrogram (hierarchical clustering) of the
#' input matrix.
#' `panel = 3` is used to visualise abundance level of samples from the input
#' matrix.
#' `panel = 4` is used to show PCA plot
#'
#' @return A graphical plot
#'
#' @importFrom dendextend labels_colors
#' @importFrom calibrate textxy
#' @importFrom pcaMethods pca
#' @importFrom graphics barplot plot boxplot par title
#'
#' @examples
#'
#' # Imputation
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
#' phospho.cells.Ins.impute[,seq_len(5)] <- ptImpute(
#'     phospho.cells.Ins.impute[,seq(6,10)],
#'     phospho.cells.Ins.impute[,seq(5)], 
#'     percent1 = 0.6, percent2 = 0, paired = FALSE)
#'
#' phospho.cells.Ins.ms <- medianScaling(phospho.cells.Ins.impute,
#'                                     scale = FALSE)
#'
#' cols <- rep(c('#ED4024', '#7FBF42', '#3F61AD', '#9B822F'), each=6)
#' par(mfrow=c(1,2))
#' plotQC(phospho.cells.Ins.filtered,
#'         labels=colnames(phospho.cells.Ins.filtered),
#'         panel = 1, cols = cols)
#' plotQC(phospho.cells.Ins.ms,
#'         labels=colnames(phospho.cells.Ins.ms),
#'         panel = 1, cols = cols)
#'
#' # Batch correction
#' data('phospho_L6_ratio')
#' data('SPSs')
#'
#' grps = gsub('_.+', '', colnames(phospho.L6.ratio))
#'
#' # Cleaning phosphosite label
#' phospho.site.names = rownames(phospho.L6.ratio)
#' L6.sites = gsub(' ', '', sapply(strsplit(rownames(phospho.L6.ratio), '~'),
#'                                 function(x){paste(toupper(x[2]), x[3], '',
#'                                                 sep=';')}))
#' phospho.L6.ratio = t(sapply(split(data.frame(phospho.L6.ratio), L6.sites),
#'                             colMeans))
#' phospho.site.names = split(phospho.site.names, L6.sites)
#'
#' # Construct a design matrix by condition
#' design = model.matrix(~ grps - 1)
#'
#' # phosphoproteomics data normalisation using RUV
#' ctl = which(rownames(phospho.L6.ratio) %in% SPSs)
#' phospho.L6.ratio.RUV = RUVphospho(phospho.L6.ratio, M = design, k = 3,
#'                                 ctl = ctl)
#'
#' cs = rainbow(length(unique(grps)))
#' colorCodes = sapply(grps, switch, AICAR=cs[1], Ins=cs[2], AICARIns=cs[3])
#'
#' # plot after batch correction
#' par(mfrow=c(1,2))
#' plotQC(phospho.L6.ratio, panel = 2, cols=colorCodes)
#' plotQC(phospho.L6.ratio.RUV, cols=colorCodes,
#'         labels = colnames(phospho.L6.ratio),
#'         panel=2, ylim=c(-20, 20), xlim=c(-30, 30))
#'
#' par(mfrow=c(1,2))
#' plotQC(phospho.L6.ratio, panel = 4, cols=colorCodes,
#'         labels = colnames(phospho.L6.ratio),
#'         main='Before Batch correction')
#' plotQC(phospho.L6.ratio.RUV, cols=colorCodes,
#'         labels = colnames(phospho.L6.ratio),
#'         panel=4, ylim=c(-20, 20), xlim=c(-30, 30),
#'         main='After Batch correction')
#'
#' @export
#'
plotQC <- function(mat, cols = NA, labels = NULL, panel = c("quantify", "dendrogram", "abundance", "pca", "all"), ...) {
    if (missing(mat))
        stop("Paramter mat is missing!")
    if (panel == "quantify") {
        quantPlot(mat, cols, ...)
    } else if (panel == "dendrogram") {
        dendPlot(mat, cols, ...)
    } else if (panel == "abundance") {
        abundPlot(mat, cols, ...)
    } else if (panel == "pca") {
        pcaPlot(mat, cols, labels, ...)
    }
    if (panel == "all") {
        graphics::par(mfrow = c(2, 2))
        quantPlot(mat, cols, ...)
        dendPlot(mat, cols, ...)
        abundPlot(mat, cols, ...)
        pcaPlot(mat, cols, labels, ...)
    }
}

quantPlot = function(mat, cols, ...) {
    graphics::barplot((1 - colSums(is.na(mat))/nrow(mat)) *
            100, las = 2, col = cols, ylab = "Quantification (%)",
        main = "Quantification per sample", ylim = c(0, 100), ...)
}

dendPlot = function(mat, cols, ...) {
    dend <- stats::as.dendrogram(stats::hclust(stats::dist(t(mat))))
    dendextend::labels_colors(dend) <- cols[stats::order.dendrogram(dend)]
    graphics::plot(dend, main = "Sample hierarchical clustering",
        ylab = "Tree height", ...)
}

abundPlot = function(mat, cols, ...) {
    graphics::boxplot(mat, las = 2, col = cols,
        ylab = "Expression/Abundance level", ...)
}

pcaPlot = function(mat, cols, labels, ...) {
    result <- pcaMethods::pca(t(mat), method = "ppca", nPcs = 2,
        seed = 123, main = "PCA")
    graphics::plot(result@scores[, 1], result@scores[, 2],
        col = cols, pch = 16, cex = 1.5, xlab = paste("PC1",
            round(result@R2[1] * 100), "%"), ylab = paste("PC2",
                round(result@R2[2] * 100), "%"), ...)
    textxy(result@scores[, 1], result@scores[, 2], labels,
        cex = 0.75)
}
