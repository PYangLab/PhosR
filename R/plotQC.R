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
#'
#' @description
#' The `panel` parameter allows different type of visualisation for output
#' object from PhosR.
#' `panel = "all"` is used to create a 2*2 panel of plots including the following.
#' `panel = "quantify"` is used to visualise percentage of quantification after
#' imputataion.
#' `panel = "dendrogram"` is used to visualise dendrogram (hierarchical clustering) of the
#' input matrix.
#' `panel = "abundance"` is used to visualise abundance level of samples from the input
#' matrix.
#' `panel = "pca"` is used to show PCA plot
#'
#' @return A graphical plot
#'
#' @importFrom dendextend labels_colors
#' @importFrom pcaMethods pca
#' @importFrom ggpubr ggarrange
#' @importFrom ggdendro ggdendrogram
#' 
#'
#' @examples
#' # Imputation
#' data('phospho.cells.Ins.sample')
#' grps = gsub('_[0-9]{1}', '', colnames(phospho.cells.Ins))
#' phospho.cells.Ins.filtered <- selectGrps(phospho.cells.Ins, grps, 0.5, n=1)
#'
#' set.seed(123)
#' phospho.cells.Ins.impute <-
#'     scImpute(
#'         phospho.cells.Ins.filtered,
#'         0.5,
#'         grps)[,colnames(phospho.cells.Ins.filtered)]
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
#' p1 = plotQC(phospho.cells.Ins.filtered,
#'         labels=colnames(phospho.cells.Ins.filtered),
#'         panel = "quantify", cols = cols)
#' p2 = plotQC(phospho.cells.Ins.ms,
#'         labels=colnames(phospho.cells.Ins.ms),
#'         panel = "quantify", cols = cols)
#' ggpubr::ggarrange(p1, p2, nrow = 1)
#' 
#' # Batch correction
#' data('phospho.L6.ratio.pe')
#' data('SPSs')
#' 
#' grps = gsub('_.+', '', phospho.L6.ratio.pe@colData@rownames)
#' 
#' # Cleaning phosphosite label
#' L6.sites = paste(sapply(phospho.L6.ratio.pe@GeneSymbol, function(x)paste(x)),
#'                  ";",
#'                  sapply(phospho.L6.ratio.pe@Residue, function(x)paste(x)),
#'                  sapply(phospho.L6.ratio.pe@Site, function(x)paste(x)),
#'                  ";", sep = "")
#' phospho.L6.ratio = t(sapply(split(data.frame(
#'     phospho.L6.ratio.pe@assays@data$Quantification), L6.sites),colMeans))
#' phospho.site.names = split(
#'     rownames(phospho.L6.ratio.pe@assays@data$Quantification), L6.sites)
#' 
#' # Construct a design matrix by condition
#' design = model.matrix(~ grps - 1)
#' 
#' # phosphoproteomics data normalisation using RUV
#' ctl = which(rownames(phospho.L6.ratio) %in% SPSs)
#' phospho.L6.ratio.RUV = RUVphospho(phospho.L6.ratio, M = design, k = 3,
#'                                   ctl = ctl)
#'                                   
#' cs = rainbow(length(unique(grps)))
#' colorCodes = sapply(grps, switch, AICAR=cs[1], Ins=cs[2], AICARIns=cs[3])
#'
#' # plot after batch correction
#' p1 = plotQC(phospho.L6.ratio, panel = "dendrogram", cols=colorCodes,
#'          labels = colnames(phospho.L6.ratio))
#' p2 = plotQC(phospho.L6.ratio.RUV, cols=colorCodes,
#'         labels = colnames(phospho.L6.ratio),
#'         panel="dendrogram")
#' ggpubr::ggarrange(p1, p2, nrow = 1)
#'
#' p1 = plotQC(phospho.L6.ratio, panel = "pca", cols=colorCodes,
#'         labels = colnames(phospho.L6.ratio)) +
#'         ggplot2::ggtitle('Before Batch correction')
#' p2 = plotQC(phospho.L6.ratio.RUV, cols=colorCodes,
#'         labels = colnames(phospho.L6.ratio),
#'         panel="pca") +
#'         ggplot2::ggtitle('After Batch correction')
#' ggpubr::ggarrange(p1, p2, nrow = 1)
#'
#' @export
#'
plotQC <- function(mat, cols, labels, panel = 
        c("quantify", "dendrogram", "abundance", "pca", "all")) {
    if (missing(mat))
        stop("Paramter mat is missing!")
    p = NULL
    if (panel == "quantify") {
        p = quantPlot(mat, cols, labels)
    } else if (panel == "dendrogram") {
        p = dendPlot(mat, cols, labels)
    } else if (panel == "abundance") {
        p = abundPlot(mat, cols, labels)
    } else if (panel == "pca") {
        p = pcaPlot(mat, cols, labels)
    }
    if (panel == "all") {
        p1 = quantPlot(mat, cols, labels)
        p2 = dendPlot(mat, cols, labels)
        p3 = abundPlot(mat, cols, labels)
        p4 = pcaPlot(mat, cols, labels) 
        p = ggpubr::ggarrange(
            p1,
            p2,
            p3,
            p4,
            nrow = 2
        )
    }
    p
}

#' @importFrom ggplot2 ggplot geom_bar coord_cartesian ggtitle labs theme aes 
#' element_text
quantPlot = function(mat, cols, labels) {
    quant = (1-colSums(is.na(mat))/nrow(mat))*100
    dat = data.frame(
        Quantification = quant,
        Sample = labels,
        Groups = cols
    )
    ggplot2::ggplot(dat, ggplot2::aes(x = Sample, y = Quantification, 
        fill = Groups)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::coord_cartesian(ylim = c(0,100)) + 
        ggplot2::ggtitle("Quantification per sample") +
        ggplot2::labs(y = "Quantification (%)") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, 
            hjust = 1))    
}

#' @importFrom ggdendro ggdendrogram
#' @importFrom ggplot2 ggtitle
dendPlot = function(mat, cols, labels) {
    dend <- stats::hclust(stats::dist(t(mat)))
    label_colors <- cols[stats::order.dendrogram(as.dendrogram(dend))]
    
    dendr = ggdendro::dendro_data(dend, type = "rectangle")
        
    ggdendro::ggdendrogram(dend) +
        ggplot2::ggtitle("Sample hierarchical clustering") +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(color = factor(label_colors))
        )
    
}


#' @importFrom dplyr %>% mutate
#' @importFrom ggplot2 ggplot geom_boxplot labs aes
abundPlot = function(mat, cols, labels) {
    rep_num = nrow(mat)
    dat = mat %>%
        as.data.frame() %>%
        dplyr::mutate(sites = rownames(.)) %>%
        tidyr::pivot_longer(-sites, names_to = "Samples", 
            values_to = "abundance") %>%
        dplyr::mutate(
            Groups = rep(cols, rep_num)
        )
    
    ggplot2::ggplot(dat, ggplot2::aes(x = Samples, y = abundance, 
        fill = Groups)) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(y = "Expression/Abundance level")
}

#' @importFrom ggplot2 ggplot geom_point geom_text labs aes
pcaPlot = function(mat, cols, labels) {
    result <- pcaMethods::pca(t(mat), method = "ppca", nPcs = 2,
        seed = 123, main = "PCA")
    
    dat = data.frame(
        PC1 = result@scores[, 1],
        PC2 = result@scores[, 2],
        col = cols,
        Samples = labels
    )
    ggplot2::ggplot(dat, ggplot2::aes(x = PC1, y = PC2, color = col, 
        label = Samples)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_text(hjust = 0, vjust = 0) +
        ggplot2::labs(
            x = paste("PC1", round(result@R2[1] * 100), "%"),
            y = paste("PC2", round(result@R2[2] * 100), "%")
        )
}
