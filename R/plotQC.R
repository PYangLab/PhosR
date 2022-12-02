#' @title A set of function for data QC plot
#'
#' @usage plotQC(mat, grps, labels, panel = 
#' c("quantify", "dendrogram", "abundance", "pca", "all"))
#'
#' @param mat A p by n matrix, where p is the number of phosphosites and n is
#' the number of samples.
#' @param grps A vector of colours to be used in the plot. The length should be
#' equal to the columns of the mat.
#' @param labels A vector of sample names. Used the label points in PCA plot
#' (panel=4)
#' @param panel A type of plot to output. See description for details.
#'
#' @description
#' The `panel` parameter allows different type of visualisation for output
#' object from PhosR.
#' `panel = "all"` is used to create a 2*2 panel of plots including the 
#' following.
#' `panel = "quantify"` is used to visualise percentage of quantification after
#' imputataion.
#' `panel = "dendrogram"` is used to visualise dendrogram (hierarchical 
#' clustering) of the input matrix.
#' `panel = "abundance"` is used to visualise abundance level of samples from 
#' the input matrix.
#' `panel = "pca"` is used to show PCA plot
#'
#' @return A graphical plot
#'
#' @importFrom dendextend labels_colors
#' @importFrom pcaMethods pca
#' @importFrom ggpubr ggarrange
#' @importFrom ggdendro ggdendrogram
#' @importFrom grDevices rainbow
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
#' p1 = plotQC(phospho.cells.Ins.filtered,
#'         labels=colnames(phospho.cells.Ins.filtered),
#'         panel = "quantify", grps = grps)
#' p2 = plotQC(phospho.cells.Ins.ms,
#'         labels=colnames(phospho.cells.Ins.ms),
#'         panel = "quantify", grps = grps)
#' ggpubr::ggarrange(p1, p2, nrow = 1)
#' 
#' # Batch correction
#' data('phospho_L6_ratio_pe')
#' data('SPSs')
#' 
#' grps = gsub('_.+', '', rownames(
#'     SummarizedExperiment::colData(phospho.L6.ratio.pe))
#' )
#' 
#' # Cleaning phosphosite label
#' L6.sites = paste(sapply(GeneSymbol(phospho.L6.ratio.pe),function(x)paste(x)),
#'                  ";",
#'                  sapply(Residue(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  sapply(Site(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  ";", sep = "")
#' phospho.L6.ratio = t(sapply(split(data.frame(
#'     SummarizedExperiment::assay(phospho.L6.ratio.pe, "Quantification")), 
#'     L6.sites),colMeans))
#' phospho.site.names = split(
#'     rownames(
#'         SummarizedExperiment::assay(phospho.L6.ratio.pe, "Quantification")
#'     ), L6.sites)
#' 
#' # Construct a design matrix by condition
#' design = model.matrix(~ grps - 1)
#' 
#' # phosphoproteomics data normalisation using RUV
#' ctl = which(rownames(phospho.L6.ratio) %in% SPSs)
#' phospho.L6.ratio.RUV = RUVphospho(phospho.L6.ratio, M = design, k = 3,
#'                                   ctl = ctl)
#'                                   
#' # plot after batch correction
#' p1 = plotQC(phospho.L6.ratio, panel = "dendrogram", grps=grps,
#'          labels = colnames(phospho.L6.ratio))
#' p2 = plotQC(phospho.L6.ratio.RUV, grps=grps,
#'         labels = colnames(phospho.L6.ratio),
#'         panel="dendrogram")
#' ggpubr::ggarrange(p1, p2, nrow = 1)
#'
#' p1 = plotQC(phospho.L6.ratio, panel = "pca", grps=grps,
#'         labels = colnames(phospho.L6.ratio)) +
#'         ggplot2::ggtitle('Before Batch correction')
#' p2 = plotQC(phospho.L6.ratio.RUV, grps=grps,
#'         labels = colnames(phospho.L6.ratio),
#'         panel="pca") +
#'         ggplot2::ggtitle('After Batch correction')
#' ggpubr::ggarrange(p1, p2, nrow = 1)
#'
#' @export
#'
plotQC <- function(mat, grps, labels, panel = 
        c("quantify", "dendrogram", "abundance", "pca", "all")) {
    if (missing(mat))
        stop("Paramter mat is missing!")
    p = NULL
    if (panel == "quantify") {
        p = quantPlot(mat, grps, labels)
    } else if (panel == "dendrogram") {
        p = dendPlot(mat, grps, labels)
    } else if (panel == "abundance") {
        p = abundPlot(mat, grps, labels)
    } else if (panel == "pca") {
        p = pcaPlot(mat, grps, labels)
    }
    if (panel == "all") {
        p1 = quantPlot(mat, grps, labels)
        p2 = dendPlot(mat, grps, labels)
        p3 = abundPlot(mat, grps, labels)
        p4 = pcaPlot(mat, grps, labels) 
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
quantPlot = function(mat, grps, labels) {
    quant = (1-colSums(is.na(mat))/nrow(mat))*100
    dat = data.frame(
        Quantification = quant,
        Sample = labels,
        Groups = grps
    )
    gNum = length(table(grps))
    ggplot2::ggplot(dat, ggplot2::aes(x = .data$Sample, 
        y = .data$Quantification, fill = .data$Groups)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::coord_cartesian(ylim = c(0,100)) + 
        ggplot2::ggtitle("Quantification per sample") +
        ggplot2::labs(y = "Quantification (%)") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                        vjust = 1, hjust = 1)) +
        ggplot2::scale_fill_manual(values = grDevices::rainbow(gNum))
}

#' @importFrom ggdendro ggdendrogram
#' @importFrom ggplot2 ggtitle
dendPlot = function(mat, grps, labels) {
    dend <- stats::hclust(stats::dist(t(mat)))
    label_grps <- grps[stats::order.dendrogram(as.dendrogram(dend))]
    
    dendr = ggdendro::dendro_data(dend, type = "rectangle")
    
    nGrps = length(table(label_grps))
    label_colors = grDevices::rainbow(nGrps)[factor(label_grps)]
    
    ggdendro::ggdendrogram(dend) +
        ggplot2::ggtitle("Sample hierarchical clustering") +
        ggplot2::theme(
            axis.text.x = ggtext::element_markdown(color = label_colors)
        )
    
}


#' @importFrom dplyr %>% mutate
#' @importFrom ggplot2 ggplot geom_boxplot labs aes
abundPlot = function(mat, grps, labels) {
    rep_num = nrow(mat)
    dat = mat %>%
        as.data.frame() %>%
        dplyr::mutate(sites = rownames(.)) %>%
        tidyr::pivot_longer(-.data$sites, names_to = "Samples", 
            values_to = "abundance") %>%
        dplyr::mutate(
            Groups = rep(grps, rep_num)
        )
    nGrps = length(table(dat$Groups))
    ggplot2::ggplot(dat, ggplot2::aes(x = .data$Samples, y = .data$abundance, 
        fill = .data$Groups)) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(y = "Expression/Abundance level") +
        ggplot2::scale_fill_manual(values = grDevices::rainbow(nGrps))
}

#' @importFrom ggplot2 ggplot geom_point geom_text labs aes
pcaPlot = function(mat, grps, labels) {
    result <- pcaMethods::pca(t(mat), method = "ppca", nPcs = 2,
        seed = 123, main = "PCA")
    
    dat = data.frame(
        PC1 = pcaMethods::scores(result)[, 1],
        PC2 = pcaMethods::scores(result)[, 2],
        grps = grps,
        Samples = labels
    )
    nGrps = length(table(grps))
    
    ggplot2::ggplot(dat, ggplot2::aes(x = .data$PC1, y = .data$PC2, 
        color = .data$grps, label = .data$Samples)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_text(hjust = 0, vjust = 0) +
        ggplot2::labs(
            x = paste("PC1", round(result@R2[1] * 100), "%"),
            y = paste("PC2", round(result@R2[2] * 100), "%")
        ) + 
        ggplot2::scale_colour_manual(values = grDevices::rainbow(nGrps))
}
