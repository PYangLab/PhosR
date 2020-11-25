#' @title PhosR Signalomes
#'
#' @description A function to generate signalomes
#'
#' @usage Signalomes(KSR, predMatrix, exprsMat, KOI, threskinaseNetwork=0.9,
#' signalomeCutoff=0.5, verbose = TRUE)
#'
#' @param KSR kinase-substrate relationship scoring results
#' @param predMatrix output of kinaseSubstratePred function
#' @param exprsMat a matrix with rows corresponding to phosphosites and columns
#' corresponding to samples
#' @param PhosphoExperiment PhosphoExperiment object containing a matrix with rows corresponding to phosphosites and columns
#' corresponding to samples
#' @param KOI a character vector that contains kinases of interest for which
#' expanded signalomes will be generated
#' @param threskinaseNetwork threshold used to select interconnected kinases for
#'  the expanded signalomes
#' @param signalomeCutoff threshold used to filter kinase-substrate
#' relationships
#' @param verbose Default to \code{TRUE} to show messages during the progress.
#' All messages will be suppressed if set to \code{FALSE}
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr count
#' @importFrom stats hclust
#' @importFrom rlang .data
#'
#' @return A list of 3 elements.
#'  \code{Signalomes}, \code{proteinModules} and \code{kinaseSubstrates}
#'
#' @examples
#'
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
#' phosphoL6 = phospho.L6.ratio.RUV
#' rownames(phosphoL6) = phospho.site.names
#'
#' # filter for up-regulated phosphosites
#' phosphoL6.mean <- meanAbundance(phosphoL6, grps = gsub('_.+', '',
#'                                 colnames(phosphoL6)))
#' aov <- matANOVA(mat=phosphoL6, grps=gsub('_.+', '', colnames(phosphoL6)))
#' phosphoL6.reg <- phosphoL6[(aov < 0.05) &
#'                         (rowSums(phosphoL6.mean > 0.5) > 0),,drop = FALSE]
#' L6.phos.std <- standardise(phosphoL6.reg)
#' rownames(L6.phos.std) <- sapply(strsplit(rownames(L6.phos.std), '~'),
#'     function(x){gsub(' ', '', paste(toupper(x[2]), x[3], '', sep=';'))})
#'
#' L6.phos.seq <- sapply(strsplit(rownames(phosphoL6.reg), '~'),
#'                     function(x)x[4])
#'
#' L6.matrices <- kinaseSubstrateScore(PhosphoSite.mouse, L6.phos.std,
#'     L6.phos.seq, numMotif = 5, numSub = 1)
#' set.seed(1)
#' L6.predMat <- kinaseSubstratePred(L6.matrices, top=30)
#'
#' kinaseOI = c('PRKAA1', 'AKT1')
#'
#' Signalomes_results <- Signalomes(KSR=L6.matrices,
#'                                 predMatrix=L6.predMat,
#'                                 exprsMat=L6.phos.std,
#'                                 KOI=kinaseOI)
#' @export

Signalomes <- function(KSR, predMatrix, exprsMat, KOI, threskinaseNetwork = 0.9,
    signalomeCutoff = 0.5, verbose = TRUE) {
    ############## generate objects required for signalome function
    protein_assignment = mapply("[[",
                                strsplit(rownames(KSR$combinedScoreMatrix),";"),
                                MoreArgs = list(1))
    KinaseFamily = PhosR::KinaseFamily
    kinaseGroup <- KinaseFamily[, "kinase_group"]
    names(kinaseGroup) <- KinaseFamily[, "gene_symbol"]
    ############## set color palette
    my_color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,
        "Accent"))
    kinase_all_color <- my_color_palette(ncol(KSR$combinedScoreMatrix))
    names(kinase_all_color) <- colnames(KSR$combinedScoreMatrix)
    kinase_signalome_color <- kinase_all_color[colnames(predMatrix)]
    my_color_palette_kinaseGroup <-
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "Set2"))
    kinaseGroup_color <-
        my_color_palette_kinaseGroup(length(unique(kinaseGroup)))
    names(kinaseGroup_color) <- unique(kinaseGroup)
    ############## interconnected kinases
    resKinaseNetwork <- .kinaseNetwork(predMatrix, KSR, threskinaseNetwork,
        kinase_signalome_color)
    ############## cluster phosphosites
    substrate_clusters <- .phosphositeClusters(KSR, verbose)
    ############## generate coassignment
    cluster_assignment <- as.factor(substrate_clusters)
    dat.long <- data.frame(table(cluster_assignment, protein_assignment))
    dftoHeatmap <- tidyr::pivot_wider(dat.long, names_from = protein_assignment, 
        values_from = Freq)[,-1]
    dftoHeatmap[is.na(dftoHeatmap)] <- 0
    dftoHeatmap[dftoHeatmap > 0] <- 1
    hclust_res <- stats::hclust(stats::dist(t(dftoHeatmap)),
        method = "ward.D")
    tree_height <- as.numeric(names(table(hclust_res$height)))
    branching <- as.numeric(table(hclust_res$height))
    hcutree <- min(tree_height[tree_height > 0])
    modules <- stats::cutree(hclust_res, h = hcutree)
    ############## generate signalomes
    signalomeSubstrates <- .phosRsignalome(predMatrix, signalomeCutoff,
        kinase_signalome_color, modules)
    ############## generate kinase-specific signalomes (extended signalome)
    signalomes_of_KOI <- .getSignalomes(predMatrix, exprsMat,
        KOI, resKinaseNetwork, signalomeSubstrates, modules,
        kinaseGroup, verbose = verbose)
    signalome_res <- list(Signalomes = signalomes_of_KOI,
                        proteinModules = modules,
                        kinaseSubstrates = signalomeSubstrates)
    return(signalome_res)
}

.kinaseNetwork <- function(predMatrix, KSR, threskinaseNetwork,
    kinase_signalome_color) {


    signalomeKinase <- colnames(predMatrix)
    kinase_cor <- stats::cor(KSR$combinedScoreMatrix)

    cor_kinase_mat <- kinase_cor
    diag(cor_kinase_mat) <- 0
    kinase_network <- lapply(seq_len(ncol(cor_kinase_mat)), function(x)
        names(which(cor_kinase_mat[, x] > threskinaseNetwork)))
    names(kinase_network) <- colnames(cor_kinase_mat)

    cor_kinase_mat <- apply(cor_kinase_mat, 2,
                            function(x) x > threskinaseNetwork)
    cor_kinase_mat[cor_kinase_mat == FALSE] <- 0
    cor_kinase_mat[cor_kinase_mat == TRUE] <- 1

    network <- igraph::graph_from_adjacency_matrix(cor_kinase_mat,
        mode = "undirected", diag = FALSE)


    kinaseNetwork.res <- list(kinaseNetwork = kinase_network,
        kinaseCor = cor_kinase_mat)

    return(kinaseNetwork.res)

}

.phosphositeClusters <- function(KSR, verbose = TRUE) {
    substrate_cor <- stats::cor(t(KSR$combinedScoreMatrix))
    substrate_hclust <- stats::hclust(stats::dist(KSR$combinedScoreMatrix),
        method = "ward.D")
    if (verbose)
        message("calculating optimal number of clusters...")
    res <- lapply(seq(2,10,1), function(x) {
        substrate_clusters <- stats::cutree(substrate_hclust, k = x)
        cor.res <- lapply(seq_len(x), function(y) {
            substrate_cor = substrate_cor[substrate_clusters == y,
                substrate_clusters == y]
            diag(substrate_cor) <- 0
            return(substrate_cor)
        })
        cor.res <- lapply(cor.res, function(x) median(x))
        cor.res <- unlist(cor.res)
        cor.logic <- sum(cor.res >= 0.5) == x
        if (cor.logic) {
            cluster = mean(cor.res)
            names(cluster) = x
            return(cluster)
        }
    })
    if (isTRUE(is.null(unlist(res)))) {
        res <- lapply(seq(2,10,1), function(x) {
            substrate_clusters <- cutree(substrate_hclust, k = x)
            cor.res <- lapply(seq_len(x), function(y) {
                substrate_cor = substrate_cor[substrate_clusters == y,
                    substrate_clusters == y]
                diag(substrate_cor) <- 0
                return(substrate_cor)
            })
            cor.res <- unlist(lapply(cor.res, function(x) median(x)))
            cor.logic <- sum(cor.res >= 0.1) == x
            if (cor.logic) {
                cluster = median(cor.res)
                names(cluster) = x
                return(cluster)
            }
        })
        res <- as.numeric(names(which(unlist(res) == max(unlist(res)))[1]))
    } else {
        res <- as.numeric(names(which(unlist(res) == max(unlist(res)))[1]))
    }
    if (verbose)
        message(paste0("optimal number of clusters = ", res))
    substrate_clusters <- stats::cutree(substrate_hclust, k = res)
    return(substrate_clusters)
}

#' @import circlize
#' @importFrom utils stack
#' @importFrom rlang .data
.phosRsignalome <- function(predMatrix, signalomeCutoff, kinase_signalome_color,
    modules) {

    signalomeKinase <- colnames(predMatrix)

    signalomeSubstrates <- list()
    for (i in seq_len(length(signalomeKinase))) {
        signalomeSubstrates[[i]] =
            mapply(function(x) x[1],
                strsplit(names(which(
                    predMatrix[,signalomeKinase[[i]]] > signalomeCutoff)),
                    ";"))
    }
    names(signalomeSubstrates) <- signalomeKinase

    ############## generate circlize plot
    dftoPlot_signalome <- stack(signalomeSubstrates)
    dftoPlot_signalome$modules <- modules[dftoPlot_signalome$values]
    # adjacencyData <- with(dftoPlot_signalome, table(ind, modules))
    adjacencyData <- table(dftoPlot_signalome$ind, dftoPlot_signalome$modules)

    grid.col <- c(kinase_signalome_color, rep("grey", length(unique(modules))))
    names(grid.col) <- c(rownames(adjacencyData), as.character(unique(modules)))

    n = length(grid.col)
    circos.par(start.degree = 180)
    circos.initialize(factors = "a", xlim = c(0, n))
    chordDiagram(adjacencyData, transparency = 0.2,
                order = c(rownames(adjacencyData),
                            rev(as.character(unique(modules)))),
                grid.col = grid.col, big.gap = 15,
                annotationTrack = c("name", "grid"), scale = TRUE)
    title("Signalomes")
    circos.clear()

    return(signalomeSubstrates)
}

#' @importFrom dplyr count %>%
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
.getSignalomes <- function(predMatrix, exprsMat, KOI, resKinaseNetwork,
    signalomeSubstrates, modules, kinaseGroup, verbose = TRUE) {
    protein = mapply("[[", strsplit(rownames(exprsMat), ";"),
        MoreArgs = list(1))
    KinaseSubstrateList <- resKinaseNetwork$kinaseNetwork

    ############## annotate kinase-substrate relationship
    kinaseAnnot = annoKinaseSubstrateRelation(predMatrix)

    ############## proportion of regulation by kinase
    signalomeMatrix <- stack(signalomeSubstrates)
    signalomeMatrix$cluster <- modules[signalomeMatrix$values]

    balloon_bycluster <- signalomeMatrix
    balloon_bycluster <- na.omit(balloon_bycluster) %>%
        dplyr::count(.data$cluster, .data$ind)
    balloon_bycluster$ind <- as.factor(balloon_bycluster$ind)
    balloon_bycluster$cluster <- as.factor(balloon_bycluster$cluster)
    balloon_bycluster <- tidyr::pivot_wider(balloon_bycluster, names_from = .data$ind,
        values_from = .data$n)[, -1]
    
    balloon_bycluster[is.na(balloon_bycluster)] <- 0
    balloon_bycluster <- do.call(rbind, lapply(seq_len(nrow(balloon_bycluster)),
        function(x) {
            res = mapply(function(y, balloon_bycluster, x) {
                y/sum(balloon_bycluster[x, ]) * 100
            }, balloon_bycluster[x, ],
            MoreArgs = list(balloon_bycluster = balloon_bycluster, x = x))
        }))
    kinaseProportions <- round(balloon_bycluster, 3)

    ############## generate kinase specific signalomes
    res = generateSignalome(kinaseAnnot, kinaseGroup, predMatrix, KOI,
                            KinaseSubstrateList, kinaseProportions,
                            signalomeSubstrates, exprsMat, protein, modules, 
                            verbose = verbose)
    return(res)
}

annoKinaseSubstrateRelation = function(predMatrix) {
    kinaseAnnot <- lapply(seq_len(nrow(predMatrix)), function(x) {
        scores <- predMatrix[x, ]
        siteName <- rownames(predMatrix)[[x]]
        kinaseTop <- names(which(scores == max(scores)))
        score <- as.numeric(max(scores))
        res <- c(siteName, kinaseTop, score)
        return(res)
    })
    kinaseAnnot <- do.call(rbind, kinaseAnnot)
    rownames(kinaseAnnot) <- kinaseAnnot[, 1]
    kinaseAnnot <- data.frame(kinaseAnnot[, -c(1)], stringsAsFactors = FALSE)
    colnames(kinaseAnnot) <- c("kinase", "score")
    kinaseAnnot$score <- as.numeric(kinaseAnnot$score)
    kinaseAnnot
}

generateSignalome = function(kinaseAnnot, kinaseGroup, predMatrix, KOI,
                            KinaseSubstrateList, kinaseProportions,
                            signalomeSubstrates, exprsMat, protein, modules,
                            verbose = TRUE) {
    annotation <- data.frame(kinase = as.factor(kinaseAnnot$kinase),
        kinaseFamily = as.factor(kinaseGroup[kinaseAnnot$kinase]),
        score = as.numeric(kinaseAnnot$score))
    rownames(annotation) <- rownames(predMatrix)
    res <- lapply(KOI, function(x) {
        if (x %in% names(KinaseSubstrateList)) {
            regModule <- which(kinaseProportions[, x] > 1)
            kinases <- unique(c(x, KinaseSubstrateList[[x]]))
            substrates <- unique(unlist(
                lapply(kinases, function(x) signalomeSubstrates[[x]])))
            exprs_dat <- lapply(regModule, function(y) {
                exprsMat[protein %in%
                        names(modules[modules == y]) &
                        protein %in% substrates, ]
            })
            exprs_dat <- do.call(rbind, exprs_dat)
            annotation_dat <- annotation[rownames(annotation) %in%
                    rownames(exprs_dat), ]
            kinaseSignalome <- list(exprs = exprs_dat,
                                    annotation = annotation_dat)
            return(kinaseSignalome)
        } else {
            if (verbose)
                message(paste0(x, " is not found"))
        }
    })
    names(res) <- KOI
    res
}
