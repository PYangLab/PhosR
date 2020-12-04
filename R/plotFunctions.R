#' Plot signalome map
#'
#' @usage plotSignalomeMap(signalomes, color)
#'
#' @param signalomes output from `Signalomes` function
#' @param color a string specifying the color vector for kinases.
#'
#' @return a ggplot object
#'
#' @examples
#' data('phospho.L6.ratio.pe')
#' data('SPSs')
#' data('PhosphoSitePlus')
#'
#' grps = gsub('_.+', '', colnames(phospho.L6.ratio.pe))
#' 
#' # Construct a design matrix by condition
#' design = model.matrix(~ grps - 1)
#' 
#' # phosphoproteomics data normalisation using RUV
#' L6.sites = paste(sapply(phospho.L6.ratio.pe@GeneSymbol, function(x)paste(x)),
#'                  ";",
#'                  sapply(phospho.L6.ratio.pe@Residue, function(x)paste(x)),
#'                  sapply(phospho.L6.ratio.pe@Site, function(x)paste(x)),
#'                  ";", sep = "")
#' ctl = which(L6.sites %in% SPSs)
#' phospho.L6.ratio.pe = RUVphospho(phospho.L6.ratio.pe,
#'                                  M = design, k = 3,ctl = ctl)
#' 
#' phosphoL6 = phospho.L6.ratio.pe@assays@data$normalised
#' phosphoL6.mean <- meanAbundance(phosphoL6, grps = grps)
#'
#' @importFrom ggplot2 ggplot
#' @importFrom tidyr spread
#' @importFrom reshape2 melt
#' @importFrom dplyr count
#' 
#' 
#'
#' @export
plotSignalomeMap <- function(signalomes, color) {
    
    dftoPlot_signalome <- stack(signalomes$kinaseSubstrates)
    modules <- signalomes$proteinModule
    names(modules) <- sapply(strsplit(as.character(names(signalomes$proteinModules)), ";"), "[[", 1)
    dftoPlot_signalome$cluster <- modules[dftoPlot_signalome$values]
    
    dftoPlot_balloon_bycluster <- dftoPlot_signalome
    dftoPlot_balloon_bycluster <- na.omit(dftoPlot_balloon_bycluster) %>% dplyr::count(cluster, ind)
    dftoPlot_balloon_bycluster$ind <- as.factor(dftoPlot_balloon_bycluster$ind)
    dftoPlot_balloon_bycluster$cluster <- as.factor(dftoPlot_balloon_bycluster$cluster)
    dftoPlot_balloon_bycluster <- tidyr::spread(dftoPlot_balloon_bycluster, ind, n)[,-1]
    dftoPlot_balloon_bycluster[is.na(dftoPlot_balloon_bycluster)] <- 0
    
    dftoPlot_balloon_bycluster <- do.call(rbind, lapply(1:nrow(dftoPlot_balloon_bycluster), function(x) {
        
        res <- sapply(dftoPlot_balloon_bycluster[x,], function(y) y/sum(dftoPlot_balloon_bycluster[x,])*100)
        
    }))
    
    dftoPlot_balloon_bycluster <- reshape2::melt(as.matrix(dftoPlot_balloon_bycluster))
    colnames(dftoPlot_balloon_bycluster) <- c("cluster", "ind", "n")
    
    ggplot2::ggplot(dftoPlot_balloon_bycluster, aes(x = ind, y = cluster)) + 
        geom_point(aes(col=ind, size=n)) + 
        scale_color_manual(values=color) + 
        scale_size_continuous(range = c(2, 17)) + 
        theme_classic() + 
        theme(
            aspect.ratio=0.25, 
            legend.position = "bottom",
            axis.line = element_blank(),            
            axis.title = element_blank(),   
            panel.grid.major.x = element_blank(),  
            panel.grid.minor.x = element_blank()) 
    
}


#' Plot kinase network map
#'
#' @usage plotKinaseNetwork(KSR, predMatrix, threshold=0.9, color)
#'
#' @param KSR kinase-substrate relationship scoring results
#' @param predMatrix output of kinaseSubstratePred function
#' @param threshold threshold used to select interconnected kinases for
#'  the expanded signalomes
#'  
#' @return a graphical plot
#' 
#' @importFrom network network
#' @importFrom reshape2 melt
#' @importFrom stats cor
#' @importFrom GGally ggnet2
#' 
#' @examples
#' data('phospho.L6.ratio.pe')
#' data('SPSs')
#' data('PhosphoSitePlus')
#'
#' grps = gsub('_.+', '', colnames(phospho.L6.ratio.pe))
#' 
#' # Construct a design matrix by condition
#' design = model.matrix(~ grps - 1)
#' 
#' # phosphoproteomics data normalisation using RUV
#' L6.sites = paste(sapply(phospho.L6.ratio.pe@GeneSymbol, function(x)paste(x)),
#'                  ";",
#'                  sapply(phospho.L6.ratio.pe@Residue, function(x)paste(x)),
#'                  sapply(phospho.L6.ratio.pe@Site, function(x)paste(x)),
#'                  ";", sep = "")
#' ctl = which(L6.sites %in% SPSs)
#' phospho.L6.ratio.pe = RUVphospho(phospho.L6.ratio.pe,
#'                                  M = design, k = 3,ctl = ctl)
#' 
#' phosphoL6 = phospho.L6.ratio.pe@assays@data$normalised
#' phosphoL6.mean <- meanAbundance(phosphoL6, grps = grps)
#'
#' @export
plotKinaseNetwork <- function(KSR, predMatrix, threshold = 0.9, color) {
    
    threskinaseNetwork = threshold
    signalomeKinase <- colnames(predMatrix)
    kinase_cor <- stats::cor(KSR$combinedScoreMatrix)
    
    cor_kinase_mat <- kinase_cor
    diag(cor_kinase_mat) <- 0
    kinase_network <- lapply(1:ncol(cor_kinase_mat), function(x) names(which(cor_kinase_mat[,x] > threskinaseNetwork))) 
    names(kinase_network) <- colnames(cor_kinase_mat)
    
    cor_kinase_mat <- apply(cor_kinase_mat, 2, function(x) x  > threskinaseNetwork)
    cor_kinase_mat[cor_kinase_mat == F] <- 0
    cor_kinase_mat[cor_kinase_mat == T] <- 1
    
    links <- reshape2::melt(cor_kinase_mat)
    links <- links[links$value == 1,]
    res <- sapply(1:length(links$Var1), function(x) { kinase_cor[rownames(kinase_cor) == links$Var1[x], colnames(kinase_cor) == links$Var2[x]]
    })
    links$cor <- res
    colnames(links) <- c("source", "target", "binary", "cor")
    
    network <- network::network(cor_kinase_mat, directed=F)
    GGally::ggnet2(network, 
                   node.size=10, 
                   node.color=color, 
                   edge.size = 0.5, 
                   size = "degree",
                   size.cut=3,
                   label=colnames(cor_kinase_mat),
                   label.size=2,
                   mode="circle",
                   label.color="black")
}