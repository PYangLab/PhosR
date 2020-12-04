#' @usage plotSignalomeMap(signalomes, color)
#'
#' @param signalomes output from `Signalomes` function
#' @param color a string specifying the color vector for kinases
#'
#' @return a ggplot object
#'
#' @examples
#'
#' @importFrom ggplot2 ggplot
#' @importFrom tidyr spread
#' @importFrom reshape2 melt
#' @importFrom dplyr count
#' 
#' @export
plotSignalomeMap <- function(signalomes, color) {
    
    df <- stack(signalomes$kinaseSubstrates)
    modules <- signalomes$proteinModule
    names(modules) <- sapply(strsplit(as.character(names(signalomes$proteinModules)), ";"), "[[", 1)
    df$cluster <- modules[df$values]
    
    df_balloon <- df
    df_balloon <- na.omit(df_balloon) %>% dplyr::count(cluster, ind)
    df_balloon$ind <- as.factor(df_balloon$ind)
    df_balloon$cluster <- as.factor(df_balloon$cluster)
    df_balloon <- tidyr::spread(df_balloon, ind, n)[,-1]
    df_balloon[is.na(df_balloon)] <- 0

    df_balloon <- do.call(rbind, lapply(1:nrow(df_balloon), function(x) {
        
        res <- sapply(df_balloon[x,], function(y) y/sum(df_balloon[x,])*100)
        
    }))
    
    df_balloon <- reshape2::melt(as.matrix(df_balloon))
    colnames(df_balloon) <- c("cluster", "ind", "n")
    
    g <- ggplot2::ggplot(df_balloon, aes(x = ind, y = cluster)) + 
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
    g
    
}

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
#' 
#' @export
plotKinaseNetwork <- function(KSR, predMatrix, threshold = 0.9, color) {
    
    cor_mat <- stats::cor(KSR$combinedScoreMatrix)
    
    diag(cor_mat) <- 0
    kinase_network <- lapply(1:ncol(cor_mat), function(x) names(which(cor_mat[,x] > threshold))) 
    names(kinase_network) <- colnames(cor_mat)
    
    cor_mat <- apply(cor_mat, 2, function(x) x  > threshold)
    cor_mat[cor_mat == F] <- 0
    cor_mat[cor_mat == T] <- 1
    
    links <- reshape2::melt(cor_mat)
    links <- links[links$value == 1,]
    res <- sapply(1:length(links$Var1), function(x) { cor_mat[rownames(cor_mat) == links$Var1[x], colnames(cor_mat) == links$Var2[x]]
    })
    links$cor <- res
    colnames(links) <- c("source", "target", "binary", "cor")
    
    g <- GGally::ggnet2(network::network(cor_mat, directed=F), 
                   node.size=10, 
                   node.color=color, 
                   edge.size = 0.5, 
                   size = "degree",
                   size.cut=3,
                   label=colnames(cor_mat),
                   label.size=2,
                   mode="circle",
                   label.color="black")
    g
}