#' Plot signalome map
#'
#' @usage plotSignalomeMap(signalomes, color)
#'
#' @param signalomes output from `Signalomes` function
#' @param color a string specifying the color vector for kinases
#'
#' @return a ggplot object
#'
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
    names(modules) <- unlist(lapply(
        strsplit(
            as.character(names(signalomes$proteinModules)), 
            ";"
        ), "[[", 1))
    df$cluster <- modules[df$values]
    
    df_balloon <- df
    df_balloon <- na.omit(df_balloon) %>% dplyr::count(.data$cluster, .data$ind)
    df_balloon$ind <- as.factor(df_balloon$ind)
    df_balloon$cluster <- as.factor(df_balloon$cluster)
    df_balloon <- tidyr::spread(df_balloon, .data$ind, .data$n)[,-1]
    df_balloon[is.na(df_balloon)] <- 0

    df_balloon <- do.call(rbind, lapply(seq(nrow(df_balloon)), function(x) {
        
        res <- unlist(lapply(df_balloon[x,], 
            function(y) y/sum(df_balloon[x,])*100))
        
    }))
    
    df_balloon <- reshape2::melt(as.matrix(df_balloon))
    colnames(df_balloon) <- c("cluster", "ind", "n")
    
    g <- ggplot2::ggplot(df_balloon, aes(x = .data$ind, y = .data$cluster)) + 
        geom_point(aes(col=.data$ind, size=.data$n)) + 
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

#' Plot kinase network
#'
#' @usage plotKinaseNetwork(KSR, predMatrix, threshold = 0.9, color, 
#' type = NULL, verbose = FALSE)
#'
#' @param KSR Kinase-substrate relationship scoring results
#' @param predMatrix Output of kinaseSubstratePred function
#' @param threshold Threshold used to select interconnected kinases for
#'  the expanded signalomes
#' @param color A string specifying the color vector for nodes
#' @param type A type (\code{graph} or \code{chord}) of plot. If NULL, network 
#' graph is plotted
#' @param verbose Default to \code{TRUE} to show messages during the progress.
#' All messages will be suppressed if set to \code{FALSE}
#'  
#' @return a graphical plot
#' 
#' @importFrom network network
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 melt
#' @importFrom stats cor
#' @importFrom GGally ggnet2
#' @import circlize
#' 
#' 
#' @export
plotKinaseNetwork <- function(KSR, predMatrix, threshold = 0.9, color=NULL, 
    type = NULL, verbose = FALSE) {
    
    if (is.null(type)) { type = "chord" } else {
        type <- match.arg(type, c("graph", "chord"),
                          several.ok = FALSE)
    }

    cor_mat <- stats::cor(KSR$combinedScoreMatrix)
    diag(cor_mat) <- 0

    cor_mat_values = cor_mat

    cor_mat <- apply(cor_mat, 2, function(x) x  > threshold)
    cor_mat[cor_mat == FALSE] <- 0
    cor_mat[cor_mat == TRUE] <- 1
    
    if (type == "graph") {
        if (is.null(color)) {
            stop("Parameter color cannot be NULL.")
        }
        if (verbose) {
            print("Generating network graph...")
        }
        
        n = network::network(cor_mat, directed=FALSE)
        network::set.edge.value(n, "cor", cor_mat_values^2)
        my_col = RColorBrewer::brewer.pal(11, "Spectral")
        
        g <- GGally::ggnet2(n, 
                       node.size=10, 
                       node.color=color, 
                       edge.size = "cor", 
                       size = "degree",
                       size.cut=3,
                       label=colnames(cor_mat),
                       label.size=2,
                       mode="circle",
                       label.color="black")
        print(g)
    
    }
    
    if (type == "chord") {
        
        if (verbose) {
            print("Generating circular graph...")
        }
        if (is.null(color)) {
            stop("Parameter color cannot be NULL.")
        }
        grid.col <- color
        names(grid.col) <- rownames(cor_mat)
        
        n = length(grid.col)
        
        circos.clear()
        circos.par(start.degree = 180)
        circos.initialize(factors = "a", xlim = c(0, n))
        chordDiagram(cor_mat, transparency = 0.2,
                     order = rownames(cor_mat),
                     grid.col = grid.col,
                     annotationTrack = c("name", "grid"), scale = TRUE)
        title("Kinase network")

    }
}
