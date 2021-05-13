#' Generate set of stable phosphoporylated sites
#'
#' @usage getSPS(phosData, assays, conds, num)
#'
#' @param phosData a list of users' PhosphoExperiment objects from which 
#' generate SPSs
#' @param assays an assay to use for each dataset in phosData
#' @param conds a list of vector contains the conditions labels for each sample 
#' in the phosphoExperiment objects
#' @param num the number of identified SPSs, by default is 100
#'
#' @return A vectors of stably phosphorylated sites
#'
#' @examples
#'
#' library(stringr)
#' 
#' data("phospho_L6_ratio_pe")
#' data("phospho.liver.Ins.TC.ratio.RUV.pe")
#' data("phospho.cells.Ins.pe")
#' 
#' ppe1 <- phospho.L6.ratio.pe
#' ppe2 <- phospho.liver.Ins.TC.ratio.RUV.pe
#' ppe3 <- phospho.cells.Ins.pe
#' grp3 = gsub('_[0-9]{1}', '', colnames(ppe3))
#' 
#' cond.list <- list(grp1 = gsub("_.+", "", colnames(ppe1)),
#'                   grp2 = stringr::str_sub(colnames(ppe2), end=-5),
#'                   grp3 = grp3)
#' 
#' ppe3 <- selectGrps(ppe3, grps = grp3, 0.5, n=1)
#' ppe3 <- tImpute(ppe3)
#' 
#' # convert matrix to ratio
#' FL83B.ratio <- SummarizedExperiment::assay(ppe3,"imputed")[, seq(12)] - 
#'     rowMeans(
#'         SummarizedExperiment::assay(ppe3,"imputed")[,grep("FL83B_Control", 
#'         colnames(ppe3))])
#' Hepa.ratio <- SummarizedExperiment::assay(ppe3,"imputed")[, seq(13,24,1)] - 
#'     rowMeans(
#'         SummarizedExperiment::assay(ppe3, "imputed")[,grep("Hepa1.6_Control", 
#'         colnames(ppe3))])
#' SummarizedExperiment::assay(ppe3, "Quantification") <- 
#'     cbind(FL83B.ratio, Hepa.ratio)
#' 
#' ppe.list <- list(ppe1, ppe2, ppe3)
#' 
#' inhouse_SPSs <- getSPS(ppe.list, conds = cond.list)
#' 
#' @export
 
getSPS <-function (phosData, assays="Quantification", conds, num = 100) {
        if (missing(phosData)) 
            stop("phosData is missing")
        if (missing(conds)) 
                stop("conds is missing")
        
        
        sites <- sites.unique <- mat.max <- list()
        n <- length(phosData)
        m <- length(conds)
        if (n < 2) {
            stop("Please use more than one dataset")
        }
        if (n != m) {
            stop("Please use the same number of datasets and conditions")
        }
        for (i in seq(n)) {
            if (!"PhosphoExperiment" %in% is(phosData[[i]])) {
                stop("Wrong phosData, need to be a PhosphoExperiment object")
            }
        }
        for (i in seq(n)) {
            sites[[i]] <- paste(toupper(GeneSymbol(phosData[[i]])), 
                                paste(Residue(phosData[[i]]), 
                                        Site(phosData[[i]]), sep = ""), 
                                sep = ";")
            
            sites.unique[[i]] <- unique(sites[[i]])
            nrep <- ncol(phosData[[i]])/length(unique(conds[[i]]))
            if (nrep == 1) {
                mat.mean <- SummarizedExperiment::assay(phosData[[i]], assays)
            } else {
                grps <- conds[[i]]
                mat.mean <- PhosR::meanAbundance(
                        SummarizedExperiment::assay(phosData[[i]], assays),grps)
            }
            fun_val = (rep(0, ncol(mat.mean)))
            names(fun_val) = colnames(mat.mean)
            sites.mean <- 
                t(
                        vapply(split(as.data.frame(mat.mean), sites[[i]]),
                               colMeans, FUN.VALUE = fun_val
                        )
                )
            sites.max <- apply(sites.mean, 1, function(x) {
                x[which.max(abs(x))]
            })
            mat.max[[i]] <- sort(abs(sites.max), decreasing = TRUE)
        }
        o <- as.data.frame(table(unlist(sites.unique)))
        if (length(which(o$Freq > 1)) < 200) {
            stop("Fewer than 200 overlapped sites")
        }
        if (length(which(o$Freq == m)) > 1000) {
            top <- as.character(o[which(o$Freq == m), 1])
        }
        else {
            message("Warning: there aren't enough overlappling sites")
            top <- as.character(o[which(o$Freq > 1), 1])
        }
        Ts <- data.frame(mat.max[[1]][top])
        for (i in seq(2,n,1)) {
            Ts <- cbind(Ts, mat.max[[i]][top])
        }
        Tc <- (apply(-abs(Ts), 2, rank) - 0.5)/nrow(Ts)
        Tt4 <- pchisq(-2 * rowSums(log(Tc)), (n - 1) * 2, lower.tail = FALSE)
        names(Tt4) <- top
        sites.sorted <- names(sort(Tt4, decreasing = TRUE))
        sites.sorted = paste0(sites.sorted, ";")
        
        return(sites.sorted[seq(num)])
}

