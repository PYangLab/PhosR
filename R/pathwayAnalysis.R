#' phosphosite/Gene set over-representation analysis
#'
#' @description This function performes phosphosite (or gene) set 
#' over-representation analysis using Fisher's exact test.
#'
#' @param geneSet an array of gene or phosphosite IDs (IDs are gene symbols etc
#' that match to your pathway annotation list).
#' @param annotation a list of pathways with each element containing an array of
#' gene or phosphosite IDs.
#' @param universe the universe/backgrond of all genes or phosphosites in your
#' profiled dataset.
#' @param alter test for enrichment ('greater', default), depletion ('less'), or
#' 'two.sided'.
#'
#' @return A matrix of pathways and their associated substrates and p-values.
#'
#'
#' @examples
#' \donttest{
#' library(limma)
#' library(org.Rn.eg.db)
#' library(reactome.db)
#' library(annotate)
#' 
#' 
#' data('phospho_L6_ratio_pe')
#' data('SPSs')
#' 
#' ppe <- phospho.L6.ratio.pe
#' sites = paste(sapply(GeneSymbol(ppe), function(x)x),";",
#'     sapply(Residue(ppe), function(x)x),
#'     sapply(Site(ppe), function(x)x),
#'     ";", sep = "")
#' grps = gsub("_.+", "", colnames(ppe))
#' design = model.matrix(~ grps - 1)
#' ctl = which(sites %in% SPSs)
#' ppe = RUVphospho(ppe, M = design, k = 3, ctl = ctl)
#' 
#' phosphoL6 = SummarizedExperiment::assay(ppe, "normalised")
#'                                   
#' # fit linear model for each phosphosite
#' f <- grps
#' X <- model.matrix(~ f - 1)
#' fit <- lmFit(phosphoL6, X)
#'
#' # extract top-ranked phosphosites for each condition compared to basal
#' table.AICAR <- topTable(eBayes(fit), number=Inf, coef = 1)
#' table.Ins <- topTable(eBayes(fit), number=Inf, coef = 3)
#' table.AICARIns <- topTable(eBayes(fit), number=Inf, coef = 2)
#'
#' DE1.RUV <- c(sum(table.AICAR[,'adj.P.Val'] < 0.05),
#'     sum(table.Ins[,'adj.P.Val'] < 0.05),
#'     sum(table.AICARIns[,'adj.P.Val'] < 0.05))
#'
#' # extract top-ranked phosphosites for each group comparison
#' contrast.matrix1 <- makeContrasts(fAICARIns-fIns, levels=X)
#' contrast.matrix2 <- makeContrasts(fAICARIns-fAICAR, levels=X)
#' fit1 <- contrasts.fit(fit, contrast.matrix1)
#' fit2 <- contrasts.fit(fit, contrast.matrix2)
#' table.AICARInsVSIns <- topTable(eBayes(fit1), number=Inf)
#' table.AICARInsVSAICAR <- topTable(eBayes(fit2), number=Inf)
#'
#' DE2.RUV <- c(sum(table.AICARInsVSIns[,'adj.P.Val'] < 0.05),
#'     sum(table.AICARInsVSAICAR[,'adj.P.Val'] < 0.05))
#'
#' o <- rownames(table.AICARInsVSIns)
#' Tc <- cbind(table.Ins[o,'logFC'], table.AICAR[o,'logFC'],
#'             table.AICARIns[o,'logFC'])
#' rownames(Tc) = gsub('(.*)(;[A-Z])([0-9]+)(;)', '\\1;\\3;', o)
#' colnames(Tc) <- c('Ins', 'AICAR', 'AICAR+Ins')
#'
#' # summary phosphosite-level information to proteins for performing downstream
#' # gene-centric analyses.
#' Tc.gene <- phosCollapse(Tc, id=gsub(';.+', '', rownames(Tc)),
#'     stat=apply(abs(Tc), 1, max), by = 'max')
#' geneSet <- names(sort(Tc.gene[,1],
#'                     decreasing = TRUE))[seq(round(nrow(Tc.gene) * 0.1))]
#' #lapply(PhosphoSite.rat, function(x){gsub(';[STY]', ';', x)})
#'
#' 
#' # Preparing Reactome annotation for our pathways analysis
#' pathways = as.list(reactomePATHID2EXTID)
#' 
#' path_names = as.list(reactomePATHID2NAME)
#' name_id = match(names(pathways), names(path_names))
#' names(pathways) = unlist(path_names)[name_id]
#' 
#' pathways = pathways[which(grepl("Rattus norvegicus", names(pathways), 
#'     ignore.case = TRUE))]
#' 
#' pathways = lapply(pathways, function(path) {
#'     gene_name = unname(getSYMBOL(path, data = "org.Rn.eg"))
#'     toupper(unique(gene_name))
#' })
#'
#'
#' # 1D gene-centric pathway analysis
#' path1 <- pathwayOverrepresent(geneSet, annotation=pathways,
#'     universe = rownames(Tc.gene), alter = 'greater')
#' }
#' @export
pathwayOverrepresent <- function(geneSet, annotation, universe,
    alter = "greater") {
    if (missing(geneSet)) {
        stop("Parameter geneSet is missing!")
    }
    if (missing(annotation)) {
        stop("Parameter annotation is missing!")
    }
    if (missing(universe)) {
        stop("Parameter universe is missing!")
    }

    fisherTest.mat <- matrix("NA", ncol = 3, nrow = length(annotation))
    colnames(fisherTest.mat) <- c("pvalue", "# of substrates",
        "substrates")
    for (i in seq_len(length(annotation))) {
        di <- length(intersect(geneSet, annotation[[i]]))
        dn <- length(setdiff(geneSet, annotation[[i]]))
        ndi <- length(setdiff(annotation[[i]], geneSet))
        ndn <- length(setdiff(universe, union(geneSet, annotation[[i]])))
        p <- stats::fisher.test(rbind(c(di, ndi), c(dn, ndn)),
            alternative = alter)$p.value
        substrates <- paste(intersect(geneSet, annotation[[i]]),
            collapse = "|")
        fisherTest.mat[i, ] <- c(p, di, substrates)
    }
    rownames(fisherTest.mat) <- names(annotation)
    fisherTest.mat <- fisherTest.mat[order(as.numeric(fisherTest.mat[,
        1])), ]
    return(fisherTest.mat)
}

#' Phosphosite/Gene set enrichment analysis
#'
#' This function performes phosphosite (or gene) set enrichment analysis using 
#' Wilcoxon Rank Sum test.
#'
#' @param geneStats an array of statistics (e.g. log2 FC) of all quantified
#' genes or phosphosite with names of the array as gene or phosphosite IDs.
#' @param annotation a list of pathways with each element containing an array of
#'  gene IDs.
#' @param alter test for enrichment ('greater', default), depletion ('less'), or
#'  'two.sided'.
#'
#' @return A matrix of pathways and their associated substrates and p-values.
#'
#'
#' @examples
#' \donttest{
#' library(limma)
#'
#' library(org.Rn.eg.db)
#' library(reactome.db)
#' library(annotate)
#' 
#' data('phospho_L6_ratio_pe')
#' data('SPSs')
#' 
#' ppe <- phospho.L6.ratio.pe
#' sites = paste(sapply(GeneSymbol(ppe), function(x)x),";",
#'     sapply(Residue(ppe), function(x)x),
#'     sapply(Site(ppe), function(x)x),
#'     ";", sep = "")
#' grps = gsub("_.+", "", colnames(ppe))
#' design = model.matrix(~ grps - 1)
#' ctl = which(sites %in% SPSs)
#' ppe = RUVphospho(ppe, M = design, k = 3, ctl = ctl)
#' 
#' phosphoL6 = SummarizedExperiment::assay(ppe, "normalised")
#'
#' # fit linear model for each phosphosite
#' f <- grps
#' X <- model.matrix(~ f - 1)
#' fit <- lmFit(phosphoL6, X)
#'
#' # extract top-ranked phosphosites for each condition compared to basal
#' table.AICAR <- topTable(eBayes(fit), number=Inf, coef = 1)
#' table.Ins <- topTable(eBayes(fit), number=Inf, coef = 3)
#' table.AICARIns <- topTable(eBayes(fit), number=Inf, coef = 2)
#'
#' DE1.RUV <- c(sum(table.AICAR[,'adj.P.Val'] < 0.05),
#'     sum(table.Ins[,'adj.P.Val'] < 0.05),
#'     sum(table.AICARIns[,'adj.P.Val'] < 0.05))
#'
#' # extract top-ranked phosphosites for each group comparison
#' contrast.matrix1 <- makeContrasts(fAICARIns-fIns, levels=X)
#' contrast.matrix2 <- makeContrasts(fAICARIns-fAICAR, levels=X)
#' fit1 <- contrasts.fit(fit, contrast.matrix1)
#' fit2 <- contrasts.fit(fit, contrast.matrix2)
#' table.AICARInsVSIns <- topTable(eBayes(fit1), number=Inf)
#' table.AICARInsVSAICAR <- topTable(eBayes(fit2), number=Inf)
#'
#' DE2.RUV <- c(sum(table.AICARInsVSIns[,'adj.P.Val'] < 0.05),
#'     sum(table.AICARInsVSAICAR[,'adj.P.Val'] < 0.05))
#'
#' o <- rownames(table.AICARInsVSIns)
#' Tc <- cbind(table.Ins[o,'logFC'], table.AICAR[o,'logFC'],
#'             table.AICARIns[o,'logFC'])
#' rownames(Tc) = gsub('(.*)(;[A-Z])([0-9]+)(;)', '\\1;\\3;', o)
#' colnames(Tc) <- c('Ins', 'AICAR', 'AICAR+Ins')
#'
#' # summary phosphosite-level information to proteins for performing downstream
#' #  gene-centric analyses.
#' Tc.gene <- phosCollapse(Tc, id=gsub(';.+', '', rownames(Tc)),
#'     stat=apply(abs(Tc), 1, max), by = 'max')
#' 
#' # Preparing Reactome annotation for our pathways analysis
#' pathways = as.list(reactomePATHID2EXTID)
#' 
#' path_names = as.list(reactomePATHID2NAME)
#' name_id = match(names(pathways), names(path_names))
#' names(pathways) = unlist(path_names)[name_id]
#' 
#' pathways = pathways[which(grepl("Rattus norvegicus", names(pathways), 
#'     ignore.case = TRUE))]
#' 
#' pathways = lapply(pathways, function(path) {
#'     gene_name = unname(getSYMBOL(path, data = "org.Rn.eg"))
#'     toupper(unique(gene_name))
#' })
#' 
#' # 1D gene-centric pathway analysis
#' path2 <- pathwayRankBasedEnrichment(Tc.gene[,1],
#'                                     annotation=pathways,
#'                                     alter = 'greater')
#' }
#' @export
pathwayRankBasedEnrichment <- function(geneStats, annotation,
    alter = "greater") {
    if (missing(geneStats))
        stop("Parameter geneStats is missing!")
    if (missing(annotation))
        stop("Parameter annotation is missing!")

    wilcoxTest.mat <- matrix("NA", ncol = 3, nrow = length(annotation))
    colnames(wilcoxTest.mat) <- c("pvalue", "# of substrates",
        "substrates")

    # perform wilcox sum rank test for pathway enrichment
    # analysis
    for (i in seq_len(length(annotation))) {
        p.in <- intersect(names(geneStats), annotation[[i]])
        p.out <- setdiff(names(geneStats), annotation[[i]])

        if (length(geneStats[p.in]) <= 3) {
            wilcoxTest.mat[i, seq_len(3)] <- NA
        } else {
            wilcoxTest.mat[i, 1] <- stats::wilcox.test(geneStats[p.in],
                geneStats[p.out], alternative = alter)$p.value
            wilcoxTest.mat[i, 2] <- length(p.in)
            wilcoxTest.mat[i, 3] <- paste(p.in, collapse = ";")
        }
    }

    rownames(wilcoxTest.mat) <- names(annotation)
    wilcoxTest.mat <- wilcoxTest.mat[order(as.numeric(wilcoxTest.mat[,
        1])), ]
    return(wilcoxTest.mat)
}
