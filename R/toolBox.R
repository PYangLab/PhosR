#' Median centering and scaling
#'
#' Median centering and scaling of an input numeric matrix
#'
#'
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to samples.
#' @param scale a boolean flag indicating whether to scale the samples.
#' @param grps a string or factor specifying the grouping (replciates).
#' @param reorder To reorder the columns by group (\code{grps}).
#' By default (\code{reorder=FALSE}), original column order is maintained.
#' @param assay an assay to be selected if \code{mat} is a PhosphoExperiment 
#' object.
#'
#' @return A median scaled matrix
#'
#' @importFrom limma normalizeMedianAbsValues
#'
#' @examples
#'
#' data('phospho.cells.Ins.sample')
#' grps = gsub('_[0-9]{1}', '', colnames(phospho.cells.Ins))
#' phospho.cells.Ins.filtered <- selectGrps(phospho.cells.Ins, grps, 0.5, n=1)
#'
#' set.seed(123)
#' phospho.cells.Ins.impute <-
#'     scImpute(phospho.cells.Ins.filtered,
#'             0.5,
#'             grps)[,colnames(phospho.cells.Ins.filtered)]
#'
#' set.seed(123)
#' phospho.cells.Ins.impute[,seq(5)] <- ptImpute(
#'     phospho.cells.Ins.impute[,seq(6,10)],
#'     phospho.cells.Ins.impute[,seq(5)], percent1 = 0.6,
#'     percent2 = 0, paired = FALSE)
#'
#' phospho.cells.Ins.ms <-
#'     medianScaling(phospho.cells.Ins.impute, scale = FALSE)
#'
#' @export
#'
medianScaling <- function(mat, scale = FALSE, grps = NULL, reorder = FALSE, 
                          assay = NULL) {
    pe = FALSE
    if (methods::is(mat, "PhosphoExperiment")) {
        pe = TRUE
        mat.orig = mat
        if (is.null(assay)) {
            mat = SummarizedExperiment::assay(mat)
        } else {
            mat = SummarizedExperiment::assay(mat, assay)
        }
    }
    mat.medianScaled <- NULL
    if (!is.null(grps)) {
        tmp <- lapply(split(seq_len(ncol(mat)), grps), function(i) mat[,
            i])
        mat.medianScaled <- do.call(cbind, lapply(tmp, medianScale,
            scale = scale))

        if (!reorder) {
            mat.medianScaled <- mat.medianScaled[, colnames(mat)]
        }
    } else {
        mat.medianScaled <- medianScale(mat, scale = scale)
    }
    
    if (pe) {
        SummarizedExperiment::assay(mat.orig, "scaled") = mat.medianScaled
        mat.medianScaled = mat.orig
    }

    return(mat.medianScaled)
}

medianScale <- function(mat, scale) {
    if (scale) {
        normcont <- stats::median(apply(mat, 2, stats::median,
            na.rm = TRUE))
        adjval <- apply(mat, 2, stats::median, na.rm = TRUE) -
            normcont
        mat.scaled <- normalizeMedianAbsValues(sweep(mat, 2,
            adjval, "-"))
        return(mat.scaled)
    } else {
        normcont <- stats::median(apply(mat, 2, stats::median,
            na.rm = TRUE))
        adjval <- apply(mat, 2, stats::median, na.rm = TRUE) -
            normcont
        mat.unscaled <- sweep(mat, 2, adjval, "-")
        return(mat.unscaled)
    }
}


#' Standardisation
#'
#' Standardisation by z-score transformation.
#'
#' @usage standardise(mat)
#'
#' @param mat a matrix (or a PhosphoExperiment object) with rows correspond 
#' to phosphosites and columns correspond to samples.
#'
#' @return A standardised matrix
#'
#' @examples
#' data('phospho_L6_ratio_pe')
#' data('SPSs')
#'
#' grps = gsub('_.+', '', colnames(phospho.L6.ratio.pe))
#'
#' # Construct a design matrix by condition
#' design = model.matrix(~ grps - 1)
#'
#' # phosphoproteomics data normalisation using RUV
#' L6.sites = paste(sapply(GeneSymbol(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  ";",
#'                  sapply(Residue(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  sapply(Site(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  ";", sep = "")
#' ctl = which(L6.sites %in% SPSs)
#' phospho.L6.ratio.pe = RUVphospho(phospho.L6.ratio.pe,
#'                                  M = design, k = 3,ctl = ctl)
#' 
#' phosphoL6 = SummarizedExperiment::assay(phospho.L6.ratio.pe, "normalised")
#' 
#' # filter for up-regulated phosphosites
#' phosphoL6.mean <- meanAbundance(phosphoL6, grps = grps)
#' aov <- matANOVA(mat=phosphoL6, grps = grps)
#' phosphoL6.reg <- phosphoL6[(aov < 0.05) &
#'                         (rowSums(phosphoL6.mean > 0.5) > 0),,drop = FALSE]
#' L6.phos.std <- standardise(phosphoL6.reg)
#'
#' @export
#'
standardise <- function(mat) {
    
    pe = FALSE
    if (methods::is(mat, "PhosphoExperiment")) {
        pe = TRUE
        mat.orig = mat
        if (is.null(assay)) {
            mat = SummarizedExperiment::assay(mat)
        } else {
            mat = SummarizedExperiment::assay(mat, assay)
        }
    }
    
    means <- apply(mat, 1, mean)
    sds <- apply(mat, 1, stats::sd)

    X2 <- sweep(mat, MARGIN = 1, STATS = means, FUN = "-")
    zX <- sweep(X2, MARGIN = 1, STATS = sds, FUN = "/")
    
    if (pe) {
        SummarizedExperiment::assay(mat.orig, "standardise") = zX
        zX = mat.orig
    }

    return(zX)
}

#' @title Multi-intersection, union
#'
#' @description A recusive loop for intersecting multiple sets.
#'
#' @aliases mUnion
#'
#' @usage
#' mIntersect(x, y, ...)
#' mUnion(x, y, ...)
#'
#' @param x,y,... objects to find intersection/union.
#'
#' @return An intersection/union of input parameters
#'
#' @examples
#'
#' data('phospho_liverInsTC_RUV_sample')
#' data('phospho_L6_ratio')
#'
#' site1 <- gsub('~[STY]', ';',
#'             sapply(strsplit(rownames(phospho.L6.ratio), ';'),
#'                     function(x){paste(toupper(x[2]), x[3], sep=';')}))
#' site2 <- rownames(phospho.liver.Ins.TC.ratio.RUV)
#'
#' # step 2: rank by fold changes
#' treatment.grps = split(seq(ncol(phospho.L6.ratio)), 
#'     gsub('_exp\\d+', '', colnames(phospho.L6.ratio)))
#' tmp <- do.call(
#'     cbind, 
#'     lapply(treatment.grps, function(i){
#'         rowMeans(phospho.L6.ratio[,i])
#'     })
#' )
#' site1 <- t(sapply(split(data.frame(tmp), site1), colMeans))[,-1]
#'
#' treatment.grps = split(
#'     seq(ncol(phospho.liver.Ins.TC.ratio.RUV)),
#'     gsub('(Intensity\\.)(.*)(\\_Bio\\d+)', '\\2', 
#'         colnames(phospho.liver.Ins.TC.ratio.RUV)
#'     )
#' )
#' tmp <- do.call(
#'     cbind, 
#'     lapply(
#'         treatment.grps,
#'         function(i){
#'             rowMeans(phospho.liver.Ins.TC.ratio.RUV[,i])
#'         }
#'     )
#' )
#' site2 <- t(sapply(split(data.frame(tmp), site2), colMeans))
#'
#' o <- mIntersect(site1, site2)
#'
#' @export mIntersect
#'
mIntersect <- function(x, y, ...) {
    if (missing(...))
        intersect(x, y) else intersect(x, mIntersect(y, ...))
}


#' @export mUnion
#'
mUnion <- function(x, y, ...) {
    if (missing(...))
        union(x, y) else union(x, mUnion(y, ...))
}



#' Minmax scaling
#'
#' Perform a minmax standardisation to scale data into 0 to 1 range
#'
#' @usage minmax(mat)
#'
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to condition
#'
#' @return Minmax standardised matrix
#'
#' @examples
#' data('phospho_L6_ratio_pe')
#' data('SPSs')
#' data('PhosphoSitePlus')
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
#' # filter for up-regulated phosphosites
#' phosphoL6.mean <- meanAbundance(phosphoL6, grps = grps)
#' aov <- matANOVA(mat=phosphoL6, grps = grps)
#' idx <- (aov < 0.05) & (rowSums(phosphoL6.mean > 0.5) > 0)
#' phosphoL6.reg <- phosphoL6[idx, ,drop = FALSE]
#' 
#' L6.phos.std <- standardise(phosphoL6.reg)
#'
#' ks.profile.list <- kinaseSubstrateProfile(PhosphoSite.mouse, L6.phos.std)
#'
#' data(KinaseMotifs)
#' 
#' numMotif = 5
#' numSub = 1
#'
#' motif.mouse.list.filtered <-
#'     motif.mouse.list[which(motif.mouse.list$NumInputSeq >= numMotif)]
#' ks.profile.list.filtered <-
#'     ks.profile.list[which(ks.profile.list$NumSub >= numSub)]
#'
#' # scoring all phosphosites against all motifs
#' motifScoreMatrix <-
#'     matrix(NA, nrow=nrow(L6.phos.std),
#'     ncol=length(motif.mouse.list.filtered))
#' rownames(motifScoreMatrix) <- rownames(L6.phos.std)
#' colnames(motifScoreMatrix) <- names(motif.mouse.list.filtered)
#' 
#' L6.phos.seq <- Sequence(ppe)[idx]
#' 
#' # extracting flanking sequences
#' seqWin = mapply(function(x) {
#'     mid <- (nchar(x)+1)/2
#'     substr(x, start=(mid-7), stop=(mid+7))
#' }, L6.phos.seq)
#'
#'
#' print('Scoring phosphosites against kinase motifs:')
#' for(i in seq_len(length(motif.mouse.list.filtered))) {
#'     motifScoreMatrix[,i] <-
#'         frequencyScoring(seqWin, motif.mouse.list.filtered[[i]])
#'         cat(paste(i, '.', sep=''))
#' }
#' motifScoreMatrix <- minmax(motifScoreMatrix)
#'
#' @export
#'
minmax <- function(mat) {
    # minmax normalise
    apply(mat, 2, function(x) {
        (x - min(x))/(max(x) - min(x))
    })
}

#' Summarising phosphosites to proteins
#'
#' Summarising phosphosite-level information to proteins for performing
#' downstream gene-centric analyses.
#'
#' @usage phosCollapse(mat, id, stat, by='min')
#'
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to samples.
#' @param id an array indicating the groupping of phosphosites etc.
#' @param stat an array containing statistics of phosphosite such as
#' phosphorylation levels.
#' @param by how to summarise phosphosites using their statistics. Either by
#' 'min' (default), 'max', or 'mid'.
#'
#' @return A matrix summarised to protein level
#'
#' @examples
#' library(limma)
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
#' ctl = which(L6.sites %in% SPSs)
#' phospho.L6.ratio.pe = RUVphospho(phospho.L6.ratio.pe, 
#'                                   M = design, k = 3, ctl = ctl)
#'
#' # fit linear model for each phosphosite
#' f <- grps
#' X <- model.matrix(~ f - 1)
#' fit <- lmFit(SummarizedExperiment::assay(phospho.L6.ratio.pe, "normalised"), X)
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
#'
#' @export
#'
phosCollapse <- function(mat, id, stat, by = "min") {
    listMat <- split(data.frame(mat, stat), id)

    matNew <- c()
    if (by == "min") {
        matNew <- as.matrix(do.call(rbind, lapply(listMat, function(x) {
            if (nrow(x) == 1) {
                as.numeric(x[1, seq_len(ncol(x) - 1)])
            } else {
                x[order(as.numeric(x[, ncol(x)]))[1], seq_len(ncol(x) - 1)]
            }
        })))
    } else if (by == "max") {
        matNew <- as.matrix(do.call(rbind, lapply(listMat, function(x) {
            if (nrow(x) == 1) {
                as.numeric(x[1, seq_len(ncol(x) - 1)])
            } else {
                x[order(as.numeric(x[, ncol(x)]), decreasing = TRUE)[1],
                    seq_len(ncol(x) - 1)]
            }
        })))
    } else if (by == "mid") {
        matNew <- as.matrix(do.call(rbind, lapply(listMat, function(x) {
            if (nrow(x) == 1) {
                as.numeric(x[1, seq_len(ncol(x) - 1)])
            } else {
                mid <- round(nrow(x)/2)
                x[order(as.numeric(x[, ncol(x)]))[mid], seq_len(ncol(x) - 1)]
            }
        })))
    } else {
        stop("Unrecognised way of collapsing the data!")
    }

    return(matNew)
}

#' ANOVA test
#'
#' @description Performs an ANOVA test and returns its adjusted p-value
#'
#'
#' @usage matANOVA(mat, grps)
#'
#' @param mat An p by n matrix where p is the number of phosphosites and n is
#' the number of samples
#' @param grps A vector of length n, with group or time point information of the
#'  samples
#'
#' @return A vector of multiple testing adjusted p-values
#'
#' @examples
#' data('phospho_L6_ratio_pe')
#' data('SPSs')
#' data('PhosphoSitePlus')
#'
#' grps = gsub('_.+', '', colnames(phospho.L6.ratio.pe))
#' 
#' # Construct a design matrix by condition
#' design = model.matrix(~ grps - 1)
#' 
#' # phosphoproteomics data normalisation using RUV
#' L6.sites = paste(sapply(GeneSymbol(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  ";",
#'                  sapply(Residue(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  sapply(Site(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  ";", sep = "")
#' ctl = which(L6.sites %in% SPSs)
#' phospho.L6.ratio.pe = RUVphospho(phospho.L6.ratio.pe,
#'                                  M = design, k = 3,ctl = ctl)
#' phosphoL6 = SummarizedExperiment::assay(phospho.L6.ratio.pe, "normalised")
#'
#' # filter for up-regulated phosphosites
#' phosphoL6.mean <- meanAbundance(phosphoL6, grps = grps)
#' aov <- matANOVA(mat=phosphoL6, grps = grps)
#' @export
matANOVA <- function(mat, grps) {
    ps <- apply(mat, 1, function(x) {
        summary(stats::aov(as.numeric(x) ~ grps))[[1]][["Pr(>F)"]][1]
    })

    # adjust for multiple testing
    ps.adj <- stats::p.adjust(ps, method = "fdr")
    return(ps.adj)
}

#' @title Obtain average expression from replicates
#'
#' @usage meanAbundance(mat, grps)
#'
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to samples.
#' @param grps a string specifying the grouping (replciates).
#'
#' @return a matrix with mean expression from replicates
#'
#' @examples
#' data('phospho_L6_ratio_pe')
#' data('SPSs')
#' data('PhosphoSitePlus')
#'
#' grps = gsub('_.+', '', colnames(phospho.L6.ratio.pe))
#' 
#' # Construct a design matrix by condition
#' design = model.matrix(~ grps - 1)
#' 
#' # phosphoproteomics data normalisation using RUV
#' L6.sites = paste(sapply(GeneSymbol(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  ";",
#'                  sapply(Residue(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  sapply(Site(phospho.L6.ratio.pe), function(x)paste(x)),
#'                  ";", sep = "")
#' ctl = which(L6.sites %in% SPSs)
#' phospho.L6.ratio.pe = RUVphospho(phospho.L6.ratio.pe,
#'                                  M = design, k = 3,ctl = ctl)
#' 
#' phosphoL6 = SummarizedExperiment::assay(phospho.L6.ratio.pe, "normalised")
#' phosphoL6.mean <- meanAbundance(phosphoL6, grps = grps)
#'
#' @export
meanAbundance <- function(mat, grps) {
    # meanMat <- sapply(split(seq_len(ncol(mat)), grps),
    # function(i) rowMeans(mat[,i]))[,unique(grps)]
    meanMat = mapply(function(i, mat) {
        rowMeans(mat[, i], na.rm = TRUE)
    }, split(seq_len(ncol(mat)), grps), MoreArgs = list(mat = mat))[,
        unique(grps)]
    return(meanMat)
}
