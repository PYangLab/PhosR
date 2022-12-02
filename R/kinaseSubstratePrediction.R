############################### Global kinase annotatiion ##

#' @title Kinase substrate scoring
#'
#' @description This function generates substrate scores for kinases that pass
#' filtering based on both motifs and dynamic profiles
#'
#' @param substrate.list A list of kinases with each element containing an array
#'  of substrates.
#' @param mat A matrix with rows correspond to phosphosites and columns
#' correspond to samples.
#' @param seqs An array containing aa sequences surrounding each of all
#' phosphosites.
#' Each sequence has length of 15 (-7, p, +7).
#' @param numMotif Minimum number of sequences used for compiling motif for
#' each kinase.
#' Default is 5.
#' @param numSub Minimum number of phosphosites used for compiling
#' phosphorylation
#' profile for each kinase. Default is 1.
#' @param species Motif list species to be used. Currently there are 
#' \code{mouse} (default), \code{human} and \code{rat}.
#' @param verbose Default to \code{TRUE} to show messages during the progress.
#' All messages will be suppressed if set to \code{FALSE}
#'
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom utils data
#'
#' @return A list of 4 elements.
#' \code{motifScoreMatrix}, \code{profileScoreMatrix},
#' \code{combinedScoreMatrix}, \code{ksActivityMatrix} (kinase activity matrix)
#' and their \code{weights}.
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
#' rownames(L6.phos.std) <- paste0(GeneSymbol(ppe), ";", Residue(ppe), 
#'     Site(ppe), ";")[idx]
#' 
#' L6.phos.seq <- Sequence(ppe)[idx]
#' 
#' L6.matrices <- kinaseSubstrateScore(PhosphoSite.mouse, L6.phos.std,
#'     L6.phos.seq, numMotif = 5, numSub = 1)
#' @export
kinaseSubstrateScore <- function(substrate.list, mat, seqs, numMotif = 5, 
                                numSub = 1, species = "mouse", verbose = TRUE) {
    ks.profile.list <- kinaseSubstrateProfile(substrate.list, mat)
    # motif.mouse.list = PhosR::motif.mouse.list
    utils::data("KinaseMotifs", envir = environment())
    if (!(species %in% c("mouse", "human", "rat"))) {
        stop("Parameter 'species' must be one of 'mouse', 'human' or 'rat'")
    }
    
    if (any(is.na(mat))) {
        stop("Phosphosite quantification matrix contains NAs. Please remove NAs in mat.")
    }
    
    
    if (species == "mouse") {
        motif.list <- motif.mouse.list
    } else if (species == "human") {
        motif.list <- motif.human.list
    } else if (species == "rat") {
        motif.list <- motif.rat.list
    }
    
    if (verbose) {
        message(paste("Number of kinases passed motif size filtering:", 
                      sum(motif.list$NumInputSeq >= numMotif)))
        message(paste("Number of kinases passed profile size filtering:",
            sum(ks.profile.list$NumSub >= numSub)))
    }
    
    motif.list.filtered <- motif.list[
        which(motif.list$NumInputSeq >= numMotif) ]
    ks.profile.list.filtered <- ks.profile.list[
        which(ks.profile.list$NumSub >= numSub) ]
    # scoring all phosphosites against all motifs
    motifScoreMatrix =
        scorePhosphositesMotifs(mat, motif.list.filtered, seqs, verbose)
    if (verbose) {
        message("done.")
        # scoring all phosphosites against all profiles
        message("Scoring phosphosites against kinase-substrate profiles:")
    }
    profileScoreMatrix = scorePhosphositeProfile(mat, ks.profile.list.filtered)
    if (verbose) {
        message("done.")
        ### prioritisation by integrating the two parts
        message("Generating combined scores for phosphosites
by motifs and phospho profiles:")
    }
    o <- intersect(colnames(motifScoreMatrix), colnames(profileScoreMatrix))
    combinedScoreMatrix <- matrix(NA, nrow = nrow(motifScoreMatrix),
        ncol = length(o))
    colnames(combinedScoreMatrix) <- o
    rownames(combinedScoreMatrix) <- rownames(motifScoreMatrix)
    # normalising weights for the two parts
    w1 <- log(rank(motif.list$NumInputSeq[o]) + 1)
    w2 <- log(rank(ks.profile.list$NumSub[o]) + 1)
    w3 <- w1 + w2
    for (i in seq_len(length(o))) {
        # weight the two parts by the number of
        # motifs and quantified known substrates
        combinedScoreMatrix[, i] <- (w1[i]/(w1[i] +
            w2[i]) * motifScoreMatrix[, o[i]]) +
            (w2[i]/(w1[i] + w2[i]) * profileScoreMatrix[, o[i]])
    }
    if (verbose)
        message("done.")
    # visualise
    ksActivityMatrix <- do.call(rbind, ks.profile.list.filtered)[o, ]
    phosScoringMatrices <- list(motifScoreMatrix = motifScoreMatrix,
        profileScoreMatrix = profileScoreMatrix,
        combinedScoreMatrix = combinedScoreMatrix,
        ksActivityMatrix = ksActivityMatrix, weights = w3)
    kinaseSubstrateHeatmap(phosScoringMatrices)
    return(phosScoringMatrices)
}

scorePhosphositesMotifs = function(mat, motif.mouse.list.filtered, seqs, 
                                verbose = TRUE) {
    motifScoreMatrix <- matrix(NA, nrow = nrow(mat),
                            ncol = length(motif.mouse.list.filtered))
    rownames(motifScoreMatrix) <- rownames(mat)
    colnames(motifScoreMatrix) <- names(motif.mouse.list.filtered)
    # extracting flanking sequences
    seqWin = mapply(function(x) {
        mid <- (nchar(x) + 1)/2
        substr(x, start = (mid - 7), stop = (mid + 7))
    }, seqs)
    
    if (verbose)
        message("Scoring phosphosites against kinase motifs:")
    for (i in seq_len(length(motif.mouse.list.filtered))) {
        motifScoreMatrix[, i] <- frequencyScoring(seqWin,
                                                motif.mouse.list.filtered[[i]])
        if (verbose)
            message(paste(i, ".", sep = ""))
    }
    motifScoreMatrix <- minmax(motifScoreMatrix)
    motifScoreMatrix
}

scorePhosphositeProfile = function(mat, ks.profile.list.filtered) {
    profileScoreMatrix <- (t(apply(mat, 1,
        cor, t(do.call(rbind, ks.profile.list.filtered)))) + 1)/2
    rownames(profileScoreMatrix) <- rownames(mat)
    colnames(profileScoreMatrix) <- names(ks.profile.list.filtered)
    profileScoreMatrix
}

#' @title Kinase substrate profiling
#'
#' @description This function generates substrate profiles for kinases that have
#' one or more substrates quantified in the phosphoproteome data.
#'
#'
#' @param substrate.list a list of kinases with each element containing an array
#' of substrates.
#' @param mat a matrix with rows correspond to phosphosites and columns
#' correspond to samples.
#' @return Kinase profile list.
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
#' @export
kinaseSubstrateProfile <- function(substrate.list, mat) {
    
    if (any(is.na(mat))) {
        stop("Phosphosite quantification matrix contains NAs. Please remove NAs in mat.")
    }
    
    # generate kinase substrate profile list
    ks.profile.list <- lapply(substrate.list,
        function(x) {
            ns <- intersect(x, rownames(mat))
            m <- c()
            if (length(ns) == 1) {
                m <- mat[ns, ]
            } else if (length(ns) > 1) {
                m <- apply(mat[ns, ], 2,
                    median)
            } else {
                m <- NA
            }
            return(m)
        })

    # ks.profile.list$NumSub <-
    # sapply(substrate.list, function(x){
    # sum(x %in% rownames(mat)) })
    ks.profile.list$NumSub <- mapply(function(x,
        mat) {
        sum(x %in% rownames(mat))
    }, substrate.list, MoreArgs = list(mat = mat))

    return(ks.profile.list)
}


kinaseActivityHeatmap <- function(ksProfileMatrix) {
    # KinaseFamily = PhosR::KinaseFamily
    utils::data("KinaseFamily", envir = environment())
    o <- intersect(rownames(ksProfileMatrix),
        rownames(KinaseFamily))
    annotation_row = data.frame(group = KinaseFamily[o,
        "kinase_group"], family = KinaseFamily[o,
        "kinase_family"])
    rownames(annotation_row) <- o

    pheatmap(ksProfileMatrix, annotation_row = annotation_row,
        fontsize = 7, main = paste("Kinase substrate profiles"))

}

#' @title Kinase-substrate annotation prioritisation heatmap
#'
#' @param phosScoringMatrices a matrix returned from kinaseSubstrateScore.
#' @param top the number of top ranked phosphosites for each kinase to be
#' included in the heatmap. Default is 1.
#' @param printPlot indicate whether the plot should be saved as a PDF
#' in the specified directory. Default is NULL, otherwise specify TRUE.
#' @param filePath path name to save the plot as a PDF file. 
#' Default saves in the working directory.
#' @param width width of PDF.
#' @param height height of PDF.
#'
#' @return a pheatmap object.
#'
#' @import pheatmap
#' @importFrom utils data
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#'
#' @examples
#' \donttest{
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
#' rownames(L6.phos.std) <- paste0(GeneSymbol(ppe), ";", Residue(ppe), 
#'     Site(ppe), ";")[idx]
#' 
#' L6.phos.seq <- Sequence(ppe)[idx]
#' 
#' L6.matrices <- kinaseSubstrateScore(PhosphoSite.mouse, L6.phos.std,
#'     L6.phos.seq, numMotif = 5, numSub = 1)
#'     
#' kinaseSubstrateHeatmap(L6.matrices)
#' kinaseSubstrateHeatmap(L6.matrices, printPlot=TRUE)
#' }
#' @export
kinaseSubstrateHeatmap <- function(phosScoringMatrices, top = 3, printPlot=NULL, 
    filePath="./kinaseSubstrateHeatmap.pdf", width=10, height=10) {
    # KinaseFamily = PhosR::KinaseFamily
    utils::data("KinaseFamily", envir = environment())
    ####### heatmap 1
    sites <- c()
    for (i in seq_len(ncol(phosScoringMatrices$combinedScoreMatrix))) {
        sites <- union(sites, names(
            sort(phosScoringMatrices$combinedScoreMatrix[,i],
                decreasing = TRUE)[seq_len(top)]))
    }

    o <- intersect(colnames(phosScoringMatrices$combinedScoreMatrix),
        rownames(KinaseFamily))
    annotation_col = data.frame(group = KinaseFamily[o,
        "kinase_group"], family = KinaseFamily[o,
        "kinase_family"])
    rownames(annotation_col) <- o

    if (is.null(printPlot)==TRUE) {
        
        pheatmap(phosScoringMatrices$combinedScoreMatrix[sites,
        ], annotation_col = annotation_col,
        cluster_rows = TRUE, cluster_cols = TRUE,
        fontsize = 7, main = paste("Top",
                                   top, "phosphosite(s) for each kinase"))
        
    } else {
        
        pdf(file=filePath, width=width, height=height)
        
        pheatmap(phosScoringMatrices$combinedScoreMatrix[sites,
        ], annotation_col = annotation_col,
        cluster_rows = TRUE, cluster_cols = TRUE,
        fontsize = 7, main = paste("Top",
                                   top, "phosphosite(s) for each kinase"))
        
        dev.off()
    }
}

#' @title Phosphosite annotation
#'
#' @description This function plots the combined scores of each of all kinases
#' for a given phosphosites
#'
#'
#' @param site site the ID of a phosphosite
#' @param phosScoringMatrices output from function kinaseSubstrateScore()
#' @param predMatrix a prediction matrix from kinaseSubstratePred()
#'
#' @return A graphical plot
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
#' rownames(L6.phos.std) <- paste0(GeneSymbol(ppe), ";", Residue(ppe), 
#'     Site(ppe), ";")[idx]
#' 
#' L6.phos.seq <- Sequence(ppe)[idx]
#' 
#' L6.matrices <- kinaseSubstrateScore(PhosphoSite.mouse, L6.phos.std,
#'     L6.phos.seq, numMotif = 5, numSub = 1)
#'     
#' set.seed(1)
#' L6.predMat <- kinaseSubstratePred(L6.matrices, top=30)
#' dev.off()
#' 
#' # We will look at the phosphosite AAK1;S677 for demonstration purpose.
#' site = "AAK1;S677;"
#' siteAnnotate(site, L6.matrices, L6.predMat)
#' @export
siteAnnotate <- function(site, phosScoringMatrices,
    predMatrix) {
    od <- order(predMatrix[site, ], decreasing = TRUE)
    kinases <- colnames(predMatrix)[od]
    # cols <- rainbow(length(kinases))

    par(mfrow = c(4, 1))
    barplot(predMatrix[site, kinases], las = 2,
        ylab = "Prediction score", col = "red3",
        main = paste("Site =", site), ylim = c(0,
            1))
    barplot(phosScoringMatrices$combinedScoreMatrix[site,
        kinases], las = 2, ylab = "Combined score",
        col = "orange2", ylim = c(0, 1))
    barplot(phosScoringMatrices$motifScoreMatrix[site,
        kinases], las = 2, ylab = "Motif score",
        col = "green4", ylim = c(0, 1))
    barplot(phosScoringMatrices$profileScoreMatrix[site,
        kinases], las = 2, ylab = "Profile score",
        col = "lightblue3", ylim = c(0, 1))
}


#' kinaseSubstratePred
#'
#' A machine learning approach for predicting specific kinase for a given
#' substrate. This prediction framework utilise adaptive sampling.
#'
#' @usage
#' kinaseSubstratePred(
#'     phosScoringMatrices,
#'     ensembleSize = 10,
#'     top = 50,
#'     cs = 0.8,
#'     inclusion = 20,
#'     iter = 5,
#'     verbose = TRUE
#' )
#'
#' @param phosScoringMatrices An output of kinaseSubstrateScore.
#' @param ensembleSize An ensemble size.
#' @param top a number to select top kinase substrates.
#' @param cs Score threshold.
#' @param inclusion A minimal number of substrates required for a kinase to be
#' selected.
#' @param iter A number of iterations for adaSampling.
#' @param verbose Default to \code{TRUE} to show messages during the progress.
#' All messages will be suppressed if set to \code{FALSE}
#'
#' @return Kinase prediction matrix
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
#' rownames(L6.phos.std) <- paste0(GeneSymbol(ppe), ";", Residue(ppe), 
#'     Site(ppe), ";")[idx]
#' 
#' L6.phos.seq <- Sequence(ppe)[idx]
#' 
#' L6.matrices <- kinaseSubstrateScore(PhosphoSite.mouse, L6.phos.std,
#'     L6.phos.seq, numMotif = 5, numSub = 1)
#' set.seed(1)
#' L6.predMat <- kinaseSubstratePred(L6.matrices, top=30)
#' @export
#'
kinaseSubstratePred <- function(phosScoringMatrices,
    ensembleSize = 10, top = 50, cs = 0.8,
    inclusion = 20, iter = 5, verbose = TRUE) {
    # create the list of kinase-substrates for prediction
    substrate.list = substrateList(phosScoringMatrices, top, cs, inclusion)

    # building the positive training set
    if (verbose)
        message("Predicting kinases for phosphosites:")
    featureMat <- phosScoringMatrices$combinedScoreMatrix
    predMatrix <- matrix(0, nrow = nrow(featureMat),
        ncol = length(substrate.list))
    colnames(predMatrix) <- names(substrate.list)
    rownames(predMatrix) <- rownames(featureMat)
    
    tmp.list = lapply(seq(length(substrate.list)), function(i) {
        positive.train <- featureMat[substrate.list[[i]],]
        positive.cls <- rep(1, length(substrate.list[[i]]))
        negative.pool <- featureMat[!(rownames(featureMat) %in%
                substrate.list[[i]]), ]
        if (verbose)
            message(paste(i, ".", sep = ""))
        tmp_col = predMatrix[,i]
        for (e in seq_len(ensembleSize)) {
            negativeSize <- length(substrate.list[[i]])
            idx <- sample(seq_len(nrow(negative.pool)),
                size = negativeSize, replace = FALSE)
            negative.samples <- rownames(negative.pool)[idx]
            negative.train <- featureMat[negative.samples,]
            negative.cls <- rep(2, length(negative.samples))
            train.mat <- rbind(positive.train, negative.train)
            cls <- as.factor(c(positive.cls, negative.cls))
            names(cls) <- rownames(train.mat)
            pred <- multiAdaSampling(train.mat,
                test.mat = featureMat, label = cls,
                kernelType = "radial", iter = iter)
            tmp_col <- tmp_col[names(pred[, 1])] + pred[, 1]
        }
        tmp_col
    })
    predMatrix = matrix(unlist(tmp.list), ncol = ncol(predMatrix))
    colnames(predMatrix) = names(substrate.list)
    rownames(predMatrix) = rownames(featureMat)
    
    predMatrix <- predMatrix/ensembleSize
    if (verbose)
        message("done")
    return(predMatrix)
}

substrateList = function(phosScoringMatrices, top, cs, inclusion) {
    kinaseSel <- c()
    substrate.list <- list()
    count <- 0
    for (i in seq_len(ncol(phosScoringMatrices$combinedScoreMatrix))) {
        sel <- names(which(sort(phosScoringMatrices$combinedScoreMatrix[, i],
            decreasing = TRUE)[seq_len(top)] > cs))
        if (length(sel) >= inclusion) {
            count <- count + 1
            substrate.list[[count]] <- sel
            kinaseSel <- c(kinaseSel,
                colnames(phosScoringMatrices$combinedScoreMatrix)[i])
        }
    }
    if (length(substrate.list)) { # If substrate.list is not empty
        names(substrate.list) <- kinaseSel
    }
    substrate.list
}


#' @import stats
#' @import e1071
#' @importFrom utils tail
multiAdaSampling <- function(train.mat, test.mat,
    label, kernelType, iter = 5) {

    X <- train.mat
    Y <- label

    model <- c()
    prob.mat <- c()

    for (i in seq_len(iter)) {
        tmp <- X
        rownames(tmp) <- NULL
        model <- e1071::svm(tmp, factor(Y),
            kernel = kernelType, probability = TRUE)
        prob.mat <- attr(predict(model, train.mat,
            decision.values = FALSE, probability = TRUE),
            "probabilities")

        X <- c()
        Y <- c()
        tmp = lapply(seq(ncol(prob.mat)), function(j) {
            voteClass <- prob.mat[label == colnames(prob.mat)[j], ]
            idx <- c()
            idx <- sample(seq_len(nrow(voteClass)),
                size = nrow(voteClass), replace = TRUE,
                prob = voteClass[, j])
            X <- rbind(X, train.mat[rownames(voteClass)[idx],])
            Y <- c(Y, label[rownames(voteClass)[idx]])
            cbind(X,Y)
        })
        tmp = do.call(rbind, tmp)
        Y = tmp[,tail(seq(ncol(tmp)), 1)]
        X = tmp[,-tail(seq(ncol(tmp)),1)]
        
    }

    pred <- attr(predict(model, newdata = test.mat,
        probability = TRUE), "prob")
    return(pred)
}
