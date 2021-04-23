#' Create frequency matrix
#'
#' @usage createFrequencyMat(substrates.seq)
#'
#' @param substrates.seq A substrate sequence
#'
#' @return A frequency matrix of amino acid from substrates.seq.
#'
#' @examples
#'
#' data("phospho_L6_ratio_pe")
#' 
#' # We will create a frequency matrix of Tfg S198 phosphosite.
#' idx = which(grepl("TFG\\;S198\\;", rownames(phospho.L6.ratio.pe)))
#' substrate.seq = Sequence(phospho.L6.ratio.pe)[idx]
#' freq.mat = createFrequencyMat(substrate.seq)
#'
#' @export
createFrequencyMat <- function(substrates.seq) {
    # substrates.seq.split <-
    # sapply(substrates.seq, strsplit, '')

    substrates.seq.split = Map(function(substrates.seq) {
        unlist(strsplit(substrates.seq, split = ""))
    }, substrates.seq)
    
    frequency.mat <- matrix(data = 0, nrow = 20,
        ncol = length(substrates.seq.split[[1]]))
    rownames(frequency.mat) <- c("A", "R",
        "N", "D", "C", "E", "Q", "G", "H",
        "I", "L", "K", "M", "F", "P", "S",
        "T", "W", "Y", "V")
    colnames(frequency.mat) <- paste("p",
        seq_len(length(substrates.seq.split[[1]])),
        sep = "")

    # calculate frequency
    frequency.list = lapply(seq(ncol(frequency.mat)), function(i) {
        aa = mapply(function(x, i) x[i],
            substrates.seq.split, MoreArgs = list(i = i))
        tmp_col = frequency.mat[,i]
        for (j in seq_len(length(aa))) {
            if (aa[j] == "_") {
                next
            }
            tmp_col[aa[j]] <- tmp_col[aa[j]] + 1
        }
        tmp_col
    })
    frequency.mat = matrix(unlist(frequency.list), ncol = ncol(frequency.mat))
    colnames(frequency.mat) = paste("p",
        seq_len(length(substrates.seq.split[[1]])),
        sep = "")
    rownames(frequency.mat) <- c("A", "R",
        "N", "D", "C", "E", "Q", "G", "H",
        "I", "L", "K", "M", "F", "P", "S",
        "T", "W", "Y", "V")
    
    frequency.mat <- frequency.mat/length(substrates.seq)
    return(frequency.mat)
}

#' Frequency scoring
#'
#' @usage frequencyScoring(sequence.list, frequency.mat)
#'
#' @param sequence.list A vector list of sequences
#' @param frequency.mat A matrix output from `createFrequencyMat`
#'
#' @return A vector of frequency score
#'
#' @examples
#'
#' data('phospho_L6_ratio_pe')
#' data('KinaseMotifs')
#'
#' # Extracting first 10 sequences for demonstration purpose
#' seqs = Sequence(phospho.L6.ratio.pe)
#' seqs = seqs[seq(10)]
#'
#' # extracting flanking sequences
#' seqWin = mapply(function(x) {
#'     mid <- (nchar(x)+1)/2
#'     substr(x, start=(mid-7), stop=(mid+7))
#' }, seqs)
#'
#' # The first 10 for demonstration purpose
#' phospho.L6.ratio = SummarizedExperiment::assay(phospho.L6.ratio.pe, 
#'     "Quantification")[seq(10),]
#'
#' # minimum number of sequences used for compiling motif for each kinase.
#' numMotif=5
#'
#' motif.mouse.list.filtered <-
#'     motif.mouse.list[which(motif.mouse.list$NumInputSeq >= numMotif)]
#'
#' # scoring all phosphosites against all motifs
#' motifScoreMatrix <-
#'     matrix(NA, nrow=nrow(phospho.L6.ratio),
#'         ncol=length(motif.mouse.list.filtered))
#' rownames(motifScoreMatrix) <- rownames(phospho.L6.ratio)
#' colnames(motifScoreMatrix) <- names(motif.mouse.list.filtered)
#'
#' # Scoring phosphosites against kinase motifs
#' for(i in seq_len(length(motif.mouse.list.filtered))) {
#'     motifScoreMatrix[,i] <-
#'         frequencyScoring(seqWin, motif.mouse.list.filtered[[i]])
#'     cat(paste(i, '.', sep=''))
#' }
#'
#' @export
frequencyScoring <- function(sequence.list, frequency.mat) {

    frequency.score <- c()

    frequency.score = lapply(seq(length(sequence.list)), function(idx) {
        if (sequence.list[idx] == "") {
            sequence.list[idx] = "_"
        }
        seqs = unlist(mapply(strsplit, sequence.list[idx],
            MoreArgs = list(split = "")))
        score <- 0
        if (is.na(sequence.list[idx])) {
            frequency.score <- c(frequency.score,
                score)
            return()
        }
        for (i in seq_len(length(seqs))) {
            aa <- c("A", "R", "N", "D", "C",
                "E", "Q", "G", "H", "I",
                "L", "K", "M", "F", "P",
                "S", "T", "W", "Y", "V")
            if (!seqs[i] %in% aa) {
                next
            }
            score <- frequency.mat[seqs[i],
                i] + score
        }
        frequency.score <- c(frequency.score,
            score)
        frequency.score
    })
    frequency.score = unlist(frequency.score)
    
    if (!is.null(names(sequence.list))) {
        names(frequency.score) <- names(sequence.list)
    }

    return(frequency.score)
}
