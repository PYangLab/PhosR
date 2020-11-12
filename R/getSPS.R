
## format of input dataset
# 1. PhosphoExperiment objects
# 2. Quantification matrix in PhosphoExperiment object is log2 fold change matrix without control

## parameter
# phosData: a list of users' PhosphoExperiment objects from which generate SPSs
# conds: a numeric vector contains the numbers of conditions in each of the datasets
# num: the number of identified SPSs, by default is 100
# residuleInfo: whether the phosphosite contains amino acid information or not, by default is FALSE

getSPS <- function(phosData = ..., conds = ..., num = 100, residueInfo = FALSE) {
    sites <- list()
    sites.unique <- list()
    mat.max <- list()
    n <- length(phosData)
    m <- length(conds)
    if (n < 2) {
        stop("Please use more than one dataset")
    }
    if (n != m) {
        stop("Please use same number of datasets and conditions")
    }
    for (i in 1:n) {
        if (!"PhosphoExperiment" %in% is(phosData[[i]])) {
            stop("Wrong phosData, need to be PhosphoExperiment object") }
        }
    for (i in 1:n) {
        if (residueInfo) {
            sites[[i]] <- paste(toupper(phosData[[i]]@GeneSymbol), paste(phosData[[i]]@Residue, phosData[[i]]@Site, sep = ""), sep = ";")
            
        } else {
            sites[[i]] <- paste(toupper(phosData[[i]]@GeneSymbol), phosData[[i]]@Site, sep = ";")
        }
        sites.unique[[i]] <- unique(sites[[i]])
        nrep <- ncol(phosData[[i]]@Quantification)/conds[i]
        if (nrep == 1) {
            mat.mean <- phosData[[i]]@Quantification # if no replicates
        } else {
            grps <- rep(seq_len(conds[i]), each = nrep)
            mat.mean <- PhosR::meanAbundance(phosData[[i]]@Quantification, grps)
        }
        
        sites.mean <- t(sapply(split(data.frame(mat.mean), sites[[i]]), colMeans))
        # maximum fold change for each phosphosite
        sites.max <- apply(sites.mean, 1, function(x){x[which.max(abs(x))]})
        mat.max[[i]] <- sites.max
    }
    
    # identify the overlapped sites in all datasets, if it is less than 1000, find top 1000 overlapped phosphosites
    o <- as.data.frame(table(unlist(sites.unique)))
    if (length(o$Freq == m) > 1000) {
        top <- as.character(o[which(o$Freq == m), 1])
    } else {
        top <- as.character(o[order(o$Freq, decreasing = TRUE), 1][1:1000])
    }
    
    Ts <- data.frame(mat.max[[1]][top])
    for (i in 2:n) {
        Ts <- cbind(Ts, mat.max[[i]][top])
    }
    
    # rank phosphosites by absolute fold change
    Tc <- (apply(-abs(Ts),2,rank)-0.5)/nrow(Ts)
    
    # overlapping
    Tt4 <- pchisq(-2*rowSums(log(Tc)), (n-1)*2, lower.tail = FALSE)
    
    sites.sorted <- names(sort(Tt4, decreasing=TRUE))
    return(sites.sorted[1:num])
}

