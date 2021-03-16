#' @export
#' @rdname PhosphoExperiment
#' @importFrom S4Vectors SimpleList
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("PhosphoExperiment", 
    slots=representation(
        UniprotID="character", 
        GeneSymbol="character", 
        Site="numeric", 
        Residue="character", 
        Sequence="character",
        Localisation="numeric"
    ),
    contains = "SummarizedExperiment"
)




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters/setters for internal slots
###
# setGeneric("Quantification", function(x, type, ...) standardGeneric("Quantification"))
# 
# setGeneric("Quantification<-", function(x, type, ..., value) standardGeneric("Quantification<-"))

#' @export
setGeneric("UniprotID", function(x, ...) standardGeneric("UniprotID"))

#' @export
setGeneric("UniprotID<-", function(x, value) standardGeneric("UniprotID<-"))

#' @export
setGeneric("GeneSymbol", function(x, ...) standardGeneric("GeneSymbol"))

#' @export
setGeneric("GeneSymbol<-", function(x, value) standardGeneric("GeneSymbol<-"))

#' @export
setGeneric("Site", function(x, ...) standardGeneric("Site"))

#' @export
setGeneric("Site<-", function(x, value) standardGeneric("Site<-"))

#' @export
setGeneric("Residue", function(x, ...) standardGeneric("Residue"))

#' @export
setGeneric("Residue<-", function(x, value) standardGeneric("Residue<-"))

#' @export
setGeneric("Sequence", function(x, ...) standardGeneric("Sequence"))

#' @export
setGeneric("Sequence<-", function(x, value) standardGeneric("Sequence<-"))

#' @export
setGeneric("Localisation", function(x, ...) standardGeneric("Localisation"))

#' @export
setGeneric("Localisation<-", function(x, value) standardGeneric("Localisation<-"))


#' @export
setMethod("UniprotID", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@UniprotID
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})

#' @export
setMethod("GeneSymbol", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@GeneSymbol
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})

#' @export
setMethod("Site", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@Site
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})

#' @export
setMethod("Residue", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@Residue
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})

#' @export
setMethod("Sequence", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@Sequence
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})

#' @export
setMethod("Localisation", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@Localisation
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})



#' @export
setReplaceMethod("UniprotID", signature="PhosphoExperiment", function(x, value) {
    x@UniprotID <- as.character(value)
    return(x)
})

#' @export
setReplaceMethod("GeneSymbol", signature="PhosphoExperiment", function(x, value) {
    x@GeneSymbol<- as.character(value)
    return(x)
})

#' @export
setReplaceMethod("Site", signature="PhosphoExperiment", function(x, value) {
    x@Site <- as.numeric(value)
    return(x)
})

#' @export
setReplaceMethod("Residue", signature="PhosphoExperiment", function(x, value) {
    x@Residue <- as.character(value)
    return(x)
})

#' @export
setReplaceMethod("Sequence", signature="PhosphoExperiment", function(x, value) {
    x@Sequence<- as.character(value)
    return(x)
})

#' @export
setReplaceMethod("Localisation", signature="PhosphoExperiment", function(x, value) {
    x@Localisation<- as.numeric(value)
    return(x)
})

################################################################################
# Subsetting
################################################################################
## Getting subset
#' @importFrom methods callNextMethod
setMethod("[", "PhosphoExperiment", function(x, i, j, drop=TRUE) {
    uID=UniprotID(x, withDimnames = FALSE)
    gs=GeneSymbol(x, withDimnames = FALSE)
    site=Site(x, withDimnames = FALSE)
    residue=Residue(x, withDimnames = FALSE)
    seq=Sequence(x, withDimnames = FALSE)
    loc=Localisation(x, withDimnames = FALSE)
    
    
    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                i, rownames(x), fmt
            )
        }
        i <- as.vector(i)
        uID <- uID[i]
        gs <- gs[i]
        site <- site[i]
        residue <- residue[i]
        seq <- seq[i]
        loc <- loc[i]
    }
    
    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
            j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                j, colnames(x), fmt
            )
        }
        j <- as.vector(j)
    }
    
    out <- methods::callNextMethod()
    BiocGenerics:::replaceSlots(out, UniprotID=uID, GeneSymbol=gs,
        Site=site, Residue=residue, 
        Sequence=seq, Localisation=loc, check=FALSE)
})


## Assigning a subset
setReplaceMethod("[", c("PhosphoExperiment"),
    function(x, i, j, ..., value) {
        uID=UniprotID(x, withDimnames = FALSE)
        gs=GeneSymbol(x, withDimnames = FALSE)
        site=Site(x, withDimnames = FALSE)
        residue=Residue(x, withDimnames = FALSE)
        seq=Sequence(x, withDimnames = FALSE)
        loc=Localisation(x, withDimnames = FALSE)
        
        if (!missing(i)) {
            if (is.character(i)) {
                fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
                i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                    i, rownames(x), fmt
                )
            }
            i <- as.vector(i)
            uID[i] <- UniprotID(value, withDimnames = FALSE)
            gs[i,] <- GeneSymbol(value, withDimnames = FALSE)
            site[,i] <- Site(value, withDimnames = FALSE)
            residue[i] <- Residue(value, withDimnames = FALSE)
            seq[i,] <- Sequence(value, withDimnames = FALSE)
            loc[,i] <- Localisation(value, withDimnames = FALSE)
            
        }
        
        if (!missing(j)) {
            if (is.character(j)) {
                fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
                j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                    j, colnames(x), fmt
                )
            }
            j <- as.vector(j)
            
        }
        
        out <- callNextMethod()
        BiocGenerics:::replaceSlots(out, UniprotID=uID, GeneSymbol=gs,
            Site=site, Residue=residue, 
            Sequence=seq, Localisation=loc, check=FALSE)
    })


################################################################################
# Combining
################################################################################
# By row
#' @export
#' @importFrom BiocGenerics rbind cbind
setMethod("rbind", "PhosphoExperiment", function(..., deparse.level=1) {
    args <- list(...)
    all.uID <- lapply(args, UniprotID, withDimnames=FALSE)
    all.gs <- lapply(args, GeneSymbol, withDimnames=FALSE)
    all.site <- lapply(args, Site, withDimnames=FALSE)
    all.residue <- lapply(args, Residue, withDimnames=FALSE)
    all.seq <- lapply(args, Sequence, withDimnames=FALSE)
    all.loc <- lapply(args, Localisation, withDimnames=FALSE)
    
    all.uID <- do.call(c, all.uID)
    all.gs <- do.call(c, all.gs)
    all.site <- do.call(c, all.site)
    all.residue <- do.call(c, all.residue)
    all.seq <- do.call(c, all.seq)
    all.loc <- do.call(c, all.loc)
    
    # Checks for identical column state.
    # ref <- args[[1]]
    
    
    # No column slots
    # for (x in args[-1]) {
    #     if (!identical(ref.cv, colVec(x, withDimnames=FALSE))
    #         || !identical(ref.ccm, colToColMat(x, withDimnames=FALSE))
    #         || !identical(ref.rcm, rowToColMat(x, withDimnames=FALSE)))
    #     {
    #         stop("per-column values are not compatible")
    #     }
    # }

    old.validity <- S4Vectors:::disableValidity()
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old.validity))
    
    out <- callNextMethod()
    
    BiocGenerics:::replaceSlots(out, UniprotID=all.uID, GeneSymbol=all.gs,
        Site=all.site, Residue=all.residue, 
        Sequence=all.seq, Localisation=all.loc, check=FALSE)
   
})

# By column
#' @export
#' @importFrom BiocGenerics rbind cbind
setMethod("cbind", "PhosphoExperiment", function(..., deparse.level=1) {
    args <- list(...)
    
    
    # Checks for identical column state.
    ref <- args[[1]]
    
    ref.uID <- UniprotID(ref, withDimnames=FALSE)
    ref.gs <- GeneSymbol(ref, withDimnames=FALSE)
    ref.site <- Site(ref, withDimnames=FALSE)
    ref.residue <- Residue(ref, withDimnames=FALSE)
    ref.seq <- Sequence(ref, withDimnames=FALSE)
    ref.loc <- Localisation(ref, withDimnames=FALSE)
    
    for (x in args[-1]) {
        if (!identical(ref.uID, UniprotID(x, withDimnames=FALSE)) 
            || !identical(ref.gs, GeneSymbol(x, withDimnames=FALSE))
            || !identical(ref.site, Site(x, withDimnames=FALSE))
            || !identical(ref.residue, Residue(x, withDimnames=FALSE))
            || !identical(ref.seq, Sequence(x, withDimnames=FALSE))
            || !identical(ref.loc, Localisation(x, withDimnames=FALSE)))
        {
            stop("per-row values are not compatible")
        }
    }
    
    old.validity <- S4Vectors:::disableValidity()
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old.validity))
    
    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, 
        check=FALSE)
})

#'The PhosphoExperiment class
#'
#'
#'
#' @docType class
#' @aliases 
#' coerce, SummarizedExperiment, PhosphoExperiment-method
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods is as
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' 
#' @examples 
#' library(PhosR)
#' quant <- as.matrix(PhosR::phospho.L6.ratio)
#' uniprot <- as.character(sapply(strsplit(rownames(quant),";"), function(x) x[[2]]))
#' symbol <- as.character(sapply(strsplit(rownames(quant),";"), function(x) x[[2]]))
#' site <- as.numeric(gsub("[STY]","",sapply(strsplit(rownames(quant),";"), function(x) x[[3]])))
#' res <- as.character(gsub("[0-9]","",sapply(strsplit(rownames(quant),";"), function(x) x[[3]])))
#' seq <- as.character(sapply(strsplit(rownames(quant),";"), function(x) x[[4]]))
#' phosData <- PhosphoExperiment(assays = list(Quantification = quant), 
#'     UniprotID = uniprot, Site = site, GeneSymbol = symbol, Residue = res, 
#'     Sequence = seq)
#'
PhosphoExperiment <- function(..., UniprotID=c(), GeneSymbol=c(), Site=c(), 
    Residue=c(), Sequence=c(), Localisation=c()) {
    se <- SummarizedExperiment(...)
    .se_to_pe(
        se,
        UniprotID,
        GeneSymbol,
        Site,
        Residue,
        Sequence,
        Localisation
    )
}

#' @importFrom S4Vectors DataFrame SimpleList
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods new
.se_to_pe = function(se, UniprotID=c(), GeneSymbol=c(), Site=c(), Residue=c(), Sequence=c(), Localisation = c()) {
    out <- new("PhosphoExperiment", se)
    UniprotID(out) <- UniprotID
    GeneSymbol(out) <- GeneSymbol
    Site(out) <- Site
    Residue(out) <- Residue
    Sequence(out) <- Sequence
    Localisation(out) <- Localisation
    
    # Warning messages for missing inputs
    if (!length(GeneSymbol))
        warning("GeneSymbol is not specified. This may affect subsequent analysis steps.\n")
    if (!length(Site))
        warning("Site is not specified. This may affect subsequent analysis steps.\n")
    if (!length(Residue))
        warning("Residue is not specified. This may affect subsequent analysis steps.\n")
    if (!length(Sequence))
        warning("Sequence is not specified. This may affect subsequent analysis steps.\n")
    if (!length(Localisation))
        warning("Localisation is not specified. This may affect subsequent analysis steps.\n")
    
    out
}



