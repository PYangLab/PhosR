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
#' PhosphoExperiment object accessors
#' 
#' @description These are methods for getting for setting accessors of 
#' \code{PhosphoExperiment} object.
#' This provides some convenience for users.
#' 
#' @param value A vector of values to set to respective accessor. See section 
#' \code{Available methods} for more details.
#' @param withDimnames A \code{logical(1)}, indicating whether the names of the 
#' vector should be applied. 
#' @param x A \code{PhosphoExperiment} object to be assigned to.
#' @param ... Ignored for accessors.
#' 
#' @section Available methods:
#'  In the following code snippets, \code{ppe} is a 
#'  \linkS4class{PhosphoExperiment} object.
#' 
#' \describe{
#' \item{\code{UniprotID(ppe)}, \code{UniprotID(ppe) <- value}:}{Get or set a Uniprot ID, where \code{value} is a character 
#' vector}
#' \item{\code{GeneSymbol(ppe)}, \code{GeneSymbol(ppe) <- value}:}{Get or set a gene symbol , where \code{value} is a character 
#' vector}
#' \item{\code{Site(ppe)}, \code{Site(ppe) <- value}:}{Get or set a phosphorylation site, where \code{value} is a 
#' numeric vector}
#' \item{\code{Residue(ppe)}, \code{Residue(ppe) <- value}:}{Get or set a residue of phosphorylation site, where \code{value} is a 
#' character}
#' \item{\code{Sequence(ppe)}, \code{Sequence(ppe) <- value}:}{Get or set a sequence, where \code{value} is a character vector}
#' \item{\code{Localisation(ppe)}, \code{Localisation(ppe) <- localisation}:}{Get or set a localisation score, where \code{localisation} is a numeric 
#' vector}
#' }
#' 
#' @author Taiyun Kim
#' 
#' @examples
#' example(PhosphoExperiment, echo = FALSE)
#' 
#' UniprotID(phosData) <- uniprot
#' head(UniprotID(phosData))
#' 
#' GeneSymbol(phosData) <- symbol
#' head(GeneSymbol(phosData))
#' 
#' Site(phosData) <- site
#' head(Site(phosData))
#' 
#' Residue(phosData) <- res
#' head(Residue(phosData))
#' 
#' Sequence(phosData) <- seq
#' head(Sequence(phosData))
#' 
#' Localisation(phosData) <- rnorm(nrow(phosData))
#' head(Localisation(phosData))
#' 
#' @docType methods
#' @name PPE-accessors
#' @rdname PhosphoExperiment-methods
NULL


#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("UniprotID", function(x, ...) standardGeneric("UniprotID"))

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("UniprotID<-", function(x, value) standardGeneric("UniprotID<-"))

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("GeneSymbol", function(x, ...) standardGeneric("GeneSymbol"))

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("GeneSymbol<-", function(x, value) standardGeneric("GeneSymbol<-"))

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("Site", function(x, ...) standardGeneric("Site"))

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("Site<-", function(x, value) standardGeneric("Site<-"))

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("Residue", function(x, ...) standardGeneric("Residue"))

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("Residue<-", function(x, value) standardGeneric("Residue<-"))

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("Sequence", function(x, ...) standardGeneric("Sequence"))

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("Sequence<-", function(x, value) standardGeneric("Sequence<-"))

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("Localisation", function(x, ...) standardGeneric("Localisation"))

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setGeneric("Localisation<-", function(x, value) standardGeneric("Localisation<-"))


#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setMethod("UniprotID", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@UniprotID
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setMethod("GeneSymbol", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@GeneSymbol
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setMethod("Site", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@Site
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setMethod("Residue", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@Residue
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setMethod("Sequence", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@Sequence
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setMethod("Localisation", "PhosphoExperiment", function(x, withDimnames=TRUE) {
    out <- x@Localisation
    if (withDimnames) {
        names(out) = rownames(x)
    }
    out
})



#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setReplaceMethod("UniprotID", signature="PhosphoExperiment", function(x, value) {
    x@UniprotID <- as.character(value)
    return(x)
})

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setReplaceMethod("GeneSymbol", signature="PhosphoExperiment", 
    function(x, value) {
    x@GeneSymbol<- as.character(value)
    return(x)
})

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setReplaceMethod("Site", signature="PhosphoExperiment", function(x, value) {
    x@Site <- as.numeric(value)
    return(x)
})

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setReplaceMethod("Residue", signature="PhosphoExperiment", function(x, value) {
    x@Residue <- as.character(value)
    return(x)
})

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setReplaceMethod("Sequence", signature="PhosphoExperiment", function(x, value) {
    x@Sequence<- as.character(value)
    return(x)
})

#' @rdname PhosphoExperiment-methods
#' @aliases  
#' UniprotID
#' UniprotID<-
#' GeneSymbol
#' GeneSymbol<-
#' Site
#' Site<-
#' Residue
#' Residue<-
#' Sequence
#' Sequence<-
#' Localisation
#' Localisation<-
#' @export
setReplaceMethod("Localisation", signature="PhosphoExperiment", 
    function(x, value) {
    x@Localisation<- as.numeric(value)
    return(x)
})



#' PhosphoExperiment object subset, combine methods
#' 
#' @description These are methods for combining or subsetting for  
#' \code{PhosphoExperiment} object. This provides some convenience for users.
#' 
#' @section Available methods:
#'  In the following code snippets, \code{ppe1} and \code{ppe2} is a 
#'  \code{PhosphoExperiment} object with matching \code{colData}. 
#'  \code{ppe3} and \code{ppe4} is a \code{PhosphoExperiment} object with 
#'  matching \code{rowData}.
#' 
#' \describe{
#' \item{\code{rbind(ppe1, ppe2)}:}{Combine row-wise}
#' \item{\code{cbind(ppe3, ppe4)}:}{Combine column-wise}
#' }
#' 
#' @param ... In \code{cbind} or \code{rbind}, a \code{PhosphoExperiment} 
#' objects
#' @param deparse.level {See \code{?base::\link[base]{cbind}} for a description 
#' of this argument.}
#' @param x A  \code{PhosphoExperiment} object
#' @param i For \code{[,PhosphoExperiment}, \code{[,PhosphoExperiment<-, i, j} 
#' are subscripts that can act to subset the rows of x
#' @param j For \code{[,PhosphoExperiment}, \code{[,PhosphoExperiment<-, i, j} 
#' are subscripts that can act to subset the columns of x
#' @param drop A \code{logical(1)}, ignored by these methods
#' @param value An object of a class specified in the S4 method signature.
#' 
#' @author Taiyun Kim
#' 
#' @seealso method \code{rbind}, \code{cbind} from 
#' \linkS4class{SummarizedExperiment} object.
#' 
#' @examples
#' example(PhosphoExperiment, echo = FALSE)
#' 
#' n = ncol(phosData)
#' ppe1 = phosData[,seq(round(n/2))]
#' ppe2 = phosData[,-seq(round(n/2))]
#' 
#' ppe = cbind(ppe1, ppe2)
#' identical(ppe, phosData)
#' 
#' ppe[,seq(round(n/2))] = ppe1
#' identical(ppe, phosData)
#' 
#' p = nrow(phosData)
#' ppe1 = phosData[seq(round(p/2)),]
#' ppe2 = phosData[-seq(round(p/2)),]
#' 
#' ppe = rbind(ppe1, ppe2)
#' identical(ppe, phosData)
#' 
#' ppe[seq(round(p/2)),] = ppe1
#' identical(ppe, phosData)
#' 
#' @docType methods
#' @name PPE-operate
#' @rdname PhosphoExperiment-operate
NULL


################################################################################
# Subsetting
################################################################################
## Getting subset
#' @rdname PhosphoExperiment-operate
#' @aliases  
#' rbind
#' cbind
#' [
#' [<-
#' @export
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
#' @rdname PhosphoExperiment-operate
#' @aliases  
#' rbind
#' cbind
#' [
#' [<-
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
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
#' @rdname PhosphoExperiment-operate
#' @aliases  
#' rbind
#' cbind
#' [
#' [<-
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
# By row
#' @rdname PhosphoExperiment-operate
#' @aliases  
#' rbind
#' cbind
#' [
#' [<-
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
#' @param ... Arguments parsed, identical to those used to create 
#' \linkS4class{SummarizedExperiment}.
#' @param UniprotID A character vector of Uniprot ID
#' @param GeneSymbol A character vector of gene symbol
#' @param Site A numeric vector of phosphorylation site
#' @param Residue A character vector of site residue
#' @param Sequence A character vector of sequences
#' @param Localisation A localisation score.
#' 
#' @docType class
#' @aliases 
#' coerce, SummarizedExperiment, PhosphoExperiment-method
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods is as
#' 
#' 
#' @examples 
#' data(phospho_L6_ratio)
#' quant <- as.matrix(phospho.L6.ratio)
#' uniprot <- as.character(sapply(strsplit(rownames(quant),";"), 
#'     function(x) x[[2]]))
#' symbol <- as.character(sapply(strsplit(rownames(quant),";"), 
#'     function(x) x[[2]]))
#' site <- as.numeric(gsub("[STY]","",sapply(strsplit(rownames(quant),";"), 
#'     function(x) x[[3]])))
#' res <- as.character(gsub("[0-9]","",sapply(strsplit(rownames(quant),";"), 
#'     function(x) x[[3]])))
#' seq <- as.character(sapply(strsplit(rownames(quant),";"), 
#'     function(x) x[[4]]))
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


getResidue = function(Seqs) {
    splitSeqs = strsplit(Seqs, "")
    pos = round(nchar(Seqs)/2)
    res = unlist(lapply(seq_along(splitSeqs), function(x) 
        splitSeqs[[x]][pos[[x]]]))
    if (sum(!names(table(res)) %in% c("S", "T", "Y", "X"))) {
        stop("One or more windows are not centered around the phosphosite")
    }
    return(res)
}

#' @importFrom methods new
.se_to_pe = function(se, UniprotID=c(), GeneSymbol=c(), Site=c(), Residue=c(), 
                    Sequence=c(), Localisation = c()) {
    out <- new("PhosphoExperiment", se)
    
    # Warning messages for missing inputs
    if (!length(GeneSymbol))
        warning("GeneSymbol is not specified. This may affect subsequent analysis steps.\n")
    if (!length(Site))
        warning("Site is not specified. This may affect subsequent analysis steps.\n")
    if (!length(Sequence)) {
        warning("Sequence is not specified. This may affect subsequent analysis steps.\n")
        if (!length(Residue))
            warning("Residue is not specified. This may affect subsequent analysis steps.\n")
    } else {
        if (!length(Residue))
            Residue = getResidue(Sequence)
    } 
    
    UniprotID(out) <- UniprotID
    GeneSymbol(out) <- GeneSymbol
    Site(out) <- Site
    Residue(out) <- Residue
    Sequence(out) <- Sequence
    Localisation(out) <- Localisation
    
    
    out
}



