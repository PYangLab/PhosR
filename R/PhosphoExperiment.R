### =========================================================================
### PhosphoExperiment object
### -------------------------------------------------------------------------
###
#' @export
#' @rdname PhosphoExperiment
#' @importFrom utils packageVersion
#' @importFrom S4Vectors SimpleList
#' @importFrom stats setNames
#' @importClassesFrom S4Vectors DataFrame SimpleList
#' 
setClass("PhosphoExperiment",
         slots=list(Quantification="matrix", 
                    UniprotID="character", 
                    GeneSymbol="character", 
                    Site="numeric", 
                    Residue="character", 
                    Sequence="character")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters/setters for internal slots
###
#' @export
setGeneric("Quantification", function(x, type, ...) standardGeneric("Quantification"))

#' @export
setGeneric("Quantification<-", function(x, type, ..., value) standardGeneric("Quantification<-"))

#' @export
setGeneric("UniprotID", function(x) standardGeneric("UniprotID"))

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
#' @rdname UniprotID
#' @importFrom PhosphoExperiment UniprotID <- UniprotID

setReplaceMethod("UniprotID", signature="PhosphoExperiment", function(x, value) {
        x@UniprotID <- as.character(value)
        return(x)
})

#' @export
#' @rdname GeneSymbol
#' @importFrom PhosphoExperiment GeneSymbol <- GeneSymbol

setReplaceMethod("GeneSymbol", signature="PhosphoExperiment", function(x, value) {
    x@GeneSymbol<- as.character(value)
    return(x)
})

#' @export
#' @rdname Site
#' @importFrom PhosphoExperiment Site <- Site

setReplaceMethod("Site", signature="PhosphoExperiment", function(x, value) {
    x@Site <- as.numeric(value)
    return(x)
})

#' @export
#' @rdname Residue
#' @importFrom PhosphoExperiment Residue <- Residue

setReplaceMethod("Residue", signature="PhosphoExperiment", function(x, value) {
    x@Residue <- as.character(value)
    return(x)
})

#' @export
#' @rdname Sequence
#' @importFrom PhosphoExperiment Sequence <- Sequence

setReplaceMethod("Sequence", signature="PhosphoExperiment", function(x, value) {
    x@Sequence<- as.character(value)
    return(x)
})

PhosphoExperiment <- function(..., UniprotID=c(), GeneSymbol=c(), Site=c(), Residue=c(), Sequence=c()) {
    
    out <- new("PhosphoExperiment", 
               Quantification=...
               )
    
    UniprotID(out) <- UniprotID
    GeneSymbol(out) <- GeneSymbol
    Site(out) <- Site
    Residue(out) <- Residue
    Sequence(out) <- Sequence
    
    out
}

#library(PhosR)
#quant <- as.matrix(PhosR::phospho.L6.ratio)
#uniprot <- as.character(sapply(strsplit(rownames(quant),"~"), function(x) x[[2]]))
#symbol <- as.character(sapply(strsplit(rownames(quant),"~"), function(x) x[[2]]))
#site <- as.numeric(gsub("[STY]","",sapply(strsplit(rownames(quant),"~"), function(x) x[[3]])))
#res <- as.character(gsub("[0-9]","",sapply(strsplit(rownames(quant),"~"), function(x) x[[3]])))
#seq <- as.character(sapply(strsplit(rownames(quant),"~"), function(x) x[[4]]))
#
#phosData <- PhosphoExperiment(Quantification = quant, UniprotID = uniprot)
#phosData <- PhosphoExperiment(Quantification = quant, Site= site)


