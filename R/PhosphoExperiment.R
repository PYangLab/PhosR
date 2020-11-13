#' @export
#' @rdname PhosphoExperiment
#' @importFrom S4Vectors SimpleList
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("PhosphoExperiment", 
    slots=list(
        UniprotID="character", 
        GeneSymbol="character", 
        Site="numeric", 
        Residue="character", 
        Sequence="character"
    ),
    contains = "SummarizedExperiment"
    
)




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters/setters for internal slots
###
setGeneric("Quantification", function(x, type, ...) standardGeneric("Quantification"))

setGeneric("Quantification<-", function(x, type, ..., value) standardGeneric("Quantification<-"))

setGeneric("UniprotID", function(x) standardGeneric("UniprotID"))

setGeneric("UniprotID<-", function(x, value) standardGeneric("UniprotID<-"))

setGeneric("GeneSymbol", function(x, ...) standardGeneric("GeneSymbol"))

setGeneric("GeneSymbol<-", function(x, value) standardGeneric("GeneSymbol<-"))

setGeneric("Site", function(x, ...) standardGeneric("Site"))

setGeneric("Site<-", function(x, value) standardGeneric("Site<-"))

setGeneric("Residue", function(x, ...) standardGeneric("Residue"))

setGeneric("Residue<-", function(x, value) standardGeneric("Residue<-"))

setGeneric("Sequence", function(x, ...) standardGeneric("Sequence"))

setGeneric("Sequence<-", function(x, value) standardGeneric("Sequence<-"))


setReplaceMethod("UniprotID", signature="PhosphoExperiment", function(x, value) {
    x@UniprotID <- as.character(value)
    return(x)
})


setReplaceMethod("GeneSymbol", signature="PhosphoExperiment", function(x, value) {
    x@GeneSymbol<- as.character(value)
    return(x)
})


setReplaceMethod("Site", signature="PhosphoExperiment", function(x, value) {
    x@Site <- as.numeric(value)
    return(x)
})


setReplaceMethod("Residue", signature="PhosphoExperiment", function(x, value) {
    x@Residue <- as.character(value)
    return(x)
})


setReplaceMethod("Sequence", signature="PhosphoExperiment", function(x, value) {
    x@Sequence<- as.character(value)
    return(x)
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
PhosphoExperiment <- function(..., UniprotID=c(), GeneSymbol=c(), Site=c(), Residue=c(), Sequence=c()) {
    se <- SummarizedExperiment(...)
    .se_to_pe(
        se,
        UniprotID,
        GeneSymbol,
        Site,
        Residue,
        Sequence
    )
}

#' @importFrom S4Vectors DataFrame SimpleList
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods new
.se_to_pe = function(se, UniprotID=c(), GeneSymbol=c(), Site=c(), Residue=c(), Sequence=c()) {
    
    
    
    out <- new("PhosphoExperiment", se)
    UniprotID(out) <- UniprotID
    GeneSymbol(out) <- GeneSymbol
    Site(out) <- Site
    Residue(out) <- Residue
    Sequence(out) <- Sequence
    
    
    
    out
}



