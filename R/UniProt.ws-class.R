#' @name UniProt.ws-class
#'
#' @title UniProt.ws objects and their related methods and functions
#'
#' @description `UniProt.ws` is the base class for interacting with the UniProt
#'   web services from Bioconductor.
#'
#' @details `UniProt.ws` is a class that is used to interact with the UniProt
#'   web services. It makes use of `AnnotationDbi` methods similarly to
#'   `AnnotationDb` objects.
#'
#'   The `UniProt.ws` will be loaded whenever you load the `UniProt.ws` package.
#'   This object will be set up to retrieve information from Homo sapiens by
#'   default, but this value can be changed to any of the species supported by
#'   UniProt.  The `species` and `taxId` methods allow users to see what species
#'   is currently being accessed, and `taxId<-` allows them to change this
#'   value.
#'
#'   `species` shows the genus and species label currently attached to the
#'   `UniProt.ws` objects database.
#'
#'   `taxId` shows the NCBI taxonomy ID currently attached to the `AnnotationDb`
#'   objects database.  Using the equivalently names replace method (`taxId<-`)
#'   allows the user to change the taxon ID, and the species represented along
#'   with it.
#'
#'   `availableUniprotSpecies` is a helper function to list out the available
#'   Species along with their official taxonomy IDs that are available by
#'   UniProt.  Because there are so many species represented at UniProt, there
#'   is also a pattern argument that can be used to restrict the range of things
#'   returned to be only those whose species names match the searth term. Please
#'   remember when using this argument that the Genus is always capitalized and
#'   the species never is.
#'
#'   `lookupUniprotSpeciesFromTaxId` is another helper that will look up the
#'   species of any tax ID that is supported by UniProt.
#'
#' @aliases taxId taxId<-
#'
#' @docType class
#'
#' @inheritParams base::grep
#'
#' @param x,object a `UniProt.ws` object.
#'
#' @param taxId `numeric(1)` a taxonomy identifier
#'
#' @param value `numeric(1)` the new taxId to set
#'
#' @param ... other arguments
#'
#' @return * `species` and `lookupUniprotSpeciesFromTaxId` each return a
#'   character vector of possible values
#'
#'   * `taxId` returns a numeric value that corresponds to the taxonomy ID
#'
#'   * `availableUniprotSpecies` returns a `data.frame`
#'
#' @author Marc Carlson
#'
#' @seealso [UniProt.ws-methods]
#'
#' @examples
#'
#' ## Make a UniProt.ws object
#' up <- UniProt.ws(taxId=9606)
#'
#' ## look at the object
#' up
#'
#' ## get the current species
#' species(up)
#'
#' ## look up available species with their tax ids
#' availableUniprotSpecies("musculus")
#'
#' ## get the current taxId
#' taxId(up)
#'
#' ## look up the species that goes with a tax id
#' lookupUniprotSpeciesFromTaxId(9606)
#'
#' ## set the taxId to something else
#' taxId(up) <- 10090
#' up
#'
#' @exportClass UniProt.ws
.UniProt.ws <- setClass(
    Class = "UniProt.ws",
    slots = c(
        taxId="numeric",
        db="character",
        taxIdUniprots="character",
        organism="character"
    )
)

#' @rdname UniProt.ws-class
#'
#' @export
UniProt.ws <- function(taxId=9606, ...) {
    ## pre-cache taxIdUniprots from the taxId
    taxId <- as.numeric(taxId)
    results <- queryUniProt(
        query = paste0("taxonomy_id:", taxId),
        fields = c("accession", "organism_name"),
        n = 25,
        pageSize = 25
    )
    taxIdUniprots <- results[["Entry"]]
    if (!length(taxIdUniprots))
        stop("No UniProt entries found for 'taxId = ", taxId, "'")
    organism <- unique(results[["Organism"]])
    .UniProt.ws(
        taxId = taxId, taxIdUniprots = taxIdUniprots, organism = organism, ...
    )
}

#' @describeIn UniProt.ws-class Show method for UniProt.ws objects
#'
#' @exportMethod show
setMethod("show", "UniProt.ws", function(object) {
    cat(class(object), "interface object:")
    cat("\nTaxonomy ID:", object@taxId)
    cat("\nSpecies name:", object@organism)
    cat("\nList species with 'availableUniprotSpecies()'\n")
})

## getters
#' @export
setGeneric("taxId", function(x) standardGeneric("taxId"))

#' @describeIn UniProt.ws-class Get the taxonomy ID from a UniProt.ws object
#'
#' @exportMethod taxId
setMethod("taxId", "UniProt.ws", function(x){ x@taxId })

.taxIdUniprots <- function(x){ x@taxIdUniprots }

.extractLink <- function(txt) {
    if (!is.null(txt)) {
        link <- vapply(strsplit(txt, ";"), `[[`, character(1L), 1L)
        txt <- gsub("^<(.*)>$", "\\1", link)
    }
    txt
}

#' @rdname UniProt.ws-class
#'
#' @export
availableUniprotSpecies <- function(pattern="") {
    specfile <- digestspecfile()
    specfile[grepl(pattern, specfile[["Official (scientific) name"]]), ]
}

#' @rdname UniProt.ws-class
#'
#' @export
lookupUniprotSpeciesFromTaxId <- function(taxId){
    specfile <- availableUniprotSpecies()
    res <- specfile[
        specfile[["Taxon Node"]] %in% taxId, "Official (scientific) name"
    ]
    res
}

#' @export
setGeneric(
    "taxId<-", signature="x", function(x, value) standardGeneric("taxId<-")
)

#' @describeIn UniProt.ws-class Set or chnage the taxonomy ID for a UniProt.ws
#'   object
#'
#' @importFrom BiocBaseUtils setSlots
#'
#' @exportMethod taxId<-
setReplaceMethod("taxId", "UniProt.ws", function(x, value) {
    value <- as.numeric(value)
    ## make sure that there is a record for the suggested taxId
    species <- lookupUniprotSpeciesFromTaxId(value)
    if (!length(species))
        stop("No species were found with the given 'taxId'")
    results <- queryUniProt(
        query = paste0("taxonomy_id:", value),
        fields = c("accession", "organism_name"),
        n = 25,
        pageSize = 25L
    )
    setSlots(
        object = x,
        taxId = value,
        taxIdUniprots = results[["Entry"]],
        organism = unique(results[["Organism"]])
    )
})

#' @describeIn UniProt.ws-class Get the species name from a UniProt.ws object
#'
#' @importFrom BiocGenerics species
#'
#' @exportMethod species
setMethod("species", "UniProt.ws", function(object) { object@organism })
