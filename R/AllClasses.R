## Class to allow interaction with Uniprot Web services.
setOldClass("package_version")  ## For S3

.UniProt.ws <- setClass(
    Class = "UniProt.ws",
    slots = c(
        taxId="numeric",
        db="character",
        taxIdUniprots="character",
        organism="character"
    )
)

UniProt.ws <- function(taxId=9606, ...) {
    ## pre-cache taxIdUniprots from the taxId
    taxId <- as.numeric(taxId)
    results <- .queryUniProt(
        qlist = paste0("taxonomy_id:", taxId),
        fields = "accession,organism_name"
    )
    taxIdUniprots <- results[["Entry"]]
    organism <- unique(results[["Organism"]])
    .UniProt.ws(
        taxId = taxId, taxIdUniprots = taxIdUniprots, organism = organism, ...
    )
}

setMethod("show", "UniProt.ws", function(object) {
    cat(class(object), "interface object:")
    cat("\nTaxonomy ID:", object@taxId)
    cat("\nSpecies name:", object@organism)
    cat("\nList species with 'availableUniprotSpecies()'\n")
})

## getters
setMethod("taxId", "UniProt.ws",
          function(x){x@taxId }
)

setMethod("db", "UniProt.ws",
          function(x){x@db}
)

## taxIdUniprots is not intended to be exported.
setMethod("taxIdUniprots", "UniProt.ws",
          function(x){x@taxIdUniprots }
)


## setters

.makePkgName <- function(taxId){
  ## now we want pkgs that match species name
  species <- lookupUniprotSpeciesFromTaxId(taxId)
  ## reformat to match use in names
  paste(c("UniProt.ws",unlist(strsplit(species, split=" ")),"db"),
                  collapse=".")
}

.makeDbPkgName <- function(taxId){  ## now we want pkgs that match species name
  species <- lookupUniprotSpeciesFromTaxId(taxId)
  paste(c("UniProt.ws",unlist(strsplit(species, split=" ")),
                       "sqlite"), collapse=".")
}

## General helper to query UniProt
.queryUniProt <- function(
    qlist = character(0L), fields = c("accession", "id")
) {
    stopifnot(isCharacter(qlist), isCharacter(fields))
    if (!length(qlist))
        stop("<internal> 'qlist' must be populated with queries")
    resp <- GET("https://rest.uniprot.org/uniprotkb/search",
        query = list(
            query = as.list(qlist),
            fields = paste(fields, collapse = ","),
            format = "tsv"
        )
    )
    read.delim(text = content(resp, encoding = "UTF-8"))
}

setReplaceMethod("taxId", "UniProt.ws", function(x, value) {
    value <- as.numeric(value)
    ## make sure that there is a record for the suggested taxId
    species <- lookupUniprotSpeciesFromTaxId(value)
    if (!length(species))
        stop("No species were found with the given 'taxId'")
    results <- .queryUniProt(
        qlist = paste0("taxonomy_id:", value),
        fields = "accession,organism_name"
    )
    setSlots(x,
        taxId = value,
        taxIdUniprots = results[["Entry"]],
        organism = unique(results[["Organism"]])
    )
})

setMethod("species", "UniProt.ws", function(object) { object@organism })
