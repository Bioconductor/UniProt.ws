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
    results <- queryUniProt(
        query = paste0("taxonomy_id:", taxId),
        fields = c("accession", "organism_name"),
        n = 25,
        pageSize = 25
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
queryUniProt <- function(
    query = character(0L), fields = c("accession", "id"), collapse = " OR ",
    n = Inf, pageSize = 25L
) {
    stopifnot(isCharacter(query), isCharacter(fields))
    if (!length(query))
        stop("<internal> 'qlist' must be populated with queries")
    .uniprot_pages(
        FUN = .search_paged1, query = query, fields = fields,
        collapse = collapse, n = n, pageSize = pageSize
    )
}

.uniprot_pages <- function(FUN, ..., n, pageSize) {
    url <- paste0(UNIPROT_REST_URL, "uniprotkb/search")
    response <- FUN(url = url, ..., pageSize = pageSize)
    result <- response$results
    bar <- NULL
    while(
        (!is.null(response$headerLink) &&
            grepl("\"next\"", response$headerLink, fixed = TRUE)) &&
        (NROW(result) < n)
    ) {
        response <- FUN(url = response$url, ..., pageSize = pageSize)
        result <- rbind.data.frame(result, response$results)

        if (is.null(bar)) {
          max <- max(min(n, as.numeric(response$totalResults)), 1L)
          bar <- txtProgressBar(max = max, style = 3L)
          on.exit(close(bar))
        }
        setTxtProgressBar(bar, min(NROW(result), n))
    }
    head(result, n)
}

.extract_link <- function(txt) {
    link <- vapply(strsplit(txt, ";"), `[[`, character(1L), 1L)
    gsub("^<(.*)>$", "\\1", link)
}

.search_paged1 <- function(url, query, fields, collapse, pageSize) {
    resp <- httpcache::GET(
        url = url,
        query = list(
            query = paste(query, collapse = collapse),
            fields = paste(fields, collapse = ","),
            format = "tsv",
            size = pageSize
        )
    )

    .stop_for_status(resp, "queryUniProt_paged")

    lst <- as.list(resp)

    resdata <- content(resp, encoding = "UTF-8")
    if (length(resdata))
        results <- read.delim(text = resdata)
    else
        results <- data.frame()

    list(
        url = .extract_link(resp$headers$link),
        headerLink = resp$headers$link,
        totalResults = resp$headers$`x-total-results`,
        results = results
    )
}

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
    setSlots(x,
        taxId = value,
        taxIdUniprots = results[["Entry"]],
        organism = unique(results[["Organism"]])
    )
})

setMethod("species", "UniProt.ws", function(object) { object@organism })
