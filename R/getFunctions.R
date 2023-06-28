.makeChunkVector <- function(chnkSize,query){
  ## how many chunks?
  chnks <- length(query) %/% chnkSize
  ## compute the remainder
  rem <- length(query) - (chnks * chnkSize)
  ## make the factor
  if(length(query) > chnkSize){
    ## make a vector
    v <- rep(1:chnks, each=chnkSize)
    ## and add on remainder
    if(chnks < length(query)/chnkSize ) v <- c(v, rep(chnks+1, each=rem))
  }else{## or we only need the remainder...
    v <- rep(chnks+1, each=rem)
  }
  v
}

.tryToGetAllChunks <- function(res, qs, FUN, ...) {
    ## call FUN for each
    ## res <- lapply(qs, FUN, ...)
    for (i in seq_along(qs)[is.na(res)]) {
        res[[i]] <- tryCatch({
            FUN(qs[[i]], ...) ## call the getter method
        }, error = function(err) {
            message(
                "error while trying to retrieve data in chunk ", i, ":",
                "\n    ", conditionMessage(err),
                "\ncontinuing to try"
            )
            NULL
        })
    }
    res
}

dataNibbler <- function(query, FUN, chnkSize=400, ...){
  ## make the vector for the chunks
  f <- .makeChunkVector(chnkSize, query)
  ## split by f
  qs <- split(query, as.factor(f))
  ## assign all vals in res to be NA
  res <- rep.int(list(NA),length(qs))

  while (anyNA(res)) {
      ## repeat till you get all the answers.
      res <- .tryToGetAllChunks(res, qs, FUN, ...)
  }

  fin <- do.call(rbind, res)

  ## return combined results
  fin
}

.dotter <- function(ndots, maxlength) {
  paste0(
    paste0(rep(".", times = ndots), collapse = ""),
    paste0(rep(" ", times = maxlength-ndots), collapse = ""),
    collapse = ""
  )
}

UNIPROT_REST_URL <- "https://rest.uniprot.org/"

.getResponse <- function(jobId) {
    url <- paste0(UNIPROT_REST_URL, "idmapping/status/", jobId)
    resp <- GET(url = url, accept_json())
    content(resp, as = "parsed")
}

.checkResponse <- function(response) {
    msgs <- response[["messages"]]
    if (!is.null(msgs)) {
        if (grepl("Resource not found", msgs))
            stop(msgs)
        else
            message(response[["messages"]])
    }
    if (!is.null(response[["failedIds"]]))
        warning(
            "IDs not mapped: ",
            paste(response[["failedIds"]], collapse = ", "),
            call. = FALSE
        )
    is.null(response[["results"]])
}

allFromKeys <- function() {
    results <- content(
        httpcache::GET(
            paste0(UNIPROT_REST_URL, "configure/idmapping/fields"),
            content_type("application/json")
        ), as = "text", encoding = "UTF-8"
    )
    allnames <- jmespath(
      results,
      paste0("groups[].items[?from==`true`].name[]")
    )
    sort(unlist(parse_json(allnames)))
}

allToKeys <- function(fromName = "UniProtKB_AC-ID") {
    results <- content(
        httpcache::GET(
            paste0(UNIPROT_REST_URL, "configure/idmapping/fields"),
            content_type("application/json")
        ), as = "text", encoding = "UTF-8"
    )
    from <- jmespath(
        results,
        paste0("groups[].items[?name=='", fromName, "'].from[]|[0]")
    )
    if (identical(from, "false"))
        stop(fromName, " cannot be a 'from' value")
    ruleId <- jmespath(
        results,
        paste0("groups[].items[?name=='", fromName, "'].ruleId[]|[0]")
    )
    tos <- parse_json(
        jmespath(
            results,
            paste0("rules[?ruleId == `", ruleId, "`].tos[]")
        )
    )
    sort(unlist(tos))
}

returnFields <- function() {
    results <- content(
        httpcache::GET(
            paste0(UNIPROT_REST_URL, "configure/uniprotkb/result-fields"),
            content_type("application/json")
        ), as = "text", encoding = "UTF-8"
    )
    gnames <- parse_json(
        jmespath(results, "[].groupName[]"), simplifyVector = TRUE
    )
    glengths <- parse_json(
        jmespath(results, "[].length(fields)"), simplifyVector = TRUE
    )
    labname <- parse_json(
        jmespath(
            results,
            "[].fields[].[label, name]"
        ),
        simplifyVector = TRUE
    )
    labname <- as.data.frame(labname)
    names(labname) <- c("label", "name")
    groupName <- rep(gnames, times = glengths)
    data.frame(groupName = groupName, labname)
}

.getResultsURL <- function(redurl, paginate, debug) {
    if (!paginate) {
        redurl <- gsub(
            "/idmapping/results/", "/idmapping/stream/", redurl, fixed = TRUE
        )
        redurl <- gsub("/results/", "/results/stream/", redurl, fixed = TRUE)
    }
    .messageDEBUG(redurl, debug)
}

.prepQuery <- function(columns, format = "tsv", paginate, pageSize) {
    qlist <- list(format = format)
    if (length(columns))
        qlist <- c(qlist, fields = paste(columns, collapse = ","))
    if (paginate)
        qlist <- c(qlist, size = pageSize)
    qlist
}

.messageDEBUG <- function(url, debug) {
    if (debug)
        message("Hitting: ", url)
    url
}

.handleResults <- function(results, debug) {
    rdata <- read.delim(text = content(results, encoding = "UTF-8"))
    while (length(headers(results)$link)) {
        nextlink <- headers(results)$link
        results <- GET(
            .messageDEBUG(gsub("<(.*)>.*", "\\1", nextlink), debug),
            accept_json()
        )
        result <- read.delim(text = content(results, encoding = "UTF-8"))
        rdata <- do.call(rbind.data.frame, list(rdata, result))
    }
    rdata
}

mapUniProt <- function(
    from = "UniProtKB_AC-ID", to = "UniRef90",
    columns = character(0L), query, verbose = FALSE, debug = FALSE,
    paginate = TRUE, pageSize = 500L
) {
    stopifnot(
        isScalarCharacter(from), isScalarCharacter(to),
        isCharacter(query) || is.list(query), isTRUEorFALSE(verbose)
    )
    if (is.character(query))
        query <- list(ids = paste(query, collapse = ","))
    else if (is.list(query))
        query[["ids"]] <- paste(query[["ids"]], collapse = ",")
    files <- c(query, list(from = from, to = to))
    resp <- httpcache::POST(
        url = .messageDEBUG(paste0(UNIPROT_REST_URL, "idmapping/run"), debug),
        body = files,
        encode = "multipart",
        accept_json()
    )
    submission <- content(resp, as = "parsed")
    jobId <- submission[["jobId"]]
    if (verbose)
      message("ID Mapping jobId: ", jobId)
    pb <- progress::progress_bar$new(
        format = "  (:spin) waiting for query completion:dots :elapsedfull",
        total = NA, clear = FALSE
    )

    while (.checkResponse(.getResponse(jobId))) {
      for (ndot in seq(0, 10)) {
        pb$tick(tokens = list(dots = .dotter(ndot, 10)))
        Sys.sleep(2/8)
      }
      cat("\n")
    }

    url <- paste0(UNIPROT_REST_URL, "idmapping/details/", jobId)
    resp <- GET(url = .messageDEBUG(url, debug), accept_json())
    details <- content(resp, as = "parsed")
    resurl <- .getResultsURL(details[["redirectURL"]], paginate, debug)
    results <- GET(
        url = resurl,
        query = .prepQuery(columns, pageSize = pageSize, paginate = paginate),
        accept_json()
    )
    httr::stop_for_status(results)
    .handleResults(results, debug)
}


## Try five times and give error if all attempts fail.
.tryReadResult <- function(url){
 for (i in 1:5) {
     result <- tryCatch({
         read.delim(URLencode(url), stringsAsFactors=FALSE)
     }, error=function(err) {
         message(
             "reading url",
             "\n    ", URLencode(url),
             "\nfailed on attempt ", i, " of 5"
         )
         NULL
     })
     if (!is.null(result)) return(result)
     Sys.sleep(5)
 }
 stop("no results after 5 attempts; please try again later")
}

## A function that take UniProt IDs and gets supplementary cols back
.getSomeUniProtGoodies  <- function(query, cols){
    ## query and cols start as a character vectors
    if (!all(c("accession", "id") %in% cols))
        cols <- union(c("accession", "id"), cols)
    queryUniProt(query = query, fields = cols, collapse = " OR ")
}

getUniProtGoodies <- function(query, cols){
  dataNibbler(query=query, FUN=.getSomeUniProtGoodies, chnkSize=400, cols=cols)
}

## Need method to return dataFrame of available species.
availableUniprotSpecies <- function(pattern="") {
    specfile <- digestspecfile()
    specfile[grepl(pattern, specfile[["Official (scientific) name"]]), ]
}

## and another method to look up the species name based on the tax ID.
lookupUniprotSpeciesFromTaxId <- function(taxId){
  specfile <- availableUniprotSpecies("")
  res <- specfile[
      specfile[["Taxon Node"]] %in% taxId, "Official (scientific) name"
  ]
  res
}












## important resources:
## query syntax
## http://www.uniprot.org/help/text-search

## how to use the REST site (generally):
## http://www.uniprot.org/faq/28

## query fields (what you can put for the "query=" part)
## http://www.uniprot.org/help/query-fields





## I need a function that can take a UniProt and retrieve associated data.
## something like below...

## these work:
## http://www.uniprot.org/uniprot/?query=P12345&format=tab&columns=id,sequence
## http://www.uniprot.org/uniprot/?query=accession:P30443&format=tab&columns=id,sequence


## This also works (gets ALL the IDs):
## http://www.uniprot.org/uniprot/?query=organism:9606&format=tab&columns=id,sequence



## And this is if you want multiple of the same kind
## http://www.uniprot.org/uniprot/?query="P04217"+or+"P30443"&format=tab&columns=id,sequence

## Notice though that I am getting extra entries.  Are they equivalant or cruft?
## The DEFINITELY are not perfectly equivalent...
## I would like to understand why some requests get me multiple hits, while others only get me one hit????

## http://www.uniprot.org/uniprot/?query="P04217"&format=tab&columns=id,sequence

## So lets look at this example here: "P30443" returns data for it and also some data for another species with the ID "G7ZFG0"  </puzzled>
## http://www.uniprot.org/uniprot/?query="P30443"&format=tab&columns=id,sequence

## But if I just ask for "G7ZFG0" I don't also get the human version "P30443"...
## http://www.uniprot.org/uniprot/?query="G7ZFG0"&format=tab&columns=id,sequence


## What's more: using single quotes seems to change thigns AGAIN (makes it even LESS specific somehow?)

## So here is how it seems to work: Double quotes indicate one level of non-specificity (exactly what I am unsure), and single quotes make it even LESS specific!  So for my purposes, I need to avoid quotes.  One more BUMMER: not using quotes is the SAME as using double quotes (a small amount of non-specificity).  So I get to use a post filter regardless...  :(



## SO I want to format the query like this:
## http://www.uniprot.org/uniprot/?query=P04217+or+P30443&format=tab&columns=id,sequence

## I can also do things like this (thus extracting the stuff from other DBs)
## http://www.uniprot.org/uniprot/?query=P04217+or+P30443&format=tab&columns=id,database%28interpro%29

## potential values for columns (for use by getUniProtGoodies)  =
## RETURN and QUERY FIELDS
## c("citation", "clusters", "comments", "domains", "domain", "ec", "id",
## "entry name", "existence","families","features","genes","go","go-id",
## "interpro","interactor","keywords","keyword-id","last-modified","length",
## "organism","organism-id","pathway","protein names","reviewed","score",
## "sequence","3d","subcellular locations","taxon","tools","version",
## "virus hosts","database(pfam)","database(pdb)")

## RETURN ## https://www.uniprot.org/help/return_fields
# fields
# citation, lit_pubmed_id
# clusters, name
# ec, ec
# entry name, id
# id, accession
# existence, protein_existence
# families, protein_families
# features, feature_count
# genes, gene_names
# go, go
# go-id, go_id
# database(PDB), xref_pdb

# ## QUERY ## https://www.uniprot.org/help/query-fields
# query
# accession, accession
# cluster, uniprot_id
# reviewed, reviewed:true

## NOTE: parameterized column values have to be "expanded" so database became
## database(pfam) and database(pdb)...








