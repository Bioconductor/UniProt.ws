## Some code to make the string into a data.frame...
.cleanup <- function(str, from, to){
  res <- read.delim(text = gsub("[\t]+", "\t", readLines(textConnection(str)), perl = TRUE), sep = "\t")
  res <- res[,c(1, 2)]
  colnames(res) <- c(from, to)
  res
}

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

.getResponse <- function(jobId) {
    url <- paste0("https://rest.uniprot.org/idmapping/status/", jobId)
    resp <- GET(url = url, accept_json())
    content(resp, as = "parsed")
}

.checkResponse <- function(response) {
    if (!is.null(response[["messages"]]))
        message(response[["messages"]])
    if (!is.null(response[["failedIds"]]))
        message(
            "IDs not mapped: ", paste(response[["failedIds"]], collapse = ", ")
        )
    is.null(response[["results"]])
}

allFromKeys <- function() {
    results <- content(
        httpcache::GET("https://rest.uniprot.org/configure/idmapping/fields",
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
        httpcache::GET("https://rest.uniprot.org/configure/idmapping/fields",
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

getStreamURL <- function(redirectURL, debug) {
    url <- gsub(
        "/idmapping/results/", "/idmapping/stream/", redirectURL, fixed = TRUE
    )
    url <- gsub("/results/", "/results/stream/", url, fixed = TRUE)
    .messageDEBUG(url, debug)
}

.prepQuery <- function(columns, format = "tsv") {
    qlist <- list(format = format)
    if (length(columns))
        qlist <- c(qlist, fields = paste(columns, collapse = ","))
    qlist
}

.messageDEBUG <- function(url, debug) {
    if (debug)
        message("Hitting: ", url)
    url
}

mapUniprot <- function(
    from = "UniProtKB_AC-ID", to = "UniRef90",
    columns = character(0L), query, verbose = FALSE, debug = FALSE
) {
    stopifnot(
        isScalarCharacter(from), isScalarCharacter(to), isCharacter(query),
        isTRUEorFALSE(verbose)
    )
    files <- list(ids = paste(query, collapse = ","), from = from, to = to)
    resp <- POST(
        url = .messageDEBUG("https://rest.uniprot.org/idmapping/run", debug),
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

    url <- paste0("https://rest.uniprot.org/idmapping/details/", jobId)
    resp <- GET(url = .messageDEBUG(url, debug), accept_json())
    details <- content(resp, as = "parsed")
    results <- GET(
        url = getStreamURL(details[["redirectURL"]], debug),
        query = .prepQuery(columns),
        accept_json()
    )
    read.delim(text = content(results, encoding = "UTF-8"))
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

## helper to fill back in missing cols.
backFillCols <- function(tab, cols){
  ## 1st we need to translate cols to be the expected headers for tab.
  ecols <- extraColsDat[,3][match(cols, extraColsDat[,2])]
  ## Get vector with NAs where we need replacement cols
  ind = match(ecols,colnames(tab))
  ## Make a blank col
  blank <- data.frame(val=rep(NA,times=dim(tab)[1]))
  ## then loop to place it whenever needed.
  res <- data.frame()
  for(i in seq_len(length(ind))){
    if(i==1){
      res <- tab[,1,drop=FALSE] ## 1st one is always the ids
    }
    else{
      if(!is.na(ind[i])){
        res <- cbind(res,tab[,ind[i],drop=FALSE])
      }else{
        res <- cbind(res, blank)
      }
    }
  }
  colnames(res) <- ecols
  res
}


## A function that take UniProt IDs and gets supplementary cols back
.getSomeUniprotGoodies  <- function(query, cols){
    message(
        "Getting extra data for ",
        paste(head(query, 3), collapse=", "),
        if (length(query) > 3)
            paste0("... (", length(query), " total)")
    )
  ## query and cols start as a character vectors
  qstring <- paste(query, collapse="+or+")
  cstring <- paste(cols, collapse=",")
  url <- 'https://www.uniprot.org/uniprot/?query='
  fullUrl <- paste0(url,qstring,'&format=tab&columns=id,',cstring)
  ## This step may need to repeat (in the event that it fails).
  dat <- .tryReadResult(fullUrl)
  ## read.delim will name mangle if colnames have repeats or [CC]:
  colnames(dat) <- sub("\\.\\d","",colnames(dat))
  colnames(dat) <- sub("\\.\\.CC\\.", "", colnames(dat))
  ## now remove things that were not in the specific original query...
  dat <- dat[dat[,1] %in% query,,drop=FALSE]
  if(dim(dat)[2]< (length(cols)+1)){## we have some empty cols.
    dat <- backFillCols(dat, cols=c("id",cols))
  }
  dat
}

getUniprotGoodies <- function(query, cols){
  dataNibbler(query=query, FUN=.getSomeUniprotGoodies, chnkSize=400, cols=cols)
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

## potential values for columns (for use by getUniprotGoodies)  =
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








