.dotter <- function(ndots, maxlength) {
    paste0(
        paste0(rep(".", times = ndots), collapse = ""),
        paste0(rep(" ", times = maxlength-ndots), collapse = ""),
        collapse = ""
    )
}

.UNIPROT_REST_URL <- "https://rest.uniprot.org/"

#' @importFrom httr GET accept_json content
.getResponse <- function(jobId) {
    url <- paste0(.UNIPROT_REST_URL, "idmapping/status/", jobId)
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

#' @rdname mapUniProt
#'
#' @importFrom rjsoncons jmespath
#' @importFrom jsonlite parse_json
#' @importFrom httr content content_type
#'
#' @export
allFromKeys <- function() {
    results <- content(
        httpcache::GET(
            paste0(.UNIPROT_REST_URL, "configure/idmapping/fields"),
            content_type("application/json")
        ), as = "text", encoding = "UTF-8"
    )
    allnames <- jmespath(
        results,
        paste0("groups[].items[?from==`true`].name[]")
    )
    sort(unlist(parse_json(allnames)))
}

#' @rdname mapUniProt
#' @export
allToKeys <- function(fromName = "UniProtKB_AC-ID") {
    results <- content(
        httpcache::GET(
            paste0(.UNIPROT_REST_URL, "configure/idmapping/fields"),
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

#' @rdname mapUniProt
#' @export
returnFields <- function() {
    results <- content(
        httpcache::GET(
            paste0(.UNIPROT_REST_URL, "configure/uniprotkb/result-fields"),
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

#' @importFrom httr headers
#' @importFrom utils read.delim
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

#' Mapping identifiers with the UniProt API
#'
#' These functions are the main workhorses for mapping identifiers from one
#' database to another. They make use of the latest UniProt API (seen at
#' <https://www.uniprot.org/help/api>).
#'
#' Note that `mapUniProt` is used internally by the `select` method
#' but made available for API queries with finer control. Provide values from
#' the `name` column in `returnFields` as the `columns` input in
#' either `mapUniProt` or `select` method.
#'
#' When using `from='Gene_Name'`, you may restrict the search results to a
#' specific organism by including e.g., `taxId=9606` in the query as a
#' named list element. See examples below.
#'
#' @param from `character(1)` The identifier type to map from, by default
#'   "UniProtKB_AC-ID", short for UniProt accession identifiers.  See a list of
#'   all 'from' type identifiers with `allFromKeys`.
#'
#' @param to `character(1)` The target mapping identifier, by default
#'   "UniRef90". It can be any one of those returned by `allToKeys` from the
#'   appropriate `fromName` argument.
#'
#' @param columns,fields `character()` Additional information to be retreived
#'   from UniProt service.  See a full list of possible input return fields at
#'   <https://www.uniprot.org/help/return_fields>. Example fields include,
#'   "accession", "id", "gene_names", "xref_pdb", "xref_hgnc", "sequence", etc.
#'
#' @param query `character()` or named `list()` Typically, a string that would
#'   indicate the target accession identifiers but can also be a named list
#'   based on the available query fields. See
#'   <https://www.uniprot.org/help/query-fields> for a list of query fields. The
#'   typical query might only include a character vector of UniProt accession
#'   identifiers, e.g., `c("A0A0C5B5G6", "A0A1B0GTW7", "A0JNW5", "A0JP26",
#'   "A0PK11", "A1A4S6")`
#'
#' @param collapse `character(1)` A string indicating either `" OR "` or
#'   `" AND "` for combining `query` clauses.
#'
#' @param n `numeric(1)` Maximum number of rows to return
#'
#' @param fromName `character(1)` A `from` key to use as the basis of mapping to
#'   other keys, by default, `"UniProtKB_AC-ID"`.
#'
#' @param verbose `logical(1)` Whether the operations should provide verbose
#'   updates (default `FALSE`).
#'
#' @param debug `logical(1)` Whether to display the URL API endpoints, for
#'   advanced debugging (default `FALSE`)
#'
#' @param paginate `logical(1)` Whether to use the pagination API (i.e.,
#'   "results" vs "stream") in the request responses. For performance, it is set
#'   to `TRUE` by default.
#'
#' @param pageSize `integer(1)` number of records per page. It corresponds to
#'   the `size` parameter in the API request.
#'
#' @return * `mapUniProt`: A data.frame of returned results
#' * `allToKeys`: A sorted character vector of possible "To" keytypes based
#'   on the given "From" type
#' * `allFromKeys`: A sorted character vector of
#'   possible "From" keytypes
#' * `returnFields`: A `data.frame` of entries for
#'   the columns input in `mapUniProt`; see 'name' column
#'
#' @author M. Ramos
#'
#' @importFrom progress progress_bar
#' @importFrom AnVILBase avstop_for_status
#' @importFrom BiocBaseUtils isScalarCharacter isTRUEorFALSE
#' @examples
#'
#' mapUniProt(
#'     from="UniProtKB_AC-ID",
#'     to='RefSeq_Protein',
#'     query=c('P13368','Q9UM73','P97793','Q17192')
#' )
#'
#' mapUniProt(
#'     from='GeneID', to='UniProtKB', query=c('1','2','3','9','10')
#' )
#'
#' mapUniProt(
#'     from = "UniProtKB_AC-ID",
#'     to = "UniProtKB",
#'     columns = c("accession", "id"),
#'     query = list(organism_id = 10090, ids = c('Q7TPG8', 'P63318'))
#' )
#'
#' ## restrict 'from = Gene_Name' result to taxId 9606
#' mapUniProt(
#'     from = "Gene_Name",
#'     to = "UniProtKB-Swiss-Prot",
#'     columns = c("accession", "id"),
#'     query = list(taxId = 9606, ids = 'TP53')
#' )
#'
#' mapUniProt(
#'     from = "UniProtKB_AC-ID", to = "UniProtKB",
#'     query = c("P31946", "P62258"),
#'     columns = c("accession", "id", "xref_pdb", "xref_hgnc", "sequence")
#' )
#'
#' queryUniProt(
#'     query = c("accession:A5YMT3", "organism_id:9606"),
#'     fields = c("accession", "id", "reviewed"),
#'     collapse = " AND "
#' )
#'
#' allToKeys(fromName = "UniRef100")
#'
#' head(allFromKeys())
#'
#' head(returnFields())
#'
#' @export
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
        url = .messageDEBUG(paste0(.UNIPROT_REST_URL, "idmapping/run"), debug),
        body = files,
        encode = "multipart",
        accept_json()
    )
    avstop_for_status(resp, "idmapping_run")
    submission <- content(resp, as = "parsed")
    jobId <- submission[["jobId"]]
    if (verbose)
        message("ID Mapping jobId: ", jobId)
    pb <- progress_bar$new(
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

    url <- paste0(.UNIPROT_REST_URL, "idmapping/details/", jobId)
    resp <- GET(url = .messageDEBUG(url, debug), accept_json())
    avstop_for_status(resp, "idmapping_details_query")
    details <- content(resp, as = "parsed")
    resurl <- .getResultsURL(details[["redirectURL"]], paginate, debug)
    results <- GET(
        url = resurl,
        query = .prepQuery(columns, pageSize = pageSize, paginate = paginate),
        accept_json()
    )
    avstop_for_status(results, "redirectURL_query")
    .handleResults(results, debug)
}

#' @rdname mapUniProt
#'
#' @importFrom BiocBaseUtils isCharacter
#"
#' @export
queryUniProt <- function(
    query = character(0L), fields = c("accession", "id"), collapse = " OR ",
    n = Inf, pageSize = 25L
) {
    stopifnot(isCharacter(query), isCharacter(fields))
    if (!length(query))
        stop("<internal> 'qlist' must be populated with queries")
    .uniprotPages(
        FUN = .searchPaged, query = query, fields = fields,
        collapse = collapse, n = n, pageSize = pageSize
    )
}

#' @importFrom utils txtProgressBar setTxtProgressBar head
.uniprotPages <- function(FUN, ..., n, pageSize) {
    url <- paste0(.UNIPROT_REST_URL, "uniprotkb/search")
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
