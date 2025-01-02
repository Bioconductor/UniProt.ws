.dotter <- function(ndots, maxlength) {
    paste0(
        paste0(rep(".", times = ndots), collapse = ""),
        paste0(rep(" ", times = maxlength-ndots), collapse = ""),
        collapse = ""
    )
}

.UNIPROT_REST_URL <- "https://rest.uniprot.org/"

.checkResponse <- function(jobId) {
    response <- request(.UNIPROT_REST_URL) |>
        req_template("idmapping/status/{jobId}") |>
        req_perform() |>
        resp_body_json()
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
#' @importFrom httr2 request req_template req_perform resp_body_string
#'
#' @export
allFromKeys <- function() {
    results <- request(.UNIPROT_REST_URL) |>
        req_template("configure/idmapping/fields") |>
        req_perform() |>
        resp_body_string()
    allnames <- jmespath(
        results,
        paste0("groups[].items[?from==`true`].name[]")
    )
    sort(unlist(parse_json(allnames)))
}

#' @rdname mapUniProt
#' @export
allToKeys <- function(fromName = "UniProtKB_AC-ID") {
    results <- request(.UNIPROT_REST_URL) |>
        req_template("configure/idmapping/fields") |>
        req_perform() |>
        resp_body_string()
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
    results <- request(.UNIPROT_REST_URL) |>
        req_template("configure/uniprotkb/result-fields") |>
        req_perform() |>
        resp_body_string()
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
    redurl
}

.prepQuery <- function(columns, format = "tsv", paginate, pageSize) {
    qlist <- list(format = format)
    if (length(columns))
        qlist <- c(qlist, fields = paste(columns, collapse = ","))
    if (paginate)
        qlist <- c(qlist, size = pageSize)
    qlist
}

.extract_link <- function(linkElement) {
    gsub("^<(.*)>.*", "\\1", linkElement)
}

#' @importFrom httr2 resp_header
#' @importFrom utils read.delim head
.resp_bind_pages <- function(response, n = Inf) {
    rdata <- resp_body_string(response) |>
        read.delim(text = _)
    pb <- progress_bar$new(
        format = "  (:spin) binding paginated requests:dots :elapsedfull",
        total = NA, clear = FALSE
    )
    on.exit(pb$terminate())
    step <- 1L
    while (length(resp_header(response, "link")) && NROW(rdata) < n) {
        response <- .extract_link(
                resp_header(response, "link")
            ) |>
            request() |>
            req_perform()
        result <- resp_body_string(response) |>
            read.delim(text = _)
        rdata <- do.call(rbind.data.frame, list(rdata, result))
        step <- step + 1
        ndots <- step %% 11
        pb$tick(tokens = list(dots = .dotter(ndots, 10)))
    }
    head(rdata, n)
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
#' @importFrom httr2 req_body_multipart resp_body_json req_url_query
#'
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
#' queryUniProt(
#'     query = c("organism_id:9606", "gene_exact:A2M"),
#'     fields = c(
#'         "id", "accession", "gene_primary",
#'         "organism_name", "protein_name", "reviewed"
#'     ),
#'     collapse = " OR ", n = 3, pageSize = 3
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
    resp <- request(.UNIPROT_REST_URL) |>
        req_template("idmapping/run") |>
        req_body_multipart(
            ids = query[["ids"]], from = from, to = to
        ) |>
        req_perform() |>
        resp_body_json()
    jobId <- resp[["jobId"]]
    if (verbose)
        message("ID Mapping jobId: ", jobId)
    pb <- progress_bar$new(
        format = "  (:spin) waiting for query completion:dots :elapsedfull",
        total = NA, clear = FALSE
    )
    on.exit(pb$terminate())

    while (.checkResponse(jobId)) {
        for (ndot in seq(0, 10)) {
            pb$tick(tokens = list(dots = .dotter(ndot, 10)))
            Sys.sleep(2/8)
        }
    }

    details <- request(.UNIPROT_REST_URL) |>
        req_template("idmapping/details/{jobId}") |>
        req_perform() |>
        resp_body_json()

    if (length(columns))
        columns <- paste(columns, collapse = ",")
    if (!paginate)
        pageSize <- NULL

    .getResultsURL(details[["redirectURL"]], paginate) |>
        request() |>
        req_url_query(
            format = "tsv",
            fields = columns,
            pageSize = pageSize
        ) |>
        req_perform() |>
        .resp_bind_pages(n = Inf)
}

#' @rdname mapUniProt
#'
#' @importFrom BiocBaseUtils isCharacter
#"
#' @export
queryUniProt <- function(
    query = character(0L), fields = c("accession", "id"),
    collapse = c(" OR ", " AND "),
    n = Inf, pageSize = 25L
) {
    stopifnot(isCharacter(query), isCharacter(fields))
    if (!length(query))
        stop("<internal> 'qlist' must be populated with queries")

    collapse <- match.arg(collapse)
    request(.UNIPROT_REST_URL) |>
        req_template("uniprotkb/search") |>
        req_url_query(
            query = paste(query, collapse = collapse),
            fields = paste(fields, collapse = ","),
            format = "tsv",
            size = pageSize
        ) |>
        req_perform() |>
        .resp_bind_pages(n = n)
}
