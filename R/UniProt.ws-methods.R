#' @name UniProt.ws-methods
#'
#' @title UniProt.ws methods from AnnotationDbi
#'
#' @description Various methods from AnnotationDbi such as `select`, `columns`,
#'   `keys`, `keytypes`, and `species` are made available for UniProt.ws
#'   objects.
#'
#' @details In much the same way as an `AnnotationDb` object allows acces to
#'   select for many other annotation packages, `UniProt.ws` is meant to allow
#'   usage of `select` methods and other supporting methods to enable the easy
#'   extraction of data from the UniProt web services.
#'
#'   `select`, `columns` and `keys` are used together to extract data via an
#'   `UniProt.ws` object.
#'
#'   `columns` shows which kinds of data can be returned for the `UniProt.ws`
#'   object.
#'
#'   `keytypes` allows the user to discover which keytypes can be passed in to
#'   `select` or `keys` via the `keytype` argument.
#'
#'   `keys` returns keys for the database contained in the `UniProt.ws` object .
#'   By default it will return the primary keys for the database, which are
#'   UniProtKB keys, but if used with the `keytype` argument, it will return the
#'   keys from that keytype.
#'
#'   `select` will retrieve the data as a data.frame based on parameters for
#'   selected `keys` and `columns` and `keytype` arguments.
#'
#' @param x a `UniProt.ws` object.
#'
#' @param keys `character()` the keys to select records for from the database.
#'   All possible keys are returned by using the `keys` method.
#'
#' @param columns `character()` The columns or kinds of things that can be
#'   retrieved from the database.  As with `keys`, all possible columns are
#'   returned by using the `columns` method.
#'
#' @param keytype `character(1)` The keytype that matches the keys used. For the
#'   `select` methods, this is used to indicate the kind of ID being used with
#'   the keys argument. For the `keys` method this is used to indicate which
#'   kind of keys are desired from `keys`
#'
#' @param ... Additional arguments passed to lower level functions, mainly used
#'   for the `to` input to `mapUniProt`.
#'
#' @returns * `keys`,`columns`,`keytypes`, return a character vector
#'   of possible values
#'
#'   * `select` returns a `data.frame`
#'
#' @importFrom AnnotationDbi keytypes
#'
#' @seealso [UniProt.ws-class]
#'
#' @examples
#' ## Make a UniProt.ws object
#' up <- UniProt.ws(taxId=9606)
#'
#' ## list the possible key types
#' head(keytypes(up))
#'
#' ## list of possible columns
#' head(columns(up))
#'
#' ## list all possible keys of type entrez gene ID
#' egs <- keys(up, "GeneID")
#'
#' ## use select to extract some data
#' res <- select(
#'     x = up,
#'     keys = c("22627","22629"),
#'     columns = c("xref_pdb","xref_hgnc","sequence"),
#'     keytype = "GeneID"
#' )
#' res
#'
#' univals <- c("A0A0C5B5G6", "A0A1B0GTW7", "A0JNW5", "A0JP26", "A0PK11")
#' res <- select(
#'     x = up,
#'     keys = univals,
#'     to = "Ensembl"
#' )
#' res
#'
NULL

#' @describeIn UniProt.ws-methods Get keytypes for a UniProt.ws object
#' @exportMethod keytypes
setMethod("keytypes", "UniProt.ws", function(x) {
    allToKeys(fromName = "UniProtKB_AC-ID")
})

#' @describeIn UniProt.ws-methods Get columns for a UniProt.ws object
#'
#' @importFrom AnnotationDbi columns
#'
#' @exportMethod columns
setMethod("columns", "UniProt.ws", function(x) {
    returnFields()[["name"]]
})

## http://www.UniProt.org/UniProt/?query=organism:9606&format=tab&columns=id,sequence

#' @describeIn UniProt.ws-methods Get keys for a UniProt.ws object
#'
#' @importFrom AnnotationDbi keys
#'
#' @exportMethod keys
setMethod("keys", "UniProt.ws", function(x, keytype) {
    if (missing(keytype))
        stop("Please supply a keytype argument.")
    if (!any(keytypes(x) %in% keytype))
        stop("keytype argument MUST match a value returned by keytypes method")
    dat <- .taxIdUniprots(x) ## pre-cached
    if (!identical(keytype, "UniProtKB")) {
        ## then convert this to be the keytype requested...
        dat <- mapUniProt(from="UniProtKB_AC-ID", to=keytype, query=dat)
        dat <- unique(dat[["To"]])
    }
    dat
})

.OLD_IDS <- c("ACC+ID", "ENTREZ_GENE", "GeneID")

#' @describeIn UniProt.ws-methods Select columns from keys
#'
#' @importFrom AnnotationDbi select
#'
#' @exportMethod select
setMethod("select", "UniProt.ws", function(x, keys, columns, keytype, ...) {
    if (missing(keytype))
        keytype <- "UniProtKB"
    args <- list(...)
    tokey <- args[["to"]]
    if (is.null(tokey))
        tokey <- "UniProtKB"
    if (!keytype %in% keytypes(x)) {
        stop("'keytype' must be one of 'keytypes(x)'")
    }
    if (missing(columns))
        columns <- "accession"
    else if (!"accession" %in% columns)
        columns <- c("accession", columns)
    columns <- columns[!columns %in% keytype]
    if (!length(columns))
        stop("'columns' should be different from 'keytype'")
    hasOLDID <- columns %in% .OLD_IDS
    oldids <- paste0(columns[hasOLDID], collapse = ", ")
    if (any(hasOLDID))
        stop("Unsupported identifiers -\n ", oldids,
            "\n See https://www.uniprot.org/help",
            call. = FALSE
        )
    if (identical(keytype, "UniProtKB"))
        keytype <- "UniProtKB_AC-ID"
    if (!keytype %in% allFromKeys())
        stop("'", keytype, "' is not a valid 'from' key")
    keys <- list(organism_id = taxId(x), ids = paste(keys, collapse = ","))
    if ("clusters" %in% tolower(columns) && !identical(keytype, "UniProtKB")) {
        columns <- columns[tolower(columns) != "clusters"]
        dat <- mapUniProt(
            from = keytype, to = tokey, query = keys, columns = columns
        )
        dat2 <- mapUniProt(
            from = "UniProtKB_AC-ID", to = "UniRef100", query = dat[["Entry"]]
        )
        dat <- merge(dat, dat2, by.x = "Entry", by.y = "From")
    } else {
        dat <- mapUniProt(
            from = keytype, to = tokey, query = keys, columns = columns
        )
    }
    dat[] <- lapply(dat, function(col) { gsub("^$", NA_character_, col) })
    dat
})
