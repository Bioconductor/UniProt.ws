#' @name utilities
#'
#' @title Translate UniProt taxon names to scientific names, taxids, or domain
#'   codes
#'
#' @description UniProt uses custom coding of organism names from which protein
#'   sequences they store. These taxon names are used also in the protein names
#'   (not in the UniProt IDs!). These functions help to translate those names to
#'   standard scientific (Latin) taxon names and other useful identifiers.
#'
#' * taxname2species(): converts UniProt taxonomy names to
#'   scientific species names
#' * taxname2taxid(): converts UniProt taxonomy names to NCBI Taxonomy IDs
#' * taxname2domain(): converts UniProt taxonomy names to the following
#'   taxonomical domains:
#'     * 'A' for archaea (=archaebacteria)
#'     * 'B' for bacteria (=prokaryota or eubacteria)
#'     * 'E' for eukaryota (=eukarya)
#'     * 'V' for viruses and phages (=viridae)
#'     * 'O' for others (such as artificial sequences)
#'
#' @param taxname Character string up to 6 uppercase characters, like HUMAN,
#'   MOUSE, or AERPX. Also works for a vector of such taxon names.
#'
#' @param specfile An optional local file where speclist.RData is saved from
#'   UniProt.org.  When `specfile` is missing, a cached file from the extdata/
#'   package directory is used.
#'
#' @returns * `taxname2species`: a character vector of scientific taxon names
#'   matching to the UniProt taxon names supplied as `taxname`.
#'   * `taxname2taxid`: a numeric vector of Taxonomy IDs matching to the
#'   UniProt taxon names supplied as `taxname`.
#'   * `taxname2domain`: a character vector of one letter domain
#'   symbols matching to the UniProt taxon names supplied as `taxname`.
#'
#' @author Csaba Ortutay
#'
#' @seealso [UniProt controlled vocabulary of
#'   species](https://www.uniprot.org/docs/speclist.txt), which defines the
#'   taxon names.
#'
#' @examples
#'
#' taxname2species("PIG")
#' taxname2species(c("PIG","HUMAN","TRIHA"))
#'
#' taxname2taxid("PIG")
#' taxname2taxid(c("PIG","HUMAN","TRIHA"))
#'
#' taxname2domain("PIG")
#' taxname2domain(c("PIG","HUMAN","TRIHA"))
#'
NULL

#' @importFrom BiocFileCache BiocFileCache bfcneedsupdate bfcrpath bfcdownload
.getSpecfile <-
    function(url)
{
    cache <- tools::R_user_dir("UniProt.ws", "cache")
    bfc <- BiocFileCache(cache, ask=FALSE)
    rpath <- bfcrpath(
        bfc, rnames = url, exact = TRUE, download = TRUE, rtype = "web"
    )
    update <- bfcneedsupdate(bfc, names(rpath))
    if (update)
        bfcdownload(bfc, names(rpath), ask = FALSE)
    rpath
}

.parseSpecfile <-
    function(specfile)
{
    rlines <- readLines(specfile)
    pattern <- "^([[:alnum:]]+) +([[:alnum:]]) +([[:digit:]]+): N=(.*)"
    codetable <- rlines[grepl(pattern, rlines)]
    os_name <- sub(pattern, "\\4", codetable)

    data.frame(
        row.names = sub(pattern, "\\1", codetable),
        kingdom = factor(sub(pattern, "\\2", codetable)),
        `Taxon Node` = as.integer(sub(pattern, "\\3", codetable)),
        ## removing strain/isolate information in parentheses
        `Official (scientific) name` = sub(" +\\(.*", "\\1", os_name),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

digestspecfile <- local({
    db <- new.env(parent=emptyenv())
    function(specfile) {
        if (missing(specfile)) {
            specfile <- "https://www.UniProt.org/docs/speclist.txt"
            if (is.null(db[[specfile]])) {
                rsrc <- .getSpecfile(specfile)
                db[[specfile]] <- .parseSpecfile(rsrc)
            }
            specfile <- db[[specfile]]
        } else if (is.character(specfile)) {
            if (is.null(db[[specfile]]))
                db[[specfile]] <- .parseSpecfile(specfile)
            specfile <- db[[specfile]]
        }
        if (!is(specfile, "data.frame"))
            stop(
                "'specfile' must be the name of a local file or ",
                "(advanced use) a 'data.frame' of appropriate format"
            )
        specfile
    }
})

#' @rdname utilities
#' @export
taxname2species <- function(taxname, specfile) {
    codetable <- digestspecfile(specfile)
    specnames <- codetable[taxname, "Official (scientific) name" ]
    specnames
}

#' @rdname utilities
#' @export
taxname2taxid  <- function(taxname, specfile) {
    codetable <- digestspecfile(specfile)
    taxids <- codetable[taxname, "Taxon Node"]
    taxids
}

#' @rdname utilities
#' @export
taxname2domain <- function(taxname, specfile) {
    codetable <- digestspecfile(specfile)
    domains <- codetable[taxname, "kingdom"]
    domains
}
