#########################################################################
## Helper functions for UniProt.ws
## These functions translate UniProt taxon names to scientific names,
## taxids, or domain codes.
## Contributed by Csaba Ortutay; csaba.ortutay@gmail.com, 17.10.2016
#########################################################################

.getSpecfile <-
    function(url)
{
    cache <- tools::R_user_dir("UniProt.ws", "cache")
    bfc <- BiocFileCache(cache, ask=FALSE)
    rpath <- BiocFileCache::bfcrpath(
        bfc, rnames = url, exact = TRUE, download = TRUE, rtype = "web"
    )
    update <- bfcneedsupdate(bfc, names(rpath))
    if (update)
        bfcdownload(bfc, names(rpath), ask = FALSE)
    rpath
}

## digest specfile
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
            specfile <- paste0(
                "https://ftp.uniprot.org/pub/databases/uniprot/",
                "current_release/knowledgebase/complete/docs/speclist"
            )
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
            stop("'specfile' must be the name of a local file or ",
                 "(advanced use) a 'data.frame' of appropriate format")
        specfile
    }
})

# Utility functions
#
# UniProt uses custom coding of organism names from which protein sequences
# they store. These taxon names are used also in the protein names (not in the
# UniProt IDs!). These functions help to translate those names to standard
# scientific (Latin) taxon names and other useful identifiers.
#
# Converting UniProt taxonomy names to scientific species names:
# taxname2species()
#
taxname2species <- function(taxname, specfile) {
    codetable <- digestspecfile(specfile)
    specnames <- codetable[taxname, "Official (scientific) name" ]
    specnames
}

# Converting UniProt taxonomy names to NCBI Taxonomy IDs: taxname2taxid()
#
taxname2taxid  <- function(taxname, specfile) {
    codetable <- digestspecfile(specfile)
    taxids <- codetable[taxname, "Taxon Node"]
    taxids
}

# Converting UniProt taxonomy names to taxonomical domains: taxname2domain().
# This function helps to map those taxon names to these domains:
#   'A' for archaea (=archaebacteria)
#   'B' for bacteria (=prokaryota or eubacteria)
#   'E' for eukaryota (=eukarya)
#   'V' for viruses and phages (=viridae)
#   'O' for others (such as artificial sequences)

taxname2domain <- function(taxname, specfile) {
    codetable <- digestspecfile(specfile)
    domains <- codetable[taxname, "kingdom"]
    domains
}
