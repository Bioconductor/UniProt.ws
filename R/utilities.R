#########################################################################
## Helper functions for UniProt.ws
## These functions translate UniProt taxon names to scientific names, 
## taxids, or domain codes.
## Contributed by Csaba Ortutay; csaba.ortutay@gmail.com, 17.10.2016
#########################################################################

.getSpecfile <-
    function(url)
{
    cache <- rappdirs::user_cache_dir(appname="UniProt.ws")
    bfc <- BiocFileCache(cache, ask=FALSE)

    rid <- tryCatch({
        rid0 <- bfcquery(bfc, url, "fpath")$rid
        if (length(rid0)) {
            update <- bfcneedsupdate(bfc, rid0)
            if (!isFALSE(update))
                bfcdownload(bfc, rid0, ask = FALSE)
        } else
            rid0 <- names(bfcadd(bfc, url))
        rid0
    }, error = function(...) integer())

    if (length(rid))
        bfcrpath(bfc, rids=rid)
    else
        system.file("extdata", "speclist.txt", package="UniProt.ws")
}

## digest specfile
.parseSpecfile <-
    function(specfile)
{
    my.text <- readLines(specfile)

    pattern <- "^([[:alnum:]]+) +([[:alnum:]]) +([[:digit:]]+): N=(.*)"
    codetable <- my.text[grepl(pattern, my.text)]

    codes <- data.frame(
        domain = factor(sub(pattern, "\\2", codetable)),
        taxId = as.integer(sub(pattern, "\\3", codetable)),
        taxname = sub(pattern, "\\4", codetable),
        row.names = sub(pattern, "\\1", codetable),
        stringsAsFactors = FALSE
    )
    codes$species <- sub(" +\\(.*", "\\1", codes$taxname)
    codes
}

digestspecfile <- local({
    db <- new.env(parent=emptyenv())
    function(specfile) {
        if (missing(specfile)) {
            specfile <- "http://www.uniprot.org/docs/speclist.txt"
            if (is.null(db[[specfile]])) {
                rsrc <- .getSpecfile()
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
    specnames <- codetable[taxname,"species"]
    specnames
}

# Converting UniProt taxonomy names to NCBI Taxonomy IDs: taxname2taxid()
# 
taxname2taxid  <- function(taxname, specfile) {
    codetable <- digestspecfile(specfile)
    taxids <- codetable[taxname,"taxId"]
    taxids
}

# Converting UniProt taxonomy names to taxonomical domains: taxname2domain(). This function helps
# to map those taxon names to these domains: 
#   'A' for archaea (=archaebacteria)
#   'B' for bacteria (=prokaryota or eubacteria)
#   'E' for eukaryota (=eukarya)
#   'V' for viruses and phages (=viridae)
#   'O' for others (such as artificial sequences)

taxname2domain <- function(taxname, specfile) {
    codetable <- digestspecfile(specfile)
    domains <- codetable[taxname,"domain"]
    domains
}

# Download the latest version of species definition file from UnoProt:
# updatespecfile().  UniProt.ws package has an archived version of the species
# definition file downloaded from http://www.uniprot.org/docs/speclist.  This
# updatespecfile helper function attempts to download the current version of
# conversion table, but if it fails to do that still you can use the archived
# copy saved to the extdata directory.

updatespecfile <- function() {
    .getSpecfile()
}
