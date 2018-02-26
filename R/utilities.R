#########################################################################
## Helper functions for UniProt.ws
## These functions translate UniProt taxon names to scientific names, 
## taxids, or domain codes.
## Contributed by Csaba Ortutay; csaba.ortutay@gmail.com, 17.10.2016
#########################################################################

.getSpecfile <- function() {
    load(system.file("extdata", "speclist.RData", package="UniProt.ws"))
    codes
}

# digest specfile
digestspecfile <- function(specfile) {
    if (is.character(specfile)) {
        my.text <- readLines(specfile)
        codetable <- my.text[grep("^[A-Z0-9]+ +[ABEVO]", my.text, 
                                  perl=TRUE, value=FALSE)]
        codes <- data.frame(
            domain = substr(codetable, 7, 7),
            taxId = as.numeric(substr(codetable, 8, 15)),
            taxname = sapply(strsplit(codetable, ": N="), '[', 2),
            row.names = 
                sapply(strsplit(codetable, " +[ABEVO] +",perl=TRUE), '[', 1))
        codes$species <- 
            gsub("^([^ ]* [^ ]*) .*$","\\1",codes$taxname, perl=TRUE)
        codes
    } else {
        if (!is(specfile, "data.frame"))
            stop(paste0("'specfile' must be the name of a local file or ",
                        "the cached 'data.frame' in extdata/"))
        specfile
    }
}


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
    if (missing(specfile))
        specfile <- .getSpecfile()
    codetable <- digestspecfile(specfile)
    specnames <- codetable[taxname,"species"]
    specnames
}

# Converting UniProt taxonomy names to NCBI Taxonomy IDs: taxname2taxid()
# 
taxname2taxid  <- function(taxname, specfile) {
    if (missing(specfile))
        specfile <- .getSpecfile()
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
    if (missing(specfile))
        specfile <- .getSpecfile()
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
    specfile <- .getSpecfile() 
    # Check if speclist.txt is available online, if yes, download, if not, 
    # fall back to local copy
    localfile <- tempfile("specfile.")
    specurl <- 'http://www.uniprot.org/docs/speclist.txt'
    test <- getBinaryURL(specurl)
    if(length(test)>1) {
        writeBin(test,localfile)
        specfile <- localfile
    }
    specfile
}

## Update specfile and save it as a data.frame in
## inst/extdata/specfile.RData
processSpecFile <-
    function()
{
    file <- updatespecfile()
    codes <- digestspecfile(file)
    rda <- system.file("extdata",
                       "speclist.RData",
                       package="UniProt.ws")
    save(codes, file=rda)
}

.cacheNeedsUpdate <- function(url) {
    message("updating resource from ", url)
    needsUpdate <- TRUE
    cache <- rappdirs::user_cache_dir(appname="UniProt.ws")

    tryCatch({
        bfc <- BiocFileCache::BiocFileCache(cache, ask=FALSE)
        query <- bfcquery(bfc, url, "rname")

        if (nrow(query) == 0L) {
            file <- bfcnew(bfc, url)
            needsUpdate <- TRUE
        } else {
            file <- query$rpath
            id <- query$rid
            mtime <- file.mtime(query$rpath) 
            expires <- httr::cache_info(httr::HEAD(url))$expires
            needsUpdate <- expires < Sys.Date()
        }
    }, error = function(err) {
        stop(
            "could not connect or cache url ", url,
            "\n  reason: ", conditionMessage(err)
        )
    })

    setNames(needsUpdate, file)
}
