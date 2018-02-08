#########################################################################
## Helper functions for UniProt.ws
## These functions translate UniProt taxon names to scientific names,
## taxids, or domain codes.
## Contributed by Csaba Ortutay; csaba.ortutay@gmail.com, 17.10.2016
#########################################################################

.getSpecfile <- local({
    codes <- NULL
    function() {
        if (is.null(codes)) {
            fpath <- .updatespecfile()
            fpath <- .digestspecfile(fpath)
            codes <- read.delim(fpath, header=TRUE, stringsAsFactors=FALSE)
            codes$domain <- as.factor(codes$domain)
            codes$taxname <- as.factor(codes$taxname)
            codes$taxId <- as.numeric(codes$taxId)
            codes <<- codes
        }
        codes
    }
})

# digest specfile
.digestspecfile <- function(specfile) {
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

        write.table(codes, specfile,
                    col.names=TRUE, sep="\t", quote=FALSE)
        specfile
    }
}

# Download the latest version of species definition file from UnoProt:
# updatespecfile().  UniProt.ws package has an archived version of the species
# definition file downloaded from http://www.uniprot.org/docs/speclist.  This
# updatespecfile helper function attempts to download the current version of
# conversion table, but if it fails to do that still you can use the archived
# copy saved to the extdata directory.

.updatespecfile <-
    function(url="http://www.uniprot.org/docs/speclist.txt",
             file= tempfile())
{
    ## Check if speclist.txt is available online, if yes, download, if not,
    ## fall back to local copy
    test <- RCurl::getBinaryURL(url)
    if (length(test) > 1) {
        writeBin(test,file)
        specfile <- file
    }
    specfile
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
taxname2species <- function(taxname, specfile) {
    if (missing(specfile))
        specfile <- .getSpecfile()
    specnames <- specfile[taxname,"species"]
    specnames
}

# Converting UniProt taxonomy names to NCBI Taxonomy IDs: taxname2taxid()
taxname2taxid  <- function(taxname, specfile) {
    if (missing(specfile))
        specfile <- .getSpecfile()
    taxids <- specfile[taxname,"taxId"]
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
    domains <- specfile[taxname,"domain"]
    domains
}
