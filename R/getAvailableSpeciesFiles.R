## Update species list for taxID and formal name and
## write to availSpecies.txt in extdata.

.updateAvailableSpeciesFiles <-
    function(url="http://www.UniProt.org/docs/speclist.txt",
             file=tempfile())
{
    message("updating available species from ", url)
    pattern <- "^[[:alnum:]]+ +[[:alnum:]] +([[:digit:]]+): N=(.*)"

    dat <- readLines(url)
    dat <- dat[grep(pattern, dat)]       # keep matching lines
    dat <- sub(pattern, "\\1\t\\2", dat) # extract id & species name
    dat <- dat[dat != "1\troot"]         # drop the root node

    writeLines(dat, file)

    file
}


## A couple of constants (referred to by many functions)
## We will maintain a list of supported key types in file in extdata..
.getAvailableSpeciesData <- local({
    species <- NULL
    function() {
        if (is.null(species)) {
            fpath <- .updateAvailableSpeciesFiles()
            species <<- read.delim(fpath, header=FALSE, stringsAsFactors=FALSE)
        }
        species
    }
})
