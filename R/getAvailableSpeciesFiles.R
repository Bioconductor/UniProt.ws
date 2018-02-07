## Update species list for taxID and formal name and
## write to availSpecies.txt in extdata. 

.updateAvailableSpeciesFiles <-
    function(url="http://www.UniProt.org/docs/speclist.txt",
             file=tempfile())
{
    message("updating available species from ", url)
    ## Read Lines of URL
    dat <- readLines(url)
    ## Remove header and footer
    header_idx <-
        grep("_____ _ _______  ____________________________________________________________",
             dat)
    footer_idx <- grep("9ZZZZ X       1: N=root", dat)
    dat <- dat[seq(header_idx+1, footer_idx-1)]
    ## Remove C= and S= lines
    xx1 <- dat[grep("C=", dat, fixed=TRUE, invert=TRUE)]
    xx2 <- xx1[grep("S=",xx1, fixed=TRUE, invert=TRUE)]
    ## Remove lines which split the file
    ## ==================================================================
    ## (2) "Virtual" codes that regroup organisms at a certain
    ## taxonomic level
    ## ==================================================================
    idx <- grep("===============================", xx2)
    xx3 <- xx2[-(idx[1]:idx[2])]
    ## Remove empty lines
    xx4 <- xx3[xx3!=""]
    
    ## Grab only taxID and formal name
    xx5 <- gsub("^(.+)\\W(.+)\\W(.+): N=(.*)", "\\3\\\t\\4", xx4)

    ## Write file to availSpecies.txt
    fileConn <- file(file)
    writeLines(xx5,fileConn)
    close(fileConn)

    ## Return file
    file
}


## A couple of constants (referred to by many functions)
## We will maintain a list of supported key types in file in extdata..
.getAvailableSpeciesData <- local({
    species <- NULL
    function() {
        if (is.null(species)) {
            fpath <- .updateAvailableSpeciesFiles()
            species <<- read.delim(
                fpath, header=FALSE, stringsAsFactors=FALSE
            )
        }
        species
    }
})
