## Species List URL
## Update species list for taxID and formal name
## and write to availSpecies.txt in extdata
processAvailableSpeciesFiles <-
    function(URL="http://www.UniProt.org/docs/speclist.txt")
{
    ## Read Lines of URL
    dat <- readLines(URL)
    
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
    ## =======================================================================
    ## (2) "Virtual" codes that regroup organisms at a certain taxonomic level
    ## =======================================================================
    idx <- grep("===============================", xx2)
    xx3 <- xx2[-(idx[1]:idx[2])]
    ## Remove empty lines
    xx4 <- xx3[xx3!=""]
    
    ## Grab only taxID and formal name
    xx5 <- gsub("^(.+)\\W(.+)\\W(.+): N=(.*)", "\\3\\\t\\4", xx4)

    ## Write file to availSpecies.txt
    fileConn <- file(
        system.file("extdata",
                    "availSpecies.txt",
                    package="UniProt.ws")
    )
    writeLines(xx5,fileConn)
    close(fileConn)
}
