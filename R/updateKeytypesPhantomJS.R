library(RSelenium)
library(wdman)
library(rvest)
library(xml2)

## The function is used to update the keytypes.
## This function assumes you have PhantomJS installed,
## and RSelenium knows how to find it in the path.
## https://cran.r-project.org/web/packages/RSelenium/vignettes/RSelenium-headless.html
## Mac OS: brew install phantomjs
updateKeytypes <-
    function(URL="https://www.uniprot.org/help/api_idmapping",
             verbose=FALSE,
             save=TRUE)
{
    ## Start phantomJS server
    ## Always start on port 4568 instead of default, so
    ## it does not interfere with the commonly used port 4567
    pJS <- wdman::phantomjs(verbose=verbose, port=4568L)
    ## give the binary a moment
    Sys.sleep(5)
    ## open connection with a remote driver
    remDr <- RSelenium::remoteDriver(browserName = "phantomjs",
                                     port = 4568L)
   
    ## Open and navigate to URL
    remDr$open()
    remDr$navigate(URL)
    ## Check if title is similar
    title <- remDr$getTitle(URL)[[1]]
    if (title != "Programmatic access - Mapping database identifiers")
        stop("The webpage has changed, please check the URL")
    
    ## Get table
    doc <- xml2::read_html(remDr$getPageSource()[[1]])
    keytypes <- rvest::html_table(doc, header=TRUE)[[1]]

    ## Omit category names
    keytypes <- na.omit(keytypes)

    ## Stop remDr and phantomJS 
    remDr$close()
    pJS$stop()     

    file <- system.file('extdata', 'keytypes.txt', package='UniProt.ws')
    if (save) {
        write.table(keytypes, file=file, row.names=FALSE, col.names=TRUE)
    }
    
    keytypes
}
