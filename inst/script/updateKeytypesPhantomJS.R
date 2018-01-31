library(RSelenium)
library(XML)

## URL with js-rendered content to be scraped
URL <- "https://www.uniprot.org/help/api_idmapping"

## The function is used to update the keytypes.
## This function assumes you have PhantomJS installed,
## and RSelenium knows how to find it in the path.
## https://cran.r-project.org/web/packages/RSelenium/vignettes/RSelenium-headless.html
## Mac OS: brew install phantomjs
updateKeytypes <-
    function(URL, silent=TRUE)
{
    pJS <- phantom()
    Sys.sleep(5) # give the binary a moment
    remDr <- remoteDriver(browserName = 'phantomjs')
    
    ## Open and navigate to URL
    remDr$open()
    remDr$navigate(URL)
    
    title <- remDr$getTitle(URL)[[1]]

    
    if (title != "Programmatic access - Mapping database identifiers")
        stop("The webpage has changed, please check the URL")
    
    ## Get table
    doc <- htmlParse(remDr$getPageSource()[[1]])

    keytypes <- readHTMLTable(doc, header=TRUE)[[1]]

    ## Omit category names
    keytypes <- na.omit(keytypes)
    
    remDr$close()
    pJS$stop()     
    
    keytypes
}
