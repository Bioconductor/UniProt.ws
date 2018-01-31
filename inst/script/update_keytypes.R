library(RSelenium)
library(XML)

## URL with js-rendered content to be scraped
URL <- "https://www.uniprot.org/help/api_idmapping"

## The function is used to update the keytypes.
## This function assumes you have docker installed,
## and this docker container is running,
## docker run -d -p 4445:4444 selenium/standalone-firefox:2.53.1
## You can use this function only if the docker container is running.
updateKeytypes <-
    function(URL, silent=TRUE)
{

    remDr <- remoteDriver(remoteServerAddr = "localhost" 
                      , port = 4445L
                      , browserName = "firefox"
                        )
    remDr$open(silent=silent) 
    ## ## Get status
    ## remDr$getStatus()
    ## ## Navigate URL
    remDr$navigate(URL)
    ## Check if webpage is the same
    title <- remDr$getTitle(URL)[[1]]
    if (title != "Programmatic access - Mapping database identifiers")
        stop("The webpage has changed, please check the URL")
    
    ## Get table
    doc <- htmlParse(remDr$getPageSource()[[1]])

    keytypes <- readHTMLTable(doc, header=TRUE)[[1]]

    ## Omit category names
    keytypes <- na.omit(keytypes)

    ## Close and stop server
    remDr$close()
    rD$server$stop()
    
    keytypes
}
