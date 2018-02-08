## The function is used to update the keytypes.  This function assumes
## you have PhantomJS installed, and RSelenium knows how to find it in
## the path.
## 
## https://cran.r-project.org/web/packages/RSelenium/vignettes/RSelenium-headless.html
## 
## Optional: Mac OS: brew install phantomjs

.updateKeytypes <-
    function(url="https://www.uniprot.org/help/api_idmapping",
             file = tempfile(), verbose=FALSE)
{
    message("updating keytype map from ", url)
    loadNamespace("wdman")
    loadNamespace("RSelenium")
    loadNamespace("rvest")

    ## Start phantomJS server
    ## Always start on port 4568 instead of default, so
    ## it does not interfere with the commonly used port 4567
    pJS <- wdman::phantomjs(verbose=verbose, port=4568L)
    ## open connection with a remote driver
    remDr <- RSelenium::remoteDriver(browserName = "phantomjs", port = 4568L)
    ## Open and navigate to url
    remDr$open(silent=TRUE)
    remDr$navigate(url)
    ## Check if title is similar
    title <- remDr$getTitle(url)[[1]]
    if (!identical(title, "Programmatic access - Mapping database identifiers"))
        stop("the webpage has changed, please check the url")
    ## Get table
    doc <- xml2::read_html(remDr$getPageSource()[[1]])
    ## Stop remDr and phantomJS
    remDr$close()
    pJS$stop()

    keytypes <- rvest::html_table(doc, header=TRUE)[[1]]
    ## Remove observations(rows) with "Category" in them,
    keytypes <- keytypes[-grep("^Category:", keytypes$Name),, drop=FALSE]

    ## Manually curate a few entries, because UniProt
    ## doesn't support them.
    keytypes$Name[keytypes$Name == "UniProtKB AC/ID"] = "UniProtKB"
    keytypes$Name[keytypes$Name == "Entrez Gene (GeneID)"] = "Entrez_Gene"
    keytypes$Name <- toupper(keytypes$Name)

    ## replace spaces with underscores
    keytypes$Name <- gsub(" ", "_", x=keytypes$Name, fixed=TRUE)

    ## Direction seems irrelavent
    keytypes$Direction <- keytypes$Abbreviation

    write.table(keytypes, file=file, row.names=FALSE, col.names=TRUE, sep="\t")

    file
}

## A couple of constants (referred to by many functions)
## We will maintain a list of supported key types in file in extdata..
.getKeytypeKeysData <- local({
    keytypeKeysData <- NULL
    function() {
        if (is.null(keytypeKeysData)) {
            fpath <- .updateKeytypes()
            keytypeKeysData <<- read.delim(
                fpath, header=TRUE, stringsAsFactors=FALSE
            )
        }
        keytypeKeysData
    }
})
