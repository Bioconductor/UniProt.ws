## Here we will define the select methods.
setMethod("keytypes", "UniProt.ws", function(x) {
    allToKeys(fromName = "UniProtKB_AC-ID")
})

setMethod("columns", "UniProt.ws", function(x) {
    returnFields()
})

## http://www.UniProt.org/UniProt/?query=organism:9606&format=tab&columns=id,sequence

setMethod("keys", "UniProt.ws", function(x, keytype) {
    if (missing(keytype))
        stop("Please supply a keytype argument.")
    if (!any(keytypes(x) %in% keytype))
        stop("keytype argument MUST match a value returned by keytypes method")
    dat <- taxIdUniprots(x) ## pre-cached
    if (identical(keytype, "UniProtKB")) {
        dat
    } else {
        ## then convert this to be the keytype requested...
        dat2 <- mapUniProt(from="UniProtKB_AC-ID", to=keytype, query=dat)
        unique(dat2[["To"]])
    }
})

OLD_IDS <- c("ACC+ID", "ENTREZ_GENE", "GeneID")

## Here is the business end of my select method.
## The big plan is to call mapUniProt() and getUniProtGoodies()
## (merging when necessary)
.select <- function(x, keys, cols, keytype){
  if (!keytype %in% keytypes(x)) {
      stop("'keytype' must be one of 'keytypes(x)'")
  }
  cols <- cols[!cols %in% keytype]  ## remove keytype from cols
  if (!length(cols))
      stop("'columns' should be different from 'keytype'")
  hasOLDID <- cols %in% OLD_IDS
  oldids <- paste0(cols[hasOLDID], collapse = ", ")
  if (any(hasOLDID))
      stop("Unsupported identifiers -\n ", oldids,
          "\n See https://www.uniprot.org/help",
          call. = FALSE
      )
  if (!"accession" %in% cols)
      cols <- c("accession", cols)
  if (identical(keytype, "UniProtKB"))
      keytype <- "UniProtKB_AC-ID"
  if (!keytype %in% allFromKeys())
      stop("'", keytype, "' is not a valid 'from' key")
  keys <- list(organism_id = x@taxId, ids = paste(keys, collapse = ","))
  if ("clusters" %in% tolower(cols) && !identical(keytype, "UniProtKB")) {
      cols <- cols[tolower(cols) != "clusters"]
      dat <- mapUniProt(
          from = keytype, to = "UniProtKB", query = keys, columns = cols
      )
      dat2 <- mapUniProt(
          from = "UniProtKB_AC-ID", to = "UniRef100", query = dat[["Entry"]]
      )
      dat <- merge(dat, dat2, by.x = "Entry", by.y = "From")
  } else {
      dat <- mapUniProt(
          from = keytype, to = "UniProtKB", query = keys, columns = cols
      )
  }
  .blankToNA <- function(col) {
      gsub(pattern="^$",replacement=NA_character_, col)
  }
  dat[] <- lapply(dat, .blankToNA)
  dat
}

setMethod("select", "UniProt.ws",
    function(x, keys, columns, keytype, ...){
          if (missing(keytype))
              keytype <- "UniProtKB"
          .select(x, keys, columns, keytype)
        }
)
