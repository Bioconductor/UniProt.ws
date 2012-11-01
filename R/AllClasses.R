## Class to allow interaction with Uniprot Web services.


setClass("UniProt.ws", representation(taxId="numeric", db="character"),
         prototype(db=character(0)))
## even though it's a singleton, I CANNOT initialize the tempfile() in the
## prototype.  R is somehow creating the stuff that is done in prototype
## during package INSTALLATION.  Which is interesting...

setMethod("show", "UniProt.ws",
          function(object)
          {
            cat("\"", class(object), "\" object:\n", sep="")
            cat("An interface object for UniProt web services")
            cat("\nCurrent Taxonomy ID:\n")
            cat(object@taxId)
            cat("\nCurrent Species name:\n")
            cat(lookupUniprotSpeciesFromTaxId(object@taxId))
            cat("\nTo change Species see: help('availableUniprotSpecies')\n")
          }
)

## getters
setMethod("taxId", "UniProt.ws",
          function(x){x@taxId }
)


setMethod("db", "UniProt.ws",
          function(x){x@db}
)



## setters

.makePkgName <- function(taxId){
  ## now we want pkgs that match species name
  species <- lookupUniprotSpeciesFromTaxId(taxId)
  ## reformat to match use in names
  paste(c("UniProt.ws",unlist(strsplit(species, split=" ")),"db"),
                  collapse=".")
}

.makeDbPkgName <- function(taxId){  ## now we want pkgs that match species name
  species <- lookupUniprotSpeciesFromTaxId(taxId)
  paste(c("UniProt.ws",unlist(strsplit(species, split=" ")),
                       "sqlite"), collapse=".")
}
  
## helper for taxId replace Method.
.getMatchingDbPkg <- function(taxId){
  ## get pkgs and deps
  pkgs <- data.frame(installed.packages())[,c("Package","Depends")]
  ## prefilter/remove pkgs that don't depend on UniProt.ws
  pkgs <- pkgs[grep(pattern="UniProt.ws",x=pkgs$Depends),]
  ## Now make expected PkgName and pkgDbName
  pkgName <- .makePkgName(taxId)
  pkgDbName <- .makeDbPkgName(taxId)
  ## and if that package is in our list then we have it, if not, we don't...
  if(any(pkgName %in% pkgs$Package)){
    return(pkgDbName)
  }else{
    return(character(0))
  }
}

## This method also has to be responsible for checking if there are any dbs
## available.  If there are, then we also need to set the value for db.
setReplaceMethod("taxId", "UniProt.ws",
    function(x, value)
    {
      value <- as.numeric(value)
      ## make sure that there is a record for the suggested taxId
      res <- lookupUniprotSpeciesFromTaxId(value)
      if(is.character(res) && length(res)==1){
        x@taxId <- value
        ## Here: do some checking and set x@db (if possible)
        ## otherwise set x@db = character(0)
        x@db <- .getMatchingDbPkg(value)
      }else(stop("Please verify that this is a legitimate taxId"))
      x 
    }
)

setMethod("species", "UniProt.ws",
          function(x){lookupUniprotSpeciesFromTaxId(x@taxId)}
)


