## This code will generate a DB and a package for one Organism in UniProt.  It
## will do so by just calling select etc. (the locked-down/API part of my code).

## The plan is that I will make a DB with tables for all the major cols,
## (named the same as with fields & also named the same), and one other table
## that will be called Goodies, (with fields named the same as the cols that
## could be expected to come from it).

## Then the helper functions will be modified (slightly) so that they will
## check if a DB exists 1st (helper function).  If it exists, then another
## helper function will get the data from the DB instead and pass it along to
## the regular machinery.  I think that there should also be a message for
## when you are using the DB instead of the web service.  (there should be a
## message REGARDLESS so that the user can see what is happening)

## All packages MUST depend on UniProt.ws (in template), and I will need to
## make UniProt.ws smart enough to: 1) check for DB packages of a particular
## kind and then 2) if one is available of the correct kind, it should fill a
## slot (normally set to NULL) to have the name of that package.  All methods
## should look at this slot when deciding what to do.  3) whenever the user
## changes the taxId, this slot needs to be reset (either to NULL or to
## another other package that is installed).


## There will be an argument to indicate which IDs I want to use (default will
## be all of them) and so I will have the option of making little test DBs
## (and unit tests)


## TODO:
## If there are versions for this data I should store that value in metadata
## I should also add a name for the schema type (UniProt)?

.createTable <- function(con, tableName){

}

.createTables <- function(con){
  ## 1st set up metadata
  sqliteQuickSQL(con, paste("CREATE TABLE IF NOT EXISTS metadata (name",
                            "VARCHAR(80) PRIMARY KEY,value VARCHAR(255))"))

  
  ## then loop through the fields and make a table for each
}

.populateTables <- function(con){
  ## loop though the tables, and for each, call a helper that gets
  ## that data and then populates the associated table
}

## make the DB (not exported)
makeUniProtDbFromUniProtWS <- function(taxId, keys){
##   require(RSQLite)
##   #require(GO.db)  ## TODO: filter the GO IDs.
##   dbName <- .makeDbPkgName(taxId)
##   con <- dbConnect(SQLite(), dbName)
  
##   ## Store taxId and species in table metadata (for posterity)
##   .createTables(con)  ## just makes the tables

##   ## Then make, download and store all the tables.
##   .populateTables(con)  ## gets data and populates the tables 
}



## make the package (exported)
makeUniProtPackageFromUniProtWS <- function(version,
                                            maintainer,
                                            author,
                                            outputDir = getwd(),
                                            taxId = 9606,
                                            keys = keys(UniProt.ws)){
  

}
