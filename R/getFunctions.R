## A couple of constants (referred to by many functions)
## We will maintain a list of supported key types in file in extdata..
keytypeKeysDat <- read.delim(system.file('extdata','keytypes.txt',
                                         package='UniProt.ws')
                             , header=FALSE, stringsAsFactors=FALSE)

## write.table(keytypeKeysDat, file="keytypes2.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

## We also keep a list of supported additional cols (things that can be
## retrieved but not used as keys in a file called extraCols.txt
extraColsDat <- read.delim(system.file('extdata','extraCols.txt',
                                         package='UniProt.ws')
                             , header=FALSE, stringsAsFactors=FALSE)

## FOR NOW: we are not supporting the following 4 cols (they give us the 505)
## Also remember: adjust/comment these in man page...
## keytypeKeysDat <- keytypeKeysDat[-c(37L,38L),]


## Some code to make the string into a data.frame...
.cleanup <- function(str, from, to){
  vec <- unlist(strsplit(str,"\n"))
  vec2 <- unlist(strsplit(vec,"\t"))
  res <- matrix(vec2, ncol=2, byrow=TRUE)
  res <- as.data.frame(res, stringsAsFactors=FALSE)
  colnames(res) <- c(from, to)
  res[-1,] 
}

## Try five times and give error if all attempts fail.
.tryGetResult <- function(url, params){
 for (i in 1:11) {
     result <- tryCatch({
           getForm(url,
                   .params=params,
                   .opts=list(FOLLOWLOCATION=TRUE))
     }, error=function(err) NULL)
     if (!is.null(result)) return(result)
     Sys.sleep(10)
 }
 stop("no results after 5 attempts; please try again later")
}


.mapUni <- function(query, from, to){
  ## query starts as a character vector...
  ## But the URL expects it to be a space separated string.
  query <- paste(query, collapse=" ")
  ## url is constant here
  url <- 'http://www.uniprot.org/mapping/'
  params <- c('from'=from,'to'= to,'format'='tab',
    'query'=query)
##   res <- getForm(url,
##                 .params=params,
##                 .opts=list(FOLLOWLOCATION=TRUE))
  res <- .tryGetResult(url, params)
  .cleanup(res, from, to)
}

.makeChunkVector <- function(chnkSize,query){
  ## how many chunks?
  chnks <- length(query) %/% chnkSize
  ## compute the remainder
  rem <- length(query) - (chnks * chnkSize)
  ## make the factor
  if(length(query) > chnkSize){
    ## make a vector
    v <- rep(1:chnks, each=chnkSize)
    ## and add on remainder
    if(chnks < length(query)/chnkSize ) v <- c(v, rep(chnks+1, each=rem))
  }else{## or we only need the remainder...
    v <- rep(chnks+1, each=rem)
  }
  v
}

.tryToGetAllChunks <- function(res, qs, FUN, ...){
  ## call FUN for each
  ## res <- lapply(qs, FUN, ...)
  for (i in seq_len(length(qs))){
    if(length(res[[i]])==1 && is.na(res[[i]])){
      res[[i]] <- tryCatch({
        FUN(qs[[i]], ...) ## call the getter method
      }, error=function(err) NA)
    }
  }
  res
}

dataNibbler <- function(query, FUN, chnkSize=400, ...){
  ## make the vector for the chunks
  f <- .makeChunkVector(chnkSize, query)  
  ## split by f
  qs <- split(query, as.factor(f))
  ## assign all vals in res to be NA
  res <- rep.int(list(NA),length(qs))
  
  while(any(is.na(res))){
    ## repeat till you get all the answers.
    res <- .tryToGetAllChunks(res, qs, FUN, ...)
  }
  
  fin <- do.call(rbind, res)
  
  ## return combined results
  fin
}


mapUniprot <- function(from, to, query){
  ## 1st we look at that query.  Is is longer than 400 long, then we tend to
  ## get a "bad request" response, so I am simplifying it here to do small
  ## bites and then reassemble them
  message("Getting mapping data for ",query[1]," ... and ",to)
  dataNibbler(query=query, FUN=.mapUni, chnkSize=400,
              from=from, to=to) ## not a typo that from is used twice here.
}


## Try five times and give error if all attempts fail.
.tryReadResult <- function(url){
 for (i in 1:6) {
     result <- tryCatch({
         read.delim(URLencode(url), stringsAsFactors=FALSE)
     }, error=function(err) NULL)
     if (!is.null(result)) return(result)
     Sys.sleep(5)
 }
 stop("no results after 5 attempts; please try again later")
}


## FOR NOW: There doesn't seem to be a need for getOneToMany to be altered to
## use dataNibbler.  So it's a stand-alone thing for now. (only used for
## making the annots at release time)
## What about something simpler?
## http://www.uniprot.org/uniprot/?query=P12345&format=tab&columns=id,sequence
## What about PFAM and Prosite (need special functions for special services).

getOneToMany <- function(taxId, type=c("PFAM","prosite","SMART")){ 
  type <- match.arg(type)
  url <- paste0("http://www.uniprot.org/uniprot/?query=organism:",taxId,"&format=tab&columns=id,database(")
  fullUrl <- paste0(url,type,")")
  message("Reading in data from UniProt web services.")
##   dat <- read.delim(fullUrl, stringsAsFactors=FALSE)
  dat <- .tryReadResult(fullUrl)
  colnames(dat) <- c('ids', 'ids2')
  ## split up the strings
  dat[[2]] <- strsplit( as.character(dat[[2]]), split=";")
  ## get number of things matched to each ID in col 1
  lens <- unlist(lapply(dat[[2]],length))
  ## make factor based on dat[[1]], repeated lens times
  ids <- rep.int(dat[[1]],lens) ## this excludes ones where lens==0
  ids2 <- unlist(dat[[2]])
  if(length(ids)==length(ids2)){
    res <- cbind(ids,ids2)
  }else{
    stop("getOneToMany: ids != ids2")
  }
  ## recover dat[[1]] where lens==0
  rem <- dat[lens==0,]
  rem[[2]] <- NA ## these values are all NA
  rbind(res, rem)
}



## helper to fill back in missing cols.
backFillCols <- function(tab, cols){
  ## 1st we need to translate cols to be the expected headers for tab.  
  ecols <- extraColsDat[,3][match(cols, extraColsDat[,2])]
  ## Get vector with NAs where we need replacement cols
  ind = match(ecols,colnames(tab))
  ## Make a blank col
  blank <- data.frame(val=rep(NA,times=dim(tab)[1]))
  ## then loop to place it whenever needed.
  res <- data.frame()
  for(i in seq_len(length(ind))){
    if(i==1){
      res <- tab[,1,drop=FALSE] ## 1st one is always the ids
    }
    else{
      if(!is.na(ind[i])){
        res <- cbind(res,tab[,ind[i],drop=FALSE])
      }else{
        res <- cbind(res, blank)
      }
    }
  }
  colnames(res) <- ecols
  res
}


## A function that take UniProt IDs and gets supplementary cols back
.getSomeUniprotGoodies  <- function(query, cols){
  message(paste0("Getting extra data for ",
          paste(head(query), collapse=",")))
  ## query and cols start as a character vectors
  qstring <- paste(query, collapse="+or+")  
  cstring <- paste(cols, collapse=",")
  url <- 'http://www.uniprot.org/uniprot/?query='
  fullUrl <- paste0(url,qstring,'&format=tab&columns=id,',cstring)
  ## This step may need to repeat (in the event that it fails).
  dat <- .tryReadResult(fullUrl)
  ## read.delim will name mangle if colnames have repeats or [CC]:
  colnames(dat) <- sub("\\.\\d","",colnames(dat)) 
  colnames(dat) <- sub("\\.\\.CC\\.", "", colnames(dat))
  ## now remove things that were not in the specific original query...
  dat <- dat[dat[,1] %in% query,,drop=FALSE]
  if(dim(dat)[2]< (length(cols)+1)){## we have some empty cols.
    dat <- backFillCols(dat, cols=c("id",cols))
  }
  dat
}

getUniprotGoodies <- function(query, cols){
  dataNibbler(query=query, FUN=.getSomeUniprotGoodies, 
              chnkSize=400, cols=cols)
}



## Need method to return dataFrame of available species.
availableUniprotSpecies <- function(pattern="", n=Inf){
  species <- read.delim(system.file('extdata','availSpecies.txt',
                                    package='UniProt.ws')
                        , header=FALSE, stringsAsFactors=FALSE)
  g <- grepl(pattern, species[,2])
  res <- species[g,]
  colnames(res) <- c("taxon ID","Species name")
  rownames(res)<- NULL
  head(res, n)
}




## and another method to look up the species name based on the tax ID.
lookupUniprotSpeciesFromTaxId <- function(taxId){
  species <- read.delim(system.file('extdata','availSpecies.txt',
                                    package='UniProt.ws')
                        , header=FALSE, stringsAsFactors=FALSE)
  g <- species[,1] %in% taxId
  res <- species[g,2]
  if(length(res)<1) stop("No species match the requested Tax Id.")
  if(length(res)>1) stop("There may be a problem with the Tax Id data file.")
  if(length(res)==1) return(res)
}












## important resources:
## query syntax
## http://www.uniprot.org/help/text-search

## how to use the REST site (generally):
## http://www.uniprot.org/faq/28

## query fields (what you can put for the "query=" part)
## http://www.uniprot.org/help/query-fields





## I need a function that can take a UniProt and retrieve associated data.
## something like below...

## these work:
## http://www.uniprot.org/uniprot/?query=P12345&format=tab&columns=id,sequence
## http://www.uniprot.org/uniprot/?query=accession:P30443&format=tab&columns=id,sequence


## This also works (gets ALL the IDs):
## http://www.uniprot.org/uniprot/?query=organism:9606&format=tab&columns=id,sequence



## And this is if you want multiple of the same kind
## http://www.uniprot.org/uniprot/?query="P04217"+or+"P30443"&format=tab&columns=id,sequence

## Notice though that I am getting extra entries.  Are they equivalant or cruft?
## The DEFINITELY are not perfectly equivalent...  
## I would like to understand why some requests get me multiple hits, while others only get me one hit????

## http://www.uniprot.org/uniprot/?query="P04217"&format=tab&columns=id,sequence

## So lets look at this example here: "P30443" returns data for it and also some data for another species with the ID "G7ZFG0"  </puzzled>
## http://www.uniprot.org/uniprot/?query="P30443"&format=tab&columns=id,sequence

## But if I just ask for "G7ZFG0" I don't also get the human version "P30443"...
## http://www.uniprot.org/uniprot/?query="G7ZFG0"&format=tab&columns=id,sequence


## What's more: using single quotes seems to change thigns AGAIN (makes it even LESS specific somehow?)

## So here is how it seems to work: Double quotes indicate one level of non-specificity (exactly what I am unsure), and single quotes make it even LESS specific!  So for my purposes, I need to avoid quotes.  One more BUMMER: not using quotes is the SAME as using double quotes (a small amount of non-specificity).  So I get to use a post filter regardless...  :(



## SO I want to format the query like this:
## http://www.uniprot.org/uniprot/?query=P04217+or+P30443&format=tab&columns=id,sequence


## I can also do things like this (thus extracting the stuff from other DBs)
## http://www.uniprot.org/uniprot/?query=P04217+or+P30443&format=tab&columns=id,database%28interpro%29




## potential values for columns (for use by getUniprotGoodies)  = c("citation","clusters","comments","domains","domain","ec","id","entry name","existence","families","features","genes","go","go-id","interpro","interactor","keywords","keyword-id","last-modified","length","organism","organism-id","pathway","protein names","reviewed","score","sequence","3d","subcellular locations","taxon","tools","version","virus hosts","database(pfam)","database(pdb)")

## NOTE: parameterized column values have to be "expanded" so database became
## database(pfam) and database(pdb)...








