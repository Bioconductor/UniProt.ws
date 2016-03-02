## Here we will define the select methods.

## no way to 'discover' the keytypes either (and they are subset of cols).
.keytypes <- function(){
  keytypeKeysDat[,1]
}

setMethod("keytypes", "UniProt.ws", function(x){.keytypes()})


## no way to 'discover' the cols, so I hard code them here.
.cols <- function(){
  c(keytypeKeysDat[,1], extraColsDat[,1])
}

setMethod("columns", "UniProt.ws", function(x){.cols()})

## http://www.UniProt.org/UniProt/?query=organism:9606&format=tab&columns=id,sequence

## To make keys work I just want to return what was asked for...
.keys <- function(x, keytype){
  if(!any(keytypes(x) %in% keytype)){
    stop("keytype argument MUST match a value returned by keytypes method")
  }
  dat <- taxIdUniprots(x) ## pre-cached
  if(keytype == "UNIPROTKB"){
    return(dat)
  }else{
    ## then convert this to be the keytype requested...
    tkt = keytypeKeysDat[keytypeKeysDat[,1] %in% keytype,2]
    dat2 <- mapUniprot(from="ACC+ID", to=tkt, query=dat)
    return(unique(as.character(dat2[,2])))
  }
}

setMethod("keys", "UniProt.ws",
    function(x, keytype){
      if(missing(keytype)){stop("Please supply a keytype argument.")}
      .keys(x, keytype)
    }
)

.mergeList <- function(list, joinType="left"){
  for(i in seq_len(length(list))){
    if(i==1){
      fin <- list[[1]]
    }else{
      if(joinType=="left"){
        fin <- merge(fin, list[[i]], by="ACC+ID", all.x=TRUE) ## left outer join
      }
      else if(joinType=="all"){
        fin <- merge(fin, list[[i]], by="ACC+ID", all=TRUE) ## full outer join
      }
    }
  }
  fin
}

.getUPMappdata <- function(colMappers, keys){
  ## get a list of mapping results (as data.frames)
  res <- lapply(colMappers, FUN=mapUniprot, from="ACC+ID", query=keys)
  ## Them merge all these mappings together based on UniProt.
  .mergeList(res, joinType="all")
}

## Here is the business end of my select method.
## The big plan is to call mapUniprot() and getUniprotGoodies()
## (merging when necessary)
.select <- function(x, keys, cols, keytype){
  if(!any(keytypes(x) %in% keytype)){
    stop("keytype argument MUST match a value returned by keytypes method")
  }
  if(!any(columns(x) %in% cols)){
    stop("columns argument MUST match a value returned by columns method")
  }

  ## process columns
  oriTabCols <- unique(c(keytype,cols))
  cols <- cols[!(cols %in% keytype)]  ## remove keytype from cols 
  if (!length(cols))
      stop("'columns' should be different from 'keytype'")
  trueKeys <- keys ## may change depending on keytype.
  ## split into 2 groups: cols in keytypeKeys and cols in extraCols 
  colMappers <- cols[cols %in% keytypeKeysDat[,1]]
  colUPGoodies <- cols[cols %in% extraColsDat[,1]]
  ## then convert those into the internally used IDs
  colMappers <- keytypeKeysDat[keytypeKeysDat[,1] %in% colMappers, 2]
  colUPGoodies <- extraColsDat[extraColsDat[,1] %in% colUPGoodies, 2]
  res <- list()
  if(keytype!="UNIPROTKB" ){
    kt <- keytypeKeysDat[keytypeKeysDat[,1] %in% keytype,2]
    dat <- mapUniprot(from=kt, to="ACC", query=keys)
    colnames(dat)[2] <-  "ACC+ID" ## always the 2nd one...
    ## capture UniProts as keys from this point on
    keys <- unique(dat[["ACC+ID"]])
    res <- c(res, list(dat))
  }

  ## All the (UNIPROTKB) possible keys for this organism
  orgSpecificKeys <- keys(x, keytype="UNIPROTKB")
  ## Now filter keys with orgSpecificKeys (uniprots intersected with uniprots)
  keys <- intersect(keys, orgSpecificKeys)
  if(length(keys)==0) stop("No data is available for the keys provided.")

  ## now get the other data (depending what was asked for)
  if (length(colMappers) > 0 && colMappers!="ACC+ID")
    res <- c(res, list(.getUPMappdata(colMappers, keys)))
  if(length(colUPGoodies) > 0){
    dat <- getUniprotGoodies(keys, colUPGoodies)
    colnames(dat)[1] <- "ACC+ID" ## always the 1st
    res <- c(res, list(dat))
  }
  ## At this point I have some results, Now I just need to merge them base on
  ## UniProt IDs (and upon whether or not they are real)
  tab <- .mergeList(res, joinType="all")
  ## rename cols:
  rosetta <- rbind(keytypeKeysDat, extraColsDat)
  ## We need the third col of rosetta to tell us what the cols will come back
  ## from the service as
  idx <- match(colnames(tab), rosetta[,3])
  colnames(tab) <- rosetta[idx,1]
  ## unique to this web service is the fact that I sometimes will have an
  ## extra UNIPROTKB col.  Regardless, we ONLY want the cols we asked for..
  ## BUT: we also can't try this if the above code has failed to rename anything
  if (all(!is.na(colnames(tab))))
    tab <- tab[,colnames(tab) %in% oriTabCols]
  ## resort
  tab <- AnnotationDbi:::.resort(tab, trueKeys, keytype, oriTabCols)
  ## Now one last cast to make NAs (and all cols) and make things "uniform"
  cnames <- colnames(tab)
  .blankToNA <- function(row){
      gsub(pattern="^$",replacement=as.character(NA),row)}
  tab <- data.frame( t(apply(tab,MARGIN=1,.blankToNA)), stringsAsFactors=FALSE)
  colnames(tab) <- cnames
  ## then return
  tab
}

setMethod("select", "UniProt.ws",
    function(x, keys, columns, keytype, ...){
          if (missing(keytype)) 
              keytype <- "UNIPROTKB"
          .select(x, keys, columns, keytype)
        }
)
