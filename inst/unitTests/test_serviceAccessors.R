##  require(RUnit)
##  require(UniProt.ws)

.check_rect_result <- function(res){
  all(checkTrue(dim(res)[1] >1),
      checkTrue(dim(res)[2] ==2))
}

## problem?
test_mapUniprot <- function(){
    mapUniprot <- UniProt.ws:::mapUniprot
    res <- mapUniprot(
        from='ACC',to='P_REFSEQ_AC',
        query=c('P13368','P20806','Q9UM73','P97793','Q17192')
    )
    .check_rect_result(res)
    checkTrue(res[1,1]=='P13368')
    checkTrue(res[1,2]=='NP_511114.2')

    ## what if I have entrezGene IDs and I want UniProts?
    res <- mapUniprot(
        from='P_ENTREZGENEID', to='ACC', query=c('1','2','3','9','10')
    )
    .check_rect_result(res)
    checkTrue(res[1,1]=='1')
    checkTrue(res[1,2]=='P04217')

    ## I can then map UniProt accessions to Unigene IDs
    res <- mapUniprot(
        from='ACC',to='UNIGENE_ID',
        query=c('P04217','P01023','F5H5R8','P18440','Q400J6')
    )
    .check_rect_result(res)
}


test_getUniprotGoodies <- function(){
  query = c('P04217','P30443')
  cols = 'sequence'
  res <- UniProt.ws:::getUniprotGoodies(query, cols)
  checkTrue(class(res) == "data.frame")
  checkTrue(dim(res)[1] == 2)
  checkTrue(dim(res)[2] == 2)

  ## can also be used to extract interpro IDs
  query = c('P13368','P20806','Q9UM73','P97793','Q17192')
  cols = 'database(interpro)'
  res <- UniProt.ws:::getUniprotGoodies(query, cols)
  checkTrue(class(res) == "data.frame")
  checkTrue(dim(res)[1] == 5)
  checkTrue(dim(res)[2] == 2)

  ## OR extract a number of other things...
  cols = c('3d','go-id','taxon')
  res <- UniProt.ws:::getUniprotGoodies(query, cols)
  checkTrue(class(res) == "data.frame")
  checkTrue(dim(res)[1] == 5)
  checkTrue(dim(res)[2] == 4)
   
}


## Faster testing:
## BiocGenerics:::testPackage(pattern="^test_serviceAccessors.*\\.R$")

