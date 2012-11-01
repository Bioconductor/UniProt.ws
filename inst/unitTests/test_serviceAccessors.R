##  require(RUnit)
##  require(UniProt.ws)

.check_rect_result <- function(res){
  all(checkTrue(dim(res)[1] >1),
      checkTrue(dim(res)[2] ==2))
}


test_mapUniprot <- function(){
  res <- UniProt.ws:::mapUniprot(from='ACC',to='P_REFSEQ_AC',
                    query=c('P13368','P20806','Q9UM73','P97793','Q17192'))
  .check_rect_result(res)
  checkTrue(res[1,1]=='P13368')
  checkTrue(res[1,2]=='NP_511114.2')

  ## what if I have entrezGene IDs and I want UniProts?
  res <- UniProt.ws:::mapUniprot(from='P_ENTREZGENEID',to='ACC',
                                query=c('1','2','3','9','10'))
  .check_rect_result(res)
  checkTrue(res[1,1]=='1')
  checkTrue(res[1,2]=='P04217')
  
  ## I can then map UniProt accessions to IPI IDs
  res <- UniProt.ws:::mapUniprot(from='ACC',to='P_IPI',
                    query=c('P04217','P01023','F5H5R8','P18440','Q400J6'))
  .check_rect_result(res)
}


## Test to be commented (but not removed).
## This is a useful test, but it takes a very long time to run. And since the
## function is only needed twice a year...
test_getOneToMany <- function(){
  taxId <- "9606" ## human taxonomy Id
  
  ## Now I can get pfams like this:
  pfs <- UniProt.ws:::getOneToMany(taxId, "PFAM")
  .check_rect_result(pfs)
  checkTrue(dim(pfs[pfs[[1]]=="P30443",])[1] > 1)

  ## And I can get prosite IDs like this:
  ps <- UniProt.ws:::getOneToMany(taxId, "prosite")
  .check_rect_result(ps)
  checkTrue(dim(ps[ps[[1]]=="P62258",])[1] > 1)
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



