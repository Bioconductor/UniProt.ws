##  require(RUnit)
##  require(UniProt.ws)



test_newTaxId <- function(){
  taxId(UniProt.ws)
  taxId(UniProt.ws) <- 10090
  checkTrue(taxId(UniProt.ws) == 10090)
}


test_availableUniprotSpecies <-function(){
  res <- availableUniprotSpecies(pattern="Homo")
  checkTrue(dim(res)[1]>1)
  checkTrue(dim(res)[2]==2)
}


test_lookupUniprotSpeciesFromTaxId<- function(){
  res <- lookupUniprotSpeciesFromTaxId(9606)
  checkTrue(res=="Homo sapiens")
}


test_species <- function(){
  res <- species(UniProt.ws)
  checkIdentical(res, "Homo sapiens")
}


test_cols_and_keytypes <- function(){
 res <- keytypes(UniProt.ws)
 checkTrue(length(res) >1)
 res2 <- cols(UniProt.ws)
 checkTrue(length(res2) >1)
 checkTrue(length(res) != length(res2))
 checkTrue(all(res %in% res2))
}




test_select_1 <- function(){
  ## 1st working select example...  
  keys <- c("P31946","P62258","Q04917")
  kt <- "UNIPROTKB"
  cols <- c("PDB","UNIGENE","SEQUENCE")
  res <- select(UniProt.ws, keys, cols, kt)
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==4)
  checkIdentical(c("UNIPROTKB","PDB","UNIGENE","SEQUENCE"), colnames(res))
}

test_select_2 <- function(){
  ## with an alternate keytype (need to think carefully about merge keys here)
  keys <- c('1','2','3','9','10')
  kt <- "ENTREZ_GENE"
  cols <- c("PDB","UNIGENE","SEQUENCE")
  res <- select(UniProt.ws, keys, cols, kt)
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==4)
  checkIdentical(c("ENTREZ_GENE" ,"PDB","UNIGENE","SEQUENCE"),
                 colnames(res))
}

test_select_3 <- function(){
  ## Unit tests I have to make pass:
  keys <- c("P31946","P62258","Q04917")
  kt <- "UNIPROTKB"
  cols <- c("VIRUS_HOSTS") ## this is not allowed
  checkException(res <- select(UniProt.ws, keys, cols, kt))
}

test_select_4 <- function(){
  keys <- c("P31946","P62258","Q04917")
  kt <- "UNIPROTKB"
  cols <- c("DATABASE(PDB)")
  res <- select(UniProt.ws, keys, cols, kt)
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==2)
  checkIdentical(c("UNIPROTKB" ,"DATABASE(PDB)"),
                 colnames(res))
}

test_select_5 <- function(){
  keys <- c('1','2','3','9','10')
  kt <- "ENTREZ_GENE"
  cols <- c("PDB","CLUSTERS")
  res <- select(UniProt.ws, keys, cols, kt)
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==3)
  checkIdentical(c("ENTREZ_GENE","PDB","CLUSTERS"),
                 colnames(res))
}

## system.time(res <- select(UniProt.ws, keys, cols="EC", kt))

test_select_6 <- function(){
  ## now lets just get a bunch of the sequences.
  keys <- keys(UniProt.ws,keytype="UNIPROTKB")
  keys <- head(keys, n=1000)
  kt <- "UNIPROTKB"
  cols <- c("SEQUENCE")
  ## cols <- cols(UniProt.ws)[96:126]
  ## debug(UniProt.ws:::.getSomeUniprotGoodies)
  res <- select(UniProt.ws, keys, cols, kt)
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==2)
  checkIdentical(c("UNIPROTKB" ,"SEQUENCE"),
                 colnames(res))
}


test_select_7 <- function(){
  ## test that we fail when the pass in bad keytype
  kt = "UNIPROT"
  checkException(select(UniProt.ws, keys, cols, kt))
}


## keys with ecs:
## keys = c("Q06278","Q9BRR6","Q86V24")
## What I learned from ecs is that if there are TRULY no hits, then you don't get the column you requested AT ALL, so code had to be added to handle this case


## keys is extremely slow, so I will comment this test to speed up the tests.
## but it might be useful to have this in the future.
test_keys <- function(){

  ## keys is VERY slow :(
  ##   keytype = "UNIPROTKB"
  ##   k = keys(UniProt.ws, keytype)
  
  ## so switch to a critter with fewer egs?
  taxId(UniProt.ws) <- 9913
  
  keytype = "ENTREZ_GENE"
  egs = keys(UniProt.ws, keytype)
  checkTrue(any("282126" %in% egs))
  checkTrue(is.character(egs))
  checkTrue(length(egs)>1)
}




## problem getting headers for: tools, keyword-id, clusters, score
