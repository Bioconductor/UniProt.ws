up <- UniProt.ws(taxId=9606)

test_replaceTaxIdMethod <- function(){
  up <- UniProt.ws(taxId=9606)
  taxId(up) <- 10090
  checkIdentical(taxId(up), 10090)
}

test_availableUniprotSpecies <-function(){
  res <- availableUniprotSpecies(pattern="Homo")
  checkTrue(nrow(res) > 1)
  checkIdentical(ncol(res), 3L)
}

test_lookupUniprotSpeciesFromTaxId<- function(){
  res <- lookupUniprotSpeciesFromTaxId(9606)
  checkIdentical(res, "Homo sapiens")
}

test_species <- function(){
  up <- UniProt.ws(taxId=9606)
  res <- species(up)
  checkIdentical(res, "Homo sapiens (Human)")
}

test_keytypes <- function(){
  res <- keytypes(up)
  checkTrue(length(res) >1)
}

test_select_1 <- function(){
    ## 1st working select example...
    keys <- c("P31946","P62258","Q04917")
    kt <- "UniProtKB"
    cols <- c("xref_pdb","xref_hgnc","sequence")
    res <- select(x = up, keys = keys, columns = cols, keytype = kt)
    checkTrue(is.data.frame(res))
    checkTrue(nrow(res) > 0L)
    checkIdentical(ncol(res), 5L)
    checkIdentical(
        c("From", "Entry", "PDB", "HGNC", "Sequence"), colnames(res)
    )
}

test_select_2 <- function(){
  ## with an alternate keytype (need to think carefully about merge keys here)
  keys <- c('1','2','3','9','10')
  kt <- "GeneID"
  cols <- c("xref_pdb", "xref_hgnc", "sequence")
  res <- select(x = up, keys = keys, columns = cols, keytype = kt)
  checkTrue(nrow(res) > 0)
  checkIdentical(ncol(res), 5L)
  checkIdentical(
      c("From", "Entry", "PDB", "HGNC", "Sequence"), colnames(res)
  )
}

test_select_3 <- function(){
  keys <- c("P31946","P62258","Q04917")
  kt <- "UniProtKB"
  cols <- "virus_hosts" ## this is not allowed
  res <- select(x = up, keys = keys, columns = cols, keytype = kt)
  checkTrue(is.data.frame(res))
  checkTrue(nrow(res) > 0L)
  checkIdentical(ncol(res), 3L)
  checkIdentical(
      c("From", "Entry", "Virus.hosts"), colnames(res)
  )
}

test_select_4 <- function(){
  keys <- c("P31946","P62258","Q04917")
  kt <- "UniProtKB"
  cols <- "xref_pdb"
  res <- select(x = up, keys = keys, columns = cols, keytype = kt)
  checkTrue(nrow(res)>0)
  checkIdentical(ncol(res), 3L)
  checkIdentical(
    c("From", "Entry", "PDB"), colnames(res)
  )
}

test_select_5 <- function(){
  keys <- c('1','2','3','9','10')
  kt <- "GeneID"
  ## CLUSTERS
  cols <- c("xref_pdb", "CLUSTERS")
  res <- select(x = up, keys = keys, columns = cols, keytype = kt)
  checkTrue(nrow(res) > 0)
  checkTrue(ncol(res) > 4)
  checkTrue(all(c("Entry", "From", "PDB", "Cluster.ID") %in% colnames(res)))
}

test_select_6 <- function(){
  ## now lets just get a bunch of the sequences.
  keys <- keys(up,keytype="UniProtKB")
  kt <- "UniProtKB"
  cols <- "sequence"
  res <- select(x = up, keys = keys, columns = cols, keytype = kt)
  checkTrue(nrow(res) > 0)
  checkIdentical(ncol(res), 3L)
  checkIdentical(
    c("From", "Entry" ,"Sequence"), colnames(res)
  )
}

test_select_7 <- function(){
  ## test that we fail when the pass in bad keytype
  kt = "UNIPROT"
  checkException(select(up, keys, cols, kt))
}

test_select_8 <- function(){
  ## GENEID was problematic
  ## create test cases to fix
  keys <- c("P31946","P62258","Q04917")
  kt <- "UniProtKB"
  cols <- "xref_geneid"
  res <- select(x = up, keys = keys, columns = cols, keytype = kt)
  checkIdentical(
    c("From","Entry", "GeneID"), colnames(res)
  )

  cols <- c("xref_geneid", "sequence")
  res <- select(x = up, keys = keys, columns = cols, keytype = kt)
  checkIdentical(
    c("From","Entry", "GeneID", "Sequence"), colnames(res)
  )

  cols <- c("GeneID", "sequence")
  checkException(select(x = up, keys = keys, columns = cols, keytype = kt))

  cols <- c("ENTREZ_GENE", "sequence")
  checkException(select(x = up, keys = keys, columns = cols, keytype = kt))

}

## keys with ecs:
## keys = c("Q06278","Q9BRR6","Q86V24")
## What I learned from ecs is that if there are TRULY no hits, then you don't get the column you requested AT ALL, so code had to be added to handle this case


## keys is extremely slow, so I will comment this test to speed up the tests.
## but it might be useful to have this in the future.
## test_keys <- function(){

##   ## keys is VERY slow :(
##   ##   keytype = "UniProtKB"
##   ##   k = keys(up, keytype)

##   ## so switch to a critter with fewer egs?
##   taxId(up) <- 9913

##   keytype = "ENTREZ_GENE"
##   egs = keys(up, keytype)
##   checkTrue(any("282126" %in% egs))
##   checkTrue(is.character(egs))
##   checkTrue(length(egs)>1)
## }




## problem getting headers for: tools, keyword-id, clusters, score




### Problems:

## The following does not work?:

## refseq <- c("YP_139402", "YP_141320", "YP_820357", "YP_006002448", "YP_006040776", "YP_006340074", "YP_006340075", "YP_138830", "YP_819836", "YP_006001858")

## select(up, keys=refseq, keytype="REFSEQ_PROTEIN", columns=c('UniProtKB'))

## ANSWER: Bug?  - Not actually (in this case).  It's just that NONE of the refseq protein IDs (at least for human) are in the list this person posted...

## These IDs appear to be from taxid = "1308" (Streptococcus thermophilus)
## taxId(up) = "1308"
## BUT even after setting to use that taxid, there doesn't seem to be any IDs inthe UniProt web service that match it.
## rs = keys(up, "REFSEQ_PROTEIN")
## any(refseq %in% rs)


### BUT: I am unsatisfied with the idea that this is not a bug.  The initial step to map these refseq IDs to a uniprot ID seems to proceed without trouble, and then there seems to be no data coming back from .mapUni().  So why is that happening????  I think it's because the .mapUni() method should not be called at all IF the final end point is a uniprot and the start was not one either...
## So this should work OK (for example):

## select(up, keys=refseq, keytype="REFSEQ_PROTEIN", columns=c('ENTREZ_GENE'))

## And it does...  ;)

## That means that I just need to fix the fact that when you ask for Uniprots, .select needs to stop sooner.

## STATUS: FIXED






## And then there is this one:

## Hi Marc!

## I am using your Uniprot.ws package and really appreciate it!! However, I have some questions about mapping back and forth from/to Uniprot itself.
## Below are my questions with R commands, hopefully self explanatory...
## Thanks a lot!
## Nik

## # why does this give me many organisms, after i selected only human?
## taxId(up) <- 9606 # set to human
## select(up, keys=c("I39.001"), columns=c("ID", "ORGANISM"), "MEROPS")

## ANSWER: Bug?


## # I don't manage to map FROM uniprot ACs to other databases, i thought it would work like this:
## select(up, keys=c("P01023"), columns=c("MEROPS"), "UniProtKB")

## ANSWER: Bug? - Yes (same one as above - and it's fixed now)



## # from the returned Uniprot ACs how do I distinguish between reviewed SwissProt and unreviewed TrEMBL ACs?
## # How can i retrieve Uniprot IDs (like "A2MG_HUMAN") instead of ACs?
## select(up, keys=c(216, 3679, 55607), columns=c("ID", "ORGANISM"), "ENTREZ_GENE")

## ANSWER: I don't think you can.  The web service does not seem to want to make a distinction.




## Faster testing:
## BiocGenerics:::testPackage(pattern="^test_select.*\\.R$")

