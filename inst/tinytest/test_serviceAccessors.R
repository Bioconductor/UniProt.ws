.check_rect_result <- function(res) {
    expect_true(
        nrow(res) > 1 && ncol(res) >= 2L
    )
}

res <- mapUniProt(
    from="UniProtKB_AC-ID", to='RefSeq_Protein',
    query=c('P13368','P20806','Q9UM73','P97793','Q17192')
)
.check_rect_result(res)
expect_identical(res[1,"From"], 'P13368')
expect_identical(res[1,"To"], 'NP_511114.2')

## what if I have entrezGene IDs and I want UniProts?
res <- mapUniProt(
    from='GeneID', to='UniProtKB', query=c('1','2','3','9','10')
)
.check_rect_result(res)
expect_identical(res[1,"From"], 1L)
expect_identical(res[1,"Entry"], 'P04217')

## I can then map UniProt accessions to Unigene IDs
res <- mapUniProt(
    from='UniProtKB_AC-ID',
    to='GeneID',
    query=c('P04217','P01023','F5H5R8','P18440','Q400J6')
)
.check_rect_result(res)
expect_identical(res[1,"To"], 1L)
expect_identical(res[1,"From"], 'P04217')

## from = "Gene_Name" can restrict by taxId in query
keys <- list(ids = "TP53", taxId = 9606)
res <- mapUniProt(
    from='Gene_Name',
    to = "UniProtKB-Swiss-Prot",
    query = keys,
    columns = c("accession", "id", "organism_id")
)
expect_true(is(res, "data.frame"))
expect_true(nrow(res) >= 1L)
expect_identical(ncol(res), 4L)
expect_true(
    all(names(res) %in%  c("From", "Entry", "Entry.Name", "Organism..ID."))
)

res <- mapUniProt(
    from='Gene_Name',
    to = "UniProtKB",
    query = keys,
    columns = c("accession", "id", "organism_id")
)
expect_true(is(res, "data.frame"))
expect_identical(ncol(res), 4L)
expect_true(
    all(names(res) %in%  c("From", "Entry", "Entry.Name", "Organism..ID."))
)
