.check_rect_result <- function(res) {
    checkTrue(
        nrow(res) > 1 && ncol(res) >= 2L
    )
}

getUniProtGoodies <- UniProt.ws:::getUniProtGoodies
mapUniProt <- UniProt.ws:::mapUniProt

test_mapUniProt <- function(){
    res <- mapUniProt(
        from="UniProtKB_AC-ID", to='RefSeq_Protein',
        query=c('P13368','P20806','Q9UM73','P97793','Q17192')
    )
    .check_rect_result(res)
    checkIdentical(res[1,"From"], 'P13368')
    checkIdentical(res[1,"To"], 'NP_511114.2')

    ## what if I have entrezGene IDs and I want UniProts?
    res <- mapUniProt(
        from='GeneID', to='UniProtKB', query=c('1','2','3','9','10')
    )
    .check_rect_result(res)
    checkIdentical(res[1,"From"], 1L)
    checkIdentical(res[1,"Entry"], 'P04217')

    ## I can then map UniProt accessions to Unigene IDs
    res <- mapUniProt(
        from='UniProtKB_AC-ID',
        to='GeneID',
        query=c('P04217','P01023','F5H5R8','P18440','Q400J6')
    )
    .check_rect_result(res)
    checkIdentical(res[1,"To"], 1L)
    checkIdentical(res[1,"From"], 'P04217')

    ## from = "Gene_Name" can restrict by taxId in query
    res <- mapUniProt(
        from='Gene_Name',
        to = "UniProtKB-Swiss-Prot",
        query=list(ids = "TP53", taxId = 9606),
        columns = c("accession", "id", "organism_id")
    )
    checkTrue(is(res, "data.frame"))
    checkIdentical(nrow(res), 1L)
    checkIdentical(ncol(res), 4L)
    checkTrue(
        all(names(res) %in%  c("From", "Entry", "Entry.Name", "Organism..ID."))
    )

    res <- mapUniProt(
        from='Gene_Name',
        to = "UniProtKB",
        query=list(ids = "TP53", taxId = 9606),
        columns = c("accession", "id", "organism_id")
    )
    checkTrue(is(res, "data.frame"))
    checkIdentical(ncol(res), 4L)
    checkTrue(
        all(names(res) %in%  c("From", "Entry", "Entry.Name", "Organism..ID."))
    )
}

test_getUniProtGoodies <- function(){
    query <- c('P04217','P30443')
    cols <- 'sequence'
    res <- getUniProtGoodies(query, cols)
    checkTrue(is(res, "data.frame"))
    checkTrue(nrow(res) >= 2L)
    checkIdentical(ncol(res), 3L)

    ## can also be used to extract interpro IDs
    query <- c('P13368','P20806','Q9UM73','P97793','Q17192')
    cols <- 'xref_interpro'
    res <- getUniProtGoodies(query, cols)
    checkTrue(is(res, "data.frame"))
    checkTrue(nrow(res) >= 5L)
    checkIdentical(ncol(res), 3L)

    ## OR extract a number of other things... ## taxon (?)
    query <- c('P13368','P20806','Q9UM73','P97793','Q17192')
    cols <- c('structure_3d','go_id')
    res <- getUniProtGoodies(query, cols)
    checkTrue(is(res, "data.frame"))
    checkTrue(nrow(res) >= 5L)
    checkIdentical(ncol(res), 4L)
}
