.check_rect_result <- function(res) {
    checkTrue(
        nrow(res) > 1 && ncol(res) >= 2L
    )
}

test_mapUniprot <- function(){
    mapUniprot <- UniProt.ws:::mapUniprot
    res <- mapUniprot(
        from="UniProtKB_AC-ID", to='RefSeq_Protein',
        query=c('P13368','P20806','Q9UM73','P97793','Q17192')
    )
    .check_rect_result(res)
    checkIdentical(res[1,"From"], 'P13368')
    checkIdentical(res[1,"To"], 'NP_511114.2')

    ## what if I have entrezGene IDs and I want UniProts?
    res <- mapUniprot(
        from='GeneID', to='UniProtKB', query=c('1','2','3','9','10')
    )
    .check_rect_result(res)
    checkIdentical(res[1,"From"], 1L)
    checkIdentical(res[1,"Entry"], 'P04217')

    ## I can then map UniProt accessions to Unigene IDs
    res <- mapUniprot(
        from='UniProtKB_AC-ID',
        to='GeneID',
        query=c('P04217','P01023','F5H5R8','P18440','Q400J6')
    )
    .check_rect_result(res)
    checkIdentical(res[1,"To"], 1L)
    checkIdentical(res[1,"From"], 'P04217')
}

test_getUniprotGoodies <- function(){
    query <- c('P04217','P30443')
    cols <- 'sequence'
    res <- UniProt.ws:::getUniprotGoodies(query, cols)
    checkTrue(is(res, "data.frame"))
    checkIdentical(nrow(res), 2L)
    checkIdentical(ncol(res), 3L)

    ## can also be used to extract interpro IDs
    query <- c('P13368','P20806','Q9UM73','P97793','Q17192')
    cols <- 'xref_interpro'
    res <- UniProt.ws:::getUniprotGoodies(query, cols)
    checkTrue(is(res, "data.frame"))
    checkIdentical(nrow(res), 5L)
    checkIdentical(ncol(res), 3L)

    ## OR extract a number of other things... ## taxon (?)
    query <- c('P13368','P20806','Q9UM73','P97793','Q17192')
    cols <- c('structure_3d','go_id')
    res <- UniProt.ws:::getUniprotGoodies(query, cols)
    checkTrue(is(res, "data.frame"))
    checkIdentical(nrow(res), 5L)
    checkIdentical(ncol(res), 4L)
}
