test_mapUniProt <- function() {
    keys <- list(organism_id = 10090, ids = c('Q7TPG8', 'P63318'))
    res <- mapUniProt(
        from = "UniProtKB_AC-ID",
        to = "UniProtKB",
        columns = c("accession", "id"),
        query = keys
    )
    checkTrue(is.data.frame(res))
    checkIdentical(nrow(res), 2L)
    checkIdentical(colnames(res), c("From", "Entry", "Entry.Name"))

    keys <- c('Q7TPG8', 'P63318')
    res <- mapUniProt(
        from = "UniProtKB_AC-ID",
        to = "UniProtKB",
        columns = c("accession", "id"),
        query = keys
    )
    checkTrue(is.data.frame(res))
    checkIdentical(nrow(res), 2L)
    checkIdentical(colnames(res), c("From", "Entry", "Entry.Name"))
}

test_returnFields <- function() {
    rf <- returnFields()
    checkIdentical(names(rf), c("groupName", "label", "name"))
    checkTrue(is.data.frame(rf))
}
