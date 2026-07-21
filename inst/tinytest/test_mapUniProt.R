keys <- list(organism_id = 10090, ids = c('Q7TPG8', 'P63318'))
res <- mapUniProt(
    from = "UniProtKB_AC-ID",
    to = "UniProtKB",
    columns = c("accession", "id"),
    query = keys
)
expect_true(is.data.frame(res))
expect_identical(nrow(res), 2L)
expect_identical(colnames(res), c("From", "Entry", "Entry.Name"))

keys <- c('Q7TPG8', 'P63318')
res <- mapUniProt(
    from = "UniProtKB_AC-ID",
    to = "UniProtKB",
    columns = c("accession", "id"),
    query = keys
)
expect_true(is.data.frame(res))
expect_identical(nrow(res), 2L)
expect_identical(colnames(res), c("From", "Entry", "Entry.Name"))

rf <- returnFields()
expect_identical(names(rf), c("groupName", "label", "name"))
expect_true(is.data.frame(rf))

expect_error(
    mapUniProt(
        from = "Gene_Name",
        to = "UniProtKB-Swiss-Prot",
        columns = c("accession", "id"),
        query = list(taxId = 9606)
    ),
    "'ids' must be a non-zero character vector"
)

expect_error(
    mapUniProt(
        from = "Gene_Name",
        to = "UniProtKB-Swiss-Prot",
        columns = c("accession", "id"),
        query = list(taxId = 9606, ids = "")
    ),
    "'ids' must be a non-zero character vector"
)

keys <- list(taxId = 9606, ids = "TP53")
## filters out nzchars
res <- mapUniProt(
    from = "Gene_Name",
    to = "UniProtKB-Swiss-Prot",
    columns = c("accession", "id"),
    query = keys
)
expect_true(is.data.frame(res))
expect_identical(nrow(res), 1L)
expect_identical(colnames(res), c("From", "Entry", "Entry.Name"))
