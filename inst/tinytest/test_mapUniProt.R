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
