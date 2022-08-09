\name{UniProtKB}

\title{Descriptions of available values for \code{columns} and \code{keytypes}.}

\description{
  The types of data represented and transmitted by the UniProt project grows
  constantly. Please refer to the \url{https://www.uniprot.org/help/return_fields}
  for a list of all possible return fields. These translate to \code{cols}
  inputs in functions. Keytypes can be discovered using \code{allToKeys}.
}

\author{Marc Carlson}

\examples{
    up <- UniProt.ws(taxId=9606)
    ## List the possible values for keytypes
    keytypes(up)
    ## default keys - UniProtKB
    keys(up, keytype="UniProtKB")
    select(
        x = up, keys = c("P31946","P62258"),
        columns = c("xref_pdb","sequence"), keytype = "UniProtKB"
    )
}
