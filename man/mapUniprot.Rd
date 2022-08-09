\name{mapping-and-querying}
\alias{mapUniProt}
\alias{queryUniProt}
\alias{allToKeys}
\alias{allFromKeys}
\alias{keytypes}
\alias{cols}
\title{
Mapping identifiers with the UniProt API
}
\description{
These functions are the main workhorses for mapping identifiers from one
database to another. They make use of the latest UniProt API (seen at
\url{https://www.uniprot.org/help/api}).
}
\usage{
mapUniProt(
    from = "UniProtKB_AC-ID",
    to = "UniRef90",
    columns = character(0L),
    query,
    verbose = FALSE,
    debug = FALSE
)
queryUniprot(
    query = character(0L),
    fields = c("accession", "id"),
    collapse = " OR "
)
allToKeys(fromName = "UniProtKB_AC-ID")
allFromKeys()
}
\arguments{
  \item{from}{
    character(1)
    The identifier type to map from, by default "UniProtKB_AC-ID", short for
    UniProt accession identifiers.
    See a list of all 'from' type identifiers with \code{allFromKeys}.
  }
  \item{to}{
    character(1)
    The target mapping identifier, by default "UniRef90". It can be any one of
    those returned by \code{allToKeys} from the appropriate \code{fromName}
    argument.
  }
  \item{columns}{
    character() Additional information to be retreived from UniProt service.
    See a full list of possible input return fields at
    \url{https://www.uniprot.org/help/return_fields}. Example fields include,
    "accession", "id", "gene_names", "xref_pdb", "xref_hgnc", "sequence", etc.
  }
  \item{query}{
    character() A string that would indicate the target query based on the
    available query fields. See \url{https://www.uniprot.org/help/query-fields}
    for a list of query fields. A typical query may only include UniProt
    accession identifiers, i.e., \code{c("A0A0C5B5G6", "A0A1B0GTW7", "A0JNW5",
    "A0JP26", "A0PK11", "A1A4S6")}
}
  \item{verbose}{
    logical(1) whether the operations should provide verbose updates (default
    \code{FALSE}).
  }
  \item{debug}{
    logical(1) whether to display the URL API endpoints, for advanced debugging
    (default \code{FALSE})
}
}
\details{
  Note that \code{mapUniProt} is used internally by the \code{select} method
  but made available for API queries with finer control.
}
\value{
\itemize{
  \item{mapUniProt}{A data.frame of returned results}
  \item{allToKeys}{A sorted character vector of possible "To" keytypes based
    on the given "From" type}
  \item{allFromKeys}{A sorted character vector of possible "From" keytypes}
}
}
\seealso{
  \code{\link{UniProtKB}}
}
\examples{
mapUniProt(
    from="UniProtKB_AC-ID",
    to='RefSeq_Protein',
    query=c('P13368','Q9UM73','P97793','Q17192')
)

mapUniProt(
    from='GeneID', to='UniProtKB', query=c('1','2','3','9','10')
)

mapUniProt(
    from = "UniProtKB_AC-ID", to = "UniProtKB",
    query = c("P31946", "P62258"),
    columns = c("accession", "id", "xref_pdb", "xref_hgnc", "sequence")
)

queryUniprot(
    query = c("accession:A5YMT3", "organism_id:9606"),
    fields = c("accession", "id", "reviewed"),
    collapse = " AND "
)

head(allToKeys())

allFromKeys(fromName = "UniRef100")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
