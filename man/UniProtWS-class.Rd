\name{UniProt.ws-objects}
\docType{class}

% Classes
\alias{UniProt.ws}
\alias{class:UniProt.ws}
\alias{UniProt.ws-class}

% Accessor-like methods
\alias{show,UniProt.ws-method}
\alias{species}
\alias{species,UniProt.ws-method}
\alias{taxId}
\alias{taxId,UniProt.ws-method}
\alias{taxId<-}
\alias{taxId<-,UniProt.ws-method}


% select and select related methods
\alias{columns}
\alias{columns,UniProt.ws-method}

\alias{keytypes}
\alias{keytypes,UniProt.ws-method}

\alias{keys}
\alias{keys,UniProt.ws-method}

\alias{select}
\alias{select,UniProt.ws-method}

% supporting functions
\alias{availableUniprotSpecies}
\alias{lookupUniprotSpeciesFromTaxId}


\title{UniProt.ws objects and their related methods and functions}

\description{
  \code{UniProt.ws} is the base class for interacting with the Uniprot
  web services from Bioconductor.

  In much the same way as an \code{AnnotationDb} object allows acces to
  select for many other annotation packages, \code{UniProt.ws} is meant
  to allow usage of \code{select} methods and other supporting methods
  to enable the easy extraction of data from the Uniprot web services.

  \code{select}, \code{columns} and \code{keys} are used together to
  extract data via an \code{UniProt.ws} object.

  \code{columns} shows which kinds of data can be returned for the
  \code{UniProt.ws} object.

  \code{keytypes} allows the user to discover which keytypes can be
  passed in to \code{select} or \code{keys} via the \code{keytype}
  argument.

  \code{keys} returns keys for the database contained in the
  \code{UniProt.ws} object .  By default it will return the primary keys
  for the database, which are UniProtKB keys, but if used with the
  \code{keytype} argument, it will return the keys from that keytype.

  \code{select} will retrieve the data as a data.frame based on
  parameters for selected \code{keys} and \code{columns} and \code{keytype}
  arguments.

  The \code{UniProt.ws} will be loaded whenever you load the
  \code{UniProt.ws} package.  This object will be set up to retrieve
  information from Homo sapiens by default, but this value can be
  changed to any of the species supported by Uniprot.  The
  \code{species} and \code{taxId} methods allow users to see what
  species is currently being accessed, and \code{taxId<-} allows them to
  change this value.

  \code{species} shows the genus and species label currently attached to
  the \code{UniProt.ws} objects database.

  \code{taxId} shows the NCBI taxonomy ID currently attached to
  the \code{AnnotationDb} objects database.  Using the equivalently
  names replace method (\code{taxId<-}) allows the user to change the
  taxon ID, and the species represented along with it.

  \code{availableUniprotSpecies} is a helper function to list out the
  available Species along with their official taxonomy IDs that are
  available by Uniprot.  Because there are so many species represented
  at UniProt, there is also a pattern argument that can be used to
  restrict the range of things returned to be only those whose species
  names match the searth term.  Please remember when using this argument
  that the Genus is always capitalized and the species never is.

  \code{lookupUniprotSpeciesFromTaxId} is another helper that will look
  up the species of any tax ID that is supported by Uniprot.

}

\usage{
  columns(x)
  keytypes(x)
  select(x, keys, columns, keytype, ...)
  species(object)
  taxId(x)

  availableUniprotSpecies(pattern)
  lookupUniprotSpeciesFromTaxId(taxId)
  UniProt.ws(taxId, ...)
}

\arguments{
  \item{x}{the \code{UniProt.ws} object. }
  \item{object}{the \code{UniProt.ws} object. }
  \item{keys}{the keys to select records for from the database.  All possible
    keys are returned by using the \code{keys} method.}
  \item{columns}{the columns or kinds of things that can be retrieved
    from the database.  As with \code{keys}, all possible columns are
    returned by using the \code{columns} method.}
  \item{keytype}{the keytype that matches the keys used.  For the
    \code{select} methods, this is used to indicate the kind of ID being used
    with the keys argument. For the \code{keys} method this is used to
    indicate which kind of keys are desired from \code{keys}
  }
  \item{pattern}{A string passed in to limit the results}
  \item{taxId}{a taxonomy id}
  \item{...}{other arguments}
}

\value{
  \code{keys},\code{columns},\code{keytypes}, \code{species} and
  \code{lookupUniprotSpeciesFromTaxId} each return a character vector of
  possible values.

  \code{taxId} returns a numeric value that corresponds to the taxonomy
  ID.

  \code{select} and \code{availableUniprotSpecies} each return a
  data.frame.
}

\author{Marc Carlson}

\seealso{
  \code{select}
}

\examples{
## Make a UniProt.ws object
up <- UniProt.ws(taxId=9606)

## look at the object
up

## get the current species
species(up)

## look up available species with their tax ids
availableUniprotSpecies("musculus")

## get the current taxId
taxId(up)

## look up the species that goes with a tax id
lookupUniprotSpeciesFromTaxId(9606)

## set the taxId to something else
taxId(up) <- 10090
up

## list the possible key types
head(keytypes(up))

## list all possible keys of type entrez gene ID
egs <- keys(up, "GeneID")

## use select to extract some data
res <- select(
    x = up,
    keys = c("22627","22629"),
    columns = c("xref_pdb","xref_hgnc","sequence"),
    keytype = "GeneID"
)
res
}


\keyword{classes}
\keyword{methods}

