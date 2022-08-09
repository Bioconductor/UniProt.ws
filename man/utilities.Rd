\name{utilities}
\alias{taxname2domain}
\alias{taxname2species}
\alias{taxname2taxid}
\alias{utilities}

\title{Utility functions}

\usage{
taxname2species(taxname, specfile)
taxname2taxid(taxname, specfile)
taxname2domain(taxname, specfile)
}

\arguments{
  \item{taxname}{
    Character string up to 6 uppercase characters, like HUMAN, MOUSE, or
    AERPX. Also works for a vector of such taxon names.
  }
  \item{specfile}{
    An optional local file where speclist.RData is saved from UniProt.org.
    When \code{specfile} is missing, a cached file from the
    extdata/ package directory is used.
  }
}

\value{
Function \code{taxname2species} returns a character vector of scientific taxon
names matching to the UniProt taxon names supplied as \code{taxname}.

Function \code{taxname2taxid} returns a numeric vector of Taxonomy IDs
matching to the UniProt taxon names supplied as \code{taxname}.

Function \code{taxname2domain} returns a character vector of one letter domain
symbols matching to the UniProt taxon names supplied as \code{taxname}.

Function \code{updatespecfile} is deprecated and no longer required as the
specfile is self updating based on \code{BiocFileCache}'s \code{bfcneedsupdate}
when necessary.
}

\description{
UniProt uses custom coding of organism names from which protein sequences they
store. These taxon names are used also in the protein names (not in the UniProt
IDs!). These functions help to translate those names to standard scientific
(Latin) taxon names and other useful identifiers.

\itemize{
  \item{taxname2species():}{
    converts UniProt taxonomy names to scientific species names
  }
  \item{taxname2taxid():}{
    converts UniProt taxonomy names to NCBI Taxonomy IDs
  }
  \item{taxname2domain():}{
    converts UniProt taxonomy names to the following taxonomical domains:
    'A' for archaea (=archaebacteria)\\
    'B' for bacteria (=prokaryota or eubacteria)\\
    'E' for eukaryota (=eukarya)\\
    'V' for viruses and phages (=viridae)\\
    'O' for others (such as artificial sequences)\\
  }
  \item{updatespecfile():}{
    The updatespecfile helper function attempts to download the
    current version of the controlled vocabulary of species table from
    \href{http://www.uniprot.org/docs/speclist}{UniProt controlled
    vocabulary of species}. If it fails to download,
    an archived version of the table in (in extdata/) will be used.
  }
}
}

\seealso{
\href{http://www.uniprot.org/docs/speclist}{UniProt controlled vocabulary of
species}, which defines the taxon names.
}

\author{Csaba Ortutay}

\examples{
taxname2species("PIG")
taxname2species(c("PIG","HUMAN","TRIHA"))

taxname2taxid("PIG")
taxname2taxid(c("PIG","HUMAN","TRIHA"))

taxname2domain("PIG")
taxname2domain(c("PIG","HUMAN","TRIHA"))
}
