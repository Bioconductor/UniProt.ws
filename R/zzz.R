.onLoad <- function(libname, pkgname)
{
  ns <- asNamespace(pkgname)
  UniProt <- UniProt.ws(taxId=9606)
  assign("UniProt.ws", UniProt, envir=ns)
  namespaceExport(ns, "UniProt.ws")
}


