keytypeKeysDat <- NULL

.onLoad <- function(...) {
    tryCatch({
        updateKeytypes()
        processAvailableSpeciesFiles()
    }, error = function(err) {
        warning(
            "failed to download:",
            "\n  ", conditionMessage(err),
            "\nusing cached version"
        )
        file <- system.file('extdata','keytypes.txt', package='UniProt.ws')=
    })
    keytypeKeysDat <<- read.delim(file, header=TRUE, stringsAsFactors=FALSE)
}
