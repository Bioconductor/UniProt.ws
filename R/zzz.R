keytypeKeysDat <- NULL

.onLoad <- function(...) {
    keytypeKeysDat <<- tryCatch({
        updateKeytypes()
        processAvailableSpeciesFiles()
    }, error = function(err) {
        warning(
            "failed to download:",
            "\n  ", conditionMessage(err),
            "\nusing cached version"
        )
        file <- system.file('extdata','keytypes.txt', package='UniProt.ws')
        read.delim(file, header=FALSE, stringsAsFactors=FALSE)
    })
}
