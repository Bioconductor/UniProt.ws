.onLoad <- function(...) {
    tryCatch({
        processAvailableSpeciesFiles()
    }, error = function(err) {
        warning(
            "failed to download:",
            "\n  ", conditionMessage(err),
            "\nusing cached version"
        )
    })
}
