## BiocBaseUtils is not available in RELEASE_3_15
setSlots <- function (object, ..., check = TRUE)
{
    if (!isTRUEorFALSE(check))
        stop("'check' must be TRUE or FALSE")
    object <- unsafe_replaceSlots(object, ...)
    if (check)
        validObject(object)
    object
}

unsafe_replaceSlots <- function (object, ..., .slotList = list())
{
    slots <- c(list(...), .slotList)
    slots_names <- names(slots)
    for (i in seq_along(slots)) {
        slot_name <- slots_names[[i]]
        if (identical(slot_name, "mcols"))
            slot_name <- "elementMetadata"
        slot_val <- slots[[i]]
        slot(object, slot_name, check = FALSE) <- slot_val
    }
    object
}

# adapted from S4Vectors
.isSingle <- function(x, na.ok, FUN)
{
    FUN(x) && identical(length(x), 1L) && (na.ok || !is.na(x))
}

isTRUEorFALSE <- function(x, na.ok = FALSE)
{
    .isSingle(x, na.ok, is.logical)
}

isScalarCharacter <- function(x, na.ok = FALSE, zchar = FALSE)
{
    identical(length(x), 1L) && isCharacter(x, na.ok, zchar)
}

isCharacter <- function(x, na.ok = FALSE, zchar = FALSE)
{
    is.character(x) && (na.ok || all(!is.na(x))) && (zchar || all(nzchar(x)))
}

