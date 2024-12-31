res <- taxname2species("PIG")
expect_identical(res, "Sus scrofa")

org <- c("PIG", "YEAST", "HUMAN", "MOUSE", "TRIHA", "THEAS", "SIVAM", "AERPX")
res <- taxname2species(org)
expect_identical(
    res,
    c(
        "Sus scrofa","Saccharomyces cerevisiae","Homo sapiens",
        "Mus musculus","Trichoderma harzianum",
        "Thermanaerovibrio acidaminovorans",
        "Simian immunodeficiency virus","Aeropyrum pernix"
    )
)

res <- taxname2taxid("PIG")
expect_true(res == "9823")

org <- c("PIG", "YEAST", "HUMAN", "MOUSE", "TRIHA", "THEAS", "SIVAM", "AERPX")
res <- taxname2taxid(org)
expect_identical(
    res,
    as.integer(c(9823, 559292, 9606, 10090, 5544, 525903, 36378, 56636))
)

res <- taxname2domain("PIG")
expect_true(res == "E")

org <- c("PIG", "YEAST", "HUMAN", "MOUSE", "TRIHA", "THEAS", "SIVAM", "AERPX")
res <- taxname2domain(c(org, "MYCGI"))
expect_identical(
    res,
    factor(
        c("E","E","E","E","E","B","V","A","B"),
        levels = c("A", "B", "E", "V", "X")
    )
)
