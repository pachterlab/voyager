test_that("Correct output of listSFEMethods", {
    out <- listSFEMethods("uni", "local")
    expect_s3_class(out, "data.frame")
    expect_named(out, c("name", "description"))
    expect_true(all(vapply(out, is.character, FUN.VALUE = logical(1))))

    out2 <- listSFEMethods("multi")
    expect_named(out2, c("name", "description"))
    expect_true(all(vapply(out2, is.character, FUN.VALUE = logical(1))))
})
