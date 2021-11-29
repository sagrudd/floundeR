

context("Infographic R6 object")
test_that("Infographic R6 object", {

    expect_error(Infographic$new(), NA)
    info <- Infographic$new()
    expect_error(info$add(123))

    item1 <- InfographicItem$new(key="Alignment", value="1234", icon="fa-sliders")
    expect_error(info$add(item1), NA)
    expect_equal(info$items, 1)

})




context("InfographicItem R6 object")
test_that("InfographicItem R6 object", {


    expect_error(InfographicItem$new(key="Alignment", value="1234", icon="fa-sliders"), NA)

})

