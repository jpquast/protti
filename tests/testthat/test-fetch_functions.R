context("test-fetch_functions")


test_that("fetch_uniprot works", {
  unis <- c("iRT", "P36578", "O43324", "Q00796")
  expect_warning(uniprot <- fetch_uniprot(unis))
  expect_is(uniprot, "data.frame")
  expect_equal(nrow(uniprot), 3)
  expect_equal(ncol(uniprot), 17)
})
