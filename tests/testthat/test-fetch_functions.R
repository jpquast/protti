context("test-fetch_functions")


test_that("fetch_uniprot works", {
  unis <- c("iRT", "P36578", "O43324", "Q00796")
  expect_warning(uniprot <- fetch_uniprot(unis))
  expect_is(uniprot, "data.frame")
  expect_equal(nrow(uniprot), 3)
  expect_equal(ncol(uniprot), 17)
})

test_that("fetch_mobidb works", {
  unis <- c("iRT", "P25437", "P30870", "P0A6P9")
  expect_warning(mobidb <- fetch_mobidb("83333", unis))
  expect_is(mobidb, "data.frame")
  expect_equal(nrow(mobidb), 2)
  expect_equal(ncol(mobidb), 7)
})

