context("test-fetch_functions")

test_that("fetch_uniprot works", {
  unis <- c("iRT", "P36578", "O43324", "Q00796")
  expect_warning(uniprot <- fetch_uniprot(unis))
  expect_is(uniprot, "data.frame")
  expect_equal(nrow(uniprot), 3)
  expect_equal(ncol(uniprot), 17)
})

test_that("fetch_uniprot_proteome works", {
  proteome <- fetch_uniprot_proteome(organism_id = "83333")
  expect_is(proteome, "data.frame")
  expect_equal(ncol(proteome), 1)
  expect_gt(nrow(proteome), 10)
})

test_that("fetch_mobidb works", {
  unis <- c("iRT", "P25437", "P30870", "P0A6P9")
  expect_warning(mobidb <- fetch_mobidb("83333", unis))
  expect_is(mobidb, "data.frame")
  expect_equal(nrow(mobidb), 2)
  expect_equal(ncol(mobidb), 7)
})

test_that("fetch_chebi works", {
  database <- fetch_chebi()
  relations <- fetch_chebi(relation = TRUE)
  expect_is(database, "data.frame")
  expect_is(relations, "data.frame")
  expect_equal(ncol(database), 13)
  expect_equal(ncol(relations), 3)
  expect_gt(nrow(database), 10)
  expect_gt(nrow(relations), 10)
})

test_that("fetch_kegg works", {
  ecoli <- fetch_kegg(species = "eco")
  expect_is(ecoli, "data.frame")
  expect_equal(ncol(ecoli), 4)
  expect_gt(nrow(ecoli), 10)
})

test_that("fetch_go works", {
  go_eco <- fetch_go(organism_id = "83333")
  go_sac <- fetch_go(organism_id = "559292")
  go_hs <- fetch_go(organism_id = "9606")
  expect_is(go_eco, "data.frame")
  expect_equal(ncol(go_eco), 17)
  expect_gt(nrow(go_eco), 10)
  expect_is(go_sac, "data.frame")
  expect_equal(ncol(go_sac), 17)
  expect_gt(nrow(go_sac), 10)
  expect_is(go_hs, "data.frame")
  expect_equal(ncol(go_hs), 17)
  expect_gt(nrow(go_hs), 10)
})


