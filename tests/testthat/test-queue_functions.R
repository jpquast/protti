context("test-queue_functions")

queue <- create_queue(
  date = c("200722"),
  instrument = c("EX1"),
  user = c("username"),
  measurement_type = c("DIA"),
  experiment_name = c("N01"),
  digestion = c("LiP", "tryptic control"),
  treatment_type_1 = c("EDTA", "H2O"),
  treatment_type_2 = c("Zeba", "unfiltered"),
  treatment_dose_1 = c(10, 30, 60),
  treatment_unit_1 = c("min"),
  n_replicates = 4,
  number_runs = FALSE,
  organism = c("E. coli"),
  exclude_combinations = list(list(
    treatment_type_1 = c("H2O"),
    treatment_type_2 = c("Zeba", "unfiltered"),
    treatment_dose_1 = c(10, 30)
  )),
  inj_vol = c(2),
  data_path = "D:\\2007_Data",
  method_path = "C:\\Xcalibur\\methods\\username\\DIA_120min_41var_AGC200",
  position_row = c("A", "B", "C", "D", "E", "F"),
  position_column = 8,
  blank_every_n = 4,
  blank_position = "1-V1",
  blank_method_path = "C:\\Xcalibur\\methods\\blank",
  export = FALSE
)

test_that("create_queue works", {
  expect_is(queue, "data.frame")
  expect_equal(ncol(queue), 21)
  expect_equal(nrow(queue), 80)
})

test_that("randomise_queue works", {
  set.seed(123)
  randomised_queue <- randomise_queue(data = queue, rows = 71:80)
  expect_is(randomised_queue, "data.frame")
  expect_equal(ncol(randomised_queue), 21)
  expect_equal(nrow(randomised_queue), 80)
  expect_equal(randomised_queue$Position[71:80], c("1-V1", "B8", "B5", "B3", "B6", "1-V1", "B7", "B4", "B1", "B2"))
})
