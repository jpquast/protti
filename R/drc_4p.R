

drc_4p <- function(data, response, dose, pb) {
  pb$tick()
  tryCatch({
    suppressWarnings(drc::drm(
      stats::as.formula(paste(ensym(response), "~", ensym(dose))),
      data = data,
      fct = drc::LL.4(names = c("hill", "min_value", "max_value", "ec_50"))
    ))
  }, error = function(error) {
    c("no_fit")
  })
}


# drc_4p <- function(data, response, dose, pb) {
#   pb$tick()
#   tryCatch({
#     suppressWarnings(dr4pl::dr4pl(
#       stats::as.formula(paste(ensym(response), "~", ensym(dose))),
#       data = data,
#       method.init = "logistic",
#       use.Hessian = TRUE
#       ))
#   }, error = function(error) {
#     c("no_fit")
#   })
# }