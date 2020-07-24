#' Replace each identified positions in protein by "x"
#'
#' Helper function for the calculation of sequence coverage, replaces identified positions with an "x"
#'
#' @param sequence Protein sequence
#' @param positions_start Start position of the identified peptide.
#' @param positions_end End position of the identified peptide.
#'
#' @return A modified protein sequence with each identified position replaced by "x"
#' @import dplyr
#' @import purrr
#' @import tidyr
#' @importFrom magrittr %>%
#' @importFrom stringr str_sub
#'
#' @examples
#' \dontrun{
#' replace_identified_by_x(
#' sequence = protein_sequence,
#' positions_start = start,
#' positions_end = end
#' )
#' }
#'
replace_identified_by_x <- function(sequence, positions_start, positions_end)
  {
  result <- purrr::map2(.x = positions_start, .y = positions_end,
                        function(x, y){
    times <- y-x + 1
    stringr::str_sub(sequence, start = x, end = y) <- paste(rep("x", times = times), collapse = "")
    sequence <<- sequence
  })
  result[[length(result)]]
}
