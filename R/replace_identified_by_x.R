#' Replace identified positions in protein sequence by "x"
#'
#' Helper function for the calculation of sequence coverage, replaces identified positions with an "x" within the protein sequence.
#'
#' @param sequence Sequence a character vector that contains the protein sequence.
#' @param positions_start A vector of start positions of the identified peptides.
#' @param positions_end A vector of end positions of the identified peptides.
#'
#' @return A modified protein sequence with each identified position replaced by "x"
#' @importFrom purrr map2
#' @importFrom stringr str_sub
#'
#' @examples
#' \dontrun{
#' replace_identified_by_x(
#' sequence = c("AEFGPEEAAVS"),
#' positions_start = c(1, 5, 8),
#' positions_end = c(3, 6, 11)
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
