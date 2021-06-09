#' Convert metal names to search pattern
#'
#' Converts a vector of metal names extracted from the \code{feature_metal_binding} column obtained with \code{fetch_uniprot} to a pattern that can be used to search for corresponding ChEBI IDs. This is used as a helper function for other functions.
#'
#' @param metal_names A character vector containing names of metals and metal containing molecules.
#'
#' @return A character vector with metal name search patterns.
#' @importFrom purrr map
#' @importFrom stringr str_detect str_split str_replace_all str_extract_all str_squish
split_metal_name <- function(metal_names) {
  purrr::map_chr(metal_names, function(x) {
    if (is.na(x)) {
      result <- NA
      return(result)
    }
    # Add full name to pattern
    result <- x
    # Add split name to pattern
    x <- stringr::str_squish(stringr::str_replace_all(x, pattern = "[//(//)]", replacement = " ")) # replace bracket with one whitespace
    if (stringr::str_detect(x, pattern = " ")) {
      x1 <- as.vector(stringr::str_split(x, pattern = " ", simplify = TRUE))
      # Remove specific words from split name
      x1 <- x1[!x1 %in% c("or", "Divalent", "Monovalent", "metal", "cation", "Iron-sulfur", "low-spin", "high-spin", "A", "A1", "A2", "A3", "C", "1", "1+", "2+", "2", "3", "3+", "4", "5", "6", "B", "b", "o", "distal", "axial", "proximal", "ligand", "heme", "siroheme", "methylcob", "Cob", "II", "III", "alamin", "b595", "b558", "b562", "b566")]
      result <- append(result, x1)
    }
    # Add split "-"-containing name to pattern
    if (stringr::str_detect(x, pattern = "-")) {
      x2 <- as.vector(
        stringr::str_extract_all(x, pattern = "\\b\\w+\\-\\w+\\b", simplify = TRUE)
      )
      # Remove specific split "-"-containing names
      x2_clean <- x2[!x2 %in% c("Iron-sulfur", "S-AdoMet", "Iron-oxo", "4Fe-2O", "low-spin", "high-spin", result)] # remove these names and also names that are already part of the result
      result <- append(result, x2_clean)
      x3 <-
        as.vector(stringr::str_split(
          x2,
          pattern = "-",
          simplify = TRUE
        ))
      # Remove specific words from split "-"-containing name
      x3 <- x3[!x3 %in% c("sulfur", "4Fe", "4S", "3Fe", "3S", "2Fe", "2S", "S", "AdoMet", "oxo", "2O", "low", "high", "spin")]
      result <- append(result, x3)
    }
    paste(result, collapse = "|")
  })
}
