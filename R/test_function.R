test_function <- function(data, sample){
  data %>% 
    dplyr::distinct({{sample}})
}
