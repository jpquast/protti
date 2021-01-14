#' Creates a synthetic limited proteolysis proteomics dataset
#'
#' This function creates a synthetic limited proteolysis proteomics dataset that can be used to test functions while knowing the ground truth.
#'
#' @param n_proteins Numeric, specifies the number of proteins in the synthetic dataset.
#' @param frac_change Numeric, the fraction of proteins that has a peptide changing in abundance. So far only one peptide per protein
#' is changing.
#' @param n_replicates Numeric, the number of replicates per condition.
#' @param n_conditions Numeric, the number of conditions.
#' @param method Character, specifies the method type for the random sampling of significantly changing peptides. If \code{method = "random_effect"},
#' the effect for each condition is randomly sampled and conditions do not depend on eachother. If \code{method = "dose_response"},
#' the effect is sampled based on a dose response curve and conditions are related to eachother depending on the curve shape. In this case
#' the concentrations argument needs to be specified.
#' @param concentrations Numeric vector of the length of number of conditions, only needs to be specified if \code{method = "dose_response"}.
#' This allows equal sampeling of peptide intensities. It ensures that the same positions of dose response curves are sampled for each peptide
#' based on the provided concentrations.
#' @param median_offset_sd Numeric, standard deviation of normal distribution that is used for sampling of inter-sample-differences. 
#' Default is 0.05.
#' @param mean_protein_intensity Numeric, mean of the protein intensity distribution. Default: 16.8.
#' @param sd_protein_intensity Numeric, standard deviation of the protein intensity distribution. Default: 1.4.
#' @param mean_n_peptides Numeric, mean number of peptides per protein. Default: 12.75.
#' @param size_n_peptides Numeric, dispersion parameter (the shape parameter of the gamma mixing distribution). Can be theoretically calculated
#' as \code{mean + mean^2/variance}, however, it should be rather obtained by fitting the negative binomial distribution to real data. This
#' can be done by using the \code{optim} function (see Example section). Default: 0.9.
#' @param mean_sd_peptides Numeric, mean of peptide intensity standard deviations within a protein. Default: 1.7.
#' @param sd_sd_peptides Numeric, standard deviation of peptide intensity standard deviation within a protein. Default: 0.75.
#' @param mean_log_replicates,sd_log_replicates Numeric, \code{meanlog} and \code{sdlog} value of the log normal
#' distribution of replicate standard deviations. Can be obtained by fitting a log normal distribution to the distribution
#' of replicate standard deviations from a real dataset. This can be done using
#' the \code{optim} function (see Example section). Default: -2.2 and 1.05.
#' @param effect_sd Numeric, standard deviation of a normal distribution around \code{mean = 0} that is used to sample the effect of
#' significantly changeing peptides. Default: 2.
#' @param dropout_curve_inflection Numeric, intensity inflection point of a probabilistic dropout curve that is used to sample intensity
#' dependent missing values. This argument determines how many missing values there are in the dataset. Default: 14.
#' @param dropout_curve_sd Numeric, standard deviation of the probabilistic dropout curve. Needs to be negative to sample a droupout towards
#' low intensities. Default: -1.2.
#' @param additional_metadata Logical, determines if metadata such as protein coverage, missed cleavages and charge state should be 
#' sampled and added to the list.
#'
#' @return A data frame that contains complete peptide intensities and peptide intensities with values that were created based on a
#' probabilistic dropout curve.
#' @importFrom stats rnorm rnbinom rlnorm pnorm runif
#' @importFrom tibble tibble
#' @importFrom dplyr group_by mutate select n
#' @importFrom tidyr unnest
#' @importFrom stringr str_extract
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' create_synthetic_data(
#'   n_proteins = 3000,
#'   frac_change = 0.01,
#'   n_replicates = 3,
#'   n_conditions = 2
#' )
#'
#' # determination of mean_n_peptides and size_n_peptides parameters based on real data (count)
#'
#' theta <- c(mu = 1, k = 1)
#' negbinom <- function(theta) {
#'   -sum(stats::dnbinom(count, mu = theta[1], size = theta[2], log = TRUE))
#' }
#' fit <- stats::optim(theta, negbinom)
#' fit
#'
#' # determination of mean_log_replicates and sd_log_replicates parameters
#' # based on real data (standard_deviations)
#'
#' theta2 <- c(meanlog = 1, sdlog = 1)
#' lognorm <- function(theta2) {
#'   -sum(stats::dlnorm(standard_deviations, meanlog = theta2[1], sdlog = theta2[2], log = TRUE))
#' }
#' fit2 <- stats::optim(theta2, lognorm)
#' fit2
#' }
create_synthetic_data <- function(
                                  n_proteins,
                                  frac_change,
                                  n_replicates,
                                  n_conditions,
                                  method = "random_effect",
                                  concentrations = NULL,
                                  median_offset_sd = 0.05,
                                  mean_protein_intensity = 16.88,
                                  sd_protein_intensity = 1.4,
                                  mean_n_peptides = 12.75,
                                  size_n_peptides = 0.9,
                                  mean_sd_peptides = 1.7,
                                  sd_sd_peptides = 0.75,
                                  mean_log_replicates = -2.2,
                                  sd_log_replicates = 1.05,
                                  effect_sd = 2,
                                  dropout_curve_inflection = 14,
                                  dropout_curve_sd = -1.2,
                                  additional_metadata = TRUE) {
  # the amount of proteins that change
  n_change <- round(n_proteins * frac_change)

  # sample protein intensities
  sampled_protein_intensities <- stats::rnorm(n = n_proteins, mean = mean_protein_intensity, sd = sd_protein_intensity)

  # sample the amount of peptides per protein. There is no relationship between amount of peptides per protein and protein intensity
  sampled_n_peptides <- stats::rnbinom(n = n_proteins * 3, size = size_n_peptides, mu = mean_n_peptides)
  sampled_n_peptides <- sampled_n_peptides[sampled_n_peptides > 0][1:n_proteins]

  # sample the standard deviations for the distribution of peptide intensities in a protein
  sampled_peptide_sd <- stats::rnorm(n = n_proteins * 2, mean = mean_sd_peptides, sd = sd_sd_peptides)
  # remove negative values for sd
  sampled_peptide_sd <- sampled_peptide_sd[sampled_peptide_sd > 0][1:n_proteins]

  # sample peptide intensities and standard deviations
  proteins <- tibble::tibble(protein = paste0("protein_", 1:n_proteins), n = sampled_n_peptides, mean = sampled_protein_intensities, sd = sampled_peptide_sd) %>%
    dplyr::group_by(.data$protein) %>%
    dplyr::mutate(peptide_intensity_mean = list(round(stats::rnorm(n = .data$n, mean = .data$mean, sd = .data$sd), digits = 4))) %>%
    tidyr::unnest(.data$peptide_intensity_mean) %>%
    dplyr::mutate(peptide = paste0("peptide_", stringr::str_extract(.data$protein, pattern = "\\d+"), "_", 1:dplyr::n())) %>%
    dplyr::mutate(replicate_sd = round(stats::rlnorm(n = 1, meanlog = mean_log_replicates, sdlog = sd_log_replicates), digits = 4)) %>%
    dplyr::select(-c(.data$mean, .data$sd))

  # sample peptide intensities for replicates and conditions
  proteins_replicates <- proteins %>%
    dplyr::group_by(.data$peptide) %>%
    dplyr::mutate(condition = list(sort(rep(paste0("condition_", 1:n_conditions), n_replicates)))) %>%
    tidyr::unnest(c(.data$condition)) %>%
    dplyr::mutate(sample = paste0("sample_", 1:(n_conditions * n_replicates))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(peptide_intensity = stats::rnorm(n = dplyr::n(), mean = .data$peptide_intensity_mean, sd = .data$replicate_sd))

  # apply median offset for each sample to simulate measurement differences between samples. for example caused by differnt concentrations.
  n_peptides <- length(unique(proteins_replicates$peptide))
  offset <- rep(stats::rnorm(n_conditions * n_replicates, mean = 0, sd = median_offset_sd), n_peptides)
  
  # sample significantly changing peptides

  if (method == "random_effect") {
    # sample the direction of the effect randomly
    proteins_replicates_change <- proteins_replicates %>%
      dplyr::mutate(change = stringr::str_extract(.data$protein, pattern = "\\d+") %in% 1:n_change) %>%
      dplyr::group_by(.data$protein) %>%
      dplyr::mutate(n_change_peptide = ifelse(.data$change == TRUE, ceiling(stats::rgamma(1, shape = 0.7, rate = 0.4)), 0)) %>% # sample number of significant peptides based on gamma distribution
      dplyr::mutate(n_change_peptide = ifelse(.data$n_change_peptide >= .data$n, .data$n, .data$n_change_peptide)) %>%
      dplyr::mutate(change_peptide = .data$change & stringr::str_extract(.data$peptide, pattern = "\\d+$") %in% 1:unique(.data$n_change_peptide)) %>%
      dplyr::mutate(effect = stats::rnorm(dplyr::n(), mean = 0, sd = effect_sd)) %>%
      dplyr::mutate(effect = ifelse(.data$change_peptide == TRUE, .data$effect, 0)) %>%
      dplyr::group_by(.data$condition, .data$peptide) %>%
      dplyr::mutate(effect = rep(.data$effect[1], n_replicates)) %>%
      dplyr::mutate(peptide_intensity = .data$peptide_intensity + .data$effect) %>%
      dplyr::bind_cols(offset = offset) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(peptide_intensity = .data$peptide_intensity + .data$offset) %>% 
      dplyr::select(-c(.data$peptide_intensity_mean, .data$replicate_sd, .data$effect, .data$n, .data$n_change_peptide, .data$offset))
  }

  if (method == "dose_response") {
    # sample the direction of the effect based on a dose response curve

    condition_concentration <- tibble::tibble(condition = paste0("condition_", 1:n_conditions), concentration = concentrations)

    proteins_replicates_change <- proteins_replicates %>%
      dplyr::mutate(change = stringr::str_extract(.data$protein, pattern = "\\d+") %in% 1:n_change) %>%
      dplyr::group_by(.data$protein) %>%
      dplyr::mutate(n_change_peptide = ifelse(.data$change == TRUE, ceiling(stats::rgamma(1, shape = 0.7, rate = 0.4)), 0)) %>% # sample number of significant peptides based on gamma distribution
      dplyr::mutate(n_change_peptide = ifelse(.data$n_change_peptide >= .data$n, .data$n, .data$n_change_peptide)) %>%
      dplyr::mutate(change_peptide = .data$change & stringr::str_extract(.data$peptide, pattern = "\\d+$") %in% 1:unique(.data$n_change_peptide)) %>%
      dplyr::left_join(condition_concentration, by = "condition") %>%
      dplyr::mutate(effect_total = stats::rnorm(dplyr::n(), mean = 0, sd = effect_sd)) %>%
      dplyr::mutate(effect_total = ifelse(.data$change_peptide == TRUE, .data$effect_total, 0)) %>%
      dplyr::group_by(.data$peptide) %>%
      dplyr::mutate(effect_total = rep(.data$effect_total[1], (n_replicates * n_conditions))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(b = sample(c(stats::rlnorm(dplyr::n(), meanlog = 0.6, sdlog = 0.4), -rlnorm(dplyr::n(), meanlog = 0.6, sdlog = 0.4)), size = dplyr::n())) %>%
      dplyr::mutate(c = stats::rlnorm(dplyr::n(), meanlog = log(mean(concentrations) / 2), sdlog = abs(log(mean(concentrations) / 2) / 3))) %>%
      dplyr::mutate(b = ifelse(.data$change_peptide == TRUE, .data$b, 0)) %>%
      dplyr::mutate(c = ifelse(.data$change_peptide == TRUE, .data$c, 0)) %>%
      dplyr::group_by(.data$peptide) %>%
      dplyr::mutate(b = rep(.data$b[1], (n_replicates * n_conditions))) %>%
      dplyr::mutate(c = rep(.data$c[1], (n_replicates * n_conditions))) %>%
      dplyr::mutate(effect = ifelse(.data$change_peptide == TRUE, .data$effect_total * (1 + (-1 / (1 + (.data$concentration / .data$c)^.data$b))), 0)) %>%
      dplyr::mutate(peptide_intensity = .data$peptide_intensity + .data$effect) %>%
      dplyr::bind_cols(offset = offset) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(peptide_intensity = .data$peptide_intensity + .data$offset) %>% 
      dplyr::select(-c(.data$peptide_intensity_mean, .data$replicate_sd, .data$n, .data$n_change_peptide, .data$effect, .data$effect_total, .data$b, .data$c, .data$offset))
    # formula for inflection point and slope sampling roughly simulates the behaviour of real data. They have been figured out by trial and error.
  }

  # sample NA values based on probabilistic dropout curve
  proteins_replicates_change_missing <- proteins_replicates_change %>%
    dplyr::mutate(dropout_probability = stats::pnorm(.data$peptide_intensity, mean = dropout_curve_inflection, sd = -dropout_curve_sd, lower.tail = FALSE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(peptide_intensity_missing = ifelse(stats::runif(dplyr::n()) > .data$dropout_probability, .data$peptide_intensity, NA)) %>%
    dplyr::select(-.data$dropout_probability) %>% 
    dplyr::group_by(.data$peptide) %>% 
    dplyr::mutate(isna = sum(!is.na(.data$peptide_intensity_missing))) %>% # remove peptides for which every intensity is NA after dropout
    dplyr::filter(.data$isna > 0) %>% 
    dplyr::select(-.data$isna) %>% 
    dplyr::ungroup()
  
  if (additional_metadata == FALSE) return(proteins_replicates_change_missing)
  if (additional_metadata == TRUE) {
    # adding coverage estimates based on gamma distribution
    coverage_sampled <- stats::rgamma(nrow(proteins_replicates_change_missing) * 2, shape = 1.2, rate = 0.05) 
    coverage_sampled <- coverage_sampled[coverage_sampled <= 100] # remove coverage over 100% because not possible
    
    coverage_data <- proteins_replicates_change_missing %>% 
      dplyr::mutate(coverage = coverage_sampled[1:dplyr::n()]) %>% 
      dplyr::group_by(.data$protein) %>%
      dplyr::mutate(coverage = rep(.data$coverage[1], dplyr::n())) %>% 
      dplyr::mutate(coverage_peptide = .data$coverage / dplyr::n_distinct(.data$peptide)) %>% 
      dplyr::group_by(.data$sample, .data$protein) %>% 
      dplyr::mutate(coverage = sum(!is.na(.data$peptide_intensity_missing)) * .data$coverage_peptide) %>% 
      dplyr::select(-.data$coverage_peptide) %>% 
      dplyr::ungroup()
    
    # adding missed cleavage estimates based on poisson distribution
    
    missed_cleavage_sampled <- stats::rpois(nrow(proteins_replicates_change_missing) * 2, 0.28)
    missed_cleavage_sampled <- missed_cleavage_sampled[missed_cleavage_sampled < 3] # remove missed cleavages over 3 because they should not occure in data
    
    missed_cleavages_data <- coverage_data %>% 
      dplyr::mutate(n_missed_cleavage = missed_cleavage_sampled[1:dplyr::n()]) %>% 
      dplyr::group_by(.data$peptide) %>% 
      dplyr::mutate(n_missed_cleavage = rep(.data$n_missed_cleavage[1], dplyr::n())) %>% 
      dplyr::ungroup()
    
    # add charge state estimates based on rounded gamma distribution
    
    charge_sampled <- round(stats::rgamma(nrow(proteins_replicates_change_missing) * 2, shape = 13.06, rate = 5.63))
    charge_sampled <- charge_sampled[charge_sampled > 0 & charge_sampled < 7] # remove chage state of 0 and higher than 6

    charge_data <- missed_cleavages_data %>% 
      dplyr::mutate(charge = charge_sampled[1:dplyr::n()]) %>% 
      dplyr::group_by(.data$peptide) %>% 
      dplyr::mutate(charge = rep(.data$charge[1], dplyr::n())) %>% 
      dplyr::ungroup()
    
    # add peptide type sampling based on the following probabilities: fully tryptic = 55%, semi tryptic = 40%, non tryptic = 5%
    peptide_types <- c(rep("fully-tryptic", 11), rep("semi-tryptic", 8), rep("non-tryptic", 1))
    
    peptide_type_data <- charge_data %>% 
      dplyr::mutate(pep_type = sample(peptide_types, size = dplyr::n(), replace = TRUE)) %>% 
      dplyr::group_by(.data$peptide) %>% 
      dplyr::mutate(pep_type = rep(.data$pep_type[1], dplyr::n())) %>% 
      dplyr::ungroup()
    
    # add peak width estimates based on gamma distribution, an associated retention time is sampled with a uniform
    # distribution form 0 to 120. This is not how peak width is actually associated with retention time, but a simple way of
    # obtaining values. The real relationship is very complex and not easy to sample.
    
    peak_width_sampled <- stats::rgamma(nrow(proteins_replicates_change_missing), shape = 10.4, rate = 36.21)
    
    peak_width_data <- peptide_type_data %>% 
      dplyr::mutate(peak_width = peak_width_sampled) %>% 
      dplyr::mutate(retention_time = stats::runif(n = dplyr::n(), min = 0, max = 120)) %>% 
      dplyr::group_by(.data$peptide) %>% 
      dplyr::mutate(peak_width = rep(.data$peak_width[1], dplyr::n())) %>% 
      dplyr::mutate(retention_time = rep(.data$retention_time[1], dplyr::n())) %>% 
      dplyr::ungroup()
    
    peak_width_data
  }
  
}
