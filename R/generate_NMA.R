#' Generate a simulated Network Meta-Analysis Dataset
#'
#' This function creates a simulated dataset for network meta-analysis with binary outcomes.
#' Users can specify the number of studies, treatments, treatment effects, and sample size range.
#'
#' @param num_studies Integer. Number of studies to simulate.
#' @param num_treatments Integer. Number of treatments in the network.
#' @param treatment_effects Numeric vector. Probabilities of success for each treatment.
#' @param sample_size_range Integer vector of length 2. Range of sample sizes (default: c(50, 200)).
#' @return A data frame with columns: study.id, treatment.id, sample.size, and response.
#' @examples
#' set.seed(12345)
#' data <- generate_NMA(
#'   num_studies = 10,
#'   num_treatments = 3,
#'   treatment_effects = c(0.3, 0.5, 0.7),
#'   sample_size_range = c(50, 200)
#' )
#' print(data)
#' @export


generate_NMA <- function(num_studies, num_treatments, treatment_effects, sample_size_range = c(50, 200)) {
  # Check inputs
  if (length(treatment_effects) != num_treatments) {
    stop("The length of treatment_effects must equal the number of treatments.")
  }

  if (any(treatment_effects < 0 | treatment_effects > 1)) {
    stop("Treatment effects must be probabilities between 0 and 1.")
  }

  # Initialize an empty data frame
  dat <- data.frame(study.id = integer(),
                     treatment.id = integer(),
                     sample.size = integer(),
                     response = integer())

  # Generate data for each study
  for (study.id in 1:num_studies) {
    # Randomly assign treatments to the study (at least 2 treatments per study)
    num_treatments_in_study <- sample(2:num_treatments, 1)
    treatments_in_study <- sample(1:num_treatments, num_treatments_in_study)

    # Generate data for each treatment in the study
    for (treatment.id in treatments_in_study) {
      # Randomly generate sample size
      sample.size <- sample(sample_size_range[1]:sample_size_range[2], 1)

      # Generate the number of successes based on the treatment effect
      response <- rbinom(1, size = sample.size, prob = treatment_effects[treatment.id])

      dat <- rbind(dat, data.frame(study.id = study.id,
                                     treatment.id = treatment.id,
                                     sample.size = sample.size,
                                     response = response))
    }
  }

  return(dat)
}

