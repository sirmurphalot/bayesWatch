# Script file for fitting a Hidden Markov Model with Dirichlet Prior
## Author: Alexander Murph
require(parallel)
require(Rcpp)
require(Matrix)
require(CholWishart)
require(Hotelling)
source("R/helpers.R")


#' Title
#'
#' @param data_woTimeValues
#' @param time_of_observations
#' @param time_points
#' @param not.cont
#' @param iterations
#' @param burnin
#' @param lower_bounds
#' @param upper_bounds
#' @param previous_states
#' @param previous_model_fits
#' @param sparsity_parameter
#' @param linger_parameter
#' @param move_parameter
#' @param g.prior
#' @param set_G
#' @param wishart_df_inital
#' @param lambda
#' @param g_sampling_distribution
#' @param n.cores
#' @param scaleMatrix
#' @param allow_for_mixture_models
#' @param dirichlet_prior
#' @param component_truncation
#' @param regime_truncation
#' @param hyperprior_b
#' @param model_params_save_every
#' @param simulation_iter
#' @param T2_window_size
#' @param determining_p_cutoff
#'
#' @return
#' @import parallel
#' @export
#'
#' @examples
fit_regime_vector = function(data_woTimeValues,
                             time_of_observations,
                             time_points,
                             not.cont = NULL,
                             iterations = 100,
                             burnin = floor(iterations / 2),
                             lower_bounds = NULL,
                             upper_bounds = NULL,
                             previous_states = NULL,
                             previous_model_fits = NULL,
                             sparsity_parameter = NULL,
                             linger_parameter = 10,
                             move_parameter = 1,
                             g.prior = 0.2,
                             set_G = NULL,
                             wishart_df_inital = 20,
                             lambda = 20,
                             g_sampling_distribution = NULL,
                             n.cores = 1,
                             scaleMatrix = NULL,
                             allow_for_mixture_models = FALSE,
                             dirichlet_prior = 0.004,
                             component_truncation = 7,
                             regime_truncation = 15,
                             hyperprior_b = 20,
                             model_params_save_every = 5,
                             simulation_iter = NULL,
                             T2_window_size = 3,
                             determining_p_cutoff = FALSE) {
  
  model_saves_list = list()
  
  # As a sanity check, component truncation should be set to 1 when we do not allow for mixture models:
  if (!allow_for_mixture_models) {
    component_truncation = 1
  }
  # | --------------------- Establish Initial Hyperparameters -------------------------- |
  print("inside the fit function")
  # Sort days and observations
  data_woTimeValues    = data_woTimeValues[order(time_of_observations),]
  time_of_observations = time_of_observations[order(time_of_observations)]
  time_points          = time_points[order(time_points)]
  data_set_list        = list()
  Z_timepoint_indices  = list(list(
    timepoint_first_index = NA,
    timepoint_last_index = NA
  ))
  negative_count       = 0
  for (i in 1:(length(time_points) - 1)) {
    if (i == 1) {
      if (length(which(time_of_observations <= time_points[i + 1 - negative_count])) < 100) {
        print(paste(
          "There are not enough observations for time point,",
          i - negative_count
        ))
        index_to_remove = i - negative_count + 1
        time_points     = time_points[-index_to_remove]
        negative_count  = negative_count + 1
        next
      }
      print(paste("Number of timepoints in day", 1, "is:", length(
        which(time_of_observations <= time_points[i + 1 - negative_count])
      )))
      
      Z_timepoint_indices[[i - negative_count]]$timepoint_first_index = 1
      Z_timepoint_indices[[i - negative_count]]$timepoint_last_index  = max(which(time_of_observations <= time_points[i +
                                                                                                                        1 - negative_count]))
    } else {
      indices_of_timepoint                           = which((time_of_observations > time_points[i -
                                                                                                   negative_count]) &
                                                               (time_of_observations <= time_points[i +
                                                                                                      1 - negative_count]))
      if (length(indices_of_timepoint) < 100) {
        print(paste(
          "There are not enough observations for time point,",
          i - negative_count
        ))
        index_to_remove = i - negative_count + 1
        time_points         = time_points[-index_to_remove]
        negative_count  = negative_count + 1
        next
      }
      print(paste(
        "Number of timepoints in day",
        i - negative_count,
        "is:",
        length(indices_of_timepoint)
      ))
      Z_timepoint_indices[[i - negative_count]]                       = list(timepoint_first_index = NA,
                                                                             timepoint_last_index = NA)
      Z_timepoint_indices[[i - negative_count]]$timepoint_first_index = min(indices_of_timepoint)
      Z_timepoint_indices[[i - negative_count]]$timepoint_last_index  = max(indices_of_timepoint)
    }
    data_set_list[[i - negative_count]] = data_woTimeValues[Z_timepoint_indices[[i -
                                                                                   negative_count]]$timepoint_first_index:Z_timepoint_indices[[i - negative_count]]$timepoint_last_index,]
  }
  
  
  if (is.null(previous_states)) {
    print("using T2 to create initial guess on state vector")
    time_since_last_change   = 1
    current_regime           = 1
    previous_states          = rep(1, times = (length(time_points) - 1))
    
    for (t in 3:(length(time_points) - 2)) {
      data_before = NULL
      data_after = NULL
      for (index in 1:(length(time_points) - 1)) {
        if ((index <= t) & ((index > (t - 3)))) {
          data_before = rbind(data_before, data_set_list[[index]])
        } else if (index == (t + 1)) {
          data_after = rbind(data_after, data_set_list[[index]])
        }
      }
      
      # Calculate the Hotelling T^2 Scan Statistic
      cov1  = empirical_cov_w_missing(data_before)
      cov2  = empirical_cov_w_missing(data_after)
      xbar1 = as.matrix(colMeans(data_before, na.rm = TRUE))
      xbar2 = as.matrix(colMeans(data_after, na.rm = TRUE))
      n1    = nrow(data_before)
      n2    = nrow(data_after)
      
      x_stuff = list(mean = xbar1,
                     cov = cov1,
                     n = n1)
      y_stuff = list(mean = xbar2,
                     cov = cov2,
                     n = n2)
      
      out = tryCatch({
        x       = hotelling.test(x_stuff, y_stuff, var.equal = FALSE)
        if ((det(x_stuff$cov) == 0) | (det(y_stuff$cov) == 0)) {
          print("one of the cov matrices is precisely singular.")
          p_value = 1
        } else {
          p_value = x$pval
        }
        p_value
      }, error = function(cond) {
        print("-----------------------------------")
        print("-----------------------------------")
        message(paste("Hotelling Calculation failed!"))
        print("-----------------------------------")
        print("-----------------------------------")
        return(1)
      })
      
      if (is.na(out)) {
        time_since_last_change     = time_since_last_change + 1
        previous_states[t + 1]       = current_regime
      } else if (out <= 0.05) {
        # Put in a change-point.
        time_since_last_change     = time_since_last_change + 1
        if (time_since_last_change > T2_window_size) {
          time_since_last_change   = 1
          current_regime           = current_regime + 1
        }
        previous_states[t + 1]       = current_regime
      } else {
        time_since_last_change     = time_since_last_change + 1
        previous_states[t + 1]       = current_regime
      }
    }
  }
  
  my_states = previous_states
  states_before_split1 = my_states[1:14]
  states_before_split2 = my_states[15:34]
  states_after_split2  = my_states[35:length(my_states)]
  cat(c(
    "[1]  -----> current states are:",
    paste(sapply(states_before_split1, paste, collapse = ' ')),
    "|",
    paste(sapply(states_before_split2, paste, collapse = ' ')),
    "|",
    paste(sapply(states_after_split2, paste, collapse = ' ')),
    "\n"
  ))
  bdgraph_object = BDgraph::bdgraph(
    data_woTimeValues[1:Z_timepoint_indices[[max(which(previous_states == 1))]]$timepoint_last_index,],
    g.prior = g.prior,
    iter = 500,
    method = "gcgm"
  )
  p              = ncol(data_woTimeValues)
  
  if (is.null(scaleMatrix)) {
    # I'm following the formulation of Lenkowski and Dobra 2011 for the scale matrix.  In contrast
    # to the Wishart Distribution from Wikipedia, the average of their distribution is df * solve(scaleMatrix).
    # Thus, to make the average of this Wishart distribution to be equal to the observed precision matrix, I'll
    # require the following scale matrix:
    sum_empirical_prec    = matrix(0, p, p)
    for (regime_index in 1:max(previous_states)) {
      # IMPORTANT NOTE: in the formulation of Lenkowski and Dobra, the average of the distribution is D^{-1} * (df.prior + p - 1).  Thus, when the
      #                 following scale matrix is inverted, and multiplied by (df.prior + p - 1), what should remain is the empirical average precision.
      min_index           = Z_timepoint_indices[[min(which(previous_states ==
                                                             regime_index))]]$timepoint_first_index
      max_index           = Z_timepoint_indices[[max(which(previous_states ==
                                                             regime_index))]]$timepoint_last_index
      temp_data           = data_woTimeValues[min_index:max_index,]
      
      temp_cov_matrix     = empirical_cov_w_missing(temp_data)
      if (any(is.null(scaleMatrix) |
              is.na(scaleMatrix)) | (det(temp_cov_matrix) == 0)) {
        print("NOTE: We were not able to get a good initial guess for the precision matrix.")
        precMatrix_temp   = bdgraph_object$K_hat
      } else {
        precMatrix_temp     = solve(temp_cov_matrix)
      }
      sum_empirical_prec  = sum_empirical_prec + precMatrix_temp
    }
    average_prec_estimate = sum_empirical_prec / max(previous_states)
    covMatrix             = solve(average_prec_estimate)
  }
  print("Empirical Covariance Estimate:")
  print(covMatrix)
  hyperparameters        = list(
    hyperprior_b = hyperprior_b,
    alpha = linger_parameter,
    beta = move_parameter,
    mu_0 = NA,
    lambda = lambda,
    p = NA,
    hyperprior_scale_matrix = 0.2 * diag(p),
    wishart_scale_matrix = covMatrix * (wishart_df_inital +
                                          p - 1),
    wishart_df = wishart_df_inital,
    log_wishart_prior_term = NA,
    dirichlet_prior = dirichlet_prior,
    component_truncation = component_truncation,
    regime_truncation = regime_truncation,
    g.prior = g.prior
  ) #, pseudoprior_df=pseudoprior_df)
  # | ------------------------ Fill in Missing Inputs ---------------------------------- |
  data_woTimeValues[is.infinite(data_woTimeValues)] = NA
  model_states                                      = list()
  transition_probabilities                          = rep(NA, times = nrow(data_woTimeValues))
  
  # Here I'm assuming that the first values in time_of_observations is the earliest
  # available time, and the last is the last available time.
  iter = iterations
  if (is.null(lower_bounds)) {
    min_value    = min(data_woTimeValues, na.rm = TRUE) - 1e100
    lower_bounds = rep(min_value, times = ncol(data_woTimeValues))
  }
  if (is.null(upper_bounds)) {
    max_value    = max(data_woTimeValues, na.rm = TRUE) + 1e100
    upper_bounds = rep(max_value, times = ncol(data_woTimeValues))
  }
  if (is.null(not.cont)) {
    not.cont     = rep(0, times = ncol(data_woTimeValues))
  }
  
  if (length(not.cont) < ncol(data_woTimeValues)) {
    not.cont = c(not.cont, rep(1, times = (
      ncol(data_woTimeValues) - length(not.cont)
    )))
  }
  
  is_continuous            = !not.cont
  hyperparameters$not.cont = not.cont
  
  if (sum(apply(data_woTimeValues, 2, class) == "numeric") != ncol(data_woTimeValues)) {
    stop("Data must be all numeric.  Makes sure that you've seperated the POSIXct objects.")
  }
  if (!all((
    class(time_of_observations) %in% c("POSIXct", "POSIXlt" , "POSIXt")
  ))) {
    stop("time_of_observations must contain only POSIXt objects.")
  }
  
  # | --------------------------- Gathering Data Values -------------------------------------- |
  if (wishart_df_inital < 3)
    stop(" 'hyperprior.df' must be >= 3.")
  if (iter < burnin)
    stop(" Number of iteration must be more than number of burn-in.")
  burnin = floor(burnin)
  
  # Gather indicator information used for the latent variable draws
  n                              = nrow(data_woTimeValues)
  p                              = ncol(data_woTimeValues)
  hyperparameters$p              = p
  lower_bound_is_equal           = data_woTimeValues
  upper_bound_is_equal           = data_woTimeValues
  is_missing                     = data_woTimeValues
  mean_wo_missing                = function(x) {
    mean(x, na.rm = T)
  }
  mu_0                           = rep(0, times = p)
  
  # If I don't have the sampling distribution over the graph structure set, start it with a
  # uniform distribution.
  print("establishing the starting graph structure.")
  if (is.null(set_G)) {
    # If G is not set, I'll draw from the distribution over all possible graphs determined by g.prior,
    # the probability of a given edge being present.
    previous_G              = BDgraph::select(bdgraph_object)
    
  } else {
    previous_G    = set_G
  }
  hyperparameters$G = previous_G
  
  if (is.null(g_sampling_distribution)) {
    g_sampling_distribution = bdgraph_object$p_links
    g_sampling_distribution = (1 - g_sampling_distribution) * (previous_G) + g_sampling_distribution *
      (1 - previous_G)
  }
  
  
  for (index in 1:p) {
    lower_bound_is_equal[, index] = (data_woTimeValues[, index] == lower_bounds[index])
    upper_bound_is_equal[, index] = (data_woTimeValues[, index] == upper_bounds[index])
    is_missing[, index]           = is.na(data_woTimeValues[, index])
    if (sum(is.na(data_woTimeValues[, index])) == nrow(data_woTimeValues)) {
      stop("Drop the data column with all missing values!")
    }
    mu_0[index]                  = mean(data_woTimeValues[, index], na.rm = TRUE)
  }
  lower_bound_is_equal[is.na(lower_bound_is_equal)] = 0
  upper_bound_is_equal[is.na(upper_bound_is_equal)] = 0
  
  hyperparameters$mu_0                   = mu_0
  hyperparameters$is_missing             = is_missing
  hyperparameters$lower_bound_is_equal   = lower_bound_is_equal
  hyperparameters$upper_bound_is_equal   = upper_bound_is_equal
  hyperparameters$lower_bounds           = lower_bounds
  hyperparameters$upper_bounds           = upper_bounds
  my_states                              = previous_states
  
  # | --------------------------- Initial Model Fits -------------------------------------- |
  transition_probabilities_temp = redraw_transition_probs(previous_states,
                                                          linger_parameter,
                                                          move_parameter,
                                                          n.cores - 1)
  transition_probabilities      = c(transition_probabilities_temp,
                                    rep(NA, times = (
                                      length(my_states)  - length(transition_probabilities_temp)
                                    )))
  print("starting to fit the models for each state")
  
  print("data has been subset -- setting up parallelization")
  n.cores_eachparent = max(min(regime_truncation, floor(n.cores / 2)), 1)
  n.cores_eachchild  = max(floor((n.cores - n.cores_eachparent) / n.cores_eachparent), 1)
  print(
    paste(
      "Setting up",
      n.cores_eachparent,
      "children, each can fork",
      n.cores_eachchild,
      "additional children."
    )
  )
  clust              = parallel::makeCluster(n.cores_eachparent) #, port = port_number, outfile = 'childrencheck.txt'
  if (is.null(previous_model_fits)) {
    # If we don't have model fits from a previous run, we need to fit models
    # based off of the states.
    df_prior_on_wish_start = wishart_df_inital
    previous_model_fits    = lapply(1:(regime_truncation + 1),
                                    function(x) {
                                      list(
                                        precision = list(),
                                        mu = list(),
                                        cluster_assignments = NA,
                                        component_log_probs = NA,
                                        component_sticks = NA
                                      )
                                    })
    parallel::clusterExport(clust, "previous_G", envir = environment())
    parallel::clusterExport(clust, "my_states", envir = environment())
    parallel::clusterExport(clust, "hyperparameters", envir = environment())
    parallel::clusterExport(clust, "Z_timepoint_indices", envir = environment())
    parallel::clusterExport(clust, "data_woTimeValues", envir = environment())
    parallel::clusterExport(clust, "regime_truncation", envir = environment())
    parallel::clusterExport(clust, "component_truncation", envir = environment())
    parallel::clusterExport(clust, "wishart_df_inital", envir = environment())
    parallel::clusterExport(clust, "scaleMatrix", envir = environment())
    models_temp = parallel::parSapply(clust, 1:(regime_truncation + 1), function(x) {
      # source("R/helpers.R")
      linger_parameter         = hyperparameters$alpha
      move_parameter           = hyperparameters$beta
      not.cont                 = hyperparameters$not.cont
      result                   = rgwish_Rcpp(
        as.double(previous_G),
        as.double(hyperparameters$hyperprior_scale_matrix),
        as.integer(hyperparameters$hyperprior_b),
        as.integer(hyperparameters$p),
        as.double(1e-8)
      )
      new_scale                = matrix(result, hyperparameters$p, hyperparameters$p)
      
      previous_model_fits      = lapply(1:(regime_truncation + 1),
                                        function(x) {
                                          list(
                                            precision = list(),
                                            mu = list(),
                                            cluster_assignments = NA,
                                            component_log_probs = NA,
                                            component_sticks = NA
                                          )
                                        })
      
      if (sum(my_states == x) > 0) {
        first_index = Z_timepoint_indices[[min(which(my_states == x))]]$timepoint_first_index
        last_index  = Z_timepoint_indices[[max(which(my_states == x))]]$timepoint_last_index
        indicies    = c(first_index:last_index)
        temp_data   = data_woTimeValues[first_index:last_index, ]
      }
      # Start by drawing a prior for everything.
      previous_model_fits = redraw_mixture_parameters(
        my_states,
        x,
        previous_model_fits,
        data_woTimeValues,
        Z_timepoint_indices,
        linger_parameter,
        move_parameter,
        not.cont,
        previous_G,
        hyperparameters
      )
      if (sum(my_states == x) > 0) {
        previous_model_fits[[x]]$cluster_assignments = rep(1, times = nrow(temp_data))
      }
      list(
        precision = previous_model_fits[[x]]$precision,
        mu = previous_model_fits[[x]]$mu,
        cluster_assignments = previous_model_fits[[x]]$cluster_assignments
      )
    })
    print("initial model fits complete.")
    
    # Now we fill in the results from our setup into previous_model_fits,
    # and do our data subset now to save on time:
    for (i in 1:regime_truncation) {
      previous_model_fits[[i]][["precision"]]           = models_temp[["precision", i]]
      previous_model_fits[[i]][["mu"]]                  = models_temp[["mu", i]]
      previous_model_fits[[i]][["cluster_assignments"]] = models_temp[["cluster_assignments", i]]
      
      # Here, I am going to fill in each of the component assignment 'sticks' from the beta (dirichlet process) prior.
      temp_component_sticks = rep(0, times = component_truncation)
      if (anyNA(models_temp[["cluster_assignments", i]])) {
        # This means that there are no data assigned to this regime.
        previous_model_fits[[i]][["component_log_probs"]]     = rbeta(component_truncation,
                                                                      1,
                                                                      hyperparameters$dirichlet_prior)
      } else {
        for (stick_index in 1:component_truncation) {
          if (stick_index %in% models_temp[["cluster_assignments", i]]) {
            temp_component_sticks[stick_index] = rbeta(
              1,
              1 + sum(models_temp[["cluster_assignments", i]] == stick_index),
              hyperparameters$dirichlet_prior +
                sum(models_temp[["cluster_assignments", i]] > stick_index)
            )
          } else {
            temp_component_sticks[stick_index] = rbeta(1, 1, hyperparameters$dirichlet_prior)
          }
        }
      }
      previous_model_fits[[i]][["component_log_probs"]] = sticks_to_log_probs(temp_component_sticks)
      previous_model_fits[[i]][["component_sticks"]]    = temp_component_sticks
    }
    
  }
  
  full_data_Z              = data_woTimeValues
  hyperparameters$raw_data = data_woTimeValues
  hyperparameters$G        = previous_G
  
  # | -------------------------------Starting MCMC-------------------------------------- |
  # If the models were not built, now they have been.  We can begin the mcmc
  # updates.
  model_states                    = list()
  transition_probabilities_states = list()
  mergesplit_accepts              = c()
  DRJ_accepts                     = c()
  gibbswap_accepts                = c()
  alphabeta_accepts               = c()
  cholesky_failed_prec            = c()
  G_Wish_portion_trace            = c()
  D_prior_dens_trace              = c()
  my_flag = 0
  
  mergesplit_accepted   = 0
  DRJ_accepted          = 0
  gibbswap_accepted     = 0
  alphabeta_accepted    = 0
  record_count          = 1
  
  print("Beginning MCMC Chain")
  print(
    paste(
      "-----> alpha:",
      hyperparameters$alpha,
      ", beta:",
      hyperparameters$beta,
      ", lambda:",
      hyperparameters$lambda
    )
  )
  for (iter in 1:iterations) {
    print(paste("Starting Iteration", iter))
    # | --------------------------- Redraw the Latent Data ------------------------------------- |
    print("-----> setting up first-level parallelization to draw from latent data")
    print(paste("-----> Setting up", n.cores, "cores."))
    parallel::clusterExport(clust, "my_states", envir = environment())
    parallel::clusterExport(clust, "full_data_Z", envir = environment())
    parallel::clusterExport(clust, "data_woTimeValues", envir = environment())
    parallel::clusterExport(clust, "Z_timepoint_indices", envir = environment())
    parallel::clusterExport(clust, "previous_model_fits", envir = environment())
    parallel::clusterExport(clust, "lower_bounds", envir = environment())
    parallel::clusterExport(clust, "upper_bounds", envir = environment())
    parallel::clusterExport(clust, "lower_bound_is_equal", envir = environment())
    parallel::clusterExport(clust, "upper_bound_is_equal", envir = environment())
    parallel::clusterExport(clust, "is_missing", envir = environment())
    parallel::clusterExport(clust, "not.cont", envir = environment())
    data_subsets_z = parallel::parLapply(clust, 1:max(my_states), function(x) {
      library(palliative.changepoints)
      # source("R/helpers.R")
      print(x)
      first_index               = Z_timepoint_indices[[min(which(my_states == x))]]$timepoint_first_index
      last_index                = Z_timepoint_indices[[max(which(my_states == x))]]$timepoint_last_index
      temp_data                 = full_data_Z[first_index:last_index,]
      temp_raw_data             = data_woTimeValues[first_index:last_index,]
      is_missing_temp           = is_missing[first_index:last_index,]
      upper_bound_is_equal_temp = upper_bound_is_equal[first_index:last_index,]
      lower_bound_is_equal_temp = lower_bound_is_equal[first_index:last_index,]
      new_latent_data           = redraw_latent_data(
        x,
        previous_model_fits,
        hyperparameters,
        temp_data,
        temp_raw_data,
        is_missing_temp,
        upper_bound_is_equal_temp,
        lower_bound_is_equal_temp
      )
      if (anyNA(new_latent_data)) {
        stop("latent data redraw didn't work on C level -- NAs persist.")
      }
      return(new_latent_data)
    })
    
    full_data_Z = NULL
    for (i in 1:max(my_states)) {
      full_data_Z = rbind(full_data_Z, data_subsets_z[[i]])
    }
    print(paste("-----> Redraw latent data complete."))
    
    # | ------------------------------ Update States: Merge-Split --------------------------------- |
    print(paste("-----> Redraw hyperparameters before Merge-Split."))
    new_hyperparameters                         = redraw_hyperparameters(hyperparameters,
                                                                         transition_probabilities,
                                                                         previous_model_fits,
                                                                         my_states)
    hyperparameters                             = new_hyperparameters$hyperparameters_item
    previous_model_fits                         = new_hyperparameters$previous_model_fits_item
    symmetric_scale                             = hyperparameters$wishart_scale_matrix
    symmetric_scale[lower.tri(symmetric_scale)] = 0
    symmetric_scale                             = symmetric_scale + t(symmetric_scale) - diag(diag(symmetric_scale))
    D_prior_dens_trace                          = c(
      D_prior_dens_trace,
      CholWishart::dWishart(
        symmetric_scale,
        hyperparameters$wishart_df,
        diag(hyperparameters$p),
        log = TRUE
      )
    )
    print(paste("-----> Hyperparameters redraw completed."))
    
    if (iter < floor(iterations / 5)) {
      launching = TRUE
    } else {
      launching = FALSE
    }
    
    print(paste("-----> Starting the Merge-Split Algorithm."))
    states_update            = update_states_mergesplit(
      my_states,
      previous_model_fits,
      full_data_Z,
      Z_timepoint_indices,
      linger_parameter,
      move_parameter,
      not.cont,
      transition_probabilities,
      previous_G,
      hyperparameters,
      n.cores,
      launching,
      allow_mixtures = TRUE,
      min_regime_length = 1
    )
    print(paste("-----> Merge-Split Algorithm complete."))
    previous_model_fits      = states_update$previous_model_fits_item
    my_states                = states_update$my_states_item
    transition_probabilities = states_update$transition_probabilities_item
    mergesplit_accepted      = states_update$accepted_item
    full_data_Z              = states_update$full_data_Z_item
    G_Wish_portion_trace     = c(G_Wish_portion_trace,
                                 states_update$G_Wish_portion_of_MH)
    
    # | --------------------------- Update States: Gibbs Sampler ------------------------------ |
    print(paste("-----> Starting the Gibbs Sweep."))
    
    states_update            = update_states_gibbs(
      my_states,
      full_data_Z,
      Z_timepoint_indices,
      previous_model_fits,
      transition_probabilities,
      hyperparameters,
      n.cores - 1,
      min_regime_length = 2
    )
    
    previous_model_fits      = states_update$previous_model_fits_item
    my_states                = states_update$my_states_item
    transition_probabilities = redraw_transition_probs(my_states,
                                                       hyperparameters$alpha,
                                                       hyperparameters$beta,
                                                       n.cores - 1)
    gibbswap_accepted        = states_update$accepted_item
    print(paste("-----> Gibbs Sweep complete."))
    
    my_states_all_equal = unique(my_states)
    values_between      = 1:max(my_states)
    if (all.equal(my_states_all_equal, values_between) != TRUE) {
      print(my_states)
      stop("Gibbs update caused an error with state relabeling.")
    }
    
    parallel::clusterExport(clust, "full_data_Z", envir = environment())
    parallel::clusterExport(clust, "n.cores_eachchild", envir = environment())
    parallel::clusterExport(clust, "Z_timepoint_indices", envir = environment())
    parallel::clusterExport(clust, "my_states", envir = environment())
    parallel::clusterExport(clust, "previous_G", envir = environment())
    parallel::clusterExport(clust, "previous_model_fits", envir = environment())
    parallel::clusterExport(clust, "transition_probabilities", envir = environment())
    parallel::clusterExport(clust, "hyperparameters", envir = environment())
    parallel::clusterExport(clust, "p", envir = environment())
    models_temp = parallel::parSapply(clust, 1:regime_truncation, function(x) {
      #source("R/helpers.R")
      linger_parameter          = hyperparameters$alpha
      move_parameter            = hyperparameters$beta
      not.cont                  = hyperparameters$not.cont
      
      
      new_parameter_values_temp = redraw_mixture_parameters(
        my_states,
        x,
        previous_model_fits,
        full_data_Z,
        Z_timepoint_indices,
        linger_parameter,
        move_parameter,
        not.cont,
        previous_G,
        hyperparameters
      )
      
      list(new_lambda_item = new_parameter_values_temp[[x]]$precision,
           new_mu_item     = new_parameter_values_temp[[x]]$mu)
    })
    
    print("-----> Finished Model draw post merge-split.  Accepted always.")
    
    if (length(models_temp[[1]][1]) != 0) {
      for (index in 1:max(my_states)) {
        previous_model_fits[[index]]$precision = models_temp[["new_lambda_item", index]]
        previous_model_fits[[index]]$mu        = models_temp[["new_mu_item", index]]
        previous_model_fits                    = shift_components(index, previous_model_fits)
      }
    } else {
      print("-----> Conjugate update failed!!")
    }
    
    print("Finished model draw after gibbs sweep.")
    # | ------------------ Merge-Split-Gibbs the Mixture Components ----------------------- |
    # I need to make sure that the component probabilities are well-behaved here,
    # so I'll update the component probabilities using a conjugate update.
    new_hyperparameters                    = redraw_hyperparameters(hyperparameters,
                                                                    transition_probabilities,
                                                                    previous_model_fits,
                                                                    my_states)
    hyperparameters                        = new_hyperparameters$hyperparameters_item
    previous_model_fits                    = new_hyperparameters$previous_model_fits_item
    
    symmetric_scale                             = hyperparameters$wishart_scale_matrix
    symmetric_scale[lower.tri(symmetric_scale)] = 0
    symmetric_scale                             = symmetric_scale + t(symmetric_scale) - diag(diag(symmetric_scale))
    D_prior_dens_trace                          = c(
      D_prior_dens_trace,
      CholWishart::dWishart(
        symmetric_scale,
        hyperparameters$wishart_df,
        diag(hyperparameters$p),
        log = TRUE
      )
    )
    
    if (((allow_for_mixture_models) &
         (iter %% 5 == 1)) |
        ((allow_for_mixture_models) &
         (iter < floor(iterations / 4)))) {
      print("about to merge-split-gibbs the mixture components")
      parallel::clusterExport(clust, "my_states", envir = environment())
      parallel::clusterExport(clust, "full_data_Z", envir = environment())
      parallel::clusterExport(clust, "Z_timepoint_indices", envir = environment())
      parallel::clusterExport(clust, "previous_model_fits", envir = environment())
      parallel::clusterExport(clust, "lower_bounds", envir = environment())
      parallel::clusterExport(clust, "upper_bounds", envir = environment())
      parallel::clusterExport(clust, "lower_bound_is_equal", envir = environment())
      parallel::clusterExport(clust, "upper_bound_is_equal", envir = environment())
      parallel::clusterExport(clust, "is_missing", envir = environment())
      parallel::clusterExport(clust, "not.cont", envir = environment())
      parallel::clusterExport(clust, "previous_G", envir = environment())
      new_mixture_components = parallel::parLapply(clust, 1:max(my_states), function(x) {
        print(paste("starting", x))
        #source("R/helpers.R")
        print("redrawing the components for state")
        first_index               = Z_timepoint_indices[[min(which(my_states == x))]]$timepoint_first_index
        last_index                = Z_timepoint_indices[[max(which(my_states == x))]]$timepoint_last_index
        temp_data                 = full_data_Z[first_index:last_index,]
        
        new_components_for_state  = splitmerge_gibbs_comps(
          my_states,
          previous_model_fits,
          temp_data,
          x,
          hyperparameters,
          Z_timepoint_indices,
          full_data_Z
        )
        
        list(
          precisions        = new_components_for_state$precisions,
          mus               = new_components_for_state$mus,
          assigns           = new_components_for_state$assigns,
          component_probs   = new_components_for_state$comp_probs,
          component_sticks  = new_components_for_state$comp_sticks
        )
      })
      for (i in 1:max(my_states)) {
        previous_model_fits[[i]]$precision           = new_mixture_components[[i]]$precisions
        previous_model_fits[[i]]$mu                  = new_mixture_components[[i]]$mus
        previous_model_fits[[i]]$cluster_assignments = new_mixture_components[[i]]$assigns
        previous_model_fits[[i]]$component_log_probs = new_mixture_components[[i]]$component_probs
        previous_model_fits[[i]]$component_sticks    = new_mixture_components[[i]]$component_sticks
        previous_model_fits                          = shift_components(i, previous_model_fits)
        print(paste("cluster assignments for regime", i, "are:"))
        print(previous_model_fits[[i]]$cluster_assignments)
        if (min(previous_model_fits[[i]]$cluster_assignments) > 1) {
          stop("Component shift failed!!")
        }
      }
    } else if (allow_for_mixture_models) {
      # This is the case where (in between the regular split merge steps that i do every 5 iterations) I
      # just do a quick gibbs swap.
      parallel::clusterExport(clust, "my_states", envir = environment())
      parallel::clusterExport(clust, "full_data_Z", envir = environment())
      parallel::clusterExport(clust, "Z_timepoint_indices", envir = environment())
      parallel::clusterExport(clust, "previous_model_fits", envir = environment())
      parallel::clusterExport(clust, "lower_bounds", envir = environment())
      parallel::clusterExport(clust, "upper_bounds", envir = environment())
      parallel::clusterExport(clust, "lower_bound_is_equal", envir = environment())
      parallel::clusterExport(clust, "upper_bound_is_equal", envir = environment())
      parallel::clusterExport(clust, "is_missing", envir = environment())
      parallel::clusterExport(clust, "not.cont", envir = environment())
      parallel::clusterExport(clust, "previous_G", envir = environment())
      
      new_mixture_components = parallel::parLapply(clust, 1:max(my_states), function(x) {
        # print(paste("starting", x))
        #source("R/helpers.R")
        # print("redrawing the components for state")
        first_index               = Z_timepoint_indices[[min(which(my_states == x))]]$timepoint_first_index
        last_index                = Z_timepoint_indices[[max(which(my_states == x))]]$timepoint_last_index
        temp_data                 = full_data_Z[first_index:last_index,]
        new_components_for_state  = gibbs_swap_comps(
          temp_data,
          previous_model_fits[[x]]$cluster_assignments,
          previous_model_fits[[x]]$component_log_probs,
          previous_model_fits[[x]]$precision,
          previous_model_fits[[x]]$mu,
          hyperparameters$component_truncation,
          2
        )
        list(assigns           = as.vector(new_components_for_state))
      })
      for (i in 1:max(my_states)) {
        previous_model_fits[[i]]$cluster_assignments = new_mixture_components[[i]]$assigns
        previous_model_fits                          = shift_components(i, previous_model_fits)
        # print(paste("cluster assignments for regime", i, "are:"))
        # print(previous_model_fits[[i]]$cluster_assignments)
      }
    }
    
    # | -------------------------------Double Reversible Jump Graph Resampling---------------------------------------- |
    # First, resample everything based on the posteriors.
    print(paste("-----> Starting to refit the model parameters."))
    
    # Periodically when we are still in burnin, we will update the graph G sampling distribution using the BDgraph package from CRAN.
    if ((iter < burnin) && ((iter) %% 100 == 0)) {
      print(paste("-----> Redrawing the BDgraph object."))
      out <- tryCatch({
        BDgraph::bdgraph(full_data_Z, g.prior = g.prior, iter = 1000)
      },
      error = function(cond) {
        message("BDgraph failed, like due to computationally singular Ds")
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
      })
      if (!is.na(out)) {
        bdgraph_object          = out
        g_sampling_distribution = bdgraph_object$p_links
        
        g_sampling_distribution = (1 - g_sampling_distribution) * (previous_G) + g_sampling_distribution *
          (1 - previous_G)
      }
      
    }
    
    # Randomly select an edge based on the sampling distribution on the graph structure.
    total_prob_mass = 0
    rand_selection  = runif(1,
                            min = 0,
                            max = sum(g_sampling_distribution))
    if (is.na(sum(g_sampling_distribution))) {
      print("Something is wrong with the g sampling distribution!")
      print(bdgraph_object$p_links)
      stop("Check the BDgraph package access.")
    }
    done = FALSE
    for (g_index_i in 1:(p - 1)) {
      for (g_index_j in (g_index_i + 1):p) {
        total_prob_mass = total_prob_mass + g_sampling_distribution[g_index_i, g_index_j]
        if (total_prob_mass >= rand_selection) {
          selected_edge_i = g_index_i
          selected_edge_j = g_index_j
          done = TRUE
          break
        }
      }
      if (done) {
        break
      }
    }
    if (!done) {
      selected_edge_i = 2
      selected_edge_j = 3
    }
    new_G_proposed  = previous_G
    if (previous_G[selected_edge_i, selected_edge_j]) {
      new_G_proposed[selected_edge_i, selected_edge_j] = 0
      new_G_proposed[selected_edge_j, selected_edge_i] = 0
    } else {
      new_G_proposed[selected_edge_i, selected_edge_j] = 1
      new_G_proposed[selected_edge_j, selected_edge_i] = 1
    }
    
    # Next, we fix this graph G and sample new precision matrices K1,...,Km for each regime.
    n.cores_eachparent = max(min(max(my_states), n.cores - 1), 1)
    n.cores_eachchild  = max(floor((n.cores - n.cores_eachparent) / n.cores_eachparent), 1)
    print("-----> Setting up to draw new proposal precision matrices")
    
    ## First, we'll just draw new parameter values for every state, which we will accept deterministically by conjugate update.
    parallel::clusterExport(clust, "full_data_Z", envir = environment())
    parallel::clusterExport(clust, "n.cores_eachchild", envir = environment())
    parallel::clusterExport(clust, "Z_timepoint_indices", envir = environment())
    parallel::clusterExport(clust, "my_states", envir = environment())
    parallel::clusterExport(clust, "previous_G", envir = environment())
    parallel::clusterExport(clust, "previous_model_fits", envir = environment())
    parallel::clusterExport(clust, "transition_probabilities", envir = environment())
    parallel::clusterExport(clust, "hyperparameters", envir = environment())
    parallel::clusterExport(clust, "p", envir = environment())
    models_temp = parallel::parSapply(clust, 1:regime_truncation, function(x) {
      #source("R/helpers.R")
      linger_parameter          = hyperparameters$alpha
      move_parameter            = hyperparameters$beta
      not.cont                  = hyperparameters$not.cont
      
      
      new_parameter_values_temp = redraw_mixture_parameters(
        my_states,
        x,
        previous_model_fits,
        full_data_Z,
        Z_timepoint_indices,
        linger_parameter,
        move_parameter,
        not.cont,
        previous_G,
        hyperparameters
      )
      
      list(new_lambda_item = new_parameter_values_temp[[x]]$precision,
           new_mu_item     = new_parameter_values_temp[[x]]$mu)
    })
    
    print("-----> Model pre-draw complete.  Accepted always.")
    if (length(models_temp[[1]][1]) != 0) {
      for (index in 1:max(my_states)) {
        previous_model_fits[[index]]$precision = models_temp[["new_lambda_item", index]]
        previous_model_fits[[index]]$mu        = models_temp[["new_mu_item", index]]
      }
    } else {
      print("-----> Conjugate update failed!!")
    }
    ## Next, we update our graph structure G using a Double Reversible Jump
    models_temp = NULL
    tryCatch({
      parallel::clusterExport(clust, "previous_model_fits", envir = environment())
      parallel::clusterExport(clust, "full_data_Z", envir = environment())
      parallel::clusterExport(clust, "n.cores_eachchild", envir = environment())
      parallel::clusterExport(clust, "Z_timepoint_indices", envir = environment())
      parallel::clusterExport(clust, "my_states", envir = environment())
      parallel::clusterExport(clust, "new_G_proposed", envir = environment())
      parallel::clusterExport(clust, "previous_G", envir = environment())
      parallel::clusterExport(clust, "selected_edge_i", envir = environment())
      parallel::clusterExport(clust, "selected_edge_j", envir = environment())
      parallel::clusterExport(clust, "hyperparameters", envir = environment())
      parallel::clusterExport(clust, "g.prior", envir = environment())
      parallel::clusterExport(clust, "p", envir = environment())
      models_temp = parallel::parSapply(clust, 1:max(my_states), function(x) {
        library(palliative.changepoints)
        
        state_to_redraw = x
        current_graph_G = previous_G
        redraw_G_output = redraw_G_with_mixture(
          my_states,
          state_to_redraw,
          previous_model_fits,
          full_data_Z,
          Z_timepoint_indices,
          current_graph_G,
          hyperparameters,
          new_G_proposed,
          selected_edge_i,
          selected_edge_j,
          g_sampling_distribution
        )
        
        list(
          new_previous_model_fits_item = redraw_G_output$previous_model_fits_item,
          MH_value_item                = redraw_G_output$log_MH_ratio
        )
      })
    },
    error = function(cond) {
      print("-----> ERROR: There was an error in the model refitting:")
      print("")
      message(cond)
      print("")
      print("")
      models_temp = NA
    })
    # We now evaluate the metropolis-hastings ratio.
    if (length(models_temp[[1]][1]) != 0) {
      cholesky_failed_prec     = c(cholesky_failed_prec, 0)
      log_MHratio_all_models   = 0
      for (index in 1:max(my_states)) {
        log_MHratio_all_models = log_MHratio_all_models + models_temp[["MH_value_item", index]]
      }
      log_unif                             = log(runif(1))
      
      if (is.na(log_MHratio_all_models)) {
        print("-----> DRJ model fits complete (did not accept).  Missing values!")
        DRJ_accepted                             = 0
      } else if ((log_unif <= log_MHratio_all_models) &
                 (!is.infinite(log_MHratio_all_models))) {
        print("-----> DRJ model fits complete (accepted).")
        # Here, we ACCEPT
        for (index in 1:max(my_states)) {
          temp_previous_model_fits               = models_temp[["new_previous_model_fits_item", index]]
          previous_model_fits[[index]]           = temp_previous_model_fits[[index]]
        }
        previous_G                               = new_G_proposed
        DRJ_accepted                             = 1
        print(paste(
          "-----> DRJ model fits DRJ-MH ratio:",
          log_MHratio_all_models
        ))
        
      } else {
        # Note that I'm rejecting all infinite MH ratio values.
        print("-----> DRJ model fits complete (did not accept).")
        DRJ_accepted                             = 0
        print(paste(
          "-----> DRJ model fits DRJ-MH ratio:",
          log_MHratio_all_models
        ))
      }
      
    } else {
      print("-----> Cholesky Failed!  DRJ model fits incomplete. (did not accept).")
      cholesky_failed_prec     = c(cholesky_failed_prec, 0)
      log_MHratio_all_models   = NA
    }
    
    # | ----------------------------- Update Hyperparameters ---------------------------------------- |
    hyperparameters$G                      = previous_G
    # Lastly, we would like to learn the move and linger parameters.
    new_hyperparameters                    = redraw_hyperparameters(hyperparameters,
                                                                    transition_probabilities,
                                                                    previous_model_fits,
                                                                    my_states)
    hyperparameters                        = new_hyperparameters$hyperparameters_item
    previous_model_fits                    = new_hyperparameters$previous_model_fits_item
    
    symmetric_scale                             = hyperparameters$wishart_scale_matrix
    symmetric_scale[lower.tri(symmetric_scale)] = 0
    symmetric_scale                             = symmetric_scale + t(symmetric_scale) - diag(diag(symmetric_scale))
    D_prior_dens_trace                          = c(
      D_prior_dens_trace,
      CholWishart::dWishart(
        symmetric_scale,
        hyperparameters$wishart_df,
        diag(hyperparameters$p),
        log = TRUE
      )
    )
    alphabeta_accepted                     = new_hyperparameters$accepted_item
    
    print(
      paste(
        "-----> alpha:",
        hyperparameters$alpha,
        ", beta:",
        hyperparameters$beta,
        ", lambda:",
        hyperparameters$lambda
      )
    )
    cat(c("[1]  -----> current mu_0 is:", hyperparameters$mu_0), "\n")
    if (length(my_states) > 30) {
      states_before_split1 = my_states[1:14]
      states_before_split2 = my_states[15:34]
      states_after_split2  = my_states[35:length(my_states)]
      cat(c(
        "[1]  -----> current states are:",
        paste(sapply(
          states_before_split1, paste, collapse = ' '
        )),
        "|",
        paste(sapply(
          states_before_split2, paste, collapse = ' '
        )),
        "|",
        paste(sapply(
          states_after_split2, paste, collapse = ' '
        )),
        "\n"
      ))
    } else {
      states_before_split = my_states[1:8]
      states_after_split  = my_states[9:length(my_states)]
      cat(c(
        "[1]  -----> current states are:",
        paste(sapply(
          states_before_split, paste, collapse = ' '
        )),
        "|",
        paste(sapply(
          states_after_split, paste, collapse = ' '
        )),
        "\n"
      ))
    }
    # | ----------------------------- Record MCMC Information ---------------------------------------- |
    mergesplit_accepts                      = c(mergesplit_accepts, mergesplit_accepted)
    print("MERGSPLIT ACCEPTS")
    print(mergesplit_accepts)
    DRJ_accepts                             = c(DRJ_accepts, DRJ_accepted)
    print("NUMBER OF DRJ ACCEPTS IS:")
    print(sum(DRJ_accepts))
    gibbswap_accepts                        = c(gibbswap_accepts, gibbswap_accepted)
    alphabeta_accepts                       = c(alphabeta_accepts, alphabeta_accepted)
    transition_probabilities_states[[iter]] = transition_probabilities
    
    print(paste("saving data at location", iter))
    if ((iter %% model_params_save_every == 0) & (iter > burnin)) {
      # I want to save model parameters associated with different state vectors for later analysis.
      # I am only going to save the model parameters for the final two regimes.
      if (length(unique(my_states)) == 1) {
        temp_list = previous_model_fits
      } else {
        temp_list = previous_model_fits[(max(my_states) - 1):max(my_states)]
      }
      if (is.null(simulation_iter)) {
        data_name = paste("Model_Fits", "/",
                          paste(my_states, collapse = ""),
                          "_",
                          iter,
                          '.RData',
                          sep = "")
      } else {
        # I think that the only information that I really want is the location of the LAST changepoint.
        # I want every model that puts a changepoint here so that I can see if they are all putting one
        # there for the same reason. I don't need an exact match.
        temp_states = rep(0, times = (length(my_states) - 1))
        for (j in 1:(length(my_states) - 1)) {
          temp_states[j] = my_states[j + 1] - my_states[j]
        }
        if (sum(temp_states) > 1) {
          indicies_of_changepoint                                               = which(temp_states ==
                                                                                          1)
          temp_states                                                           = rep(0, times = (length(my_states) -
                                                                                                    1))
          temp_states[indicies_of_changepoint[length(indicies_of_changepoint)]] = 1
        }
        data_name = paste(
          "Model_Fits_Simulation/",
          paste(temp_states, collapse = ""),
          "_",
          simulation_iter,
          "_",
          iter,
          '.RData',
          sep = ""
        )
      }
      
      # Usually, I only need to save the model information for the last two regimes (fault detection)
      # However, when I am searching for a reasonable probability cutoff on the posterior, I need them all.
      model_saves_list$data_name = temp_list
      # saveRDS(temp_list, file = data_name)
    }
    
    if (determining_p_cutoff) {
      data_name = paste(
        "Cuttoff_Model_Logs/",
        paste(my_states, collapse = "n"),
        "_",
        simulation_iter,
        "_",
        iter,
        '.RData',
        sep = ""
      )
      saveRDS(previous_model_fits, file = data_name)
    }
    
    model_states[[iter]]         = my_states
    record_count                 = record_count + 1
    
  }
  parallel::stopCluster(clust)
  hyperparameters$current_G = previous_G
  
  # Grab the changepoint probabilities.
  my_states = model_states
  states_df = data.frame(matrix(my_states[[1]], nrow = 1))
  for (i in 2:length(my_states)) {
    states_df = rbind(states_df, data.frame(matrix(my_states[[i]], nrow = 1)))
  }
  states_df_burnt = states_df[floor(2 * nrow(states_df) / 4):nrow(states_df), ]
  prop_mat = matrix(0, nrow = 1, ncol = (ncol(states_df_burnt) - 1))
  number_of_segments = nrow(states_df)
  for (j in 1:(ncol(states_df_burnt) - 1)) {
    temp_difference = states_df_burnt[, j + 1] - states_df_burnt[, j]
    prop_mat[1, j] = mean(temp_difference == 1)
  }
  time_points     = matrix(1:ncol(prop_mat),
                           ncol = ncol(prop_mat),
                           nrow = 1)
  prop_mat        = t(rbind(time_points, prop_mat))
  names(prop_mat) = c("Time-point", "Change-point probability")
  
  mylist = list(
    changepoint_probabilities = prop_mat,
    models = previous_model_fits,
    states = model_states,
    transition_probs = transition_probabilities_states,
    mergesplit_accepts = mergesplit_accepts,
    DRJ_accepts = DRJ_accepts,
    gibbswap_accepts = gibbswap_accepts,
    alphabeta_accepts = alphabeta_accepts,
    G_Wish_portion_trace = G_Wish_portion_trace,
    D_prior_dens_trace = D_prior_dens_trace,
    cholesky_failed_prec = cholesky_failed_prec,
    hyperparameters = hyperparameters,
    latent_data = full_data_Z,
    raw_data = data_woTimeValues,
    Z_timepoint_indices = Z_timepoint_indices,
    fault_detection_logs = model_saves_list
  )
  class(mylist) = "bayesWatch_object"
  
  print.bayesWatch_object <- function(x, ...) {
    cat("bayesWatch object\n")
    cat("-----------------\n")
    print(x$prop_mat)
  }
  return(mylist)
}


#' Title
#'
#' @param regime_fit_output
#' @param prob_cutoff
#'
#' @return
#' @export
#'
get_point_estimate = function(regime_fit_object, prob_cutoff) {
  ### Grab the change-point probabilities
  yy = regime_fit_object
  MVN_model_mergesplit_accepts = yy$mergesplit_accepts
  my_states = yy$states
  states_df = data.frame(matrix(my_states[[1]], nrow = 1))
  for (i in 2:length(my_states)) {
    states_df = rbind(states_df, data.frame(matrix(my_states[[i]], nrow = 1)))
  }
  states_df_burnt = states_df[floor(2 * nrow(states_df) / 4):nrow(states_df), ]
  prop_mat = matrix(0, nrow = 1, ncol = (ncol(states_df_burnt) - 1))
  number_of_segments = nrow(states_df)
  for (j in 1:(ncol(states_df_burnt) - 1)) {
    temp_difference = states_df_burnt[, j + 1] - states_df_burnt[, j]
    prop_mat[1, j] = mean(temp_difference == 1)
  }
  changepoint_probs   = prop_mat
  changepoints        = changepoint_probs >= prob_cutoff
  changepoints        = data.frame(matrix(
    changepoints,
    nrow = 1,
    ncol = length(changepoints)
  ))
  return(changepoints)
}


#' Title
#'
#' @param regime_fit_output
#' @param prob_cutoff
#'
#' @return
#' @export
#'
#' @examples
graph_changepoints = function(regime_fit_object, prob_cutoff) {
  ### Grab the change-point probabilities
  yy = regime_fit_object
  MVN_model_mergesplit_accepts = yy$mergesplit_accepts
  my_states = yy$states
  states_df = data.frame(matrix(my_states[[1]], nrow = 1))
  for (i in 2:length(my_states)) {
    states_df = rbind(states_df, data.frame(matrix(my_states[[i]], nrow = 1)))
  }
  states_df_burnt = states_df[floor(2 * nrow(states_df) / 4):nrow(states_df), ]
  prop_mat = matrix(0, nrow = 1, ncol = (ncol(states_df_burnt) - 1))
  number_of_segments = nrow(states_df)
  for (j in 1:(ncol(states_df_burnt) - 1)) {
    temp_difference = states_df_burnt[, j + 1] - states_df_burnt[, j]
    prop_mat[1, j] = mean(temp_difference == 1)
  }
  changepoint_probs   = prop_mat
  changepoints        = changepoint_probs >= prob_cutoff
  changepoints_data   = data.frame(prob_value = changepoints,
                                   time_point = 1:length(changepoints))
  
  my_plot = ggplot(changepoints_data, aes(x = time_point, y = prob_value)) + geom_bar(stat = "identity") +
    geom_hline(yintercept = prob_cutoff, color = "red") +
    annotate("text", label = "Probability Cutoff Value")
  
  return(my_plot)
}


#' Title
#'
#' @param regime_fit_object 
#'
#' @return
#' @export
#'
#' @examples
get_fault_detection_graph = function(regime_fit_object){
  if(is.null(regime_fit_object$fault_graph)){
    stop("This regime fit was run without fault detection.")
  } else {
    return(regime_fit_object$fault_graph)
  }
}
