## Alexander Murph
## Date: 01/17/2022
source("helpers.R")
library(latex2exp)
library(ggplot2)



#' Title
#'
#' @param changepoint_probs 
#' @param prob_cutoff 
#' @param variable_names 
#' @param num_cores 
#' @param simulation_number 
#'
#' @return
#'
#' @examples
determine_marginals = function(changepoint_probs, 
                               prob_cutoff, variable_names, num_cores = 30,
                               simulation_number = NULL){
  ###############################################################################
  # Start by making a state vector that mirrors the set of expected change-points.
  p            = length(variable_names)
  num_states   = length(changepoint_probs)+1
  changepoints = changepoint_probs >= prob_cutoff
  my_states    = rep(0, times = length(changepoints))
  
  if(sum(changepoints)==0){
    print("was unable to perform fault detection because no changepoints were detected!")
    return(NA)
  }
  
  my_states[max(which(changepoints==1))] = 1
  string_of_my_states = paste(my_states,collapse="")
  # We will use regular expressions to find all saved model parameters through the course
  # of the algorithm that matched with this state vector.
  if(is.null(simulation_number)){
    file_names_of_models = list.files(path="Model_Fits")
  }else{
    file_names_of_models = list.files(path="Model_Fits_Simulation")
  }
  
  location_of_files_wanted = c()
  for(file_index in 1:length(file_names_of_models)){
    temp_file_name     = strsplit(file_names_of_models[file_index], split = "_")
    if(!is.null(simulation_number)){
      # We only want to gather the models associated with this simulation.
      temp_file_sim_number = as.integer(temp_file_name[[1]][2])
      if(temp_file_sim_number!=simulation_number){
        location_of_files_wanted = c(location_of_files_wanted,0)
        next
      }
      string_name = temp_file_name[[1]][1]
      print(string_name)
      if(string_name == string_of_my_states){
        location_of_files_wanted = c(location_of_files_wanted,1)
        next
      } else {
        location_of_files_wanted = c(location_of_files_wanted,0)
        next
      }
      
    } else {
      # If I'm not working with a simulation, do it this way.
      string_name = temp_file_name[[1]][1]
      if(string_name == string_of_my_states){
        location_of_files_wanted = c(location_of_files_wanted,1)
        next
      } else {
        location_of_files_wanted = c(location_of_files_wanted,0)
        next
      }
    }
  }

  if(sum(location_of_files_wanted)==0){
    print("no model was saved that has this changepoint!!")
    print(string_of_my_states)
    return(NA)
  }
  
  file_names_of_my_states_models = file_names_of_models[as.logical(location_of_files_wanted)]
  
  # My temporary solution is to replace all the existing cluster assignments with 1's.  This seems a good 
  ## choice anyway because I am forcing there to only be a single precision matrix/mu vector.
  
  ###############################################################################
  # Now, gather latent data+parameters that mimic this state vector.
  print("--> Fitting the models and latent data to the state vector")
  learning_repetitions = 10
  count                = 0
  number_of_regimes    = max(my_states)
  monte_carlo_iters    = 10
  
  #################################################################################################
  # Now previous_model_fits and latent_data should accurately mirror the changepoint probabilities.
  # We can use these values to look at the marginal distributions for each of the p variables.
  # With these marginals, we approximate the JSD with Monte Carlo sampling.
  print("--> Calculating metrics")
  print(length(file_names_of_my_states_models))
  print(variable_names)
  differences_of_marginals = data.frame(Variable         = rep(variable_names,times = length(file_names_of_my_states_models)), 
                                        JSD              = rep(0, times = p*(length(file_names_of_my_states_models))),
                                        GeomJSD          = rep(0, times = p*(length(file_names_of_my_states_models))),
                                        KL               = rep(0, times = p*(length(file_names_of_my_states_models))),
                                        Jeff             = rep(0, times = p*(length(file_names_of_my_states_models))),
                                        KL_single        = rep(0, times = p*(length(file_names_of_my_states_models))),
                                        Hellinger        = rep(0, times = p*(length(file_names_of_my_states_models))),
                                        Hellinger_single = rep(0, times = p*(length(file_names_of_my_states_models))),
                                        Jeff_single      = rep(0, times = p*(length(file_names_of_my_states_models))),
                                        Iteration        = rep(1:(length(file_names_of_my_states_models)), each=p) )
  all_indices              = 1:p
  
  # Change this for actual applications.  This is just to save time (for the time being).
  # file_names_of_my_states_models = file_names_of_my_states_models[(length(file_names_of_my_states_models)-60):(length(file_names_of_my_states_models)-50)]
  print(paste("There are a total of", length(file_names_of_my_states_models), "instances of model data."))
  for(model_index in 1:length(file_names_of_my_states_models)){
    # Grab the MVN parameters for in-between each regime.
    regime_index = 1
    if(is.null(simulation_number)){
      print(paste("Model_Fits/",file_names_of_my_states_models[model_index],sep=""))
      previous_model_fits = readRDS(paste("Model_Fits/",file_names_of_my_states_models[model_index],sep=""))
    }else {
      print(paste("Model_Fits_Simulation/",file_names_of_my_states_models[model_index],sep=""))
      previous_model_fits = readRDS(paste("Model_Fits_Simulation/",file_names_of_my_states_models[model_index],sep=""))
    }
    comp_probabilities  = exp(previous_model_fits[[regime_index]]$component_log_probs)

    probs_of_model_before = exp(previous_model_fits[[regime_index]]$component_log_probs)
    probs_of_model_after  = exp(previous_model_fits[[regime_index+1]]$component_log_probs)
    precisions_before     = previous_model_fits[[regime_index]]$precision
    precisions_after      = previous_model_fits[[regime_index+1]]$precision
    mus_before            = previous_model_fits[[regime_index]]$mu
    mus_after             = previous_model_fits[[regime_index+1]]$mu
    comps_before          = previous_model_fits[[regime_index]]$cluster_assignments
    comps_after           = previous_model_fits[[regime_index+1]]$cluster_assignments
    comps                 = c(comps_before, comps_after)
  
    cov_before            = 0
    cov_after             = 0
    mu_before             = 0
    mu_after              = 0
    for(comp_index in 1:max(comps)){
      comp_prob_before  = probs_of_model_before[comp_index]
      comp_prob_after   = probs_of_model_after[comp_index]
      cov_before        = cov_before + comp_prob_before*solve(precisions_before[[comp_index]])
      cov_after         = cov_after + comp_prob_after*solve(precisions_after[[comp_index]])
      mu_before         = mu_before + comp_prob_before*mus_before[[comp_index]]
      mu_after          = mu_after + comp_prob_after*mus_after[[comp_index]]
    }
    prec_before           = solve(cov_before)
    prec_after            = solve(cov_after)
    
    # Calculate the full Hellinger Distance
    Hellinger_full = sqrt(1 - (((det(cov_before)**(0.25))*(det(cov_after)**(0.25)))/(det(0.5*cov_before+0.5*cov_after)**(0.5)))*
                        exp(-(1/8)*t(mu_after - mu_before)%*%solve(((0.5*cov_before+0.5*cov_after)))%*%as.matrix(mu_after - mu_before) ))
      
    # Calculate the full KL Divergence
    KL_full     = 0.5*( sum(diag(prec_after%*%cov_before))-p +
                   t(mu_after - mu_before)%*%prec_after%*%as.matrix(mu_after - mu_before) + 
                   log(det(cov_after)/det(cov_before)) )
    
    # Calculate the full Geometric JSD according to F. Nielsen 2019
    harmonic_barycenter_prec = 0.5*prec_before + 0.5*prec_after
    harmonic_barycenter_cov  = solve(0.5*prec_before + 0.5*prec_after)
    harmonic_barycenter_mu   = harmonic_barycenter_cov%*%(0.5*prec_before%*%mu_before + 
                                                            0.5*prec_after%*%mu_after)
    log_det_barycenter_cov   = log(det(harmonic_barycenter_cov))
    log_det_cov_before       = log(det(cov_before))
    log_det_cov_after        = log(det(cov_after))
    
    Geom_JSD_full = sum(diag(harmonic_barycenter_prec*(0.5*cov_before + 0.5*cov_after)))
    Geom_JSD_full = Geom_JSD_full + log_det_barycenter_cov - 0.5*log_det_cov_before - 
      0.5*log_det_cov_after
    Geom_JSD_full = Geom_JSD_full + 0.5*t(harmonic_barycenter_mu - mu_before)%*%
      harmonic_barycenter_prec%*%
      (harmonic_barycenter_mu - mu_before)
    Geom_JSD_full = Geom_JSD_full + 0.5*t(harmonic_barycenter_mu - mu_after)%*%
      harmonic_barycenter_prec%*%
      (harmonic_barycenter_mu - mu_after)
    Geom_JSD_full = Geom_JSD_full - p
    Geom_JSD_full = 0.5*Geom_JSD_full 
    full_JSD             = 1
    
    ## Swap the covariances and mu's to calculate Jeffrey's Divergence
    cov_before_wo_p  = solve(prec_after)
    cov_after_wo_p   = solve(prec_before)
    mu_after_wo_p    = mu_before
    mu_before_wo_p   = mu_after
    prec_after_wo_p  = prec_before
    full_Jeff        = 0.5*( sum(diag(prec_after_wo_p%*%cov_before_wo_p))-p +
                               t(mu_after_wo_p - mu_before_wo_p)%*%prec_after_wo_p%*%as.matrix(mu_after_wo_p - mu_before_wo_p) + 
                               log(det(cov_after_wo_p)/det(cov_before_wo_p)) ) + KL_full
    
    for(variable_index in 1:p){
      # I'll start by calculating the marginal KL divergence for the variable at this index.
      cov_before_of_p  = cov_before[variable_index,variable_index]
      cov_after_of_p   = cov_after[variable_index, variable_index]
      prec_after_of_p  = 1/cov_after_of_p
      prec_before_of_p = 1/cov_before_of_p
      mu_before_of_p   = mu_before[variable_index]
      mu_after_of_p    = mu_after[variable_index]
      KL_of_p          = 0.5*( (prec_after_of_p*cov_before_of_p)-1 +
                           ((mu_after_of_p - mu_before_of_p)**2)%*%prec_after_of_p + 
                           log(cov_after_of_p/cov_before_of_p) )
      
      differences_of_marginals[(model_index-1)*p + variable_index,6] = (KL_of_p / KL_full)
      Hellinger_of_p = sqrt(1 - ((((cov_before_of_p)**(0.25))*((cov_after_of_p)**(0.25)))/((0.5*cov_before_of_p+0.5*cov_after_of_p)**(0.5)))*
                               exp(-(1/8)*t(mu_after_of_p-mu_before_of_p)%*%(((0.5*cov_before_of_p+0.5*cov_after_of_p)**(-1)))%*%as.matrix(mu_after_of_p-mu_before_of_p) ))
      differences_of_marginals[(model_index-1)*p + variable_index,8] = (Hellinger_of_p / Hellinger_full)
      
      # Now get the KL in the other direction so to use Jeffrey's
      cov_before_of_p  = cov_after[variable_index, variable_index]
      cov_after_of_p   = cov_before[variable_index,variable_index]
      prec_after_of_p  = 1/cov_after_of_p
      prec_before_of_p = 1/cov_before_of_p
      mu_before_of_p   = mu_after[variable_index]
      mu_after_of_p    = mu_before[variable_index]
      KL_reverse_of_p  = 0.5*( (prec_after_of_p*cov_before_of_p)-1 +
                                 ((mu_after_of_p - mu_before_of_p)**2)%*%prec_after_of_p + 
                                 log(cov_after_of_p/cov_before_of_p) )
      Jeff_of_p        = KL_of_p + KL_reverse_of_p
      
      differences_of_marginals[(model_index-1)*p + variable_index,9] = (Jeff_of_p / full_Jeff)
      
      
      ## Get the parameters for the marginal distributions WITHOUT the variable from variable_index.
      indices_wo_p     = all_indices[-variable_index]
      prec_before_wo_p = prec_before[indices_wo_p,indices_wo_p]-
                          (as.matrix(prec_before[indices_wo_p,variable_index])%*%prec_before[variable_index,indices_wo_p])/prec_before[variable_index,variable_index]
      mu_before_wo_p   = mu_before[indices_wo_p]
      prec_after_wo_p  = prec_after[indices_wo_p,indices_wo_p]-
                          (as.matrix(prec_after[indices_wo_p,variable_index])%*%prec_after[variable_index,indices_wo_p])/prec_after[variable_index,variable_index]
      mu_after_wo_p    = mu_after[indices_wo_p]
      
      cov_before_wo_p  = solve(prec_before_wo_p)
      cov_after_wo_p   = solve(prec_after_wo_p)
      
      # Hellinger Distance
      Hellinger = sqrt(1 - (((det(cov_before_wo_p)**(0.25))*(det(cov_after_wo_p)**(0.25)))/(det(0.5*cov_before_wo_p+0.5*cov_after_wo_p)**(0.5)))*
                         exp(-(1/8)*t(mu_after_wo_p-mu_before_wo_p)%*%solve(((0.5*cov_before_wo_p+0.5*cov_after_wo_p)))%*%as.matrix(mu_after_wo_p-mu_before_wo_p) ))
      
      # Calculate KL Divergence
      KL   = 0.5*( sum(diag(prec_after_wo_p%*%cov_before_wo_p))-p+1 +
                  (mu_after_wo_p - mu_before_wo_p)%*%prec_after_wo_p%*%as.matrix(mu_after_wo_p - mu_before_wo_p) + 
                  log(det(cov_after_wo_p)/det(cov_before_wo_p)) )
      
      # Calculate Geometric JSD according to F. Nielsen 2019
      harmonic_barycenter_prec = 0.5*prec_before_wo_p + 0.5*prec_after_wo_p
      harmonic_barycenter_cov  = solve(0.5*prec_before_wo_p + 0.5*prec_after_wo_p)
      harmonic_barycenter_mu   = harmonic_barycenter_cov%*%(0.5*prec_before_wo_p%*%mu_before_wo_p + 
                                                             0.5*prec_after_wo_p%*%mu_after_wo_p)
      log_det_barycenter_cov   = log(det(harmonic_barycenter_cov))
      log_det_cov_before_wo_p  = log(det(cov_before_wo_p))
      log_det_cov_after_wo_p   = log(det(cov_after_wo_p))
      
      Geom_JSD = sum(diag(harmonic_barycenter_prec*(0.5*cov_before_wo_p + 0.5*cov_after_wo_p)))
      Geom_JSD = Geom_JSD + log_det_barycenter_cov - 0.5*log_det_cov_before_wo_p - 
                                                      0.5*log_det_cov_after_wo_p
      Geom_JSD = Geom_JSD + 0.5*t(harmonic_barycenter_mu - mu_before_wo_p)%*%
                                  harmonic_barycenter_prec%*%
                                  (harmonic_barycenter_mu - mu_before_wo_p)
      Geom_JSD = Geom_JSD + 0.5*t(harmonic_barycenter_mu - mu_after_wo_p)%*%
                                  harmonic_barycenter_prec%*%
                                  (harmonic_barycenter_mu - mu_after_wo_p)
      Geom_JSD = Geom_JSD - p + 1
      Geom_JSD = 0.5*Geom_JSD 

      differences_of_marginals[(model_index-1)*p + variable_index,2] = 1

      differences_of_marginals[(model_index-1)*p + variable_index,4] = (1 - KL / KL_full)
      differences_of_marginals[(model_index-1)*p + variable_index,7] = (1 - Hellinger / Hellinger_full)
      differences_of_marginals[(model_index-1)*p + variable_index,3] = (1 - Geom_JSD / Geom_JSD_full)

      ## Swap the covariances and mu's to calculate Jeffrey's Divergence
      cov_before_wo_p  = solve(prec_after_wo_p)
      cov_after_wo_p   = solve(prec_before_wo_p)
      mu_after_wo_p    = mu_before[indices_wo_p]
      mu_before_wo_p   = mu_after[indices_wo_p]
      prec_after_wo_p  = prec_before_wo_p
      Jeff             = 0.5*( sum(diag(prec_after_wo_p%*%cov_before_wo_p))-p+1 +
                         (mu_after_wo_p - mu_before_wo_p)%*%prec_after_wo_p%*%as.matrix(mu_after_wo_p - mu_before_wo_p) + 
                         log(det(cov_after_wo_p)/det(cov_before_wo_p)) ) + KL
      
      differences_of_marginals[(model_index-1)*p + variable_index,5] = (1 - Jeff / full_Jeff)
      if(is.na(differences_of_marginals[(model_index-1)*p + variable_index,4])){
        browser()
      }
    }

  }

  #################################################################################################
  # Print out the results.
  cat(c("[1]  -----> the states are:", paste(sapply(my_states, paste, collapse=' ')), "\n"))
  regime_index = 1
  cat(c("[1]  ", as.character(max(my_states)-1), " --> ", as.character(max(my_states)), "\n", "\n"))
  
  data_for_this_regime           = differences_of_marginals
  colnames(data_for_this_regime) = c("Variable Names", "Approx JSD w/o Variable", "Geometric JSD (1-norm)", 
                                     "Exact KL: Before || After (1-norm)", "Jeffrey's Divergence (1-norm)",
                                     "Iteration Number")
  
  return(differences_of_marginals)
}


#' Title
#'
#' @param differences_of_marginals 
#' @param num_of_dimensions 
#' @param changepoint_probs 
#' @param prob_cutoff 
#' @param true_parameter_values 
#' @param simulation_number 
#'
#' @return
#' @export
#'
#' @examples
graph_posterior_distributions = function(differences_of_marginals, num_of_dimensions, 
                                         changepoint_probs, prob_cutoff, true_parameter_values = NULL,
                                         simulation_number = NULL){
  print("inside of the posterior graph function")
  ###############################################################################
  # Start by making a state vector that mirrors thet set of expected changepoints.
  p                    = num_of_dimensions
  
  regime_index = 1
  num_states   = length(changepoint_probs)+1
  changepoints = changepoint_probs >= prob_cutoff
  my_states    = rep(0, times = num_states)
  count        = 1
  my_states[1] = count
  
  for(index in 1:length(changepoint_probs)){
    if(changepoints[index]){
      count            = count + 1
    }
    my_states[index+1] = count
  }
  graph_title          = paste("Regime Change", max(my_states)-1, "-->", max(my_states))
  
  changepoints = changepoint_probs >= prob_cutoff
  my_states    = rep(0, times = length(changepoints))
  if(sum(changepoints)==0){
    stop("was unable to perform fault detection because no changepoints were detected!")
  }
  
  my_states[max(which(changepoints==1))] = 1
  string_of_my_states = paste(my_states,collapse="")

  ###############################################################################
  # KL DIVERGENCE GRAPHS
  average_values = c()
  for(var_index in 1:length(unique(differences_of_marginals$Variable))){
    temp_data = differences_of_marginals[which(differences_of_marginals$Variable == var_index),]
    average_values = c(average_values, mean(temp_data$KL))
  }
  possible_variables = 1:length(unique(differences_of_marginals$Variable))
  ordered_variables  = possible_variables[order(average_values, decreasing=TRUE)]
  top_5_variables    = ordered_variables[1:5]
  
  differences_of_marginals_top5          = differences_of_marginals[which(differences_of_marginals$Variable%in%top_5_variables),]
  differences_of_marginals_top5$Variable = factor(differences_of_marginals_top5$Variable, levels=top_5_variables)
  
  
  if(is.null(true_parameter_values)){
    # If the true parameter values aren't known (i.e. this is NOT a simulation test):
    p1 <- ggplot(data = differences_of_marginals_top5, aes(x=KL)) + geom_density(alpha=.4, fill="#FF6666")
    p1 <- p1 + facet_grid(rows = vars(Variable), scales="free")
    p1 <- p1 + ggtitle(paste("Total-Effect KL Loss,", graph_title)) + theme(text = element_text(size = 13))
    p1 = p1 + ylab("Density") + xlab("Total-Effect KL") + theme(text = element_text(size = 13))
    ggsave(
      filename = "Posterior_Graphs/posterior_KL.pdf",
      plot = p1
    )
    
    p3 <- ggplot(data = differences_of_marginals_top5, aes(x=KL_single)) + geom_density(alpha=.4, fill="#FF6666")
    p3 <- p3 + facet_grid(rows = vars(Variable), scales="free")
    p3 <- p3 + ggtitle(paste("First-Order KL Loss,", graph_title)) + theme(text = element_text(size = 13))
    p3 = p3 + ylab("Density") + xlab("First-Order KL") + theme(text = element_text(size = 13))
    ggsave(
      filename = "Posterior_Graphs/posterior_KL_single.pdf",
      plot = p3
    )
  } else {
    # If the true parameter values are known:
    prec_before     = true_parameter_values$sigma1
    prec_after      = true_parameter_values$sigma2
    mu_before       = true_parameter_values$mu1
    mu_after        = true_parameter_values$mu2
    
    cov_before      = solve(prec_before)
    cov_after       = solve(prec_after)
    
    Hellinger_full = sqrt(1 - (((det(cov_before)**(0.25))*(det(cov_after)**(0.25)))/(det(0.5*cov_before+0.5*cov_after)**(0.5)))*
                            exp(-(1/8)*t(mu_after - mu_before)%*%solve(((0.5*cov_before+0.5*cov_after)))%*%as.matrix(mu_after - mu_before) ))
    KL_full     = 0.5*( sum(diag(prec_after%*%cov_before))-p +
                          t(mu_after - mu_before)%*%prec_after%*%as.matrix(mu_after - mu_before) + 
                          log(det(cov_after)/det(cov_before)) )
    cov_before_wo_p  = solve(prec_after)
    cov_after_wo_p   = solve(prec_before)
    mu_after_wo_p    = mu_before
    mu_before_wo_p   = mu_after
    prec_after_wo_p  = prec_before
    full_Jeff        = 0.5*( sum(diag(prec_after_wo_p%*%cov_before_wo_p))-p +
                               t(mu_after_wo_p - mu_before_wo_p)%*%prec_after_wo_p%*%as.matrix(mu_after_wo_p - mu_before_wo_p) + 
                               log(det(cov_after_wo_p)/det(cov_before_wo_p)) ) + KL_full
    
    # Iterate through the top 5 variables:
    KL_marginals        = c()
    Hellinger_marginals = c()
    Jeff_marginals      = c()
    KL_one_out          = c()
    Hellinger_one_out   = c()
    Jeff_one_out        = c()
    all_indices         = 1:p
    for(variable_index in top_5_variables){
      # Marginals OF p
      cov_before_of_p  = cov_before[variable_index,variable_index]
      cov_after_of_p   = cov_after[variable_index, variable_index]
      prec_after_of_p  = 1/cov_after_of_p
      prec_before_of_p = 1/cov_before_of_p
      mu_before_of_p   = mu_before[variable_index]
      mu_after_of_p    = mu_after[variable_index]
      KL_of_p          = 0.5*( (prec_after_of_p*cov_before_of_p)-1 +
                                 ((mu_after_of_p - mu_before_of_p)**2)%*%prec_after_of_p + 
                                 log(cov_after_of_p/cov_before_of_p) )
      KL_marginals = c(KL_marginals, (KL_of_p / KL_full))
      
      Hellinger_of_p = sqrt(1 - ((((cov_before_of_p)**(0.25))*((cov_after_of_p)**(0.25)))/((0.5*cov_before_of_p+0.5*cov_after_of_p)**(0.5)))*
                              exp(-(1/8)*t(mu_after_of_p-mu_before_of_p)%*%(((0.5*cov_before_of_p+0.5*cov_after_of_p)**(-1)))%*%as.matrix(mu_after_of_p-mu_before_of_p) ))
      Hellinger_marginals = c(Hellinger_marginals, Hellinger_of_p)
      
      # Now get the KL in the other direction so to use Jeffrey's
      cov_before_of_p  = cov_after[variable_index, variable_index]
      cov_after_of_p   = cov_before[variable_index, variable_index]
      prec_after_of_p  = 1/cov_after_of_p
      prec_before_of_p = 1/cov_before_of_p
      mu_before_of_p   = mu_after[variable_index]
      mu_after_of_p    = mu_before[variable_index]
      KL_reverse_of_p  = 0.5*( (prec_after_of_p*cov_before_of_p)-1 +
                                 ((mu_after_of_p - mu_before_of_p)**2)%*%prec_after_of_p + 
                                 log(cov_after_of_p/cov_before_of_p) )
      Jeff_of_p        = KL_of_p + KL_reverse_of_p
      Jeff_marginals   = c(Jeff_marginals, Jeff_of_p / full_Jeff)
      
      # Marginals WITHOUT p
      indices_wo_p     = all_indices[-variable_index]
      prec_before_wo_p = prec_before[indices_wo_p,indices_wo_p]-
        (as.matrix(prec_before[indices_wo_p,variable_index])%*%prec_before[variable_index,indices_wo_p])/prec_before[variable_index,variable_index]
      mu_before_wo_p   = mu_before[indices_wo_p]
      prec_after_wo_p  = prec_after[indices_wo_p,indices_wo_p]-
        (as.matrix(prec_after[indices_wo_p,variable_index])%*%prec_after[variable_index,indices_wo_p])/prec_after[variable_index,variable_index]
      mu_after_wo_p    = mu_after[indices_wo_p]
      
      cov_before_wo_p  = solve(prec_before_wo_p)
      cov_after_wo_p   = solve(prec_after_wo_p)
      
      Hellinger         = sqrt(1 - (((det(cov_before_wo_p)**(0.25))*(det(cov_after_wo_p)**(0.25)))/(det(0.5*cov_before_wo_p+0.5*cov_after_wo_p)**(0.5)))*
                              exp(-(1/8)*t(mu_after_wo_p-mu_before_wo_p)%*%solve(((0.5*cov_before_wo_p+0.5*cov_after_wo_p)))%*%as.matrix(mu_after_wo_p-mu_before_wo_p) ))
      Hellinger_one_out = c(Hellinger_one_out, (1 - Hellinger/Hellinger_full))
      
      KL         = 0.5*( sum(diag(prec_after_wo_p%*%cov_before_wo_p))-p+1 +
                     (mu_after_wo_p - mu_before_wo_p)%*%prec_after_wo_p%*%as.matrix(mu_after_wo_p - mu_before_wo_p) + 
                     log(det(cov_after_wo_p)/det(cov_before_wo_p)) )
      KL_one_out = c(KL_one_out, (1 - KL/KL_full))
      
      cov_before_wo_p  = solve(prec_after_wo_p)
      cov_after_wo_p   = solve(prec_before_wo_p)
      mu_after_wo_p    = mu_before[indices_wo_p]
      mu_before_wo_p   = mu_after[indices_wo_p]
      prec_after_wo_p  = prec_before_wo_p
      Jeff             = 0.5*( sum(diag(prec_after_wo_p%*%cov_before_wo_p))-p+1 +
                                 (mu_after_wo_p - mu_before_wo_p)%*%prec_after_wo_p%*%as.matrix(mu_after_wo_p - mu_before_wo_p) + 
                                 log(det(cov_after_wo_p)/det(cov_before_wo_p)) ) + KL
      Jeff_one_out     = c(Jeff_one_out, 1 - Jeff/full_Jeff)
    }
    
    ##### Now we actually use all the above information to make the graphs:
    top_5_variables_factor =  factor(as.character(top_5_variables))
    
    data_vline <- data.frame(Variable = top_5_variables_factor, 
                             vline    = KL_one_out)
    text_locations_x = KL_one_out - 0.01
    # To get the text y locations, I need to know the density values for each of the 5 variables
    text_locations_y = c()
    for(var_index in top_5_variables){
      temp_top = differences_of_marginals_top5[differences_of_marginals_top5$Variable==var_index,]
      d <- density(temp_top$KL)
      text_locations_y = c(text_locations_y, 0.5*max(d$y))
    }
    
    dat_text = data.frame(
                            labels    = rep("True Total-Effect", times = 5),
                            vline     = text_locations_x,
                            y_heights = text_locations_y,
                            Variable  = top_5_variables_factor
                          )
    p1 = ggplot(data = differences_of_marginals_top5, aes(x=KL)) + geom_density(alpha=.4, fill="#FF6666")
    p1 = p1 + facet_grid(rows = vars(Variable), scales="free")
    p1 = p1 + ggtitle(paste("Total-Effect KL Loss,", graph_title))
    p1 = p1 + geom_vline(data = data_vline,
                        aes(xintercept = vline), colour = "blue")
    p1 = p1 + geom_text(data = dat_text, aes(x=vline, label=labels, y=y_heights), colour="blue", angle=90,  size=3)
    p1 = p1 + ylab("Density") + xlab("Total-Effect KL") + theme(text = element_text(size = 13))
    ggsave(
      filename = "Posterior_Graphs/posterior_KL.pdf",
      plot = p1
    )
    
    data_vline <- data.frame(Variable = top_5_variables_factor, 
                             vline    = KL_marginals)
    text_locations_x = KL_marginals - 0.01
    # To get the text y locations, I need to know the density values for each of the 5 variables
    text_locations_y = c()
    for(var_index in top_5_variables){
      temp_top = differences_of_marginals_top5[differences_of_marginals_top5$Variable==var_index,]
      d <- density(temp_top$KL_single)
      text_locations_y = c(text_locations_y, 0.5*max(d$y))
    }
    
    dat_text = data.frame(
      labels    = rep("True First-Order", times = 5),
      vline     = text_locations_x,
      y_heights = text_locations_y,
      Variable  = top_5_variables_factor
    )
    p3 = ggplot(data = differences_of_marginals_top5, aes(x=KL_single)) + geom_density(alpha=.4, fill="#FF6666")
    p3 = p3 + facet_grid(rows = vars(Variable), scales="free")
    p3 = p3 + ggtitle(paste("First-Order KL Loss,", graph_title))
    p3 = p3 + geom_vline(data = data_vline,
                        aes(xintercept = vline))
    p3 = p3 + geom_text(data = dat_text, aes(x=vline, label=labels, y=y_heights), colour="purple", angle=90,  size=3)
    p3 = p3 + ylab("Density") + xlab("First-Order KL") + theme(text = element_text(size = 13))

    ggsave(
      filename = "Posterior_Graphs/posterior_KL_single.pdf",
      plot = p3
    )
    
  }
  
  # Using a representative model output, graph the marginals before and after the changepoint.
  if(is.null(simulation_number)){
    file_names_of_models = list.files(path="Model_Fits")
    my_states_of_files = c()
    for(file_index in 1:length(file_names_of_models)){
      temp_file_name     = strsplit(file_names_of_models[file_index], split = "_")
      temp_file_name     = temp_file_name[[1]][1]
      if(temp_file_name==string_of_my_states){
        my_states_of_files = c(my_states_of_files, 1)
      } else {
        my_states_of_files = c(my_states_of_files, 0)
      }
      
    }
  } else {
    file_names_of_models = list.files(path="Model_Fits_Simulation")
    my_states_of_files   = c()
    for(file_index in 1:length(file_names_of_models)){
      temp_file_name     = strsplit(file_names_of_models[file_index], split = "_")
      if(!is.null(simulation_number)){
        # We only want to gather the models associated with this simulation.
        temp_file_sim_number = as.integer(temp_file_name[[1]][2])
        if(temp_file_sim_number!=simulation_number){
          my_states_of_files = c(my_states_of_files, 0)
          next
        }
        string_name = temp_file_name[[1]][1]
        if(string_name==string_of_my_states){
          print(file_names_of_models[file_index])
          my_states_of_files = c(my_states_of_files, 1)
        } else {
          my_states_of_files = c(my_states_of_files, 0)
        }
      }
    }
  }
  
  if(length(my_states_of_files)==0){
    stop("There were no files from that simulation!")
  }
  
  
  # This is making sure that I'm grabbing a model with the correct regime vector:
  file_names_of_my_states_models = file_names_of_models[which(as.logical(my_states_of_files))]
  if(length(file_names_of_my_states_models)==0){
    stop("There were no files with that change-point!")
  }
  
  # I take the last one:
  representative_file            = file_names_of_my_states_models[length(file_names_of_my_states_models)]
  if(is.null(simulation_number)){
    previous_model_fits            = readRDS(paste("Model_Fits/",representative_file,sep=""))
  } else {
    previous_model_fits            = readRDS(paste("Model_Fits_Simulation/",representative_file,sep=""))
  }
  
  comp_probabilities             = exp(previous_model_fits[[regime_index]]$component_log_probs)
  
  probs_of_model_before = exp(previous_model_fits[[regime_index]]$component_log_probs)
  probs_of_model_after  = exp(previous_model_fits[[regime_index+1]]$component_log_probs)
  precisions_before     = previous_model_fits[[regime_index]]$precision
  precisions_after      = previous_model_fits[[regime_index+1]]$precision
  mus_before            = previous_model_fits[[regime_index]]$mu
  mus_after             = previous_model_fits[[regime_index+1]]$mu
  comps_before          = previous_model_fits[[regime_index]]$cluster_assignments
  comps_after           = previous_model_fits[[regime_index+1]]$cluster_assignments
  comps                 = c(comps_before, comps_after)
  cov_before            = 0
  cov_after             = 0
  mu_before             = 0
  mu_after              = 0
  for(comp_index in 1:max(comps)){
    comp_prob_before  = probs_of_model_before[comp_index]
    comp_prob_after   = probs_of_model_after[comp_index]
    cov_before        = cov_before + comp_prob_before*solve(precisions_before[[comp_index]])
    cov_after         = cov_after + comp_prob_after*solve(precisions_after[[comp_index]])
    mu_before         = mu_before + comp_prob_before*mus_before[[comp_index]]
    mu_after          = mu_after + comp_prob_after*mus_after[[comp_index]]
  }
  prec_before           = solve(cov_before)
  prec_after            = solve(cov_after)
  
  marginal_normals_data = data.frame(Value = NA, x_value=NA, Variable = NA, Before_or_After = NA)
  
  num_to_sample         = 1000
  max_var  = sqrt(max(c(diag(cov_before),diag(cov_after))))
  max_mu   = max(c(mu_before,mu_after))
  min_mu   = min(c(mu_before,mu_after))
  upper_x  = max_mu + 2.5*max_var
  lower_x  = min_mu - 2.5*max_var
  x_values = (1:num_to_sample/(num_to_sample))*(upper_x - lower_x) + lower_x
  
  for(var_index in unique(top_5_variables)){
    temp_rows1 = data.frame(Value = dnorm(x_values,mean=mu_before[var_index],sd=sqrt(cov_before[var_index,var_index])),
                            x_value = x_values,
                            Variable = rep(var_index,times=num_to_sample),
                            Before_or_After = rep("Before", times = num_to_sample))
    temp_rows2 = data.frame(Value = dnorm(x_values, mean=mu_after[var_index],sd=sqrt(cov_after[var_index,var_index])),
                            x_value = x_values,
                            Variable = rep(var_index,times=num_to_sample),
                            Before_or_After = rep("After", times = num_to_sample))
    marginal_normals_data = rbind(marginal_normals_data,temp_rows1)
    marginal_normals_data = rbind(marginal_normals_data,temp_rows2)
  }
  marginal_normals_data               =  marginal_normals_data[-1,]
  marginal_normals_data$Variable      = factor(marginal_normals_data$Variable, levels = as.character(top_5_variables))
  
  if(is.null(true_parameter_values)){
    p2 <- ggplot(data = marginal_normals_data) + geom_area(aes(x=x_value, y=Value, fill=Before_or_After), position = "identity", alpha=0.4)
    p2 <- p2 + facet_grid(rows = vars(Variable), scales = 'free')
    p2 <- p2 + ggtitle(paste("Representative Marginal Distributions")) + theme(text = element_text(size = 13))
    p2 = p2 + ylab("Density") + xlab("Variable Value") + labs(fill = "Change-point\nOccurs")
  } else {
    # To get the text y locations, I need to know the density values for each of the 5 variables
    text_locations_y_before        = c()
    text_locations_y_after         = c()
    vertical_line_locations_before = c()
    vertical_line_locations_after  = c()
    for(var_index in top_5_variables){
      temp_top                = marginal_normals_data[which( (marginal_normals_data$Variable==var_index)&
                                                              (marginal_normals_data$Before_or_After=="Before") ),]
      d1                      = temp_top$Value
      
      temp_top                = marginal_normals_data[which( (marginal_normals_data$Variable==var_index)&
                                                               (marginal_normals_data$Before_or_After=="After") ),]
      d2                      = temp_top$Value
      text_locations_y_before = c(text_locations_y_before, 0.5*max(c(d1,d2)))
      text_locations_y_after  = c(text_locations_y_after, 0.5*max(c(d1,d2)))

      true_average                   = true_parameter_values$mu1[var_index]
      vertical_line_locations_before = c(vertical_line_locations_before, true_average)
      true_average                   = true_parameter_values$mu2[var_index]
      vertical_line_locations_after  = c(vertical_line_locations_after, true_average)
    }
    
    data_vline_before <- data.frame(Variable = top_5_variables_factor, 
                                    vline    = vertical_line_locations_before)
    data_vline_after  <- data.frame(Variable = top_5_variables_factor, 
                                    vline    = vertical_line_locations_after)
    
    text_locations_x_before = vertical_line_locations_before - (upper_x - lower_x)/40
    text_locations_x_after  = vertical_line_locations_after + (upper_x - lower_x)/40
    dat_text_before = data.frame(
      labels    = rep("True Center Before", times = 5),
      vline     = text_locations_x_before,
      y_heights = text_locations_y_before,
      Variable  = top_5_variables_factor
    )
    dat_text_after = data.frame(
      labels    = rep("True Center After", times = 5),
      vline     = text_locations_x_after,
      y_heights = text_locations_y_after,
      Variable  = top_5_variables_factor
    )
    
    p2 = ggplot(data = marginal_normals_data) + geom_area(aes(x=x_value, y=Value, fill=Before_or_After), 
                                                          position = "identity", alpha=0.4)
    p2 = p2 + facet_grid(rows = vars(Variable), scales = 'free')
    
    p2 = p2 + geom_vline(data = data_vline_before,
                         aes(xintercept = vline), color = "purple")
    p2 = p2 + geom_text(data = dat_text_before, aes(x=vline, 
                                                    label=labels, 
                                                    y=y_heights), colour="purple", angle=90,  size=2.7)
    
    p2 = p2 + geom_vline(data = data_vline_after,
                         aes(xintercept = vline), color = "green")
    p2 = p2 + geom_text(data = dat_text_after, aes(x=vline, label= labels, 
                                                   y=y_heights), colour="green", angle=90,  size=2.7)
    
    p2 = p2 + ggtitle(paste("Representative Marginal Distributions")) + labs(fill="Change-point\nOccurs")
    p2 = p2 + ylab("Density") + xlab("Variable Value") + theme(text = element_text(size = 13))
  }
  
  prob_plot_list = list(p1,p3,p2)
  #### Create a graphic with all of the plots
  ggsave(
    filename = "Posterior_Graphs/posterior_KL_w_Marginals.pdf",
    plot = marrangeGrob(prob_plot_list, nrow=1, ncol=3, widths = c(1,1,1.5)),
    width = 12, height = 9
  )
  
  ###############################################################################
  # HELLINGER DISTANCE GRAPHS
  average_values = c()
  for(var_index in 1:length(unique(differences_of_marginals$Variable))){
    temp_data = differences_of_marginals[which(differences_of_marginals$Variable == var_index),]
    average_values = c(average_values, mean(temp_data$Hellinger))
  }
  possible_variables = 1:length(unique(differences_of_marginals$Variable))
  ordered_variables  = possible_variables[order(average_values, decreasing=TRUE)]
  top_5_variables    = ordered_variables[1:5]
  
  differences_of_marginals_top5 = differences_of_marginals[which(differences_of_marginals$Variable%in%top_5_variables),]
  differences_of_marginals_top5$Variable = factor(differences_of_marginals_top5$Variable, levels=top_5_variables)
  
  if(is.null(true_parameter_values)){
    # If the true parameter values aren't known (i.e. this is NOT a simulation test):
    p1 <- ggplot(data = differences_of_marginals_top5, aes(x=Hellinger)) + geom_density(alpha=.4, fill="#FF6666")#aes(fill=grp2), alpha = 0.4
    p1 <- p1 + facet_grid(rows = vars(Variable), scales="free")
    p1 <- p1 + ggtitle(paste("Total-Effect Hellinger Loss"))
    p1 = p1 + ylab("Density") + xlab("Total-Effect Hellinger") + theme(text = element_text(size = 13))
    ggsave(
      filename = "Posterior_Graphs/posterior_Hellinger.pdf",
      plot = p1
    )
    
    p3 <- ggplot(data = differences_of_marginals_top5, aes(x=Hellinger_single)) + geom_density(alpha=.4, fill="#FF6666")#aes(fill=grp2), alpha = 0.4
    p3 <- p3 + facet_grid(rows = vars(Variable), scales="free")
    p3 <- p3 + ggtitle(paste("First-Order Hellinger Loss"))
    p3 = p3 + ylab("Density") + xlab("First-Order Hellinger") + theme(text = element_text(size = 13))
    ggsave(
      filename = "Posterior_Graphs/posterior_Hellinger_single.pdf",
      plot = p3
    )
  } else {
    # If the true parameter values are known:
    prec_before     = true_parameter_values$sigma1
    prec_after      = true_parameter_values$sigma2
    mu_before       = true_parameter_values$mu1
    mu_after        = true_parameter_values$mu2
    
    cov_before      = solve(prec_before)
    cov_after       = solve(prec_after)
    
    Hellinger_full = sqrt(1 - (((det(cov_before)**(0.25))*(det(cov_after)**(0.25)))/(det(0.5*cov_before+0.5*cov_after)**(0.5)))*
                            exp(-(1/8)*t(mu_after - mu_before)%*%solve(((0.5*cov_before+0.5*cov_after)))%*%as.matrix(mu_after - mu_before) ))
    
    # Iterate through the top 5 variables:
    Hellinger_marginals = c()
    Hellinger_one_out   = c()
    all_indices         = 1:p
    for(variable_index in top_5_variables){
      cov_before_of_p  = cov_before[variable_index,variable_index]
      cov_after_of_p   = cov_after[variable_index, variable_index]
      prec_after_of_p  = 1/cov_after_of_p
      prec_before_of_p = 1/cov_before_of_p
      mu_before_of_p   = mu_before[variable_index]
      mu_after_of_p    = mu_after[variable_index]
      
      Hellinger_of_p = sqrt(1 - ((((cov_before_of_p)**(0.25))*((cov_after_of_p)**(0.25)))/((0.5*cov_before_of_p+0.5*cov_after_of_p)**(0.5)))*
                              exp(-(1/8)*t(mu_after_of_p-mu_before_of_p)%*%(((0.5*cov_before_of_p+0.5*cov_after_of_p)**(-1)))%*%as.matrix(mu_after_of_p-mu_before_of_p) ))
      Hellinger_marginals = c(Hellinger_marginals, Hellinger_of_p)
      
      # Marginals WITHOUT p
      indices_wo_p     = all_indices[-variable_index]
      prec_before_wo_p = prec_before[indices_wo_p,indices_wo_p]-
        (as.matrix(prec_before[indices_wo_p,variable_index])%*%prec_before[variable_index,indices_wo_p])/prec_before[variable_index,variable_index]
      mu_before_wo_p   = mu_before[indices_wo_p]
      prec_after_wo_p  = prec_after[indices_wo_p,indices_wo_p]-
        (as.matrix(prec_after[indices_wo_p,variable_index])%*%prec_after[variable_index,indices_wo_p])/prec_after[variable_index,variable_index]
      mu_after_wo_p    = mu_after[indices_wo_p]
      
      cov_before_wo_p  = solve(prec_before_wo_p)
      cov_after_wo_p   = solve(prec_after_wo_p)
      
      Hellinger         = sqrt(1 - (((det(cov_before_wo_p)**(0.25))*(det(cov_after_wo_p)**(0.25)))/(det(0.5*cov_before_wo_p+0.5*cov_after_wo_p)**(0.5)))*
                                 exp(-(1/8)*t(mu_after_wo_p-mu_before_wo_p)%*%solve(((0.5*cov_before_wo_p+0.5*cov_after_wo_p)))%*%as.matrix(mu_after_wo_p-mu_before_wo_p) ))
      Hellinger_one_out = c(Hellinger_one_out, (1 - Hellinger/Hellinger_full))
    }
    
    
    ##### Now we actually use all the above information to make the graphs:
    top_5_variables_factor =  factor(as.character(top_5_variables))
    
    data_vline <- data.frame(Variable = top_5_variables_factor, 
                             vline    = Hellinger_one_out)
    text_locations_x = Hellinger_one_out - 0.01
    # To get the text y locations, I need to know the density values for each of the 5 variables
    text_locations_y = c()
    for(var_index in top_5_variables){
      temp_top = differences_of_marginals_top5[differences_of_marginals_top5$Variable==var_index,]
      d <- density(temp_top$Hellinger)
      text_locations_y = c(text_locations_y, 0.5*max(d$y))
    }
    
    dat_text = data.frame(
      labels    = rep("True Total-Effect", times = 5),
      vline     = text_locations_x,
      y_heights = text_locations_y,
      Variable  = top_5_variables_factor
    )
    
    p1 = ggplot(data = differences_of_marginals_top5, aes(x=Hellinger)) + geom_density(alpha=.4, fill="#FF6666")
    p1 = p1 + facet_grid(rows = vars(Variable), scales="free")
    p1 = p1 + ggtitle(paste("Total-Effect Hellinger Loss"))
    p1 = p1 + geom_vline(data = data_vline,
                       aes(xintercept = vline), colour = "blue")
    p1 = p1 + geom_text(data = dat_text, aes(x=vline, label=labels, y=y_heights), colour="blue", angle=90,  size=3)
    p1 = p1 + ylab("Density") + xlab("Total-Effect Hellinger") + theme(text = element_text(size = 13))
    ggsave(
      filename = "Posterior_Graphs/posterior_Hellinger.pdf",
      plot = p1
    )
    
    data_vline <- data.frame(Variable = top_5_variables_factor, 
                             vline    = Hellinger_marginals)
    text_locations_x = Hellinger_marginals - 0.01
    # To get the text y locations, I need to know the density values for each of the 5 variables
    text_locations_y = c()
    for(var_index in top_5_variables){
      temp_top = differences_of_marginals_top5[differences_of_marginals_top5$Variable==var_index,]
      d <- density(temp_top$Hellinger_single)
      text_locations_y = c(text_locations_y, 0.5*max(d$y))
    }
    
    dat_text = data.frame(
      labels    = rep("True First-Order", times = 5),
      vline     = text_locations_x,
      y_heights = text_locations_y,
      Variable  = top_5_variables_factor
    )
    
    p3 = ggplot(data = differences_of_marginals_top5, aes(x=Hellinger_single)) + geom_density(alpha=.4, fill="#FF6666")#aes(fill=grp2), alpha = 0.4
    p3 = p3 + facet_grid(rows = vars(Variable), scales="free")
    p3 = p3 + ggtitle(paste("First-Order Hellinger Loss"))
    p3 = p3 + geom_vline(data = data_vline,
                         aes(xintercept = vline))
    p3 = p3 + geom_text(data = dat_text, aes(x=vline, label=labels, y=y_heights), colour="purple", angle=90,  size=3)
    p3 = p3 + ylab("Density") + xlab("First-Order Hellinger") + theme(text = element_text(size = 13))
    
    ggsave(
      filename = "Posterior_Graphs/posterior_Hellinger_single.pdf",
      plot = p3
    )
    
  }
  
  prob_plot_list = list(p1,p3,p2)
  #### Create a graphic with all of the plots
  ggsave(
    filename = "Posterior_Graphs/posterior_Hellinger_w_Marginals.pdf",
    plot = marrangeGrob(prob_plot_list, nrow=1, ncol=3, widths = c(1,1,1.5)),
    width = 12, height = 9
  )
  
  ###############################################################################
  # JEFFREY'S DIVERGENCE GRAPHS
  average_values = c()
  for(var_index in 1:length(unique(differences_of_marginals$Variable))){
    temp_data = differences_of_marginals[which(differences_of_marginals$Variable == var_index),]
    average_values = c(average_values, mean(temp_data$Jeff))
  }
  possible_variables = 1:length(unique(differences_of_marginals$Variable))
  ordered_variables  = possible_variables[order(average_values, decreasing=TRUE)]
  top_5_variables    = ordered_variables[1:5]
  
  differences_of_marginals_top5 = differences_of_marginals[which(differences_of_marginals$Variable%in%top_5_variables),]
  differences_of_marginals_top5$Variable = factor(differences_of_marginals_top5$Variable, levels=top_5_variables)
  
  if(is.null(true_parameter_values)){
    # If the true parameter values aren't known (i.e. this is NOT a simulation test):
    p1 <- ggplot(data = differences_of_marginals_top5, aes(x=Jeff)) + geom_density(alpha=.4, fill="#FF6666")
    p1 <- p1 + facet_grid(rows = vars(Variable), scales="free")
    p1 <- p1 + ggtitle(paste("Total-Effect Jeff Loss,", graph_title))
    p1 = p1 + ylab("Density") + xlab("Total-Effect Jeff") + theme(text = element_text(size = 13))
    ggsave(
      filename = "Posterior_Graphs/posterior_Jeff.pdf",
      plot = p1
    )
    
    p3 <- ggplot(data = differences_of_marginals_top5, aes(x=Jeff_single)) + geom_density(alpha=.4, fill="#FF6666")
    p3 <- p3 + facet_grid(rows = vars(Variable), scales="free")
    p3 <- p3 + ggtitle(paste("First-Order Jeff Loss,", graph_title))
    p3 = p3 + ylab("Density") + xlab("First-Order Jeff") + theme(text = element_text(size = 13))
    ggsave(
      filename = "Posterior_Graphs/posterior_Jeff_single.pdf",
      plot = p3
    )
  } else {
    # If the true parameter values are known:
    prec_before     = true_parameter_values$sigma1
    prec_after      = true_parameter_values$sigma2
    mu_before       = true_parameter_values$mu1
    mu_after        = true_parameter_values$mu2
    
    cov_before      = solve(prec_before)
    cov_after       = solve(prec_after)

    KL_full     = 0.5*( sum(diag(prec_after%*%cov_before))-p +
                          t(mu_after - mu_before)%*%prec_after%*%as.matrix(mu_after - mu_before) + 
                          log(det(cov_after)/det(cov_before)) )
    cov_before_wo_p  = solve(prec_after)
    cov_after_wo_p   = solve(prec_before)
    mu_after_wo_p    = mu_before
    mu_before_wo_p   = mu_after
    prec_after_wo_p  = prec_before
    full_Jeff        = 0.5*( sum(diag(prec_after_wo_p%*%cov_before_wo_p))-p +
                               t(mu_after_wo_p - mu_before_wo_p)%*%prec_after_wo_p%*%as.matrix(mu_after_wo_p - mu_before_wo_p) + 
                               log(det(cov_after_wo_p)/det(cov_before_wo_p)) ) + KL_full
    
    # Iterate through the top 5 variables:
    Jeff_marginals      = c()
    Jeff_one_out        = c()
    all_indices         = 1:p
    for(variable_index in top_5_variables){
      # Marginals OF p
      cov_before_of_p  = cov_before[variable_index,variable_index]
      cov_after_of_p   = cov_after[variable_index, variable_index]
      prec_after_of_p  = 1/cov_after_of_p
      prec_before_of_p = 1/cov_before_of_p
      mu_before_of_p   = mu_before[variable_index]
      mu_after_of_p    = mu_after[variable_index]
      KL_of_p          = 0.5*( (prec_after_of_p*cov_before_of_p)-1 +
                                 ((mu_after_of_p - mu_before_of_p)**2)%*%prec_after_of_p + 
                                 log(cov_after_of_p/cov_before_of_p) )
      
      # Now get the KL in the other direction so to use Jeffrey's
      cov_before_of_p  = cov_after[variable_index, variable_index]
      cov_after_of_p   = cov_before[variable_index, variable_index]
      prec_after_of_p  = 1/cov_after_of_p
      prec_before_of_p = 1/cov_before_of_p
      mu_before_of_p   = mu_after[variable_index]
      mu_after_of_p    = mu_before[variable_index]
      KL_reverse_of_p  = 0.5*( (prec_after_of_p*cov_before_of_p)-1 +
                                 ((mu_after_of_p - mu_before_of_p)**2)%*%prec_after_of_p + 
                                 log(cov_after_of_p/cov_before_of_p) )
      Jeff_of_p        = KL_of_p + KL_reverse_of_p
      Jeff_marginals   = c(Jeff_marginals, Jeff_of_p / full_Jeff)
      
      # Marginals WITHOUT p
      indices_wo_p     = all_indices[-variable_index]
      prec_before_wo_p = prec_before[indices_wo_p,indices_wo_p]-
        (as.matrix(prec_before[indices_wo_p,variable_index])%*%prec_before[variable_index,indices_wo_p])/prec_before[variable_index,variable_index]
      mu_before_wo_p   = mu_before[indices_wo_p]
      prec_after_wo_p  = prec_after[indices_wo_p,indices_wo_p]-
        (as.matrix(prec_after[indices_wo_p,variable_index])%*%prec_after[variable_index,indices_wo_p])/prec_after[variable_index,variable_index]
      mu_after_wo_p    = mu_after[indices_wo_p]
      
      cov_before_wo_p  = solve(prec_before_wo_p)
      cov_after_wo_p   = solve(prec_after_wo_p)
      
      KL         = 0.5*( sum(diag(prec_after_wo_p%*%cov_before_wo_p))-p+1 +
                           (mu_after_wo_p - mu_before_wo_p)%*%prec_after_wo_p%*%as.matrix(mu_after_wo_p - mu_before_wo_p) + 
                           log(det(cov_after_wo_p)/det(cov_before_wo_p)) )
      
      cov_before_wo_p  = solve(prec_after_wo_p)
      cov_after_wo_p   = solve(prec_before_wo_p)
      mu_after_wo_p    = mu_before[indices_wo_p]
      mu_before_wo_p   = mu_after[indices_wo_p]
      prec_after_wo_p  = prec_before_wo_p
      Jeff             = 0.5*( sum(diag(prec_after_wo_p%*%cov_before_wo_p))-p+1 +
                                 (mu_after_wo_p - mu_before_wo_p)%*%prec_after_wo_p%*%as.matrix(mu_after_wo_p - mu_before_wo_p) + 
                                 log(det(cov_after_wo_p)/det(cov_before_wo_p)) ) + KL
      Jeff_one_out     = c(Jeff_one_out, 1 - Jeff/full_Jeff)
    }
    
    
    ##### Now we actually use all the above information to make the graphs:
    top_5_variables_factor =  factor(as.character(top_5_variables))
    
    data_vline <- data.frame(Variable = top_5_variables_factor, 
                             vline    = Jeff_one_out)
    text_locations_x = Jeff_one_out - 0.01
    # To get the text y locations, I need to know the density values for each of the 5 variables
    text_locations_y = c()
    for(var_index in top_5_variables){
      temp_top = differences_of_marginals_top5[differences_of_marginals_top5$Variable==var_index,]
      d <- density(temp_top$Jeff)
      text_locations_y = c(text_locations_y, 0.5*max(d$y))
    }
    
    dat_text = data.frame(
      labels    = rep("True Total-Effect", times = 5),
      vline     = text_locations_x,
      y_heights = text_locations_y,
      Variable  = top_5_variables_factor
    )
    
    p1 = ggplot(data = differences_of_marginals_top5, aes(x=Jeff)) + geom_density(alpha=.4, fill="#FF6666")
    p1 = p1 + facet_grid(rows = vars(Variable), scales="free")
    p1 = p1 + ggtitle(paste("Total-Effect Jeffrey's Divergence Loss,", graph_title))
    p1 = p1 + geom_vline(data = data_vline,
                       aes(xintercept = vline), colour = "blue")
    p1 = p1 + geom_text(data = dat_text, aes(x=vline, label=labels, y=y_heights), colour="blue", angle=90,  size=3)
    p1 = p1 + ylab("Density") + xlab("Total-Effect Jeffrey's Divergence") + theme(text = element_text(size = 13))
    ggsave(
      filename = "Posterior_Graphs/posterior_Jeff.pdf",
      plot = p1
    )
    
    data_vline <- data.frame(Variable = top_5_variables_factor, 
                             vline    = Jeff_marginals)
    text_locations_x = Jeff_marginals - 0.01
    # To get the text y locations, I need to know the density values for each of the 5 variables
    text_locations_y = c()
    for(var_index in top_5_variables){
      temp_top = differences_of_marginals_top5[differences_of_marginals_top5$Variable==var_index,]
      d <- density(temp_top$Jeff_single)
      text_locations_y = c(text_locations_y, 0.5*max(d$y))
    }
    
    dat_text = data.frame(
      labels    = rep("True First-Order", times = 5),
      vline     = text_locations_x,
      y_heights = text_locations_y,
      Variable  = top_5_variables_factor
    )
    
    p3 = ggplot(data = differences_of_marginals_top5, aes(x=Jeff_single)) + geom_density(alpha=.4, fill="#FF6666")
    p3 = p3 + facet_grid(rows = vars(Variable), scales="free")
    p3 = p3 + ggtitle(paste("First-Order Jeffrey's Divergence Loss,", graph_title))
    p3 = p3 + geom_vline(data = data_vline,
                         aes(xintercept = vline))
    p3 = p3 + geom_text(data = dat_text, aes(x=vline, label=labels, y=y_heights), colour="purple", angle=90,  size=3)
    p3 = p3 + ylab("Density") + xlab("First-Order Jeffrey's Divergence") + theme(text = element_text(size = 13))
    
    ggsave(
      filename = "Posterior_Graphs/posterior_Jeff_single.pdf",
      plot = p3
    )
    
  }
  
  prob_plot_list = list(p1,p3,p2)
  #### Create a graphic with all of the plots
  ggsave(
    filename = "Posterior_Graphs/posterior_Jeff_w_Marginals.pdf",
    plot = marrangeGrob(prob_plot_list, nrow=1, ncol=3, widths = c(1,1,1.5)),
    width = 12, height = 9
  )
  
  
}


#' Title
#'
#' @param mcmc_kernel 
#' @param changepoint_probs 
#' @param prob_cutoff 
#' @param output_file_name 
#'
#' @return
#' @export
#'
#' @examples
graph_marginals = function(mcmc_kernel, changepoint_probs, 
                               prob_cutoff, output_file_name = "Posterior_Graphs/data_over_time_w_changepoints.pdf"){
  ###############################################################################
  # Start by making a state vector that mirrors thet set of expected changepoints.
  num_states   = length(changepoint_probs)+1
  changepoints = changepoint_probs >= prob_cutoff
  my_states    = rep(0, times = num_states)
  count        = 1
  my_states[1] = count
  
  for(index in 1:length(changepoint_probs)){
    if(changepoints[index]){
      count            = count + 1
    }
    my_states[index+1] = count
  }
  
  latent_data          = mcmc_kernel$latent_data
  hyperparameters      = mcmc_kernel$hyperparameters
  p                    = hyperparameters$p
  Z_timepoint_indices  = mcmc_kernel$Z_timepoint_indices
  
  ###############################################################################
  # Gather the data for the graph.
  graph_data = data.frame(AverageValue = NA, StateNumber = NA, VariableNumber = NA)
  for(state_index in 1:num_states){
    Z_index_of_regime_start   = Z_timepoint_indices[[state_index]]$timepoint_first_index
    Z_index_of_regime_end     = Z_timepoint_indices[[state_index]]$timepoint_last_index
    current_data              = latent_data[Z_index_of_regime_start:Z_index_of_regime_end,]
    current_n                 = nrow(current_data)
    temp_row                  = data.frame(AverageValue = colMeans(current_data), 
                                           StateNumber = rep(state_index, times = p), 
                                           VariableNumber = 1:p)
    graph_data                = rbind(graph_data, temp_row)
  }
  changepoint_indices       = which(changepoints)+0.5
  graph_data                = graph_data[-1,]
  graph_data$VariableNumber = as.factor(graph_data$VariableNumber)
  
  ###############################################################################
  # Graph the data.
  mygg = ggplot(graph_data, aes(StateNumber, AverageValue, linetype=VariableNumber))+
    geom_line() + geom_vline(xintercept=c(changepoint_indices), color=rep("magenta",times=length(changepoint_indices)))
  ggsave(
    filename = output_file_name,
    plot = mygg
  )
}
