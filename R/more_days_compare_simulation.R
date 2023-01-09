# Simulation study of our method verses a T^2 control chart.
## Author: Alexander Murph
library(rlist)
library(ggplot2)
library(BDgraph)
library(LaplacesDemon)
library(gridExtra)
library(Hotelling)
library(boot)
library(gbm)
library(pROC)
library(forecast)
source("find_p_cutoff_functions.R")

job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
print(job_num)
if(is.na(job_num)){
  job_num = 2
}
# if(job_num == 2){
#   stop("The p isn't ready yet.")
# }


create_S2_w_diff_missing_struct = function(p, sigma){
  ## p is the total number of columns.
  ## This method takes in current covariance, sigma, and changes the last 5 columns.  
  
  shur_complement = rinvwishart(2+p, diag(5))
  undo_complement = shur_complement + sigma[(p-4):p,1:(p-5)]%*%
                      solve(sigma[1:(p-5),1:(p-5)])%*%
                      sigma[1:(p-5),(p-4):p]
  new_sigma = sigma
  new_sigma[(p-4):p,(p-4):p] = 2*undo_complement
  return(new_sigma)
}


# I think that I'd like to run this script 2600 times.
# 200 times for each of the 13 categories.

# Parameters for the whole script
parameters_number         = job_num%%8+1

if(parameters_number == 7){
  stop("I'm not doing linear drift anymore.")
}

if(parameters_number == 1){
  p                       = 15
  orig_p = 10
  # data.sim.mixed_group1 <- bdgraph.sim( n = rpois(1,100), p = p,
  #                                       graph = "cluster")
  # graph_structure       = data.sim.mixed_group1$G
  # save(graph_structure, file="graph_structure_missStruct.Rdata")
  # 
  # data.sim.mixed_group1 <- bdgraph.sim( n = rpois(1,100), p = p,
  #                                       graph = data.sim.mixed_group1$G)
  # 
  # S1 = data.sim.mixed_group1$sigma
  # save(S1, file="S1_missStruct.Rdata")
  
  load("S1_missStruct.Rdata")
  load("graph_structure_missStruct.Rdata")
  # S1              = S1_missStruct
  # graph_structure = graph_structure_missStruct
  
} else {
  p                       = 10
  orig_p = p
  load("S1.Rdata")
  load("S2.Rdata")
  load("S3.Rdata")
  load("S4.Rdata")
  load("graph_structure.Rdata")
}



if( ((parameters_number>1)&(parameters_number<4))|(parameters_number==6) ){
  discrete_vals             = c()
  has_discretes = 0
} else {
  if(parameters_number==1){
    discrete_vals = c(11,12,13,14,15)
    has_discretes = 1
  } else {
    discrete_vals = c(1,5,6,7,8)
    has_discretes = 1
    # discrete_vals             = c()
    # has_discretes = 0
  }
}

# It is recommended to use at least 30 cores when redrawing the pcutoffs.
redraw_p_cutoffs          = TRUE
p_cutoff_iterations       = 900
num_of_iterations         = 900
# mean_shifters             = c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75)
mean_shift                = 0.24
prob_cutoff               = 0.1
percent_missing_index     = parameters_number
base_data_size            = 1500
half_dist_between_modes   = mean_shift/2
num_of_days               = 51


possible_missings_percent = c(20,0,20,20,20,0,20,20)
# missing_nonmissing        = c(1,0,1,1,1,1,1,1)
missing_nonmissing        = c(1,0,1,1,1,0,1,1)
distribution_change       = c('miss','mean','cov','cov','cov','mean','mean', 'mean')
modal_type                = c('unimodal','unimodal','unimodal','unimodal','bimodal','bimodal','unimodal','unimodal')
MAR_list                  = c(0,1,1,1,1,1,1,1)
MAR                       = MAR_list[parameters_number]

should_include_missings   = missing_nonmissing[parameters_number]

if(redraw_p_cutoffs){
  set.seed(job_num + 400)
} else {
  set.seed(job_num)
}

orig_job_num = job_num
port_number = (job_num-1)*31 + 11001

not.cont = rep(0, times = p)

# data.sim.mixed_group1 <- bdgraph.sim( n = rpois(1,100), p = p,
#                                       graph = "cluster")
# graph_structure       = data.sim.mixed_group1$G
# save(graph_structure, file="graph_structure.Rdata")
# 
# data.sim.mixed_group1 <- bdgraph.sim( n = rpois(1,100), p = p,
#                                       graph = data.sim.mixed_group1$G)
# data.sim.mixed_group2 <- bdgraph.sim( n = rpois(1,100), p = p,
#                                       graph = data.sim.mixed_group1$G)
# data.sim.mixed_group3 <- bdgraph.sim( n = rpois(1,100), p = p,
#                                       graph = data.sim.mixed_group1$G)
# data.sim.mixed_group4 <- bdgraph.sim( n = rpois(1,100), p = p,
#                                       graph = data.sim.mixed_group1$G)
# 
# S1 = data.sim.mixed_group1$sigma
# save(S1, file="S1.Rdata")
# S2 = data.sim.mixed_group2$sigma
# save(S2, file="S2.Rdata")
# S3 = data.sim.mixed_group3$sigma
# save(S3, file="S3.Rdata")
# S4 = data.sim.mixed_group4$sigma
# save(S4, file="S4.Rdata")


if(distribution_change[parameters_number] == 'mean') {
  
  S1    = S1*0.5
  S2    = S1
  S3    = S1
  mean1 = rep(0, p)
  # Note -- I've done this change when I make the data, so that it's a change that might ONLY be picked up by the mixture model.
  ## Note further: the only simulation that uses the following is the 'diaboloical bimodal shift' example.  The rest are handled individually.
  mean2 = mean1 #+ mean_shift
  mean3 = mean1
  if(parameters_number%in%c(2,8)){
    mean1 = rep(0, p)
    mean2[3:4] = mean1[3:4] + mean_shift
    mean2[4]   = mean1[4] + mean_shift
    mean3[3:4] = mean1[3:4] + 2*mean_shift
    mean3[4]   = mean1[4] + 2*mean_shift
    print("mean 1 is:")
    print(mean1)
    print("mean 2 is:")
    print(mean2)
    print("mean 3 is:")
    print(mean3)
  } else if(parameters_number==6){
    S1    = S1 * 1.2
    S2    = S1 * 0.92
    S3    = S1 * 0.81
  }
  data_type = "MeanIsDifferent"
  if(parameters_number == 7){
    data_type = "LinearDrift"
  }
  
} else if(distribution_change[parameters_number] == 'cov') {
  data.sim.mixed_group2 <- bdgraph.sim( n = rpois(1,100), p = p,
                                        graph = graph_structure)

  S1      = S1 
  S2      = S1
  S2[3,3] = S2[3,3]*2
  S2[4,4] = S2[4,4]*1.5
  S3      = S1
  S3[3,3] = S3[3,3]*3.5
  S3[4,4] = S3[4,4]*3
  
  mean1   = rep(0, p)
  mean2   = mean1 #+ mean_shift
  mean3   = mean1
  data_type = "CovIsDifferent"
  
} else {
  # I might have to reconsider how I do this slightly.
  data.sim.mixed_group1 <- bdgraph.sim( n = rpois(1,100), p = p,
                                        graph = graph_structure)
  S1              = S1
  S2              = create_S2_w_diff_missing_struct(p, S1)
  S3              = S2
  
  # graph_structure = data.sim.mixed_group1$G
  mean1           = rep(0, p)
  mean2           = mean1
  mean3           = mean1
  mean3[11:15]    = mean3[11:15] + 2*mean_shift
  data_type       = "MissStructDiff"
}

day_dts = seq(as.POSIXct("2020-03-01 01:00:00", tz="CST6CDT"),as.POSIXct("2020-04-22 01:00:00", tz="CST6CDT"),by="1 day")

data_group_1 = NULL
data_group_2 = NULL
data_group_3 = NULL

full_data            = NULL
count                = 0
day_of_observations  = NULL
data_set_list        = list()
for(day in 1:(length(day_dts)-1)){
  count=count+1
  if(is.null(full_data)){
    
    if(modal_type[parameters_number]=='unimodal'){
      data.sim.mixed_group1 <- bdgraph.sim( n = 2*base_data_size, p = p,
                                            graph = "fixed", 
                                            sigma=S1, mean = mean1)
      full_data             = data.sim.mixed_group1$data
      temp_day_log          = full_data
      day_of_observations   = rep(day_dts[day+1], times = nrow(data.sim.mixed_group1$data))
      print(paste("means for day", day, "is:"))
      print(paste(apply(temp_day_log, 2, mean)))
      
    } else {
      data.sim.mixed_group1 <- bdgraph.sim( n = base_data_size, p = p,
                                            graph = "fixed", 
                                            sigma=S1, mean = mean1-half_dist_between_modes)
      full_data = data.sim.mixed_group1$data
      day_of_observations = rep(day_dts[day+1], times = nrow(data.sim.mixed_group1$data))
      temp_day_log = full_data
      
      data.sim.mixed_group3 <- bdgraph.sim( n = base_data_size, p = p,
                                            graph = "fixed", 
                                            sigma=S1, mean = mean1+half_dist_between_modes)
      full_data = rbind(full_data, data.sim.mixed_group3$data)
      day_of_observations =  c(day_of_observations,rep(day_dts[day+1], times = nrow(data.sim.mixed_group3$data)))
      
      temp_day_log = rbind(temp_day_log, data.sim.mixed_group3$data)
      print(paste("means for day", day, "is:"))
      print(paste(apply(temp_day_log, 2, mean)))
    }
    
    
    
  }else if(count < 15){
    
    if(modal_type[parameters_number]=='unimodal'){
      data.sim.mixed_group1 <- bdgraph.sim( n = 2*base_data_size, p = p, 
                                            graph = "fixed",
                                            sigma=S1, mean = mean1)
      full_data = rbind(full_data, data.sim.mixed_group1$data)
      day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group1$data)))
      
      temp_day_log = data.sim.mixed_group1$data
      print(paste("means for day", day, "is:"))
      print(paste(apply(temp_day_log, 2, mean)))
      
    } else {
      data.sim.mixed_group1 <- bdgraph.sim( n = base_data_size, p = p, 
                                            graph = "fixed",
                                            sigma=S1, mean = mean1-half_dist_between_modes)
      full_data = rbind(full_data, data.sim.mixed_group1$data)
      temp_day_log = data.sim.mixed_group1$data
      day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group1$data)))
      data.sim.mixed_group3 <- bdgraph.sim( n = base_data_size, p = p,
                                            graph = "fixed", 
                                            sigma=S1, mean = mean1+half_dist_between_modes)
      full_data = rbind(full_data, data.sim.mixed_group3$data)
      day_of_observations =  c(day_of_observations,rep(day_dts[day+1], times = nrow(data.sim.mixed_group3$data)))
      
      temp_day_log = rbind(temp_day_log, data.sim.mixed_group3$data)
      print(paste("means for day", day, "is:"))
      print(paste(apply(temp_day_log, 2, mean)))
      
    }
    
  } else if(count < 35) {
    if(count == 15){
      data_group_1 = full_data
    }
    
    if(modal_type[parameters_number]=='unimodal'){
      
      if(data_type == "LinearDrift"){
        if(distribution_change[parameters_number] == 'mean'){
          mean1 = mean1 + mean_shift
          mean2 = mean1
        } else {
          print(paste("count is currently",count))
          alpha = (count-15)/21
          S2 = solve(solve(S1)*(1-alpha) + solve(orig_S2)*(alpha))
        }
      }
      data.sim.mixed_group2 = bdgraph.sim( n = 2*base_data_size, p = p, 
                                            graph = "fixed",
                                            sigma=S2, mean = mean2)
      full_data             = rbind(full_data, data.sim.mixed_group2$data)
      day_of_observations   = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group2$data)))
      
      temp_day_log = data.sim.mixed_group2$data
      print(paste("means for day", day, "is:"))
      print(paste(apply(temp_day_log, 2, mean)))
      
      
    } else {
      if(data_type == "MeanIsDifferent"){
        print("We DID actually shift the means away from one other for the diabolical bimodal example.")
        bimodal_shift = mean_shift/2
      } else {
        bimodal_shift = 0
      }
      
      
      data.sim.mixed_group2 <- bdgraph.sim( n = base_data_size, p = p, 
                                            graph = "fixed",
                                            sigma=S2, mean = mean2-half_dist_between_modes - bimodal_shift)
      full_data = rbind(full_data, data.sim.mixed_group2$data)
      temp_day_log = data.sim.mixed_group2$data
      day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group2$data)))
      
      data_group_2 = rbind(data_group_2, data.sim.mixed_group2$data)
      
      data.sim.mixed_group4 <- bdgraph.sim( n = base_data_size, p = p, 
                                            graph = "fixed",
                                            sigma=S2, mean = mean2+half_dist_between_modes + bimodal_shift)
      full_data = rbind(full_data, data.sim.mixed_group4$data)
      day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group4$data)))
      
      data_group_2 = rbind(data_group_2, data.sim.mixed_group4$data)
      
      temp_day_log = rbind(temp_day_log, data.sim.mixed_group4$data)
    }
    # print(paste("means for day", day, "is:"))
    # print(paste(apply(temp_day_log, 2, mean)))
  } else {
    
    
    if(modal_type[parameters_number]=='unimodal'){
      
      if(data_type == "LinearDrift"){
        if((distribution_change[parameters_number] == 'mean')&&(count==35)){
          mean1 = mean1 + 0.2*mean_shift
          mean3 = mean1
        } else if (distribution_change[parameters_number] == 'mean') {
          mean3 = mean1
        } else {
          S3 = orig_S2
        }
      }
      data.sim.mixed_group3 <- bdgraph.sim( n = 2*base_data_size, p = p, 
                                            graph = "fixed",
                                            sigma=S3, mean = mean3)
      full_data = rbind(full_data, data.sim.mixed_group3$data)
      day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group3$data)))
      
      temp_day_log = data.sim.mixed_group3$data
      print(paste("means for day", day, "is:"))
      print(paste(apply(temp_day_log, 2, mean)))
      
      
    } else {
      if(data_type == "MeanIsDifferent"){
        print("We DID actually shift the means away from one other for the diabolical bimodal example (next regime change).")
        bimodal_shift = mean_shift
      } else {
        bimodal_shift = 0
      }
      data.sim.mixed_group3 <- bdgraph.sim( n = base_data_size, p = p, 
                                            graph = "fixed",
                                            sigma=S3, mean = mean3-half_dist_between_modes - bimodal_shift)
      full_data = rbind(full_data, data.sim.mixed_group3$data)
      data_group_3 = rbind(data_group_3, data.sim.mixed_group3$data)
      temp_day_log = data.sim.mixed_group3$data
      day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group3$data)))
      
      data.sim.mixed_group4 <- bdgraph.sim( n = base_data_size, p = p, 
                                            graph = "fixed",
                                            sigma=S3, mean = mean3+half_dist_between_modes + bimodal_shift)
      full_data = rbind(full_data, data.sim.mixed_group4$data)
      data_group_3 = rbind(data_group_3, data.sim.mixed_group4$data)
      day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group4$data)))
      
      
      temp_day_log = rbind(temp_day_log, data.sim.mixed_group4$data)
      
      # print(paste("means for day", day, "is:"))
      # print(paste(apply(temp_day_log, 2, mean)))
    }
  }
}

# WLOG, let's assume that the discrete values are all in the FINAL few columns.
all_indices = 1:p
not.cont = as.numeric(1:p %in% discrete_vals)
# Transform some of the variables into binary.

for(discrete_index in all_indices[which(not.cont == 1)]){
  while_loop_count = 0
  while(TRUE){
    while_loop_count = while_loop_count + 1
    cutoff_value = runif(1, -3,3)
    if(sum((full_data[, discrete_index] >= cutoff_value)) >= (nrow(full_data)-300) ){
      print("CREATED DISCRETE DATA OF TOO MANY 1s")
    } else if(sum((full_data[, discrete_index] < cutoff_value)) >= (nrow(full_data)-300) ){
      print("CREATED DISCRETE DATA OF TOO MANY 0s")
    } else if( (mean((full_data[, discrete_index] >= cutoff_value)) > 0.3) & (discrete_index > 10) ){
      print("MNAR PROCESS CREATED TOO MANY MISSING")
    } else {
      full_data[,discrete_index] = (full_data[,discrete_index] >= cutoff_value)
      break
    }
  }
  if(while_loop_count>1000){
    error("WE COULD NOT CREATE THE DISCRETE VARIABLES")
  }
}

# Mean center non-discrete rows.  (dirty move heh)
for(discrete_index in all_indices[which(not.cont == 0)]){
  full_data[,discrete_index] = full_data[,discrete_index] - mean(full_data[,discrete_index])
}

# Add in some missing values.
count = 0
missing_nums              = floor((possible_missings_percent[percent_missing_index]/100)*nrow(full_data))
# 
if(MAR){
  ## Put in missingness, following a simple Missing At Random structure.
  ## I assume that the discrete values are always non-missing.
  continuous_col_numbers = 1:ncol(full_data)
  continuous_col_numbers = continuous_col_numbers[which(not.cont==0)]
  list_of_cols_w_missing = c()
  if(should_include_missings){
    while(count < missing_nums){
      random_i               = sample(1:nrow(full_data), size = 1)
      random_j               = sample(1:ncol(full_data), size = 1)
      list_of_cols_w_missing = c(list_of_cols_w_missing, random_j)
      if(!is.na(full_data[random_i,random_j])){
        full_data[random_i,random_j] = NA
        count = count + 1
      }
    }
    list_of_cols_w_missing = unlist(sort(unique(list_of_cols_w_missing)))
    new_full_data = matrix(0,nrow=nrow(full_data),ncol = (ncol(full_data)+length(list_of_cols_w_missing)) )
    new_full_data[,1:p] = full_data
    
    count = 1
    for(column in list_of_cols_w_missing){
      # 
      new_col = matrix(as.numeric(is.na(full_data[,column]),ncol=1))
      new_full_data[,p+count] = new_col
      count = count + 1
    }
    full_data = new_full_data
    
    old_p                                = p
    p                                    = ncol(full_data)
    new_graph_structure                  = matrix(0,p,p)
    new_graph_structure[1:old_p,1:old_p] = graph_structure
    graph_structure                      = new_graph_structure
    
    # 
  }
} else {
  ## Otherwise, let's attempt to add a noteable structure to the missingness pattern itself.
  ### I'll assume n_cols_missing are the columns for which missingness is present.
  ### Thus, the last n_cols_missing are discrete variables, and the next-to-last n_cols_missing
  ### colums are the ONLY columns in which missingness is present.
  print("Doing missing NOT at random")
  # 
  # Set up the missingness so that the last 5 columns represent the next to last 5 column's missingness.
  for(j in 11:15){
    for(i in 1:nrow(full_data)){
      if(full_data[i,j]){
        full_data[i,j-5] = NA
      }
    }
  }
  
}

### Fill the data into a list for use by the T2 Control Chart
lower_index   = 1
upper_index   = 2*base_data_size
data_set_list = list()
for(day in 1:(length(day_dts)-1)){
  data_set_list[[day]] = full_data[lower_index:upper_index,]
  lower_index   = lower_index + 2*base_data_size
  upper_index   = upper_index + 2*base_data_size
}

# Perform Box-Cox Transformations on the variables as a pre-processing step
# all_lambdas = c(0.9225489, 0.9633708, 0.9767880, 0.9522332, 1.0148248, 1.0311973, 1.0168315,
#                 1.0531853, 0.9962194, 0.9801962)
# if(parameters_number!=1){
  all_lambdas = c()
  for(transform_index in 1:length(not.cont)){
    if(!as.logical(not.cont[transform_index])){
      # Then this variable is continuous, and I should try to do an automated box-cox transform:
      print(paste("Starting column:", transform_index))
      temp_data_col                 = full_data[,transform_index]
      minimum_value                 = min(temp_data_col,na.rm=T)
      temp_data_col                 = temp_data_col - minimum_value + 1
      temp_lambda                   = BoxCox.lambda(temp_data_col, method = "guerrero")
      # temp_lambda                   = all_lambdas[transform_index]
      full_data[,transform_index]   = BoxCox(temp_data_col, temp_lambda)
      full_data[,transform_index]   = full_data[,transform_index] + minimum_value - 1
      print(paste("Finished column:", transform_index))
      all_lambdas = c(all_lambdas, temp_lambda)
    }
  }
  print(all_lambdas)
# }

detach("package:BDgraph") 
num_of_nodes = Sys.getenv("NSLOTS")
print(num_of_nodes)
if(is.na(as.numeric(num_of_nodes))){
  print("not working on cluster!")
  num_of_nodes = 30
}

p = ncol(full_data)
source("fit_w_reversible_jump.R")
data_name         = paste("MIXTUREMODEL", data_type,orig_job_num, modal_type[parameters_number],
                           distribution_change[parameters_number],"includesMissing",should_include_missings, sep="_")
cuttoff_file_name = paste("P_Cutoffs/mixture_p_for", parameters_number, ".RData", sep = "")

if((orig_job_num < 9)&redraw_p_cutoffs){
  yy = fit_hmm(full_data, day_of_observations, day_dts, not.cont = not.cont, #set_G = graph_structure,
               iterations = num_of_iterations, n.cores = as.integer(num_of_nodes),
               allow_for_mixture_models = TRUE, model_params_save_every = 5, simulation_iter = orig_job_num,
               determining_p_cutoff = TRUE)
  list.save(yy, paste("Simulation_Model_Logs/", data_name,"_MIXTUREMODEL", ".rds", sep = ""))
  # yy = list.load(paste("Simulation_Model_Logs/", data_name,"_MIXTUREMODEL", ".rds", sep = ""))
  sim_hyperparameters = yy$hyperparameters
  simulation_hyperparameters = list(not.cont = not.cont, day_dts = day_dts,
                                    day_of_observations = day_of_observations,
                                    allow_for_mixture_models = TRUE,
                                    num_of_iterations = p_cutoff_iterations,
                                    set_G = sim_hyperparameters$current_G)
  p_cutoff = find_p_cutoff(orig_job_num, simulation_hyperparameters, n.cores = as.integer(num_of_nodes))
  save(p_cutoff, file = cuttoff_file_name)
} else if( (parameters_number==3)|(parameters_number==2)) {
  yy = fit_hmm(full_data, day_of_observations, day_dts, not.cont = not.cont, #set_G = graph_structure,
               iterations = num_of_iterations, n.cores = as.integer(num_of_nodes),
               allow_for_mixture_models = TRUE, model_params_save_every = 5, simulation_iter = orig_job_num)
  list.save(yy, paste("Simulation_Model_Logs/", data_name,"_MIXTUREMODEL", ".rds", sep = ""))
} else {
  yy = fit_hmm(full_data, day_of_observations, day_dts, not.cont = not.cont, #set_G = graph_structure,
               iterations = num_of_iterations, n.cores = as.integer(num_of_nodes),
               allow_for_mixture_models = TRUE, model_params_save_every = 1e9)
  list.save(yy, paste("Simulation_Model_Logs/", data_name,"_MIXTUREMODEL", ".rds", sep = ""))
}

load(cuttoff_file_name)
if(p_cutoff == 0){
  p_cutoff = 1e-5
}
prob_cutoff = p_cutoff
# yy = list.load(paste("Simulation_Model_Logs/", data_name,"_MIXTUREMODEL", ".rds", sep = ""))

### Grab the changepoint probabilities
mixture_model_mergesplit_accepts = yy$mergesplit_accepts

my_states = yy$states
states_df = data.frame(matrix(my_states[[1]], nrow=1))
for(i in 2:length(my_states)){
  states_df = rbind(states_df, data.frame(matrix(my_states[[i]], nrow=1)))
}
states_df_burnt = states_df[floor(2*nrow(states_df)/4):nrow(states_df),]
prop_mat = matrix(0, nrow=1, ncol=(ncol(states_df_burnt)-1))
number_of_segments = nrow(states_df)
for(j in 1:(ncol(states_df_burnt)-1)){
  temp_difference = states_df_burnt[,j+1] - states_df_burnt[,j]
  prop_mat[1,j] = mean(temp_difference==1)
}
changepoint_probs   = prop_mat

if(parameters_number == 7){
  changepoints = c(0,0)
  for(j in 3:(ncol(states_df_burnt)-1)){
    changepoints = c(changepoints, any(changepoint_probs[(j-2):j] >= prob_cutoff))
  }
  changepoints        = data.frame(matrix(changepoints, nrow = 1, ncol = length(changepoints)))
} else {
  changepoints        = changepoint_probs >= prob_cutoff
  changepoints        = data.frame(matrix(changepoints, nrow = 1, ncol = length(changepoints)))
}
#

changepoints$Iter   = orig_job_num
changepoints$Method = "MIXTUREMODEL"
changepoints$ShiftType           = data_type
changepoints$has_missing         = should_include_missings
changepoints$distribution_change = distribution_change[parameters_number]
changepoints$any_discretes       = has_discretes
changepoints$modal_type          = modal_type[parameters_number]

write.csv(changepoints, paste("Simulation_Data/", data_name, ".csv", sep = ""))

if( ((parameters_number==3)|(parameters_number==2))&(!redraw_p_cutoffs) ){
  source("posterior_analysis.R")
  fault_data        = data.frame(matrix(0,ncol=7,nrow=3))
  names(fault_data) = c(as.character(1:5),"Method","Iteration")
  p                 = ncol(full_data)
  posterior_data    = determine_marginals(changepoint_probs, prob_cutoff, as.character(1:p),
                                       simulation_number = orig_job_num)

  if(!anyNA(posterior_data)){
    # From the Hellinger Distance data, find the top-5 variables for Total-Effect Hellinger Distance Loss
    average_values = c()
    differences_of_marginals = posterior_data
    for(var_index in 1:length(unique(differences_of_marginals$Variable))){
      temp_data = differences_of_marginals[which(differences_of_marginals$Variable == var_index),]
      average_values = c(average_values, mean(temp_data$Hellinger))
    }
    possible_variables      = 1:length(unique(differences_of_marginals$Variable))
    ordered_variables       = possible_variables[order(average_values, decreasing=TRUE)]
    top_5_variables         = ordered_variables[1:5]
    fault_data[1,1:5]       = top_5_variables
    fault_data$Method[1]    = "Mixture_Model"
    fault_data$Iteration[1] = orig_job_num
  } else {
    fault_data[1,1:5]       = rep(NA, times = 5)
    fault_data$Method[1]    = "Mixture_Model"
    fault_data$Iteration[1] = orig_job_num
  }
} 

data_name         = paste("SINGLEMVNMODEL",data_type,orig_job_num, modal_type[parameters_number],
                             distribution_change[parameters_number],"includesMissing",should_include_missings, sep="_")
cuttoff_file_name = paste("P_Cutoffs/mvn_p_for", parameters_number, ".RData", sep = "")

## NOTE: component_trunction MUST be set to 1.
if((orig_job_num < 9)&(redraw_p_cutoffs)){
  yy = fit_hmm(full_data, day_of_observations, day_dts, not.cont = not.cont,
               iterations = num_of_iterations, n.cores = as.integer(num_of_nodes), #set_G = graph_structure, #port_number = port_number,
               allow_for_mixture_models = FALSE, model_params_save_every = 5, simulation_iter = (orig_job_num+8),
               component_truncation = 1, determining_p_cutoff = TRUE)
  list.save(yy, paste("Simulation_Model_Logs/", data_name,"_SINGLEMODEL", ".rds", sep = ""))
  sim_hyperparameters = yy$hyperparameters
  simulation_hyperparameters = list(not.cont = not.cont, day_dts = day_dts,
                                    day_of_observations = day_of_observations,
                                    allow_for_mixture_models = FALSE,
                                    num_of_iterations = p_cutoff_iterations,
                                    set_G = sim_hyperparameters$current_G)
  p_cutoff = find_p_cutoff((orig_job_num+8), simulation_hyperparameters, n.cores = as.integer(num_of_nodes))
  save(p_cutoff, file = cuttoff_file_name)
} else if((parameters_number == 3)|(parameters_number == 2)) {
  yy = fit_hmm(full_data, day_of_observations, day_dts, not.cont = not.cont,
               iterations = num_of_iterations, n.cores = as.integer(num_of_nodes), #set_G = graph_structure, #port_number = port_number,
               allow_for_mixture_models = FALSE, model_params_save_every = 5,  simulation_iter = orig_job_num)
  list.save(yy, paste("Simulation_Model_Logs/", data_name,"_SINGLEMODEL", ".rds", sep = ""))

} else {
  yy = fit_hmm(full_data, day_of_observations, day_dts, not.cont = not.cont,
               iterations = num_of_iterations, n.cores = as.integer(num_of_nodes), #set_G = graph_structure, #port_number = port_number,
               allow_for_mixture_models = FALSE, model_params_save_every = 1e9,
               component_truncation = 1)
  list.save(yy, paste("Simulation_Model_Logs/", data_name,"_SINGLEMODEL", ".rds", sep = ""))
}
load(cuttoff_file_name)
if(p_cutoff == 0){
  p_cutoff = 1e-5
}
prob_cutoff = p_cutoff
# yy = list.load(paste("Simulation_Model_Logs/", data_name,"_SINGLEMODEL", ".rds", sep = ""))


### Grab the changepoint probabilities
MVN_model_mergesplit_accepts = yy$mergesplit_accepts
my_states = yy$states
states_df = data.frame(matrix(my_states[[1]], nrow=1))
for(i in 2:length(my_states)){
  states_df = rbind(states_df, data.frame(matrix(my_states[[i]], nrow=1)))
}
states_df_burnt = states_df[floor(2*nrow(states_df)/4):nrow(states_df),]
prop_mat = matrix(0, nrow=1, ncol=(ncol(states_df_burnt)-1))
number_of_segments = nrow(states_df)
for(j in 1:(ncol(states_df_burnt)-1)){
  temp_difference = states_df_burnt[,j+1] - states_df_burnt[,j]
  prop_mat[1,j] = mean(temp_difference==1)
}
changepoint_probs   = prop_mat

if(parameters_number == 7){
  changepoints = c(0,0)
  for(j in 3:(ncol(states_df_burnt)-1)){
    changepoints = c(changepoints, any(changepoint_probs[(j-2):j] >= prob_cutoff))
  }
  changepoints        = data.frame(matrix(changepoints, nrow = 1, ncol = length(changepoints)))
} else {
  changepoints        = changepoint_probs >= prob_cutoff
  changepoints        = data.frame(matrix(changepoints, nrow = 1, ncol = length(changepoints)))
}

changepoints$Iter                = orig_job_num
changepoints$Method              = "SINGLEMVNMODEL"
changepoints$ShiftType           = data_type
changepoints$has_missing         = should_include_missings
changepoints$distribution_change = distribution_change[parameters_number]
changepoints$any_discretes       = has_discretes
changepoints$modal_type          = modal_type[parameters_number]
write.csv(changepoints, paste("Simulation_Data/", data_name, ".csv", sep = ""))

# For simulations 3 and 8, we perform fault detection.
if( ((parameters_number==3)|(parameters_number==2))&(!redraw_p_cutoffs) ){
  posterior_data = determine_marginals(changepoint_probs, prob_cutoff, as.character(1:p),
                                       simulation_number = orig_job_num)
  if(!is.na(posterior_data)){
    # From the Hellinger Distance data, find the top-5 variables for Total-Effect Hellinger Distance Loss
    average_values = c()
    differences_of_marginals = posterior_data
    for(var_index in 1:length(unique(differences_of_marginals$Variable))){
      temp_data = differences_of_marginals[which(differences_of_marginals$Variable == var_index),]
      average_values = c(average_values, mean(temp_data$Hellinger))
    }
    possible_variables      = 1:length(unique(differences_of_marginals$Variable))
    ordered_variables       = possible_variables[order(average_values, decreasing=TRUE)]
    top_5_variables         = ordered_variables[1:5]
    fault_data[2,1:5]       = top_5_variables
    fault_data$Method[2]    = "Single_Model"
    fault_data$Iteration[2] = orig_job_num
  } else {
    fault_data[2,1:5]       = rep(NA, times = 5)
    fault_data$Method[2]    = "Single_Model"
    fault_data$Iteration[2] = orig_job_num
  }
}
  

# LASTLY, WE DO THE T^2 CONTROL CHARTS

# Now use this to make a control chart:
T2_values_twosamples     = c()
T2_UCLs                  = c()
GLR_values               = c()
chisq_values             = c()
all_p_vals               = c()
full_covariance          = cov(full_data)
n_full                   = nrow(full_data)
changepoint_detected     = c(0,0)
changepoint_detected_GBM = c(0,0)

changepoint_detected_lepage = c(0,0)
#3:13
for(t in 3:(num_of_days-1)){
  data_before = NULL
  data_after = NULL
  for(index in 1:(num_of_days)){
    if( (index <= t)&( (index > (t-3))) ){
      data_before = rbind(data_before, data_set_list[[index]])
    } else if(index == (t+1)) {
      data_after = rbind(data_after, data_set_list[[index]])
    }
  }

  # Calculate the Hotelling T^2 Scan Statistic
  cov1  = cov(data_before[,1:orig_p], use = "complete.obs")
  cov2  = cov(data_after[,1:orig_p], use = "complete.obs")
  xbar1 = as.matrix(colMeans(data_before[,1:orig_p], na.rm = TRUE))
  xbar2 = as.matrix(colMeans(data_after[,1:orig_p], na.rm = TRUE))
  n1    = nrow(data_before)
  n2    = nrow(data_after)

  x_stuff = list(mean=xbar1,cov=cov1,n=n1)
  y_stuff = list(mean=xbar2,cov=cov2,n=n2)
  
  out = tryCatch({
    x       = hotelling.test(x_stuff,y_stuff, var.equal=FALSE)
    if( (det(x_stuff$cov)==0)|(det(y_stuff$cov)==0) ){
      print("one of the cov matrices is precisely singular.")
      p_value = 1
    } else {
      p_value = x$pval
    }
    p_value
  }, error=function(cond) {
    print("-----------------------------------")
    print("-----------------------------------")
    message(paste("Hotelling Calculation failed!"))
    print("-----------------------------------")
    print("-----------------------------------")
    return(1)
  })

  print(paste("the return value is: ", out, ", for changepoint at: ", t, sep = ""))
  changepoint_detected        = c(changepoint_detected, as.numeric(out <= 0.05))
  all_p_vals                  = c(all_p_vals, out)
  changepoint_detected_lepage = c(changepoint_detected_lepage, 0)

  print("-------------------------")
  print(paste("Info about data for changepoint at:", t))
  my_func_2 = function(x){var(x,na.rm=T)}
  my_func_1 = function(x){mean(x,na.rm=T)}
  my_func_3 = function(x){mean(is.na(x))}
  mean_before = apply(data_before, 2, my_func_1)
  var_before = diag(cov1)
  prop_na_before = apply(data_before, 2, my_func_3)
  mean_after = apply(data_after, 2, my_func_1)
  var_after = diag(cov2)
  prop_na_after = apply(data_after, 2, my_func_3)
  temp_mat = matrix(c(mean_before, mean_after, prop_na_before, var_before, var_after, prop_na_after), nrow=6, ncol =p, byrow = T)
  rownames(temp_mat) = c("mean_before", "mean_after", "prop_na_before", "var_before", "var_after", "prop_na_after")
  print(temp_mat)

  ##### PATRICK THIS IS WHERE I DO THE RTC CODE!!!!!
  # Calculate the Real-Time Contrasts with a GBM
  sliding_window_data       = data.frame(rbind(data_before, data_after))
  #
  sliding_window_classifier = c(rep(0, times = nrow(data_before)),rep(1,times = nrow(data_after)))
  sliding_window_data$class = sliding_window_classifier
  cv_assignments            = rep(1:10, times = ceiling(nrow(sliding_window_data)/10))
  cv_assignments            = sample(cv_assignments)

  data_for_roc = data.frame(gbm_pred = rep(NA,times = nrow(sliding_window_data)),
                            true_class = rep(NA,times = nrow(sliding_window_data)) )

  for(cv_index in 1:10){
    indicies_in  = (cv_assignments==cv_index)
    indicies_out = !indicies_in
    train_data   = sliding_window_data[indicies_in,]
    test_data    = sliding_window_data[indicies_out,]

    ### Build the GBM and get output
    gbm_model = gbm(class~., distribution="bernoulli", data = train_data)
    gbm_predictions_for_out_data = predict(gbm_model, newdata = test_data, type = "response")
    ###

    data_for_roc$true_class[indicies_out] = sliding_window_classifier[indicies_out]
    data_for_roc$gbm_pred[indicies_out]   = gbm_predictions_for_out_data
  }


  ### PATRICK NOTICE HOW I BOOTSTRAP HERE TO GET AN EMPIRICAL DISTRIBUTION ON THE AUC!
  ### With the ROC info, bootstrap a distribution on AUC.
  n_bootstraps = 1500
  aucs         = rep(0,times = n_bootstraps)
  all_possible = 1:nrow(sliding_window_data)
  for(boot_index in 1:n_bootstraps){
    boot_indicies = sample(all_possible,nrow(sliding_window_data),replace = TRUE)
    boot_roc_data = data_for_roc[boot_indicies,]
    ### Calculate the AUC for this subset
    boot_auc = as.numeric(auc(boot_roc_data$true_class, boot_roc_data$gbm_pred, direction = "<"))
    ###
    aucs[boot_index] = boot_auc
  }


  lower025_auc = quantile(aucs,p=0.02)
  upper975_auc = quantile(aucs,p=0.98)
  changepoint_detected_GBM = c(changepoint_detected_GBM, as.numeric((lower025_auc>0.5)|(upper975_auc<0.5)) )
}

# For simulations 3 and 8, we perform fault detection.
if( ((parameters_number==3)|(parameters_number==2))&(!redraw_p_cutoffs) ){

  if(sum(changepoint_detected_GBM)==0){
    fault_data[3,1:5]       = rep(NA,times = 5)
    fault_data$Method[3]    = "GBM"
    fault_data$Iteration[3] = orig_job_num
  } else {
    # Grab the index of the last changepoint.
    t = which(changepoint_detected_GBM==1)
    t = t[length(t)]

    data_before = NULL
    data_after  = NULL
    for(index in 1:num_of_days){
      if( (index <= t)&( (index > (t-3))) ){
        data_before = rbind(data_before, data_set_list[[index]])
      } else if(index == (t+1)) {
        data_after = rbind(data_after, data_set_list[[index]])
      }
    }

    # Calculate the Real-Time Contrasts with a GBM
    sliding_window_data       = data.frame(rbind(data_before, data_after))
    sliding_window_classifier = c(rep(0, times = nrow(data_before)),rep(1,times = nrow(data_after)))
    sliding_window_data$class = sliding_window_classifier
    cv_assignments            = rep(1:10, times = ceiling(nrow(sliding_window_data)/10))
    cv_assignments            = sample(cv_assignments)

    data_for_roc = data.frame(gbm_pred = rep(NA,times = nrow(sliding_window_data)),
                              true_class = rep(NA,times = nrow(sliding_window_data)) )
    curr_max_auc = 0
    for(cv_index in 1:10){
      indicies_in  = (cv_assignments==cv_index)
      indicies_out = !indicies_in
      train_data   = sliding_window_data[indicies_in,]
      test_data    = sliding_window_data[indicies_out,]

      ### Build the GBM and get output
      gbm_model = gbm(class~., distribution="bernoulli", data = train_data)
      gbm_predictions_for_out_data = predict(gbm_model, newdata = test_data, type = "response")
      temp_auc = as.numeric(auc(sliding_window_classifier[indicies_out], gbm_predictions_for_out_data, direction = "<"))
      ###
      if(temp_auc>curr_max_auc){
        curr_max_auc = temp_auc
        curr_best_cv = cv_index
      }

      data_for_roc$true_class[indicies_out] = sliding_window_classifier[indicies_out]
      data_for_roc$gbm_pred[indicies_out]   = gbm_predictions_for_out_data
    }

    # Use the best CV auc to look at fault detection:
    indicies_in  = (cv_assignments==curr_best_cv)
    indicies_out = !indicies_in
    train_data   = sliding_window_data[indicies_in,]
    test_data    = sliding_window_data[indicies_out,]
    gbm_model = gbm(class~., distribution="bernoulli", data = train_data)
    gbm_predictions_for_out_data = predict(gbm_model, newdata = test_data, type = "response")
    xx = summary.gbm(gbm_model)
    top_5_vars = xx$var[1:5]
    top_5_variables = c()
    for(variable_num in 1:5){
      individual_characters = strsplit(top_5_vars[variable_num], split = "")[[1]]
      just_the_numbers      = paste(individual_characters[2:length(individual_characters)], collapse="")
      top_5_variables       = c(top_5_variables,as.numeric(just_the_numbers))
    }

    # Record the Real Time Monitoring Data.
    fault_data[3,1:5]       = top_5_variables
    fault_data$Method[3]    = "GBM"
    fault_data$Iteration[3] = orig_job_num
  }
  # Save all of the data for this iteration of the simulation
  data_name = paste("Fault_Detection_Results/","faults_for_iter",orig_job_num,".csv",sep="_")
  write.csv(fault_data,file = data_name)
}

# Record and save the T2 data.
# break_point = rep(3:13, times = 2)
# type        = rep(c("Observed T2", "95% UCL"), each = 11)
# value       = c(T2_values_twosamples, T2_UCLs)
# graph_data  = data.frame(break_point = break_point, type = type, value = value)
graph_data        = data.frame(matrix(changepoint_detected, nrow=1, ncol = length(changepoint_detected)))
graph_data$Iter   = orig_job_num
graph_data$Method = "T2"
graph_data$ShiftType           = data_type
graph_data$has_missing         = should_include_missings
graph_data$distribution_change = distribution_change[parameters_number]
graph_data$any_discretes       = has_discretes
graph_data$modal_type          = modal_type[parameters_number]
#
data_name = paste("T2",data_type,orig_job_num, modal_type[parameters_number],distribution_change[parameters_number],"includesMissing",should_include_missings, sep="_")
print(data_name)
write.csv(graph_data, paste("Simulation_Data/", data_name, ".csv", sep = ""))


# Record and save the Real Time Monitoring Data.
graph_data        = data.frame(matrix(changepoint_detected_GBM, nrow=1, ncol = length(changepoint_detected_GBM)))
graph_data$Iter   = orig_job_num
graph_data$Method = "Real-Time Monitoring GBM"
graph_data$ShiftType           = data_type
graph_data$has_missing         = should_include_missings
graph_data$distribution_change = distribution_change[parameters_number]
graph_data$any_discretes       = has_discretes
graph_data$modal_type          = modal_type[parameters_number]
#
data_name = paste("GBM", data_type,orig_job_num, modal_type[parameters_number],distribution_change[parameters_number],"includesMissing",should_include_missings, sep="_")
print(data_name)
write.csv(graph_data, paste("Simulation_Data/", data_name, ".csv", sep = ""))

print("Printing out all acceptances for mixture model:")
print(mixture_model_mergesplit_accepts)
print("Printing out all acceptances for MVN model:")
print(MVN_model_mergesplit_accepts)

print("The merge-split acceptance rate for the mixture model was:")
print(mean(mixture_model_mergesplit_accepts))
print("The merge-split acceptance rate for the MVN model was:")
print(mean(MVN_model_mergesplit_accepts))
