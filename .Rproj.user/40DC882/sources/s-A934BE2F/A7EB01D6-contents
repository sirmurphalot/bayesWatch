## Author: Alexander C. Murph
# Grab all the data and put into a single data frame.

file_names = list.files(path ="Fault_Detection_Results")
full_ranked_variables        = data.frame(matrix(0,ncol=7,nrow=1))
names(full_ranked_variables) = c(as.character(1:5),"Method","Iteration")
experiment_number = 3

for(name in file_names){
  split_file_name   = unlist(strsplit(name, split = "_"))
  simulation_num    = as.numeric(split_file_name[(length(split_file_name)-1)])
  parameters_number = simulation_num%%8+1
  if(parameters_number==experiment_number){
    temp_data_frame            = read.csv(paste("Fault_Detection_Results/",name,sep=""))[,-1]
    names(temp_data_frame)     = c(as.character(1:5),"Method","Iteration")
    full_ranked_variables      = rbind(full_ranked_variables,temp_data_frame)
  }
}
full_ranked_variables        = full_ranked_variables[-1,]
full_ranked_variables        = full_ranked_variables[complete.cases(full_ranked_variables),]
vars_of_interest             = c(1:max(full_ranked_variables[,1:5]))
ranked_variables_data        = data.frame(matrix(NA,nrow = 150, ncol = 4))
names(ranked_variables_data) = c("Proportion","Rank", "Method", "Variable")
count = 0
p = 5
for(var in vars_of_interest){
  ## Mixture Model
  temp_data = full_ranked_variables[which(full_ranked_variables$Method == "Mixture_Model"),]
  temp_data = colMeans(temp_data[,1:5]==var)
  ranked_variables_data[1:p + count, 1]     = temp_data
  ranked_variables_data[1:p + count, 2]     = 1:p
  ranked_variables_data[1:p + count, 3]     = rep("Mixture Model", times = p)
  ranked_variables_data[1:p + count, 4]     = rep(var, times = p)
  count = count + p
  
  ## Single Model
  temp_data = full_ranked_variables[which(full_ranked_variables$Method =="Single_Model"),]
  temp_data = colMeans(temp_data[,1:5]==var)
  ranked_variables_data[1:p + count, 1]     = temp_data
  ranked_variables_data[1:p + count, 2]     = 1:p
  ranked_variables_data[1:p + count, 3]     = rep("Single Model", times = p)
  ranked_variables_data[1:p + count, 4]     = rep(var, times = p)
  count = count + p
  
  ## GBM Model
  temp_data = full_ranked_variables[which(full_ranked_variables$Method =="GBM"),]
  temp_data = colMeans(temp_data[,1:5]==var)
  ranked_variables_data[1:p + count, 1]     = temp_data
  ranked_variables_data[1:p + count, 2]     = 1:p
  ranked_variables_data[1:p + count, 3]     = rep("GBM Sliding Window RTC", times = p)
  ranked_variables_data[1:p + count, 4]     = rep(var, times = p)
  count = count + p
}

for(rank_index in 1:5){
  temp_data             = ranked_variables_data[which( (ranked_variables_data$Rank==rank_index)&
                                                       (!(ranked_variables_data$Variable%in%c(3,4)))&
                                                       (ranked_variables_data$Method=="Mixture Model") ),]
  temp_row              = data.frame(Proportion = sum(temp_data$Proportion, na.rm = T), Rank = rank_index, 
                                     Method = "Mixture Model", Variable = "Other Vars")
  ranked_variables_data = rbind(ranked_variables_data, temp_row)
  
  temp_data             = ranked_variables_data[which( (ranked_variables_data$Rank==rank_index)&
                                                         (!(ranked_variables_data$Variable%in%c(3,4)))&
                                                         (ranked_variables_data$Method=="Single Model") ),]
  temp_row              = data.frame(Proportion = sum(temp_data$Proportion), Rank = rank_index, 
                                     Method = "Single Model", Variable = "Other Vars")
  ranked_variables_data = rbind(ranked_variables_data, temp_row)
  
  temp_data             = ranked_variables_data[which( (ranked_variables_data$Rank==rank_index)&
                                                         (!(ranked_variables_data$Variable%in%c(3,4)))&
                                                         (ranked_variables_data$Method=="GBM Sliding Window RTC") ),]
  temp_row              = data.frame(Proportion = sum(temp_data$Proportion), Rank = rank_index, 
                                     Method = "GBM Sliding Window RTC", Variable = "Other Vars")
  ranked_variables_data = rbind(ranked_variables_data, temp_row)
}
ranked_variables_data$Variable = (as.character(ranked_variables_data$Variable))
ranked_variables_data          = ranked_variables_data[which((ranked_variables_data$Variable%in%c("3","4","Other Vars"))),]
colorset                       = c('3'='deepskyblue','4'='deepskyblue4','Other Vars'='red2'  )

## Graph for GBM Method
temp_ranked_variables = ranked_variables_data[which(ranked_variables_data$Method == "GBM Sliding Window RTC"),]
p1 = ggplot(temp_ranked_variables, aes(x=Rank,y=Proportion,fill=Variable)) + geom_bar(stat="identity", position=position_dodge()) +
            xlab("Rank") + ylab("Proportion of Simulations") + ggtitle("GBM RTC, Simulation A") + ylim(0,1) +
            scale_fill_manual(values=colorset,breaks = c('3','4','Other Vars')) + theme(text = element_text(size = 13), legend.position = "none")

## Graph for Single Model Method
temp_ranked_variables = ranked_variables_data[which(ranked_variables_data$Method == "Single Model"),]
p2 = ggplot(temp_ranked_variables, aes(x=Rank,y=Proportion,fill=Variable)) + geom_bar(stat="identity", position=position_dodge()) +
  xlab("Rank") + ylab("Proportion of Simulations") + ggtitle("MVN Model, Simulation A") + ylim(0,1) +
  scale_fill_manual(values=colorset,breaks = c('3','4','Other Vars')) + theme(text = element_text(size = 13), legend.position = "none")

## Graph for Mixture Model Method
temp_ranked_variables = ranked_variables_data[which(ranked_variables_data$Method == "Mixture Model"),]
p3 = ggplot(temp_ranked_variables, aes(x=Rank,y=Proportion,fill=Variable)) + geom_bar(stat="identity", position=position_dodge()) +
  xlab("Rank") + ylab("Proportion of Simulations") + ggtitle("Mixture Model, Simulation A") + ylim(0,1) +
  scale_fill_manual(values=colorset,breaks = c('3','4','Other Vars')) + theme(text = element_text(size = 13), legend.position = "none")

prob_plot_list1 = list(p1,p2,p3)
#### Create a graphic with all of the plots 
# ggsave(
#   filename = "Posterior_Fault_Detection_Graphs.pdf",
#   plot = marrangeGrob(prob_plot_list, nrow=3, ncol=1),
#   width = 10, height = 9
# )

######## Now do the cov change experiment.
file_names = list.files(path ="Fault_Detection_Results")
full_ranked_variables        = data.frame(matrix(0,ncol=7,nrow=1))
names(full_ranked_variables) = c(as.character(1:5),"Method","Iteration")
experiment_number = 2

for(name in file_names){
  split_file_name   = unlist(strsplit(name, split = "_"))
  simulation_num    = as.numeric(split_file_name[(length(split_file_name)-1)])
  parameters_number = simulation_num%%8+1
  if(parameters_number==experiment_number){
    temp_data_frame            = read.csv(paste("Fault_Detection_Results/",name,sep=""))[,-1]
    names(temp_data_frame)     = c(as.character(1:5),"Method","Iteration")
    full_ranked_variables      = rbind(full_ranked_variables,temp_data_frame)
  }
}
full_ranked_variables        = full_ranked_variables[-1,]
full_ranked_variables        = full_ranked_variables[complete.cases(full_ranked_variables),]
vars_of_interest             = c(1:max(full_ranked_variables[,1:5]))
ranked_variables_data        = data.frame(matrix(NA,nrow = 150, ncol = 4))
names(ranked_variables_data) = c("Proportion","Rank", "Method", "Variable")
count = 0
p = 5
for(var in vars_of_interest){
  ## Mixture Model
  temp_data = full_ranked_variables[which(full_ranked_variables$Method == "Mixture_Model"),]
  temp_data = colMeans(temp_data[,1:5]==var)
  ranked_variables_data[1:p + count, 1]     = temp_data
  ranked_variables_data[1:p + count, 2]     = 1:p
  ranked_variables_data[1:p + count, 3]     = rep("Mixture Model", times = p)
  ranked_variables_data[1:p + count, 4]     = rep(var, times = p)
  count = count + p
  
  ## Single Model
  temp_data = full_ranked_variables[which(full_ranked_variables$Method =="Single_Model"),]
  temp_data = colMeans(temp_data[,1:5]==var)
  ranked_variables_data[1:p + count, 1]     = temp_data
  ranked_variables_data[1:p + count, 2]     = 1:p
  ranked_variables_data[1:p + count, 3]     = rep("Single Model", times = p)
  ranked_variables_data[1:p + count, 4]     = rep(var, times = p)
  count = count + p
  
  ## GBM Model
  temp_data = full_ranked_variables[which(full_ranked_variables$Method =="GBM"),]
  temp_data = colMeans(temp_data[,1:5]==var)
  ranked_variables_data[1:p + count, 1]     = temp_data
  ranked_variables_data[1:p + count, 2]     = 1:p
  ranked_variables_data[1:p + count, 3]     = rep("GBM Sliding Window RTC", times = p)
  ranked_variables_data[1:p + count, 4]     = rep(var, times = p)
  count = count + p
}

for(rank_index in 1:5){
  temp_data             = ranked_variables_data[which( (ranked_variables_data$Rank==rank_index)&
                                                         (!(ranked_variables_data$Variable%in%c(3,4)))&
                                                         (ranked_variables_data$Method=="Mixture Model") ),]
  temp_row              = data.frame(Proportion = sum(temp_data$Proportion, na.rm = T), Rank = rank_index, 
                                     Method = "Mixture Model", Variable = "Other Vars")
  ranked_variables_data = rbind(ranked_variables_data, temp_row)
  
  temp_data             = ranked_variables_data[which( (ranked_variables_data$Rank==rank_index)&
                                                         (!(ranked_variables_data$Variable%in%c(3,4)))&
                                                         (ranked_variables_data$Method=="Single Model") ),]
  temp_row              = data.frame(Proportion = sum(temp_data$Proportion), Rank = rank_index, 
                                     Method = "Single Model", Variable = "Other Vars")
  ranked_variables_data = rbind(ranked_variables_data, temp_row)
  
  temp_data             = ranked_variables_data[which( (ranked_variables_data$Rank==rank_index)&
                                                         (!(ranked_variables_data$Variable%in%c(3,4)))&
                                                         (ranked_variables_data$Method=="GBM Sliding Window RTC") ),]
  temp_row              = data.frame(Proportion = sum(temp_data$Proportion), Rank = rank_index, 
                                     Method = "GBM Sliding Window RTC", Variable = "Other Vars")
  ranked_variables_data = rbind(ranked_variables_data, temp_row)
}
ranked_variables_data$Variable = (as.character(ranked_variables_data$Variable))
ranked_variables_data          = ranked_variables_data[which((ranked_variables_data$Variable%in%c("3","4","Other Vars"))),]
colorset                       = c('3'='deepskyblue','4'='deepskyblue4','Other Vars'='red2'  )

## Graph for GBM Method
temp_ranked_variables = ranked_variables_data[which(ranked_variables_data$Method == "GBM Sliding Window RTC"),]
p1 = ggplot(temp_ranked_variables, aes(x=Rank,y=Proportion,fill=Variable)) + geom_bar(stat="identity", position=position_dodge()) +
  xlab("Rank") + ylab("Proportion of Simulations") + ggtitle("GBM RTC, Simulation B") + ylim(0,1) +
  scale_fill_manual(values=colorset,breaks = c('3','4','Other Vars')) + theme(text = element_text(size = 13))

## Graph for Single Model Method
temp_ranked_variables = ranked_variables_data[which(ranked_variables_data$Method == "Single Model"),]
p2 = ggplot(temp_ranked_variables, aes(x=Rank,y=Proportion,fill=Variable)) + geom_bar(stat="identity", position=position_dodge()) +
  xlab("Rank") + ylab("Proportion of Simulations") + ggtitle("MVN Model, Simulation B") + ylim(0,1) +
  scale_fill_manual(values=colorset,breaks = c('3','4','Other Vars')) + theme(text = element_text(size = 13))

## Graph for Mixture Model Method
temp_ranked_variables = ranked_variables_data[which(ranked_variables_data$Method == "Mixture Model"),]
p3 = ggplot(temp_ranked_variables, aes(x=Rank,y=Proportion,fill=Variable)) + geom_bar(stat="identity", position=position_dodge()) +
  xlab("Rank") + ylab("Proportion of Simulations") + ggtitle("Mixture Model, Simulation B") + ylim(0,1) +
  scale_fill_manual(values=colorset,breaks = c('3','4','Other Vars')) + theme(text = element_text(size = 13))

prob_plot_list2 = list(p1,p2,p3)
prob_plot_list  = c(prob_plot_list1, prob_plot_list2)
#### Create a graphic with all of the plots 
ggsave(
  filename = "Posterior_Fault_Detection_Graphs.pdf",
  plot = marrangeGrob(prob_plot_list, nrow=3, ncol=2, widths=c(1, 1.5)),
  width = 10, height = 9
)


