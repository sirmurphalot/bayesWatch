## T2 COMPARE SIMULATION GRAPHICS
library(ggplot2)
library(ggtext)
library(gridExtra)
library(stringr)
library(tidyverse)
library(knitr)
modal_type                = c('unimodal','unimodal','unimodal','unimodal','bimodal','bimodal','unimodal','unimodal')
# parameters_number       = job_num%%5+1
min_change_to_visualize   = 1
max_change_to_visualize   = 50
num_of_change_places      = 50
data_type                 ="BothAfirebrickifferent"
font_size                 = 8
base_data_size            = 400
grep_func                 = function(x){grep('X',x)}
location_of_changes       = c(14, 15, 16, 34, 35, 36)
true_changepoints         = c(14, 34)

x_axis_labels             = as.character(min_change_to_visualize:max_change_to_visualize)
x_axis_labels[1:9]        = ""
x_axis_labels[11:13]      = ""
x_axis_labels[15:19]      = ""
x_axis_labels[21:29]      = ""
x_axis_labels[31:33]      = ""
x_axis_labels[35:39]      = ""
x_axis_labels[41:49]      = ""

# x_axis_labels = paste("<span style = 'color: ", ifelse( (x_axis_labels=="14")|(x_axis_labels=="34"), 
#                                                         "firebrick", "black"), ";'>", x_axis_labels, "</span>", sep="")
# x_axis_labels = fct_reorder(x_axis_labels)

x_axis_colors             = rep("black", times = length(x_axis_labels))
x_axis_colors[14]         = "firebrick"
x_axis_colors[34]         = "firebrick"
 

if(length(list.files(path = "Simulation_Graphics/"))>0){
  stop("You've forgetted to delete the files in Simulation_Graphics")
}

job_num = 1
name                       = list.files(path = "Simulation_Data/")[1]
full_data = read.csv(paste("Simulation_Data/",name,sep=""))[,-1]
full_data$modal_type = modal_type[full_data$Iter%%8+1]

for(name in list.files(path = "Simulation_Data/")[-1]){
  temp_data            = read.csv(paste("Simulation_Data/", name, sep=""))[,-1]
  if(anyNA(temp_data) | ncol(temp_data)!=57){
    print(name)
    next
  }
  temp_data$modal_type = modal_type[temp_data$Iter%%8+1]
  
  full_data            = rbind(full_data, temp_data)
}

# full_data[which(full_data$Method == "Real")]


# #### Continuous, Non-missing, Unimodal, mean change
data_for_group = full_data[which( (full_data$any_discretes==0)&
                                  (full_data$has_missing==0)&
                                  (full_data$modal_type=='unimodal')&
                                  (full_data$distribution_change=='mean')&
                                  (full_data$ShiftType=='MeanIsDifferent') ),]
plot_name      = c("B: Continuous", "NoMissing", "Unimodal", "MeanChange")

data_for_mix   = data_for_group[which(data_for_group$Method == 'MIXTUREMODEL'),]
values         = colMeans(data_for_mix[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('Mixture Model', times = num_of_change_places)
graph_data     = data.frame(Values = values, Location = location, Method = method)

num_of_regimes                 = length(unlist(sapply(names(data_for_mix), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mixtures = mean(rowMeans(data_for_mix[,indicies_of_non_changes]))

data_for_mvn   = data_for_group[which(data_for_group$Method == 'SINGLEMVNMODEL'),]
values         = colMeans(data_for_mvn[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('MVN Model', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_mvn), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mvn      = mean(rowMeans(data_for_mvn[,indicies_of_non_changes]))

data_for_t2    = data_for_group[which(data_for_group$Method == 'T2'),]
values         = colMeans(data_for_t2[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('T2 Scan Stat', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_t2), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_t2       = mean(rowMeans(data_for_t2[,indicies_of_non_changes]))

data_for_gbm   = data_for_group[which(data_for_group$Method == 'Real-Time Monitoring GBM'),]
values         = colMeans(data_for_gbm[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('RTM GBM', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_gbm), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_gbm      = mean(rowMeans(data_for_gbm[,indicies_of_non_changes]))

# g1 = ggplot(graph_data, aes(x=Location, y=Values, fill = Method)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   ggtitle(paste(plot_name,collapse=", ")) + ylim(0,1) + xlim(min_change_to_visualize,max_change_to_visualize)+ 
#   xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") + theme(text = element_text(size = font_size))
# ggsave(
#   filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
#   plot = g1
# )

summary.table <- structure(list(Method = c("MVN Model", "Mixture Model", "RTM GBM", "T2 Scan Stat"),
                                prob.fp = c(average_false_pos_for_mvn, average_false_pos_for_mixtures, 
                                            average_false_pos_for_gbm, average_false_pos_for_t2)), 
                           row.names = c(NA, -4L), class = "data.frame")
table = summary.table %>%
  rename(`Average FP %` = prob.fp) %>%
  kable %>%
  gsub('|', ' ', ., fixed = T) %>%
  strsplit('\n') %>%
  trimws
header = table[[1]]
header = paste0(header, '\n', paste0(rep('─', nchar(header)), collapse =''))
table = table[-(1:2)]
# table = do.call(rbind, list(table))[,1]
table = data.frame(Method=summary.table$Method, lab = table)

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

for(location_index in 1:nrow(plot_data)){
  plot_data$Location[location_index] = paste("<span style = 'color: ", 
                             ifelse((plot_data$Location[location_index]==14)|(plot_data$Location[location_index]==34), 
                                     "firebrick", "black"),";'>", plot_data$Location[location_index], "</span>", sep="")
}


x_axis_labels             = unique(plot_data$Location)
x_axis_labels[1:9]        = ""
x_axis_labels[11:13]      = ""
x_axis_labels[15:19]      = ""
x_axis_labels[21:29]      = ""
x_axis_labels[31:33]      = ""
x_axis_labels[35:39]      = ""
x_axis_labels[41:49]      = ""

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

g1 = ggplot(plot_data, aes(x = Location, y = Values, group = Method, fill = lab)) +
  geom_bar(stat="identity", position = position_dodge()) +
  ggtitle(paste(plot_name,collapse=", ")) + geom_vline(xintercept = c(true_changepoints), color="firebrick",linetype = "dashed") + 
  scale_x_discrete(limits=factor(x_axis_labels), expand=c(0,0)) +  scale_y_continuous(limits=c(0,1)) +
  xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") +
  guides(fill = guide_legend(header, reverse=TRUE, 
                              label.position = "right", 
                              title.theme = element_text(size=13, family='mono'),
                              label.theme = element_text(size=11, family='mono'))) +
  theme( axis.text.x = element_markdown(),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.background = element_blank())
ggsave(
  filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
  plot = g1,
  width = 11, 
  height = 2
)
  
  

# #### Continuous, Missing, bimodal, mean change
data_for_group = full_data[which( (full_data$any_discretes==0)&
                                    (full_data$has_missing==0)&
                                    (full_data$modal_type=='bimodal')&
                                    (full_data$distribution_change=='mean')&
                                    (full_data$ShiftType=='MeanIsDifferent') ),]
plot_name      = c("C: Continuous", "No Missing Values", "Bimodal", "Mean Change")

data_for_mix   = data_for_group[which(data_for_group$Method == 'MIXTUREMODEL'),]
values         = colMeans(data_for_mix[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('Mixture Model', times = num_of_change_places)
graph_data     = data.frame(Values = values, Location = location, Method = method)

num_of_regimes                 = length(unlist(sapply(names(data_for_mix), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mixtures = mean(rowMeans(data_for_mix[,indicies_of_non_changes]))

data_for_mvn   = data_for_group[which(data_for_group$Method == 'SINGLEMVNMODEL'),]
values         = colMeans(data_for_mvn[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('MVN Model', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_mvn), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mvn      = mean(rowMeans(data_for_mvn[,indicies_of_non_changes]))

data_for_t2    = data_for_group[which(data_for_group$Method == 'T2'),]
values         = colMeans(data_for_t2[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('T2 Scan Stat', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_t2), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_t2       = mean(rowMeans(data_for_t2[,indicies_of_non_changes]))

data_for_gbm   = data_for_group[which(data_for_group$Method == 'Real-Time Monitoring GBM'),]
values         = colMeans(data_for_gbm[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('RTM GBM', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_gbm), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_gbm      = mean(rowMeans(data_for_gbm[,indicies_of_non_changes]))

# g1 = ggplot(graph_data, aes(x=Location, y=Values, fill = Method)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   ggtitle(paste(plot_name,collapse=", ")) + ylim(0,1) + xlim(min_change_to_visualize,max_change_to_visualize)+ 
#   xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") + theme(text = element_text(size = font_size))
# ggsave(
#   filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
#   plot = g1
# )

summary.table <- structure(list(Method = c("MVN Model", "Mixture Model", "RTM GBM", "T2 Scan Stat"),
                                prob.fp = c(average_false_pos_for_mvn, average_false_pos_for_mixtures, 
                                            average_false_pos_for_gbm, average_false_pos_for_t2)), 
                           row.names = c(NA, -4L), class = "data.frame")
table = summary.table %>%
  rename(`Average FP %` = prob.fp) %>%
  kable %>%
  gsub('|', ' ', ., fixed = T) %>%
  strsplit('\n') %>%
  trimws
header = table[[1]]
header = paste0(header, '\n', paste0(rep('─', nchar(header)), collapse =''))
table = table[-(1:2)]
# table = do.call(rbind, list(table))[,1]
table = data.frame(Method=summary.table$Method, lab = table)

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

for(location_index in 1:nrow(plot_data)){
  plot_data$Location[location_index] = paste("<span style = 'color: ", 
                                             ifelse((plot_data$Location[location_index]==14)|(plot_data$Location[location_index]==34), 
                                                    "firebrick", "black"),";'>", plot_data$Location[location_index], "</span>", sep="")
}


x_axis_labels             = unique(plot_data$Location)
x_axis_labels[1:9]        = ""
x_axis_labels[11:13]      = ""
x_axis_labels[15:19]      = ""
x_axis_labels[21:29]      = ""
x_axis_labels[31:33]      = ""
x_axis_labels[35:39]      = ""
x_axis_labels[41:49]      = ""

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

g4 = ggplot(plot_data, aes(x = Location, y = Values, group = Method, fill = lab)) +
  geom_bar(stat="identity", position = position_dodge()) +
  ggtitle(paste(plot_name,collapse=", ")) + geom_vline(xintercept = c(true_changepoints), color="firebrick", linetype="dashed") + 
  scale_x_discrete(limits=factor(x_axis_labels), expand=c(0,0)) +  scale_y_continuous(limits=c(0,1)) +
  xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") +
  guides(fill = guide_legend(header, reverse=TRUE, 
                             label.position = "right", 
                             title.theme = element_text(size=13, family='mono'),
                             label.theme = element_text(size=11, family='mono'))) +
  theme(axis.text.x = element_markdown(),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.background = element_blank())
ggsave(
  filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
  plot = g4,
  width = 11, 
  height = 2
)

# 
# #### Continuous, Missing, Unimodal, cov change
data_for_group = full_data[which( (full_data$any_discretes==0)&
                                    (full_data$has_missing==1)&
                                    (full_data$modal_type=='unimodal')&
                                    (full_data$distribution_change=='cov')&
                                    (full_data$ShiftType=='CovIsDifferent') ),]
plot_name      = c("A: Continuous", "Has Missing Values", "Unimodal", "Covariance Change")

data_for_mix   = data_for_group[which(data_for_group$Method == 'MIXTUREMODEL'),]
values         = colMeans(data_for_mix[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('Mixture Model', times = num_of_change_places)
graph_data     = data.frame(Values = values, Location = location, Method = method)

num_of_regimes                 = length(unlist(sapply(names(data_for_mix), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mixtures = mean(rowMeans(data_for_mix[,indicies_of_non_changes]))

data_for_mvn   = data_for_group[which(data_for_group$Method == 'SINGLEMVNMODEL'),]
values         = colMeans(data_for_mvn[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('MVN Model', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_mvn), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mvn      = mean(rowMeans(data_for_mvn[,indicies_of_non_changes]))

data_for_t2    = data_for_group[which(data_for_group$Method == 'T2'),]
values         = colMeans(data_for_t2[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('T2 Scan Stat', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_t2), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_t2       = mean(rowMeans(data_for_t2[,indicies_of_non_changes]))

data_for_gbm   = data_for_group[which(data_for_group$Method == 'Real-Time Monitoring GBM'),]
values         = colMeans(data_for_gbm[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('RTM GBM', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_gbm), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_gbm      = mean(rowMeans(data_for_gbm[,indicies_of_non_changes]))

# g1 = ggplot(graph_data, aes(x=Location, y=Values, fill = Method)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   ggtitle(paste(plot_name,collapse=", ")) + ylim(0,1) + xlim(min_change_to_visualize,max_change_to_visualize)+ 
#   xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") + theme(text = element_text(size = font_size))
# ggsave(
#   filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
#   plot = g1
# )

summary.table <- structure(list(Method = c("MVN Model", "Mixture Model", "RTM GBM", "T2 Scan Stat"),
                                prob.fp = c(average_false_pos_for_mvn, average_false_pos_for_mixtures, 
                                            average_false_pos_for_gbm, average_false_pos_for_t2)), 
                           row.names = c(NA, -4L), class = "data.frame")
table = summary.table %>%
  rename(`Average FP %` = prob.fp) %>%
  kable %>%
  gsub('|', ' ', ., fixed = T) %>%
  strsplit('\n') %>%
  trimws
header = table[[1]]
header = paste0(header, '\n', paste0(rep('─', nchar(header)), collapse =''))
table = table[-(1:2)]
# table = do.call(rbind, list(table))[,1]
table = data.frame(Method=summary.table$Method, lab = table)

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

for(location_index in 1:nrow(plot_data)){
  plot_data$Location[location_index] = paste("<span style = 'color: ", 
                                             ifelse((plot_data$Location[location_index]==14)|(plot_data$Location[location_index]==34), 
                                                    "firebrick", "black"),";'>", plot_data$Location[location_index], "</span>", sep="")
}


x_axis_labels             = unique(plot_data$Location)
x_axis_labels[1:9]        = ""
x_axis_labels[11:13]      = ""
x_axis_labels[15:19]      = ""
x_axis_labels[21:29]      = ""
x_axis_labels[31:33]      = ""
x_axis_labels[35:39]      = ""
x_axis_labels[41:49]      = ""

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

g5 = ggplot(plot_data, aes(x = Location, y = Values, group = Method, fill = lab)) +
  geom_bar(stat="identity", position = position_dodge()) +
  ggtitle(paste(plot_name,collapse=", ")) + geom_vline(xintercept = c(true_changepoints), color="firebrick", linetype="dashed") + 
  scale_x_discrete(limits=factor(x_axis_labels), expand=c(0,0)) +  scale_y_continuous(limits=c(0,1)) +
  xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") +
  guides(fill = guide_legend(header, reverse=TRUE, 
                             label.position = "right", 
                             title.theme = element_text(size=13, family='mono'),
                             label.theme = element_text(size=11, family='mono'))) +
  theme(axis.text.x = element_markdown(),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.background = element_blank())
ggsave(
  filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
  plot = g5,
  width = 11, 
  height = 2
)

# #### Mixed, Missing, Unimodal, cov change
data_for_group = full_data[which( (full_data$any_discretes==1)&
                                    (full_data$has_missing==1)&
                                    (full_data$modal_type=='unimodal')&
                                    (full_data$distribution_change=='cov')&
                                    (full_data$ShiftType=='CovIsDifferent') ),]
plot_name      = c("F: Mixed", "Has Missing Values", "Unimodal", "Covariance Change")

data_for_mix   = data_for_group[which(data_for_group$Method == 'MIXTUREMODEL'),]
values         = colMeans(data_for_mix[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('Mixture Model', times = num_of_change_places)
graph_data     = data.frame(Values = values, Location = location, Method = method)

num_of_regimes                 = length(unlist(sapply(names(data_for_mix), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mixtures = mean(rowMeans(data_for_mix[,indicies_of_non_changes]))

data_for_mvn   = data_for_group[which(data_for_group$Method == 'SINGLEMVNMODEL'),]
values         = colMeans(data_for_mvn[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('MVN Model', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_mvn), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mvn      = mean(rowMeans(data_for_mvn[,indicies_of_non_changes]))

data_for_t2    = data_for_group[which(data_for_group$Method == 'T2'),]
values         = colMeans(data_for_t2[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('T2 Scan Stat', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_t2), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_t2       = mean(rowMeans(data_for_t2[,indicies_of_non_changes]))

data_for_gbm   = data_for_group[which(data_for_group$Method == 'Real-Time Monitoring GBM'),]
values         = colMeans(data_for_gbm[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('RTM GBM', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_gbm), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_gbm      = mean(rowMeans(data_for_gbm[,indicies_of_non_changes]))

# g1 = ggplot(graph_data, aes(x=Location, y=Values, fill = Method)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   ggtitle(paste(plot_name,collapse=", ")) + ylim(0,1) + xlim(min_change_to_visualize,max_change_to_visualize)+ 
#   xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") + theme(text = element_text(size = font_size))
# ggsave(
#   filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
#   plot = g1
# )

summary.table <- structure(list(Method = c("MVN Model", "Mixture Model", "RTM GBM", "T2 Scan Stat"),
                                prob.fp = c(average_false_pos_for_mvn, average_false_pos_for_mixtures, 
                                            average_false_pos_for_gbm, average_false_pos_for_t2)), 
                           row.names = c(NA, -4L), class = "data.frame")
table = summary.table %>%
  rename(`Average FP %` = prob.fp) %>%
  kable %>%
  gsub('|', ' ', ., fixed = T) %>%
  strsplit('\n') %>%
  trimws
header = table[[1]]
header = paste0(header, '\n', paste0(rep('─', nchar(header)), collapse =''))
table = table[-(1:2)]
# table = do.call(rbind, list(table))[,1]
table = data.frame(Method=summary.table$Method, lab = table)

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

for(location_index in 1:nrow(plot_data)){
  plot_data$Location[location_index] = paste("<span style = 'color: ", 
                                             ifelse((plot_data$Location[location_index]==14)|(plot_data$Location[location_index]==34), 
                                                    "firebrick", "black"),";'>", plot_data$Location[location_index], "</span>", sep="")
}


x_axis_labels             = unique(plot_data$Location)
x_axis_labels[1:9]        = ""
x_axis_labels[11:13]      = ""
x_axis_labels[15:19]      = ""
x_axis_labels[21:29]      = ""
x_axis_labels[31:33]      = ""
x_axis_labels[35:39]      = ""
x_axis_labels[41:49]      = ""

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

g8 = ggplot(plot_data, aes(x = Location, y = Values, group = Method, fill = lab)) +
  geom_bar(stat="identity", position = position_dodge()) +
  ggtitle(paste(plot_name,collapse=", ")) + geom_vline(xintercept = c(true_changepoints), color="firebrick", linetype="dashed") + 
  scale_x_discrete(limits=factor(x_axis_labels), expand=c(0,0)) +  scale_y_continuous(limits=c(0,1)) +
  xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") +
  guides(fill = guide_legend(header, reverse=TRUE, 
                             label.position = "right", 
                             title.theme = element_text(size=13, family='mono'),
                             label.theme = element_text(size=11, family='mono'))) +
  theme(axis.text.x = element_markdown(),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.background = element_blank())
ggsave(
  filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
  plot = g8,
  width = 11, 
  height = 2
)

# #### Mixed, Missing, Bimodal, cov change
data_for_group = full_data[which( (full_data$any_discretes==1)&
                                    (full_data$has_missing==1)&
                                    (full_data$modal_type=='bimodal')&
                                    (full_data$distribution_change=='cov')&
                                    (full_data$ShiftType=='CovIsDifferent') ),]
plot_name      = c("E: Mixed", "Has Missing Values", "Bimodal", "Covariance Change")

data_for_mix   = data_for_group[which(data_for_group$Method == 'MIXTUREMODEL'),]
values         = colMeans(data_for_mix[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('Mixture Model', times = num_of_change_places)
graph_data     = data.frame(Values = values, Location = location, Method = method)

num_of_regimes                 = length(unlist(sapply(names(data_for_mix), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mixtures = mean(rowMeans(data_for_mix[,indicies_of_non_changes]))

data_for_mvn   = data_for_group[which(data_for_group$Method == 'SINGLEMVNMODEL'),]
values         = colMeans(data_for_mvn[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('MVN Model', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_mvn), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mvn      = mean(rowMeans(data_for_mvn[,indicies_of_non_changes]))

data_for_t2    = data_for_group[which(data_for_group$Method == 'T2'),]
values         = colMeans(data_for_t2[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('T2 Scan Stat', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_t2), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_t2       = mean(rowMeans(data_for_t2[,indicies_of_non_changes]))

data_for_gbm   = data_for_group[which(data_for_group$Method == 'Real-Time Monitoring GBM'),]
values         = colMeans(data_for_gbm[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('RTM GBM', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_gbm), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_gbm      = mean(rowMeans(data_for_gbm[,indicies_of_non_changes]))

# g1 = ggplot(graph_data, aes(x=Location, y=Values, fill = Method)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   ggtitle(paste(plot_name,collapse=", ")) + ylim(0,1) + xlim(min_change_to_visualize,max_change_to_visualize)+ 
#   xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") + theme(text = element_text(size = font_size))
# ggsave(
#   filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
#   plot = g1
# )

summary.table <- structure(list(Method = c("MVN Model", "Mixture Model", "RTM GBM", "T2 Scan Stat"),
                                prob.fp = c(average_false_pos_for_mvn, average_false_pos_for_mixtures, 
                                            average_false_pos_for_gbm, average_false_pos_for_t2)), 
                           row.names = c(NA, -4L), class = "data.frame")
table = summary.table %>%
  rename(`Average FP %` = prob.fp) %>%
  kable %>%
  gsub('|', ' ', ., fixed = T) %>%
  strsplit('\n') %>%
  trimws
header = table[[1]]
header = paste0(header, '\n', paste0(rep('─', nchar(header)), collapse =''))
table = table[-(1:2)]
# table = do.call(rbind, list(table))[,1]
table = data.frame(Method=summary.table$Method, lab = table)

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

g11 = ggplot(plot_data, aes(x = Location, y = Values, group = Method, fill = lab)) +
  geom_bar(stat="identity", position = position_dodge()) +
  ggtitle(paste(plot_name,collapse=", ")) + geom_vline(xintercept = c(true_changepoints), color="firebrick", linetype="dashed") + 
  scale_x_discrete(limits=factor(x_axis_labels), expand=c(0,0)) +  scale_y_continuous(limits=c(0,1)) +
  xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") +
  guides(fill = guide_legend(header, reverse=TRUE, 
                             label.position = "right", 
                             title.theme = element_text(size=13, family='mono'),
                             label.theme = element_text(size=11, family='mono'))) +
  theme(axis.text.x = element_markdown(),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.background = element_blank())
ggsave(
  filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
  plot = g11,
  width = 11, 
  height = 2
)


#### Missing Structure Change
data_for_group = full_data[which( (full_data$distribution_change=='miss')&
                                    (full_data$ShiftType=='MissStructDiff') ),]
plot_name      = c("D: Missingness Structure Change")

data_for_mix   = data_for_group[which(data_for_group$Method == 'MIXTUREMODEL'),]
values         = colMeans(data_for_mix[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('Mixture Model', times = num_of_change_places)
graph_data     = data.frame(Values = values, Location = location, Method = method)

num_of_regimes                 = length(unlist(sapply(names(data_for_mix), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mixtures = mean(rowMeans(data_for_mix[,indicies_of_non_changes]))

data_for_mvn   = data_for_group[which(data_for_group$Method == 'SINGLEMVNMODEL'),]
values         = colMeans(data_for_mvn[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('MVN Model', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_mvn), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mvn      = mean(rowMeans(data_for_mvn[,indicies_of_non_changes]))

data_for_t2    = data_for_group[which(data_for_group$Method == 'T2'),]
values         = colMeans(data_for_t2[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('T2 Scan Stat', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_t2), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_t2       = mean(rowMeans(data_for_t2[,indicies_of_non_changes]))

data_for_gbm   = data_for_group[which(data_for_group$Method == 'Real-Time Monitoring GBM'),]
values         = colMeans(data_for_gbm[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('RTM GBM', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_gbm), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_gbm      = mean(rowMeans(data_for_gbm[,indicies_of_non_changes]))

# g1 = ggplot(graph_data, aes(x=Location, y=Values, fill = Method)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   ggtitle(paste(plot_name,collapse=", ")) + ylim(0,1) + xlim(min_change_to_visualize,max_change_to_visualize)+ 
#   xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") + theme(text = element_text(size = font_size))
# ggsave(
#   filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
#   plot = g1
# )

summary.table <- structure(list(Method = c("MVN Model", "Mixture Model", "RTM GBM", "T2 Scan Stat"),
                                prob.fp = c(average_false_pos_for_mvn, average_false_pos_for_mixtures, 
                                            average_false_pos_for_gbm, average_false_pos_for_t2)), 
                           row.names = c(NA, -4L), class = "data.frame")
table = summary.table %>%
  rename(`Average FP %` = prob.fp) %>%
  kable %>%
  gsub('|', ' ', ., fixed = T) %>%
  strsplit('\n') %>%
  trimws
header = table[[1]]
header = paste0(header, '\n', paste0(rep('─', nchar(header)), collapse =''))
table = table[-(1:2)]
# table = do.call(rbind, list(table))[,1]
table = data.frame(Method=summary.table$Method, lab = table)

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

for(location_index in 1:nrow(plot_data)){
  plot_data$Location[location_index] = paste("<span style = 'color: ", 
                                             ifelse((plot_data$Location[location_index]==14)|(plot_data$Location[location_index]==34), 
                                                    "firebrick", "black"),";'>", plot_data$Location[location_index], "</span>", sep="")
}


x_axis_labels             = unique(plot_data$Location)
x_axis_labels[1:9]        = ""
x_axis_labels[11:13]      = ""
x_axis_labels[15:19]      = ""
x_axis_labels[21:29]      = ""
x_axis_labels[31:33]      = ""
x_axis_labels[35:39]      = ""
x_axis_labels[41:49]      = ""

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

g13 = ggplot(plot_data, aes(x = Location, y = Values, group = Method, fill = lab)) +
  geom_bar(stat="identity", position = position_dodge()) +
  ggtitle(paste(plot_name,collapse=", ")) + geom_vline(xintercept = c(true_changepoints), color="firebrick",linetype="dashed") + 
  scale_x_discrete(limits=factor(x_axis_labels), expand=c(0,0)) +  scale_y_continuous(limits=c(0,1)) +
  xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") +
  guides(fill = guide_legend(header, reverse=TRUE, 
                             label.position = "right", 
                             title.theme = element_text(size=13, family='mono'),
                             label.theme = element_text(size=11, family='mono'))) +
  theme(axis.text.x = element_markdown(),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.background = element_blank())
ggsave(
  filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
  plot = g13,
  width = 11, 
  height = 2
)


# #### Linear Drift Change, mean
# location_of_changes_drift = c(14:34)
# data_for_group = full_data[which( (full_data$distribution_change=='mean')&(full_data$ShiftType=='LinearDrift') ),]
# plot_name      = c("Linear Drift Change of Mean")
# 
# data_for_mix   = data_for_group[which(data_for_group$Method == 'MIXTUREMODEL'),]
# values         = colMeans(data_for_mix[,1:num_of_change_places])
# location       = 1:num_of_change_places
# method         = rep('Mixture Model', times = num_of_change_places)
# graph_data     = data.frame(Values = values, Location = location, Method = method)
# 
# num_of_regimes                 = length(unlist(sapply(names(data_for_mix), grep_func)))
# indicies_of_non_changes        = 1:num_of_regimes
# indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes_drift]
# average_false_pos_for_mixtures = mean(rowMeans(data_for_mix[,indicies_of_non_changes]))
# 
# data_for_mvn   = data_for_group[which(data_for_group$Method == 'SINGLEMVNMODEL'),]
# values         = colMeans(data_for_mvn[,1:num_of_change_places])
# location       = 1:num_of_change_places
# method         = rep('MVN Model', times = num_of_change_places)
# graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))
# 
# num_of_regimes                 = length(unlist(sapply(names(data_for_mvn), grep_func)))
# indicies_of_non_changes        = 1:num_of_regimes
# indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes_drift]
# average_false_pos_for_mvn      = mean(rowMeans(data_for_mvn[,indicies_of_non_changes]))
# 
# data_for_t2    = data_for_group[which(data_for_group$Method == 'T2'),]
# values         = colMeans(data_for_t2[,1:num_of_change_places])
# location       = 1:num_of_change_places
# method         = rep('T2 Scan Stat', times = num_of_change_places)
# graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))
# 
# num_of_regimes                 = length(unlist(sapply(names(data_for_t2), grep_func)))
# indicies_of_non_changes        = 1:num_of_regimes
# indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes_drift]
# average_false_pos_for_t2       = mean(rowMeans(data_for_t2[,indicies_of_non_changes]))
# 
# data_for_gbm   = data_for_group[which(data_for_group$Method == 'Real-Time Monitoring GBM'),]
# values         = colMeans(data_for_gbm[,1:num_of_change_places])
# location       = 1:num_of_change_places
# method         = rep('RTM GBM', times = num_of_change_places)
# graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))
# 
# num_of_regimes                 = length(unlist(sapply(names(data_for_gbm), grep_func)))
# indicies_of_non_changes        = 1:num_of_regimes
# indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes_drift]
# average_false_pos_for_gbm      = mean(rowMeans(data_for_gbm[,indicies_of_non_changes]))
# 
# # g1 = ggplot(graph_data, aes(x=Location, y=Values, fill = Method)) +
# #   geom_bar(stat="identity", position=position_dodge()) +
# #   ggtitle(paste(plot_name,collapse=", ")) + ylim(0,1) + xlim(min_change_to_visualize,max_change_to_visualize)+ 
# #   xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") + theme(text = element_text(size = font_size))
# # ggsave(
# #   filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
# #   plot = g1
# # )
# 
# summary.table <- structure(list(Method = c("MVN Model", "Mixture Model", "RTM GBM", "T2 Scan Stat"),
#                                 prob.fp = c(average_false_pos_for_mvn, average_false_pos_for_mixtures, 
#                                             average_false_pos_for_gbm, average_false_pos_for_t2)), 
#                            row.names = c(NA, -4L), class = "data.frame")
# table = summary.table %>%
#   rename(`Average FP %` = prob.fp) %>%
#   kable %>%
#   gsub('|', ' ', ., fixed = T) %>%
#   strsplit('\n') %>%
#   trimws
# header = table[[1]]
# header = paste0(header, '\n', paste0(rep('─', nchar(header)), collapse =''))
# table = table[-(1:2)]
# # table = do.call(rbind, list(table))[,1]
# table = data.frame(Method=summary.table$Method, lab = table)
# 
# plot_data = graph_data %>%
#   group_by(Method) %>%
#   left_join(table)
# 
# g14 = ggplot(plot_data, aes(x = Location, y = Values, group = Method, fill = lab)) +
#   geom_bar(stat="identity", position = position_dodge()) +
#   ggtitle(paste(plot_name,collapse=", ")) + scale_x_discrete(limits=factor(min_change_to_visualize:max_change_to_visualize)) + scale_y_continuous(limits=c(0,1)) +
#   xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") +
#   guides(fill = guide_legend(header, reverse=TRUE, 
#                              label.position = "right", 
#                              title.theme = element_text(size=13, family='mono'),
#                              label.theme = element_text(size=11, family='mono'))) +
#   theme(
#     legend.key = element_rect(fill = NA, colour = NA),
#     legend.spacing.y = unit(0, "pt"),
#     legend.key.height = unit(10, "pt"),
#     legend.background = element_blank())
# ggsave(
#   filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
#   plot = g14,
#   width = 11, 
#   height = 2
# )


#### Mean example with missing 
data_for_group = full_data[which( (full_data$any_discretes==1)&
                                    (full_data$has_missing==1)&
                                    (full_data$modal_type=='unimodal')&
                                    (full_data$distribution_change=='mean')&
                                    (full_data$ShiftType=='MeanIsDifferent') ),]
plot_name      = c("G: Mixed", "Has Missing Values", "Unimodal", "Mean Change")

data_for_mix   = data_for_group[which(data_for_group$Method == 'MIXTUREMODEL'),]
values         = colMeans(data_for_mix[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('Mixture Model', times = num_of_change_places)
graph_data     = data.frame(Values = values, Location = location, Method = method)

num_of_regimes                 = length(unlist(sapply(names(data_for_mix), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mixtures = mean(rowMeans(data_for_mix[,indicies_of_non_changes]))

data_for_mvn   = data_for_group[which(data_for_group$Method == 'SINGLEMVNMODEL'),]
values         = colMeans(data_for_mvn[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('MVN Model', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_mvn), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_mvn      = mean(rowMeans(data_for_mvn[,indicies_of_non_changes]))

data_for_t2    = data_for_group[which(data_for_group$Method == 'T2'),]
values         = colMeans(data_for_t2[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('T2 Scan Stat', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_t2), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_t2       = mean(rowMeans(data_for_t2[,indicies_of_non_changes]))

data_for_gbm   = data_for_group[which(data_for_group$Method == 'Real-Time Monitoring GBM'),]
values         = colMeans(data_for_gbm[,1:num_of_change_places])
location       = 1:num_of_change_places
method         = rep('RTM GBM', times = num_of_change_places)
graph_data     = rbind(graph_data, data.frame(Values = values, Location = location, Method = method))

num_of_regimes                 = length(unlist(sapply(names(data_for_gbm), grep_func)))
indicies_of_non_changes        = 1:num_of_regimes
indicies_of_non_changes        = indicies_of_non_changes[-location_of_changes]
average_false_pos_for_gbm      = mean(rowMeans(data_for_gbm[,indicies_of_non_changes]))

# g1 = ggplot(graph_data, aes(x=Location, y=Values, fill = Method)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   ggtitle(paste(plot_name,collapse=", ")) + ylim(0,1) + xlim(min_change_to_visualize,max_change_to_visualize)+ 
#   xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") + theme(text = element_text(size = font_size))
# ggsave(
#   filename = paste("Simulation_Graphics/", gsub(" ","_", paste(plot_name,collapse="_") ), ".png", sep = ""),
#   plot = g1
# )

summary.table <- structure(list(Method = c("MVN Model", "Mixture Model", "RTM GBM", "T2 Scan Stat"),
                                prob.fp = c(average_false_pos_for_mvn, average_false_pos_for_mixtures, 
                                            average_false_pos_for_gbm, average_false_pos_for_t2)), 
                           row.names = c(NA, -4L), class = "data.frame")
table = summary.table %>%
  rename(`Average FP %` = prob.fp) %>%
  kable %>%
  gsub('|', ' ', ., fixed = T) %>%
  strsplit('\n') %>%
  trimws
header = table[[1]]
header = paste0(header, '\n', paste0(rep('─', nchar(header)), collapse =''))
table = table[-(1:2)]
# table = do.call(rbind, list(table))[,1]
table = data.frame(Method=summary.table$Method, lab = table)

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

for(location_index in 1:nrow(plot_data)){
  plot_data$Location[location_index] = paste("<span style = 'color: ", 
                                             ifelse((plot_data$Location[location_index]==14)|(plot_data$Location[location_index]==34), 
                                                    "firebrick", "black"),";'>", plot_data$Location[location_index], "</span>", sep="")
}


x_axis_labels             = unique(plot_data$Location)
x_axis_labels[1:9]        = ""
x_axis_labels[11:13]      = ""
x_axis_labels[15:19]      = ""
x_axis_labels[21:29]      = ""
x_axis_labels[31:33]      = ""
x_axis_labels[35:39]      = ""
x_axis_labels[41:49]      = ""

plot_data = graph_data %>%
  group_by(Method) %>%
  left_join(table)

g15 = ggplot(plot_data, aes(x = Location, y = Values, group = Method, fill = lab)) +
  geom_bar(stat="identity", position = position_dodge()) +
    ggtitle(paste(plot_name,collapse=", ")) + geom_vline(xintercept = c(true_changepoints), color="firebrick", linetype="dashed") + 
  scale_x_discrete(limits=factor(x_axis_labels), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0,1)) +
    xlab("Day Before Possible Regime Change") + ylab("Proportion of \n Simulations Flagged") +
  guides(fill = guide_legend(header, reverse=TRUE,  
                             label.position = "right", 
                             title.theme = element_text(size=13, family='mono'),
                             label.theme = element_text(size=11, family='mono'))) +
  theme(axis.text.x = element_markdown(),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(10, "pt"),
    legend.background = element_blank())
ggsave(
  filename = paste("Simulation_Graphics/", gsub(" ","_", gsub(" ","_", paste(plot_name,collapse="_") ) ), ".png", sep = ""),
  plot = g15,
  width = 11, 
  height = 2
)

# # prob_plot_list = list(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13)
# prob_plot_list = list(g1, g4,g5,g8,g11,g13,g14,g15)
# #### Create a graphic with all of the plots
# ggsave(
#   filename = paste("Simulation_All_Graphics_Sep26_presentationVersion", ".pdf", sep = ''),
#   plot = marrangeGrob(prob_plot_list, nrow=4, ncol=1),
#   width = 11, height = 8
# )

# example image

library("png")
# setup plot
pdf(paste("Jan8PRELIM_lambda20_df20_hyperb20_dirichlet0_004_basedata1000", ".pdf", sep = ''), width=8, height=11)
par(mar=rep(0,4)) # no margins
layout(matrix(1:7, ncol=1,nrow=7, byrow=TRUE))
lf = list.files(path = "Simulation_Graphics/") # image filenames
name_func  = function(x){paste("Simulation_Graphics/", x, sep="")}
PNG_names = unlist(lapply(lf, name_func))
PNGs      = lapply(PNG_names, readPNG)
# do the plotting
for(i in 1:length(PNG_names)) {
  if(i == 5){
    # grid.newpage()
  }
  plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
  rasterImage(PNGs[[i]],0,0,1,1)
}
# write to PDF
dev.off()
# dev.print(pdf, "output.pdf")

