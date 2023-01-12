
simulated_data = create_simulated_data()

full_data = simulated_data[[1]]
usethis::use_data(full_data, overwrite = TRUE)

day_of_observations = simulated_data[[2]]
usethis::use_data(day_of_observations, overwrite = TRUE)

day_dts = simulated_data[[3]]
usethis::use_data(day_dts, overwrite = TRUE)
