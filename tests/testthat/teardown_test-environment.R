# remove sample data
unlink("extdata", recursive = TRUE)

# stop cluster object
parallel::stopCluster(cl)
