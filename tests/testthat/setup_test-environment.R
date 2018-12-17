# copy sample data
file.copy(system.file("extdata", package = "QuasR"), ".", recursive = TRUE)

# create cluster object
clObj <- parallel::makeCluster(2L)
