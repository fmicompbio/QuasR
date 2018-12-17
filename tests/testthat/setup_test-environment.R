# copy sample data
file.copy(system.file("extdata", package = "QuasR"), ".", recursive = TRUE)

# create cluster object
cl <- parallel::makeCluster(2L)
