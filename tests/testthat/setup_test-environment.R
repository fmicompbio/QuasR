# copy sample data
file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

# create cluster object
cl <- parallel::makeCluster(2L)
