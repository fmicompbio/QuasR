# remove sample data
unlink("extdata", recursive = TRUE)

# remove QuasR log files
unlink(list.files(path = ".", pattern = "^QuasR_log_.+\\.txt$"))

# stop cluster object
parallel::stopCluster(clObj)

# restore original values of R_LIBS
Sys.setenv(R_LIBS = oldRlibs)
