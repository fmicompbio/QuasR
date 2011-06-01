.requirePkg <- function(pkgname, strict=TRUE)
{
    ## check package installed
    if(!require(pkgname, character.only=TRUE, quietly=TRUE))
        stop("Can not find the '", pkgname, "' package.")
    ## check package version
    current <- installed.packages()[pkgname, 'Version']
    pkgs <- unlist(strsplit(installed.packages()['QuasR','Enhances'], ","))
    targetpkg <- grep(pkgname, pkgs, value=T)
    ## check if target is equal current version
    if(length(grep(paste("\\(== ", current, "\\)", sep=""), targetpkg)) == 0){
        if(strict)
            stop(gettextf("package '%s %s' was found, but '%s' is required by 'QuasR'",
                          pkgname, current, targetpkg))
        else
            warning(gettextf("package '%s' %s was found, but '%s' is required by 'QuasR'",
                             pkgname, current, targetpkg))
    }
} 

.baseFileName <- function(filename)
{
    fn <- unlist(strsplit(basename(filename),"\\."))
    if(length(fn) == 1L)
        return(fn)
    else
        return(paste(fn[-length(fn)], collapse="."))
}

.fileExtension <- function(filename)
{
    fn <- unlist(strsplit(basename(filename),"\\."))
    return(fn[length(fn)])
}


.progressReport <- function(msg, phase=0)
{
    if(phase==0)
        cat("done\n")
    if(phase<=0)
        msg <- paste(msg, "...", sep="")
    if(phase>0)
        cat("done\n")
    cat(msg)
    if(phase>0)
        cat("\n")
}
    
   
