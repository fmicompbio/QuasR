###
###
###

loadIndex <- function(){
    index <- list(
        name="@PROVIDERVERSION@-@ALIGNER@",
        path=file.path(system.file("index", package="@PKGNAME@", "@GENOMENAME@")),
        aligner="@ALIGNER@",
        alignerversion="@ALIGNERVERSION@",
        organism="@ORGANISM@",
        provider="@PROVIDER@",
        providerversion="@PROVIDERVERSION@",
        sourceurl="@SRCDATAFILES@"
    )
}

.onLoad <- function(libname, pkgname)
{
    indexDir <- system.file("alignerIndex", package=pkgname, lib.loc=libname)
    index <- list(
        name="@PROVIDERVERSION@",
        path=file.path(indexDir, "@GENOMENAME@"),
        aligner="@ALIGNER@",
        alignerversion="@ALIGNERVERSION@",
        organism="@ORGANISM@",
        provider="@PROVIDER@",
        providerversion="@PROVIDERVERSION@",
        sourceurl="@SRCDATAFILES@"
    )
#    objname <- "@ALIGNERINDEXNAME@"
    objname <- "@PROVIDERVERSION@@ALIGNER@Index"
    ns <- asNamespace(pkgname)
    assign(objname, index, envir=ns)
    namespaceExport(ns, objname)
}

