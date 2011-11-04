###
###
###

#TODO define index list only once
# loadIndex <- function(){
#     
#     index <- list(
#         name="@GENOMENAME@@",
#         shortname="@PROVIDERVERSION@",
#         path=system.file("alignerindex", package="@PKGNAME@", "@GENOMENAME@"),
#         path=file.path(indexDir, "@GENOMENAME@"),
#         aligner="@ALIGNER@",
#         alignerversion="@ALIGNERVERSION@",
#         organism="@ORGANISM@",
#         provider="@PROVIDER@",
#         providerversion="@PROVIDERVERSION@",
#         sourceurl="@SRCDATAFILES@",
#         md5sum="@MD5SUM@"
#     )
# }

.onLoad <- function(libname, pkgname)
{
    indexDir <- system.file("alignerIndex", package=pkgname, lib.loc=libname)
    index <- list(
        name="@GENOMENAME@",
        shortname="@PROVIDERVERSION@",
        path=file.path(indexDir, "@GENOMENAME@"),
        aligner="@ALIGNER@",
        alignerversion="@ALIGNERVERSION@",
        organism="@ORGANISM@",
        provider="@PROVIDER@",
        providerversion="@PROVIDERVERSION@",
        sourceurl="@SRCDATAFILES@",
        md5sum="@MD5SUM@"
    )
#    objname <- "@ALIGNERINDEXNAME@"
    objname <- "@PROVIDERVERSION@@ALIGNER@Index"
    ns <- asNamespace(pkgname)
    assign(objname, index, envir=ns)
    namespaceExport(ns, objname)
}

