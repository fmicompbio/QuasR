#include "hello_world.h"


SEXP hello_world(SEXP name)
{
  if (!Rf_isString(name) || 1 != Rf_length(name)){
    Rf_error("'name' must be character(1)");
  }

  const char * cname =  Rf_translateChar(STRING_ELT(name, 0));
  Rprintf("Hallo %s\n", cname);

  return name;
}
