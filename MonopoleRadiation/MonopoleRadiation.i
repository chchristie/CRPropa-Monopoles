/* name of the plugin: MonopoleRadiation*/
%module(directors="1", threads="1", allprotected="1") MonopoleRadiation

/* Exceptions required */
%include "exception.i"

/*  define headers to include into the wrapper. These are the plugin headers
 *  and the CRPRopa headers.
 */
%{
#include "CRPropa.h"
#include "MonopoleRadiation.h"
%}

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"

/* include plugin parts to generate wrappers for */
%include "MonopoleRadiation.h"





