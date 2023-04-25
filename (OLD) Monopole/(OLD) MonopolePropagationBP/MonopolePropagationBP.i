/* name of the plugin: MonopolePropagationBP*/
%module(directors="1", threads="1", allprotected="1") MonopolePropagationBP

/* Exceptions required */
%include "exception.i"

/*  define headers to include into the wrapper. These are the plugin headers
 *  and the CRPRopa headers.
 */
%{
#include "CRPropa.h"
#include"../Monopole.h"
#include "MonopolePropagationBP.h"
%}

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"
%import (module="Monopole") "Monopole.i"

/* include plugin parts to generate wrappers for */
%include "MonopolePropagationBP.h"





