/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Cartisen GTO to spheric or spinor GTO transformation
 */

/*************************************************
 *
 * transform matrix
 *
 *************************************************/
#include <complex.h>
#include "g1e.h"
#include "simd.h"

void c2s_sph_1e_simd(__MD *opij, __MD *gctr, int *dims, CINTEnvVars *envs, double *cache);

void c2s_cart_1e_simd(__MD *opij, __MD *gctr, int *dims, CINTEnvVars *envs, double *cache);

void c2s_mset0(__MD *out, int *dims, int *counts);
