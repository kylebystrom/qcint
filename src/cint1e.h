/*
 * Qcint is a general GTO integral library for computational chemistry
 * Copyright (C) 2014- Qiming Sun <osirpt.sun@gmail.com>
 *
 * This file is part of Qcint.
 *
 * Qcint is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <complex.h>
#include "simd.h"

int CINT1e_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
               double *cache, void (*f_c2s)());
int CINT1e_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs,
                      CINTOpt *opt, double *cache, void (*f_c2s)());
int int1e_cache_size(CINTEnvVars *envs);

int CINT3c1e_cart_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                      double *cache);
int CINT3c1e_spheric_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                         double *cache);
int CINT3c1e_spinor_drv(double complex *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)());
int int3c1e_cache_size(CINTEnvVars *envs);
