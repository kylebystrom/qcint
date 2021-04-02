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

#include "simd.h"
#if !defined HAVE_DEFINED_CINTENVVARS_R_H
#define HAVE_DEFINED_CINTENVVARS_R_H
// ref to CINTinit_int1e_EnvVars, CINTinit_int2e_EnvVars
typedef struct {
        int *atm;
        int *bas;
        double *env;
        int *shls;
        int natm;
        int nbas;

        int i_l;
        int j_l;
        int k_l;
        int l_l;
        int nfi;  // number of cartesion components
        int nfj;
        int nfk;
        int nfl;
        int nf;  // = nfi*nfj*nfk*nfl;
        int _padding;
        int x_ctr[4];

        int gbits;
        int ncomp_e1; // = 1 if spin free, = 4 when spin included, it
        int ncomp_e2; // corresponds to POSX,POSY,POSZ,POS1, see cint_const.h
        int ncomp_tensor; // e.g. = 3 for gradients

        /* values may diff based on the g0_2d4d algorithm */
        int li_ceil; // power of x, == i_l if nabla is involved, otherwise == i_l
        int lj_ceil;
        int lk_ceil;
        int ll_ceil;
        int g_stride_i; // nrys_roots * shift of (i++,k,l,j)
        int g_stride_k; // nrys_roots * shift of (i,k++,l,j)
        int g_stride_l; // nrys_roots * shift of (i,k,l++,j)
        int g_stride_j; // nrys_roots * shift of (i,k,l,j++)
        int nrys_roots;
        int g_size;  // ref to cint2e.c g = malloc(sizeof(double)*g_size)

        int g2d_ijmax;
        int g2d_klmax;
        double common_factor;
        double fac;
        double expcutoff;
        double rirj[3]; // diff by sign in different g0_2d4d algorithm
        double rkrl[3];
        double *rx_in_rijrx;
        double *rx_in_rklrx;

        double *ri;
        double *rj;
        ALIGNMM double rk[12];

        int (*f_g0_2e)();
        void (*f_g0_2d4d)();
        void (*f_gout)();

        /* values are assigned during calculation */
        int *idx;
        double ai;
        double aj;
        double ak;
        double al;
        double rij[3];
        double rijrx[3];
        double aij;
        double rkl[3];
        double rklrx[3];
        double akl;
} CINTEnvVarsR;
#endif

#define RYS_ROOTS       6
void CINTinit_int1e1r_EnvVars(CINTEnvVarsR *envs, int *ng, int *shls,
                              int *atm, int natm, int *bas, int nbas, double *env);

void CINTg2cr_index_xyz(int *idx, CINTEnvVarsR *envs);
void CINTg1e1r_index_xyz(int *idx, CINTEnvVarsR *envs);

void CINTg1e1r_rinv(double *g, CINTEnvVarsR *envs);

void CINTrprim_to_ctr(__MD *gc, int nf, __MD *gp,
                     int inc, int nprim, int nctr, double *pcoeff);

double CINTcommon_fac_sp_r(int l);
