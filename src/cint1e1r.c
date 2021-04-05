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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "misc.h"
#include "g1e1r.h"
#include "optimizer.h"
#include "cint1e1r.h"
#include "cart2sph_simd.h"
#include "c2f.h"

#define SHLTYPi       0
#define SHLTYPj       1


#define ALIAS_ADDR_IF_EQUAL(x, y) \
        if (y##_ctr == 1) { \
                gctr[SHLTYP##x] = gctr[SHLTYP##y]; \
                x##empty = y##empty; \
        } else { \
                gctr[SHLTYP##x] = g1; \
                g1 += len##x; \
        }

#define PRIM2CTR(ctrsymb, gp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        fp2c[np2c] = CINTprim_to_ctr_0; \
                } else { \
                        fp2c[np2c] = CINTprim_to_ctr_1; \
                } \
                shltyp[np2c] = SHLTYP##ctrsymb; \
                gprim[np2c] = gp; \
                iprim[np2c] = ctrsymb##p; \
                np2c++; \
        } \
        *ctrsymb##empty = 0; \

#define POP_PRIM2CTR \
        for (i = 0; i < np2c; i++) { \
                it = shltyp[i]; \
                im = iprim[i]; \
                (*(fp2c[i]))(gctr[it], gprim[i], coeff[it]+im, \
                             ngp[it], x_prim[it], x_ctr[it], \
                             non0ctr[it][im], non0idx[it]+im*x_ctr[it]); \
        } \
        cum = 0; \
        np2c = 0;

#define PUSH(RIJ, EXPIJ) \
        if (cum == SIMDD) { \
                (*envs->f_gout)(gout, g, idx, envs, cum); \
                POP_PRIM2CTR; \
        } \
        envs->ai = ai[ip]; \
        envs->aj = aj[jp]; \
        envs->rij[0] = *(RIJ+0); \
        envs->rij[1] = *(RIJ+1); \
        envs->rij[2] = *(RIJ+2); \
        fac1i = fac1j * EXPIJ; \
        envs->fac = fac1i; \
        if (*iempty) { \
                fp2c[np2c] = CINTiprim_to_ctr_0; \
                *iempty = 0; \
        } else { \
                fp2c[np2c] = CINTiprim_to_ctr_1; \
        } \
        gprim[np2c] = gout + cum * ngp[0]; \
        iprim[np2c] = ip; \
        shltyp[np2c] = 0; \
        cum++; \
        np2c++;

#define INITSIMD \
        int cum = 0; \
        int np2c = 0; \
        double *gprim[SIMDD*2]; \
        int shltyp[SIMDD*2]; \
        int iprim[SIMDD*2]; \
        void (*fp2c[SIMDD*2])(); \
        MM_STORE(envs->ai, MM_SET1(1.)); \
        MM_STORE(envs->aj, MM_SET1(1.)); \
        MM_STORE(envs->fac, MM_SET1(0.));

#define RUN_REST \
        if (cum > 0) { \
                (*envs->f_gout)(gout, g, idx, envs, cum); \
        } \
        POP_PRIM2CTR


int CINT1e1r_loop_nopt(__MD *gctr, CINTEnvVarsR *envs, double *cache)
{
        int *shls  = envs->shls;
        int *atm = envs->atm;
        int *bas = envs->bas;
        double *env = envs->env;
        int i_sh = shls[0];
        int j_sh = shls[1];
        int i_l = envs->i_l;
        int j_l = envs->j_l;
        int i_ctr = envs->x_ctr[0];
        int j_ctr = envs->x_ctr[1];
        int i_prim = bas(NPRIM_OF, i_sh);
        int j_prim = bas(NPRIM_OF, j_sh);
        int nf = envs->nf;
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *ri = envs->ri;
        double *rj = envs->rj;
        double *ai = env + bas(PTR_EXP, i_sh);
        double *aj = env + bas(PTR_EXP, j_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *cj = env + bas(PTR_COEFF, j_sh);
        int ip, jp, i, n;
        int has_value = 0;
        double tau;
        __MD *cr;
        double x, u[MXRYSROOTS], w[MXRYSROOTS];
        int *idx = malloc(sizeof(int) * nf * 3);
        double rij[3], aij, dij, eij, rrij, t2;
        double fac1j, fac1i;
        __MD *g, *gout, *gctri;
        MALLOC_INSTACK(g, envs->g_size * 3 * ((1<<envs->gbits)+1)); // +1 as buffer
        MALLOC_INSTACK(gout, nf * n_comp);
        MALLOC_INSTACK(gctri, nf * i_ctr * n_comp);
        double expcutoff = envs->expcutoff;

        CINTg1e1r_index_xyz(idx, envs);

        rrij = CINTsquare_dist(ri, rj);
        double common_factor = envs->common_factor
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l);

        for (jp = 0; jp < j_prim; jp++) {
                envs->aj = aj[jp];
                if (j_ctr == 1) {
                        fac1j = common_factor;// * cj[jp];
                } else {
                        fac1j = common_factor;
                }
                CINTdset0(nf * i_ctr * n_comp * SIMDD, (double *) gctri);
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        aij = 1 / (ai[ip] + aj[jp]);
                        eij = (ai[ip] * aj[jp] * aij) * rrij;
                        if (eij > expcutoff)
                                continue;
                        has_value = 1;

                        envs->rij[0] = (ai[ip] * ri[0] + aj[jp] * rj[0]) * aij;
                        envs->rij[1] = (ai[ip] * ri[1] + aj[jp] * rj[1]) * aij;
                        envs->rij[2] = (ai[ip] * ri[2] + aj[jp] * rj[2]) * aij;
                        fac1i = fac1j * exp(-eij);
                        envs->fac = fac1i;

                        CINTdset0(nf * n_comp * SIMDD, (double *) gout);
                        (*envs->f_gout)(gout, g, idx, envs);

                        n = nf * n_comp;
                        CINTrprim_to_ctr(gctri, n, gout, 1, i_prim, i_ctr, ci+ip);
                }
                n = nf * i_ctr;
                CINTrprim_to_ctr(gctr, n, gctri, n_comp, j_prim, j_ctr, cj+jp);
        }
        free(idx);
        return has_value;
}

int CINT1e1r_loop(__MD *out, CINTEnvVarsR *envs, CINTOpt *opt, double *cache) {
    return CINT1e1r_loop_nopt(out, envs, cache);
}

// little endian on x86
//typedef union {
//    double d;
//    unsigned short s[4];
//} type_IEEE754;
static double approx_log(double x)
{
        //type_IEEE754 y;
        //y.d = x;
        //return ((double)(y.s[3] >> 4) - 1023) * 0.7;
        //return log(x);
        return 2.5;
}

int int1e1r_cache_size(CINTEnvVarsR *envs)
{
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        int leng = envs->g_size*3*((1<<envs->gbits)+1)*SIMDD;
        int len0 = envs->nf*n_comp * SIMDD;
        int cache_size = MAX(leng+len0*2+nc*n_comp*2,
                             nc*n_comp + envs->nf*8*OF_CMPLX) + SIMDD*2;
        return cache_size;
}

int CINT1e1r_drv(__MD *out, int *dims, CINTEnvVarsR *envs, CINTOpt *opt,
               double *cache, void (*f_c2s)())
{
        if (out == NULL) {
                return int1e1r_cache_size(envs);
        }
        int *x_ctr = envs->x_ctr;
        int nc = envs->nf * x_ctr[0] * x_ctr[1];
        int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
        double *stack = NULL;
        if (cache == NULL) {
                int cache_size = int1e1r_cache_size(envs);
                stack = _mm_malloc(sizeof(double)*cache_size*SIMDD, sizeof(double)*SIMDD);
                cache = stack;
        }
        __MD *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);
        CINTdset0(nc*n_comp*SIMDD, (double *) gctr);

        int n, has_value;
        if (opt != NULL) {
                has_value = CINT1e1r_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT1e1r_loop_nopt(gctr, envs, cache);
        }

        int counts[4];
        if (dims == NULL) {
                dims = counts;
        }
        if (f_c2s == c2s_sph_1e_simd) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->j_l*2+1) * x_ctr[1];
        } else {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfj * x_ctr[1];
        }
        counts[2] = 1;
        counts[3] = 1;
        int nout = dims[0] * dims[1];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_mset0(out+nout*n, dims, counts);
                }
        }

        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}

void CINTgout1e1r_rinv(__MD *gout, __MD *g, int *idx, CINTEnvVarsR *envs)
{
        int nf = envs->nf;
        int nfc = nf;
        __MD *gtmp = gout;// + nf;
        int nrys_roots = envs->nrys_roots;
        int n, i;
        __MD *gx, *gy, *gz;
        __MD r0;
        __MD zero = MM_SET1(0.0);
        for (n = 0; n < nf; n++) {
                gtmp[n] = zero;
        }
        CINTg1e1r_rinv(g, envs);
        for (n = 0; n < nf; n++) {
                gx = g + idx[n*3+0];
                gy = g + idx[n*3+1];
                gz = g + idx[n*3+2];
                r0 = zero;//gtmp[n];
                for (i = 0; i < nrys_roots; i++) {
                        r0 += gx[i] * gy[i] * gz[i];
                }
                gtmp[n] = r0;
        }
        //CINTsort_gout(gout, gtmp, nfc, SIMDD);
}

int int1e1r_rinv_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVarsR envs;
        CINTinit_int1e1r_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e1r_rinv;
        return CINT1e1r_drv((__MD*)out, dims, &envs, opt, cache, &c2s_sph_1e_simd);
}
void int1e1r_rinv_optimizer(CINTOpt **opt, int *atm, int natm,
                         int *bas, int nbas, double *env)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTall_1e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
int int1e1r_rinv_cart(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
        int ng[] = {0, 0, 0, 0, 0, 1, 0, 1};
        CINTEnvVarsR envs;
        CINTinit_int1e1r_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout1e1r_rinv;
        return CINT1e1r_drv((__MD*)out, dims, &envs, opt, cache, &c2s_cart_1e_simd);
}


#define SOME_CINT1E(NAME) \
int c##NAME##_cart(double *out, int *shls, int *atm, int natm, \
            int *bas, int nbas, double *env) { \
        return NAME##_cart(out, NULL, shls, atm, natm, bas, nbas, env, NULL, NULL); \
} \
int c##NAME##_sph(double *out, int *shls, int *atm, int natm, \
            int *bas, int nbas, double *env) { \
        return NAME##_sph(out, NULL, shls, atm, natm, bas, nbas, env, NULL, NULL); \
} \

SOME_CINT1E(int1e1r_rinv)

