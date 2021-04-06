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

/*
 * optimizer for 2e integrals.  Note if CINT2e_drv is only called a few
 * hundred times, this optimizer cannot really speed up the integration. 
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cint_const.h"
#include "cint_bas.h"
#include "g1e.h"
#include "g2e.h"
#include "optimizer.h"
#include "rys_roots.h"
#include "misc.h"

// generate caller to CINTinit_2e_optimizer for each type of function
void CINTinit_2e_optimizer(CINTOpt **opt, int *atm, int natm,
                           int *bas, int nbas, double *env)
{
        CINTOpt *opt0 = (CINTOpt *)malloc(sizeof(CINTOpt));
        opt0->index_xyz_array = NULL;
        opt0->non0ctr = NULL;
        opt0->sortedidx = NULL;
        opt0->nbas = nbas;
        opt0->log_max_coeff = NULL;
        opt0->pairdata = NULL;
        *opt = opt0;
}
void CINTinit_optimizer(CINTOpt **opt, int *atm, int natm,
                        int *bas, int nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
}

void CINTdel_2e_optimizer(CINTOpt **opt)
{
        CINTOpt *opt0 = *opt;
        if (opt0 == NULL) { // when opt is created by CINTno_optimizer
                return;
        }

        if (opt0->index_xyz_array != NULL) {
                free(opt0->index_xyz_array[0]);
                free(opt0->index_xyz_array);
        }

        if (opt0->non0ctr != NULL) {
                free(opt0->sortedidx[0]);
                free(opt0->sortedidx);
                free(opt0->non0ctr[0]);
                free(opt0->non0ctr);
        }

        if (opt0->log_max_coeff != NULL) {
                free(opt0->log_max_coeff[0]);
                free(opt0->log_max_coeff);
        }

        CINTdel_pairdata_optimizer(opt0);

        free(opt0);
        *opt = NULL;
}
void CINTdel_optimizer(CINTOpt **opt)
{
        CINTdel_2e_optimizer(opt);
}

void CINTno_optimizer(CINTOpt **opt, int *atm, int natm,
                      int *bas, int nbas, double *env)
{
        *opt = NULL;
}

static int _make_fakebas(int *fakebas, int *bas, int nbas, double *env)
{
        int i;
        int max_l = 0;
        for (i = 0; i < nbas; i++) {
                max_l = MAX(max_l, bas(ANG_OF,i));
        }

        int fakenbas = max_l + 1;
        for (i = 0; i < BAS_SLOTS*fakenbas; i++) {
                fakebas[i] = 0;
        }
        // fakebas only initializes ANG_OF, since the others does not
        // affect index_xyz
        for (i = 0; i <= max_l; i++) {
                fakebas[BAS_SLOTS*i+ANG_OF] = i;
        }
        return max_l;
}
static int *_allocate_index_xyz(CINTOpt *opt, int max_l, int order)
{
        int i;
        int cumcart = (max_l+1) * (max_l+2) * (max_l+3) / 6;
        int ll = max_l + 1;
        int cc = cumcart;
        for (i = 1; i < order; i++) {
                ll *= LMAX1;
                cc *= cumcart;
        }
        int *buf = malloc(sizeof(int) * cc * 3);
        int **ppbuf = malloc(sizeof(int*) * ll);
        ppbuf[0] = buf;
        for (i = 1; i < ll; i++) {
                ppbuf[i] = NULL;
        }
        opt->index_xyz_array = ppbuf;
        return buf;
}
static void gen_idx(CINTOpt *opt, void (*finit)(), void (*findex_xyz)(),
                    int order, int max_l, int *ng,
                    int *atm, int natm, int *bas, int nbas, double *env)
{
        int i, j, k, l, ptr;
        int fakebas[BAS_SLOTS*LMAX1];
        int max_l1 = _make_fakebas(fakebas, bas, nbas, env);
        if (max_l == 0) {
                max_l = max_l1;
        } else {
                max_l = MIN(max_l, max_l1);
        }
        int fakenbas = max_l+1;
        int *buf = _allocate_index_xyz(opt, max_l, order);

        CINTEnvVars envs;
        int shls[4];
        if (order == 2) {
                for (i = 0; i <= max_l; i++) {
                for (j = 0; j <= max_l; j++) {
                        shls[0] = i; shls[1] = j;
                        (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                        ptr = i*LMAX1 + j;
                        opt->index_xyz_array[ptr] = buf;
                        (*findex_xyz)(opt->index_xyz_array[ptr], &envs);
                        buf += envs.nf * 3;
                } }

        } else if (order == 3) {
                for (i = 0; i <= max_l; i++) {
                for (j = 0; j <= max_l; j++) {
                for (k = 0; k <= max_l; k++) {
                        shls[0] = i; shls[1] = j; shls[2] = k;
                        (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                        ptr = i*LMAX1*LMAX1 + j*LMAX1 + k;
                        opt->index_xyz_array[ptr] = buf;
                        (*findex_xyz)(opt->index_xyz_array[ptr], &envs);
                        buf += envs.nf * 3;
                } } }

        } else {
                for (i = 0; i <= max_l; i++) {
                for (j = 0; j <= max_l; j++) {
                for (k = 0; k <= max_l; k++) {
                for (l = 0; l <= max_l; l++) {
                        shls[0] = i; shls[1] = j; shls[2] = k; shls[3] = l;
                        (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                        ptr = i*LMAX1*LMAX1*LMAX1
                            + j*LMAX1*LMAX1
                            + k*LMAX1
                            + l;
                        opt->index_xyz_array[ptr] = buf;
                        (*findex_xyz)(opt->index_xyz_array[ptr], &envs);
                        buf += envs.nf * 3;
                } } } }
        }
}

void CINTall_1e_optimizer(CINTOpt **opt, int *ng,
                          int *atm, int natm, int *bas, int nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_set_log_maxc(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int1e_EnvVars, &CINTg2c_index_xyz,
                2, 0, ng, atm, natm, bas, nbas, env);
}

void CINTall_1e1r_optimizer(CINTOpt **opt, int *ng,
                            int *atm, int natm, int *bas, int nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_log_maxc(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int1e_EnvVars, &CINTg2c_index_xyz,
                2, 0, ng, atm, natm, bas, nbas, env);
}

void CINTall_2e_optimizer(CINTOpt **opt, int *ng,
                          int *atm, int natm, int *bas, int nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int2e_EnvVars, &CINTg4c_index_xyz,
                4, 6, ng, atm, natm, bas, nbas, env);
}

void CINTall_3c2e_optimizer(CINTOpt **opt, int *ng,
                            int *atm, int natm, int *bas, int nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int3c2e_EnvVars, &CINTg4c_index_xyz,
                3, 0, ng, atm, natm, bas, nbas, env);
}

void CINTall_2c2e_optimizer(CINTOpt **opt, int *ng,
                            int *atm, int natm, int *bas, int nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_set_log_maxc(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int2c2e_EnvVars, &CINTg2c_index_xyz,
                2, 0, ng, atm, natm, bas, nbas, env);
}

void CINTall_3c1e_optimizer(CINTOpt **opt, int *ng,
                            int *atm, int natm, int *bas, int nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int3c1e_EnvVars, &CINTg4c_index_xyz,
                3, 0, ng, atm, natm, bas, nbas, env);
}

#ifdef WITH_F12
void CINTall_2e_stg_optimizer(CINTOpt **opt, int *ng,
                              int *atm, int natm, int *bas, int nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int2e_stg_EnvVars, &CINTg4c_index_xyz,
                4, 6, ng, atm, natm, bas, nbas, env);
}
#endif
 
#ifdef WITH_GTG
void CINTall_2e_gtg_optimizer(CINTOpt **opt, int *ng,
                              int *atm, int natm, int *bas, int nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int2e_gtg_EnvVars, &CINTg4c_index_xyz,
                4, 6, ng, atm, natm, bas, nbas, env);
}

void CINTall_3c2e_gtg_optimizer(CINTOpt **opt, int *ng,
                                int *atm, int natm, int *bas, int nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int3c2e_gtg_EnvVars, &CINTg4c_index_xyz,
                3, 0, ng, atm, natm, bas, nbas, env);
}

void CINTall_2c2e_gtg_optimizer(CINTOpt **opt, int *ng,
                                int *atm, int natm, int *bas, int nbas, double *env)
{
        CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
        CINTOpt_set_log_maxc(*opt, atm, natm, bas, nbas, env);
        CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
        gen_idx(*opt, &CINTinit_int2c2e_gtg_EnvVars, &CINTg2c_index_xyz,
                2, 0, ng, atm, natm, bas, nbas, env);
}
#endif


// little endian on x86
typedef union {
    double d;
    unsigned short s[4];
} type_IEEE754;
// ~4 times faster than built-in log
static inline double approx_log(double x)
{
        type_IEEE754 y;
        y.d = x;
        return ((y.s[3] >> 4) - 1023 + 1) * 0.693145751953125;
}

void CINTOpt_log_max_pgto_coeff(double *log_maxc, double *coeff, int nprim, int nctr)
{
        int i, ip;
        double maxc;
        for (ip = 0; ip < nprim; ip++) {
                maxc = 0;
                for (i = 0; i < nctr; i++) {
                        maxc = MAX(maxc, fabs(coeff[i*nprim+ip]));
                }
                log_maxc[ip] = approx_log(maxc);
        }
}


void CINTOpt_set_log_maxc(CINTOpt *opt, int *atm, int natm,
                          int *bas, int nbas, double *env)
{
        int i, iprim, ictr;
        double *ci;
        size_t tot_prim = 0;
        for (i = 0; i < nbas; i++) {
                tot_prim += bas(NPRIM_OF, i);
        }
        if (tot_prim == 0) {
                return;
        }

        opt->log_max_coeff = malloc(sizeof(double *) * MAX(nbas, 1));
        double *plog_maxc = malloc(sizeof(double) * tot_prim);
        opt->log_max_coeff[0] = plog_maxc;
        for (i = 0; i < nbas; i++) {
                iprim = bas(NPRIM_OF, i);
                ictr = bas(NCTR_OF, i);
                ci = env + bas(PTR_COEFF, i);
                opt->log_max_coeff[i] = plog_maxc;
                CINTOpt_log_max_pgto_coeff(plog_maxc, ci, iprim, ictr);
                plog_maxc += iprim;
        }
}

int CINTset_pairdata(PairData *pairdata, double *ai, double *aj, double *ri, double *rj,
                     double *log_maxci, double *log_maxcj,
                     int li_ceil, int lj_ceil, int iprim, int jprim,
                     double rr_ij, double expcutoff)
{
        int ip, jp, n;
        double aij, eij, cceij;
        double log_rr_ij = (li_ceil+lj_ceil+1) * approx_log(rr_ij+1) / 2;
        PairData *pdata;

        int empty = 1;
        for (n = 0, jp = 0; jp < jprim; jp++) {
                for (ip = 0; ip < iprim; ip++, n++) {
                        aij = 1/(ai[ip] + aj[jp]);
                        eij = rr_ij * ai[ip] * aj[jp] * aij;
                        cceij = eij - log_rr_ij - log_maxci[ip] - log_maxcj[jp];
                        pdata = pairdata + n;
                        pdata->cceij = cceij;
                        if (cceij < expcutoff) {
                                empty = 0;
                                pdata->rij[0] = (ai[ip]*ri[0] + aj[jp]*rj[0]) * aij;
                                pdata->rij[1] = (ai[ip]*ri[1] + aj[jp]*rj[1]) * aij;
                                pdata->rij[2] = (ai[ip]*ri[2] + aj[jp]*rj[2]) * aij;
                                pdata->eij = exp(-eij);
                        } else {
                                pdata->rij[0] = 0;
                                pdata->rij[1] = 0;
                                pdata->rij[2] = 0;
                                pdata->eij = 0;
                        }
                }
        }
        return empty;
}

void CINTOpt_setij(CINTOpt *opt, int *ng,
                   int *atm, int natm, int *bas, int nbas, double *env)
{
        int i, j, ip, jp;
        int iprim, jprim, li, lj;
        double *ai, *aj, *ri, *rj;
        double expcutoff;
        if (env[PTR_EXPCUTOFF] == 0) {
                expcutoff = EXPCUTOFF;
        } else {
                expcutoff = MAX(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
        }

        if (opt->log_max_coeff == NULL) {
                CINTOpt_set_log_maxc(opt, atm, natm, bas, nbas, env);
        }
        double **log_max_coeff = opt->log_max_coeff;
        double *log_maxci, *log_maxcj;

        size_t tot_prim = 0;
        for (i = 0; i < nbas; i++) {
                tot_prim += bas(NPRIM_OF, i);
        }
        if (tot_prim == 0 || tot_prim > MAX_PGTO_FOR_PAIRDATA) {
                return;
        }
        opt->pairdata = malloc(sizeof(PairData *) * MAX(nbas * nbas, 1));
        PairData *pdata = malloc(sizeof(PairData) * tot_prim * tot_prim);
        opt->pairdata[0] = pdata;

        int ijkl_inc;
        if ((ng[IINC]+ng[JINC]) > (ng[KINC]+ng[LINC])) {
                ijkl_inc = ng[IINC] + ng[JINC];
        } else {
                ijkl_inc = ng[KINC] + ng[LINC];
        }

        int empty;
        double rr;
        PairData *pdata0;
        for (i = 0; i < nbas; i++) {
                ri = env + atm(PTR_COORD,bas(ATOM_OF,i));
                ai = env + bas(PTR_EXP,i);
                iprim = bas(NPRIM_OF,i);
                li = bas(ANG_OF,i);
                log_maxci = log_max_coeff[i];

                for (j = 0; j <= i; j++) {
                        rj = env + atm(PTR_COORD,bas(ATOM_OF,j));
                        aj = env + bas(PTR_EXP,j);
                        jprim = bas(NPRIM_OF,j);
                        lj = bas(ANG_OF,j);
                        log_maxcj = log_max_coeff[j];
                        rr = (ri[0]-rj[0])*(ri[0]-rj[0])
                           + (ri[1]-rj[1])*(ri[1]-rj[1])
                           + (ri[2]-rj[2])*(ri[2]-rj[2]);

                        empty = CINTset_pairdata(pdata, ai, aj, ri, rj, log_maxci, log_maxcj,
                                                 li+ijkl_inc, lj, iprim, jprim, rr, expcutoff);
                        if (i == 0 && j == 0) {
                                opt->pairdata[0] = pdata;
                                pdata += iprim * jprim;
                        } else if (!empty) {
                                opt->pairdata[i*nbas+j] = pdata;
                                pdata += iprim * jprim;
                                if (i != j) {
                                        opt->pairdata[j*nbas+i] = pdata;
                                        pdata0 = opt->pairdata[i*nbas+j];
                                        // transpose pairdata
                                        for (ip = 0; ip < iprim; ip++) {
                                        for (jp = 0; jp < jprim; jp++, pdata++) {
                                                memcpy(pdata, pdata0+jp*iprim+ip,
                                                       sizeof(PairData));
                                        } }
                                }
                        } else {
                                opt->pairdata[i*nbas+j] = NOVALUE;
                                opt->pairdata[j*nbas+i] = NOVALUE;
                        }
                }
        }
}

void CINTdel_pairdata_optimizer(CINTOpt *cintopt)
{
        if (cintopt != NULL && cintopt->pairdata != NULL) {
                free(cintopt->pairdata[0]);
                free(cintopt->pairdata);
                cintopt->pairdata = NULL;
        }
}

void CINTOpt_non0coeff_byshell(int *sortedidx, int *non0ctr, double *ci,
                               int iprim, int ictr)
{
        int ip, j, k, kp;
        int zeroidx[ictr];
        for (ip = 0; ip < iprim; ip++) {
                for (j = 0, k = 0, kp = 0; j < ictr; j++) {
                        if (ci[iprim*j+ip] != 0) {
                                sortedidx[k] = j;
                                k++;
                        } else {
                                zeroidx[kp] = j;
                                kp++;
                        }
                }
// Append the index of zero-coeff to sortedidx for function CINTprim_to_ctr_0
                for (j = 0; j < kp; j++) {
                        sortedidx[k+j] = zeroidx[j];
                }
                non0ctr[ip] = k;
                sortedidx += ictr;
        }
}

void CINTOpt_set_non0coeff(CINTOpt *opt, int *atm, int natm,
                           int *bas, int nbas, double *env)
{
        int i, iprim, ictr;
        double *ci;
        size_t tot_prim = 0;
        size_t tot_prim_ctr = 0;
        for (i = 0; i < nbas; i++) {
                tot_prim += bas(NPRIM_OF, i);
                tot_prim_ctr += bas(NPRIM_OF, i) * bas(NCTR_OF,i);
        }
        if (tot_prim == 0) {
                return;
        }

        opt->non0ctr = malloc(sizeof(int *) * MAX(nbas, 1));
        opt->sortedidx = malloc(sizeof(int *) * MAX(nbas, 1));
        int *pnon0ctr = malloc(sizeof(int) * tot_prim);
        int *psortedidx = malloc(sizeof(int) * tot_prim_ctr);
        opt->non0ctr[0] = pnon0ctr;
        opt->sortedidx[0] = psortedidx;
        for (i = 0; i < nbas; i++) {
                iprim = bas(NPRIM_OF, i);
                ictr = bas(NCTR_OF, i);
                ci = env + bas(PTR_COEFF, i);
                opt->non0ctr[i] = pnon0ctr;
                opt->sortedidx[i] = psortedidx;
                CINTOpt_non0coeff_byshell(psortedidx, pnon0ctr, ci, iprim, ictr);
                pnon0ctr += iprim;
                psortedidx += iprim * ictr;
        }
}
