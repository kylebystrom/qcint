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

#include <string.h>
#include <math.h>
#include <assert.h>
#include "cint_bas.h"
#include "simd.h"
#include "misc.h"
#include "g1e1r.h"
#include "rys_roots.h"

#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *GX = G; \
        type *GY = G + envs->g_size     * SIMDD; \
        type *GZ = G + envs->g_size * 2 * SIMDD

void CINTinit_int1e1r_EnvVars(CINTEnvVarsR *envs, int *ng, int *shls,
                              int *atm, int natm, int *bas, int nbas, double *env)
{
        int i, j;
        double *rk;

        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        const int i_sh = shls[0];
        const int j_sh = shls[1];
        const int k_sh = shls[2]; // should have a third shell number for real-space
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nf = envs->nfi * envs->nfj;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));

        for (i = 0; (i < SIMDD) && (k_sh+i < nbas); i++) {
            rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh+i));
            envs->rk[i] = rk[0];
            envs->rk[i+SIMDD] = rk[1];
            envs->rk[i+2*SIMDD] = rk[2];
        }

        envs->common_factor = 1;
        if (env[PTR_EXPCUTOFF] == 0) {
                envs->expcutoff = EXPCUTOFF;
        } else {
                envs->expcutoff = MAX(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
        }

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        if (ng[RYS_ROOTS] > 0) {
                envs->nrys_roots = ng[RYS_ROOTS];
        } else {
                envs->nrys_roots = (envs->li_ceil + envs->lj_ceil)/2 + 1;
        }

        int dli = envs->li_ceil + envs->lj_ceil + 1;
        int dlj = envs->lj_ceil + 1;
        // TODO these differ from libcint/g1e.c
        envs->g_stride_i = envs->nrys_roots;
        envs->g_stride_j = dli * envs->nrys_roots;
        envs->g_size     = dli * dlj * envs->nrys_roots;

        envs->rirj[0] = envs->ri[0] - envs->rj[0];
        envs->rirj[1] = envs->ri[1] - envs->rj[1];
        envs->rirj[2] = envs->ri[2] - envs->rj[2];

        envs->lk_ceil = 1;
        envs->ll_ceil = 1;
        envs->g_stride_k = 0;
        envs->g_stride_l = 0;
}

void CINTg2cr_index_xyz(int *idx, CINTEnvVarsR *envs)
{
        int i_l = envs->i_l;
        int j_l = envs->j_l;
        int nfi = envs->nfi;
        int nfj = envs->nfj;
        int di = envs->g_stride_i;
        int dj = envs->g_stride_j;
        int i, j, n;
        int ofx, ofjx;
        int ofy, ofjy;
        int ofz, ofjz;
        int i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
        int j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];

        CINTcart_comp(i_nx, i_ny, i_nz, i_l);
        CINTcart_comp(j_nx, j_ny, j_nz, j_l);

        ofx = 0;
        ofy = envs->g_size;
        ofz = envs->g_size * 2;
        n = 0;
        for (j = 0; j < nfj; j++) {
                ofjx = ofx + dj * j_nx[j];
                ofjy = ofy + dj * j_ny[j];
                ofjz = ofz + dj * j_nz[j];
                switch (i_l) {
                case 0:
                        idx[n+0] = ofjx;
                        idx[n+1] = ofjy;
                        idx[n+2] = ofjz;
                        n += 3;
                        break;
                case 1:
                        idx[n+0] = ofjx + di;
                        idx[n+1] = ofjy;
                        idx[n+2] = ofjz;
                        idx[n+3] = ofjx;
                        idx[n+4] = ofjy + di;
                        idx[n+5] = ofjz;
                        idx[n+6] = ofjx;
                        idx[n+7] = ofjy;
                        idx[n+8] = ofjz + di;
                        n += 9;
                        break;
                case 2:
                        idx[n+0 ] = ofjx + di*2;
                        idx[n+1 ] = ofjy;
                        idx[n+2 ] = ofjz;
                        idx[n+3 ] = ofjx + di;
                        idx[n+4 ] = ofjy + di;
                        idx[n+5 ] = ofjz;
                        idx[n+6 ] = ofjx + di;
                        idx[n+7 ] = ofjy;
                        idx[n+8 ] = ofjz + di;
                        idx[n+9 ] = ofjx;
                        idx[n+10] = ofjy + di*2;
                        idx[n+11] = ofjz;
                        idx[n+12] = ofjx;
                        idx[n+13] = ofjy + di;
                        idx[n+14] = ofjz + di;
                        idx[n+15] = ofjx;
                        idx[n+16] = ofjy;
                        idx[n+17] = ofjz + di*2;
                        n += 18;
                        break;
                default:
                        for (i = 0; i < nfi; i++) {
                                idx[n+0] = ofjx + di * i_nx[i];
                                idx[n+1] = ofjy + di * i_ny[i];
                                idx[n+2] = ofjz + di * i_nz[i];
                                n += 3;
                        }
                }
        }
}
void CINTg1e1r_index_xyz(int *idx, CINTEnvVarsR *envs)
{
        CINTg2cr_index_xyz(idx, envs);
}

void CINTg1e1r_rinv(double *g, CINTEnvVarsR *envs)
{
        int nmax = envs->li_ceil + envs->lj_ceil;
        int lj = envs->lj_ceil;
        int di = envs->g_stride_i;
        int dj = envs->g_stride_j;
        int nrys_roots = envs->nrys_roots;
        int *atm = envs->atm;
        double *env = envs->env;
        double *RESTRICT ri = envs->ri;
        double *RESTRICT rij = envs->rij;
        double *RESTRICT rirj = envs->rirj;
        __MD tau;
        ALIGNMM double x[SIMDD];
        ALIGNMM double u[MXRYSROOTS*SIMDD];
        ALIGNMM double w[MXRYSROOTS*SIMDD];
        ALIGNMM double t2[MXRYSROOTS*SIMDD];
        ALIGNMM double r0ix[MXRYSROOTS*SIMDD];
        ALIGNMM double r0iy[MXRYSROOTS*SIMDD];
        ALIGNMM double r0iz[MXRYSROOTS*SIMDD];
        double *p0x, *p0y, *p0z, *p1x, *p1y, *p1z;
        int i, j, k, n, ptr;
        __MD crij[3];
        __MD r0, r1, r2, rt2, fac1;
        double aij;

        aij = envs->ai + envs->aj;
        tau = MM_SET1(1.0);
        // TODO maybe can just set tau to some value?

        fac1 = 2*M_PI * envs->fac * tau / aij;
        crij[0] = MM_LOAD(envs->rk + 0*SIMDD) - MM_SET1(rij[0]);
        crij[1] = MM_LOAD(envs->rk + 1*SIMDD) - MM_SET1(rij[1]);
        crij[2] = MM_LOAD(envs->rk + 2*SIMDD) - MM_SET1(rij[2]);
        MM_STORE(x, MM_SET1(aij) * tau * tau * SQUARE(crij));
        for (i = 0; i < nrys_roots*SIMDD; i++) {
                u[i] = 0;
                w[i] = 0;
        }
        CINTrys_roots(nrys_roots, x, u, w, SIMDD);

        double *gx = g;
        double *gy = g + envs->g_size     * SIMDD;
        double *gz = g + envs->g_size * 2 * SIMDD;
        r1 = MM_SET1(1.);
        for (i = 0; i < nrys_roots; i++) {
                MM_STORE(gx+i*SIMDD, r1);
                MM_STORE(gy+i*SIMDD, r1);
                MM_STORE(gz+i*SIMDD, fac1 * MM_LOAD(w+i*SIMDD));
        }
        if (envs->g_size == 1) {
                return;
        }

        r2 = tau * tau;
        r1 = MM_SET1(1.);
        r0 = MM_SET1(0.5);
        for (i = 0; i < nrys_roots; i++) {
                rt2 = MM_DIV(r2 * MM_LOAD(u+i*SIMDD), r1 + MM_LOAD(u+i*SIMDD));
                MM_STORE(r0ix+i*SIMDD, MM_SET1(rij[0]) + rt2 * crij[0] - MM_SET1(ri[0]));
                MM_STORE(r0iy+i*SIMDD, MM_SET1(rij[1]) + rt2 * crij[1] - MM_SET1(ri[1]));
                MM_STORE(r0iz+i*SIMDD, MM_SET1(rij[2]) + rt2 * crij[2] - MM_SET1(ri[2]));
                MM_STORE(t2+i*SIMDD, MM_DIV(r0 - r0 * rt2, MM_SET1(aij)));
        }

        if (nmax > 0) {
                p0x = gx + di * SIMDD;
                p0y = gy + di * SIMDD;
                p0z = gz + di * SIMDD;
                for (n = 0; n < nrys_roots; n++) {
//gx[1] = -rir0[0] * gx[0];
//gy[1] = -rir0[1] * gy[0];
//gz[1] = -rir0[2] * gz[0];
MM_STORE(p0x+n*SIMDD, MM_LOAD(r0ix+n*SIMDD) * MM_LOAD(gx+n*SIMDD));
MM_STORE(p0y+n*SIMDD, MM_LOAD(r0iy+n*SIMDD) * MM_LOAD(gy+n*SIMDD));
MM_STORE(p0z+n*SIMDD, MM_LOAD(r0iz+n*SIMDD) * MM_LOAD(gz+n*SIMDD));
                }

                p0x = gx + di * SIMDD;
                p0y = gy + di * SIMDD;
                p0z = gz + di * SIMDD;
                p1x = gx - di * SIMDD;
                p1y = gy - di * SIMDD;
                p1z = gz - di * SIMDD;
                for (i = 1; i < nmax; i++) {
                        ptr = i * di;
                        for (n = 0; n < nrys_roots; n++) {
                                r1 = MM_LOAD(t2+n*SIMDD) * MM_SET1(i);
//gx[i+1] = 0.5 * (1 - t2) / aij * i * gx[i-1] - rir0[0] * gx[i];
//gy[i+1] = 0.5 * (1 - t2) / aij * i * gy[i-1] - rir0[1] * gy[i];
//gz[i+1] = 0.5 * (1 - t2) / aij * i * gz[i-1] - rir0[2] * gz[i];
MM_STORE(p0x+(ptr+n)*SIMDD, r1 * MM_LOAD(p1x+(ptr+n)*SIMDD) + MM_LOAD(r0ix+n*SIMDD) * MM_LOAD(gx+(ptr+n)*SIMDD));
MM_STORE(p0y+(ptr+n)*SIMDD, r1 * MM_LOAD(p1y+(ptr+n)*SIMDD) + MM_LOAD(r0iy+n*SIMDD) * MM_LOAD(gy+(ptr+n)*SIMDD));
MM_STORE(p0z+(ptr+n)*SIMDD, r1 * MM_LOAD(p1z+(ptr+n)*SIMDD) + MM_LOAD(r0iz+n*SIMDD) * MM_LOAD(gz+(ptr+n)*SIMDD));
                        }
                }

                p0x = gx  - dj * SIMDD;
                p0y = gy  - dj * SIMDD;
                p0z = gz  - dj * SIMDD;
                p1x = p0x + di * SIMDD;
                p1y = p0y + di * SIMDD;
                p1z = p0z + di * SIMDD;
                for (j = 1; j <= lj; j++) {
                        for (i = 0; i <= nmax - j; i++) {
                                ptr = dj * j + i * di;
                                for (n = ptr; n < ptr+nrys_roots; n++) {
//gx[i] = gx[i+1-dj] + rirj[0] * gx[i-dj];
//gy[i] = gy[i+1-dj] + rirj[1] * gy[i-dj];
//gz[i] = gz[i+1-dj] + rirj[2] * gz[i-dj];
MM_STORE(gx+n*SIMDD, MM_LOAD(p1x+n*SIMDD) + MM_SET1(rirj[0]) * MM_LOAD(p0x+n*SIMDD));
MM_STORE(gy+n*SIMDD, MM_LOAD(p1y+n*SIMDD) + MM_SET1(rirj[1]) * MM_LOAD(p0y+n*SIMDD));
MM_STORE(gz+n*SIMDD, MM_LOAD(p1z+n*SIMDD) + MM_SET1(rirj[2]) * MM_LOAD(p0z+n*SIMDD));
                                }
                        }
                }
        }
}

/*
 * to optimize memory copy in cart2sph.c, remove the common factor for s
 * and p function in cart2sph
 */
double CINTcommon_fac_sp_r(int l)
{
        switch (l) {
                case 0: return 0.282094791773878143;
                case 1: return 0.488602511902919921;
                default: return 1;
        }
}


void CINTrprim_to_ctr(__MD *gc, int nf, __MD *gp,
                      int inc, int nprim, int nctr, double *coeff)
{
        int n, i, k;
        __MD *pgc = gc;
        double c;

        for (i = 0; i < inc; i++) {
                //dger(nf, nctr, 1.d0, gp(i+1), inc, env(ptr), nprim, gc(1,i*nctr+1), nf)
                for (n = 0; n < nctr; n++) {
                        c = coeff[nprim*n];
                        if (c != 0) {
                                for (k = 0; k < nf; k++) {
                                        pgc[k] += MM_SET1(c) * gp[k*inc+i];
                                }
                        }
                        // next cgto block
                        pgc += nf;
                }
        }
}
