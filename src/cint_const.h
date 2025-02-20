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

// global parameters in env

#define PTR_EXPCUTOFF           0
#define PTR_COMMON_ORIG         1
#define PTR_RINV_ORIG           4
#define PTR_RINV_ZETA           7
// omega parameter in range-separated coulomb operator
// LR interaction: erf(omega*r12)/r12 if omega > 0
// SR interaction: erfc(omega*r12)/r12 if omega < 0
#define PTR_RANGE_OMEGA         8
// Yukawa potential and Slater-type geminal e^{-zeta r}
#define PTR_F12_ZETA            9
#define PTR_GTG_ZETA            10
#define PTR_ENV_START           20

// slots of atm
#define CHARGE_OF       0
#define PTR_COORD       1
#define NUC_MOD_OF      2
#define PTR_ZETA        3
#define PTR_FRAC_CHARGE 3
#define RESERVE_ATMLOT1 4
#define RESERVE_ATMLOT2 5
#define ATM_SLOTS       6


// slots of bas
#define ATOM_OF         0
#define ANG_OF          1
#define NPRIM_OF        2
#define NCTR_OF         3
#define KAPPA_OF        4
#define PTR_EXP         5
#define PTR_COEFF       6
#define RESERVE_BASLOT  7
#define BAS_SLOTS       8

// slots of gout
#define POSX            0
#define POSY            1
#define POSZ            2
#define POS1            3
#define POSXX           0
#define POSYX           1
#define POSZX           2
#define POS1X           3
#define POSXY           4
#define POSYY           5
#define POSZY           6
#define POS1Y           7
#define POSXZ           8
#define POSYZ           9
#define POSZZ           10
#define POS1Z           11
#define POSX1           12
#define POSY1           13
#define POSZ1           14
#define POS11           15

// tensor
#define TSRX        0
#define TSRY        1
#define TSRZ        2
#define TSRXX       0
#define TSRXY       1
#define TSRXZ       2
#define TSRYX       3
#define TSRYY       4
#define TSRYZ       5
#define TSRZX       6
#define TSRZY       7
#define TSRZZ       8

// ng[*]
#define IINC            0
#define JINC            1
#define KINC            2
#define LINC            3
#define GSHIFT          4
#define POS_E1          5
#define POS_E2          6
#define TENSOR          7

#include "config.h"

// some other boundaries
#ifdef HAVE_QUADMATH_H
#define MXRYSROOTS      32 // > ANG_MAX*2+1 for 4c2e
#else
#define MXRYSROOTS      16 // > ANG_MAX*2+1 for 4c2e
#endif
#define ANG_MAX         12 // l = 0..12
#define LMAX1           16 // > ANG_MAX
#define CART_MAX        128 // > (ANG_MAX*(ANG_MAX+1)/2)
#define SHLS_MAX        0x7fffffff
#define NPRIM_MAX       64
#define NCTR_MAX        64

#if !defined STATIC_DOUBLE_SIZE
#define STATIC_DOUBLE_SIZE    500000
#endif

#define EXPCUTOFF       60
#ifndef MIN_EXPCUTOFF
// ~ 1e-15
#define MIN_EXPCUTOFF   40
#endif

#define OF_CMPLX        2

#define PI              3.1415926535897932384626433832795028
#ifndef M_PI
#define M_PI            PI
#endif
#define SQRTPI          1.7724538509055160272981674833411451

#define POINT_NUC       1
#define GAUSSIAN_NUC    2
#define FRAC_CHARGE_NUC 3
