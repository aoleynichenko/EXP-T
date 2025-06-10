/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2025 The EXP-T developers.
 *
 *  This file is part of EXP-T.
 *
 *  EXP-T is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EXP-T is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  E-mail:        exp-t-program@googlegroups.com
 *  Google Groups: https://groups.google.com/d/forum/exp-t-program
 */

#ifndef CC_TT_CCSD_H_INCLUDED
#define CC_TT_CCSD_H_INCLUDED

#include <stdint.h>
#include "complex.h"

extern void* s_fm_create(const float *fm, const uint32_t shape, const uint32_t occ);
extern void* d_fm_create(const double *fm, const uint64_t shape, const uint64_t occ);
extern void* c_fm_create(const float complex *fm, const uint32_t shape, const uint32_t occ);
extern void* z_fm_create(const double complex *fm, const uint64_t shape, const uint64_t occ);


extern void* s_eri_create(const float *eri, const uint32_t shape, const uint32_t occ, const float tol, const float tol_eri);
extern void* d_eri_create(const double *eri, const uint64_t shape, const uint64_t occ, const double tol, const double tol_eri);
extern void* c_eri_create(const float complex *eri, const uint32_t shape, const uint32_t occ, const float tol, const float tol_eri);
extern void* z_eri_create(const double complex *eri, const uint64_t shape, const uint64_t occ, const double tol, const double tol_eri);


/*
extern void* s_ccsd_create(const void *eri, const void *fm, const float tol, _Bool use_tt);
extern void* d_ccsd_create(const void *eri, const void *fm, const double tol, _Bool use_tt);
extern void* c_ccsd_create(const void *eri, const void *fm, const float tol, _Bool use_tt);
extern void* z_ccsd_create(const void *eri, const void *fm, const double tol, _Bool use_tt);
*/

extern void* s_ccsd_update(void *ccsd, const void *eri, const void *fm, const float tol);
extern void* d_ccsd_update(void *ccsd, const void *eri, const void *fm, const double tol);
extern void* c_ccsd_update(void *ccsd, const void *eri, const void *fm, const float tol);
extern void* z_ccsd_update(void *ccsd, const void *eri, const void *fm, const double tol);


extern void* s_ccsd_create(const float *eri, const float *fm, const uint32_t shape, const uint32_t occ, const float tol, _Bool use_tt);
extern void* d_ccsd_create(const double *eri, const double *fm, const uint32_t shape, const uint32_t occ, const double tol, _Bool use_tt);
extern void* c_ccsd_create(const float complex *eri, const float complex *fm, const uint32_t shape, const uint32_t occ, const float tol, _Bool use_tt);
extern void* z_ccsd_create(const double complex *eri, const double complex *fm, const uint32_t shape, const uint32_t occ, const double tol, _Bool use_tt);

float s_ccsd_energy(void *ccsd, const void *eri, const void *fm);
double d_ccsd_energy(void *ccsd, const void *eri, const void *fm);
float c_ccsd_energy(void *ccsd, const void *eri, const void *fm);
double z_ccsd_energy(void *ccsd, const void *eri, const void *fm);

float s_mp2_energy(void *ccsd, const void *eri, const void *fm);
double d_mp2_energy(void *ccsd, const void *eri, const void *fm);
float c_mp2_energy(void *ccsd, const void *eri, const void *fm);
double z_mp2_energy(void *ccsd, const void *eri, const void *fm);

float s_max_t2(void *ccsd);
double d_max_t2(void *ccsd);
float c_max_t2(void *ccsd);
double z_max_t2(void *ccsd);

void print_all(const void *ccsd, const void *eri, const void *fm);



void tt_ccsd();

#endif // CC_TT_CCSD_H_INCLUDED
