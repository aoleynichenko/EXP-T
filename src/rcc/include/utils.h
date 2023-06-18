/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2023 The EXP-T developers.
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

/*
 * Small utility functions.
 */

#ifndef CC_UTILS_H_INCLUDED
#define CC_UTILS_H_INCLUDED

#include <stdint.h>
#include <complex.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void decompose_time(double t, int *d, int *h, int *m, int *s, int *ms);

enum {
    CC_TYPE_FLOAT = sizeof(float),
    CC_TYPE_DOUBLE = sizeof(double),
    CC_TYPE_FLOAT_COMPLEX = sizeof(float complex),
    CC_TYPE_COMPLEX = sizeof(double complex)
};

void *cc_zero_matrix(size_t n, size_t m, int type);

double complex *cc_zero_matrix_complex(size_t n, size_t m);

double *cc_zero_matrix_real(size_t n, size_t m);

int sgn(double x);

int is_whitespace_string(const char *s);

int cmpfunc_double(const void *a, const void *b);

int intcmp(size_t n, int *a, int *b);

int max(int a, int b);

int min(int a, int b);

int imax3(int a, int b, int c);

int imax(size_t count, int *values);

int64_t int_pow(int base, int exp);

size_t size_t_max(size_t count, size_t *values);

size_t run_id();

int int_search(size_t count, int *arr, int x);

void reverse_perm(int n, int *a, int *ainv);

int in_range(double x, double a, double b);

void intcpy(int *a, int *b, size_t n);

void str_replace(char *s, char x, char y);

void int2str(int *dig, char *str, size_t n);

void isort(int n, int *src, int *dst);

double sum_doubles(size_t n, double *v);

void get_real_parts(size_t n, double complex *src, double *dst);

void get_imag_parts(size_t n, double complex *src, double *dst);

size_t strcount(char *s, char c);

void print_hyphens(int offs, int len);

void print_asctime();

void array_get_real_part(size_t n, double complex *array_complex, double *array_real);

#endif /* CC_UTILS_H_INCLUDED */
