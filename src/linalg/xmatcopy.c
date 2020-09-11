/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2020 The EXP-T developers.
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

#include "linalg.h"

#include <assert.h>
#include <ctype.h>
#include <string.h>

void dimatcopy(char op, size_t rows, size_t cols, double *A);

void zimatcopy(char op, size_t rows, size_t cols, double complex *A);

void domatcopy(char op, size_t rows, size_t cols, double *A, double *B);

void zomatcopy(char op, size_t rows, size_t cols, double complex *A, double complex *B);

void transpose_inplace_double(double *m, int w, int h);

void transpose_inplace_double_complex(double complex *m, int w, int h);


/**
 * Transposition/copying/conjugation of matrices:
 *   B := op(A)
 * operation 'op' can be:
 *   'N' or 'n'   normal matrix copy
 *   'T' or 't'   transposition
 *   'C' or 'c'   conjugate transposition
 *   'R' or 'r'   conjugation
 * Both in-place and out-place operations are allowed.
 * This routine resembles the '?[i/o]matcopy' routines from MKL and OpenBLAS:
 *   https://software.intel.com/en-us/mkl-developer-reference-c-mkl-imatcopy
 *   https://software.intel.com/en-us/mkl-developer-reference-c-mkl-omatcopy
 *   https://github.com/xianyi/OpenBLAS/wiki/OpenBLAS-Extensions
 *
 * @param type  CC_DOUBLE or CC_DOUBLE_COMPLEX
 * @param op    operation type
 * @param rows  number of rows
 * @param cols  number of columns
 * @param A     input matrix
 * @param B     output matrix (can coincide with the input one)
 *
 * @todo test all this stuff
 */
void xmatcopy(data_type_t type, char op, size_t rows, size_t cols, void *A, void *B)
{
    assert(type == CC_DOUBLE || type == CC_DOUBLE_COMPLEX);

    if (type == CC_DOUBLE) {
        if (A == B) {
            dimatcopy(op, rows, cols, (double *) A);
        }
        else {
            domatcopy(op, rows, cols, (double *) A, (double *) B);
        }
    }
    else{
        if (A == B) {
            zimatcopy(op, rows, cols, (double complex *) A);
        }
        else {
            zomatcopy(op, rows, cols, (double complex *) A, (double complex *) B);
        }
    }
}


/**
 * in-place real version
 */
void dimatcopy(char op, size_t rows, size_t cols, double *A)
{
    op = tolower(op);
    assert(op == 'n' || op == 't');

    if (op == 'n') {
        return;
    }
    else{
        transpose_inplace_double(A, rows, cols);
    }
}


/**
 * in-place complex version
 */
void zimatcopy(char op, size_t rows, size_t cols, double complex *A)
{
    op = tolower(op);
    assert(op == 'n' || op == 't' || op == 'c' || op == 'r');

    // copy
    if (op == 'n') {
        return;
    }
        // conjugation only
    else if (op == 'r') {
        for (size_t i = 0; i < rows * cols; i++){
            A[i] = conj(A[i]);
        }
    }
        // transpose
    else if (op == 't') {
        transpose_inplace_double_complex(A, rows, cols);
    }
        // conjugate transpose
    else{
        transpose_inplace_double_complex(A, rows, cols);
        for (size_t i = 0; i < rows * cols; i++){
            A[i] = conj(A[i]);
        }
    }
}


/**
 * out-place real version
 */
void domatcopy(char op, size_t rows, size_t cols, double *A, double *B)
{
    op = tolower(op);
    assert(A != B);
    assert(op == 'n' || op == 't');

    if (op == 'n') {
        memcpy(B, A, sizeof(double) * rows * cols);
    }
    else{
        for (size_t i = 0; i < rows; i++){
            for (size_t j = 0; j < cols; j++){
                // B_ji = A_ij
                B[j * rows + i] = A[i * cols + j];
            }
        }
    }
}


/**
 * out-place complex version
 */
void zomatcopy(char op, size_t rows, size_t cols, double complex *A, double complex *B)
{
    op = tolower(op);
    assert(A != B);
    assert(op == 'n' || op == 't' || op == 'c' || op == 'r');

    // copy
    if (op == 'n') {
        memcpy(B, A, sizeof(double complex) * rows * cols);
    }
        // conjugation only
    else if (op == 'r') {
        for (size_t i = 0; i < rows * cols; i++){
            B[i] = conj(A[i]);
        }
    }
        // transpose
    else if (op == 't') {
        for (size_t i = 0; i < rows; i++){
            for (size_t j = 0; j < cols; j++){
                // B_ji = A_ij
                B[j * rows + i] = A[i * cols + j];
            }
        }
    }
        // conjugate transpose
    else{
        for (size_t i = 0; i < rows; i++){
            for (size_t j = 0; j < cols; j++){
                // B_ji = (A_ij)^{*}
                B[j * rows + i] = conj(A[i * cols + j]);
            }
        }
    }
}


/**
 * in-place matrix transpositions (for real and complex cases).
 * From RosettaCode:
 * https://rosettacode.org/wiki/Matrix_transposition#C
 */
void transpose_inplace_double(double *m, int w, int h)
{
    int start, next, i;
    double tmp;

    for (start = 0; start <= w * h - 1; start++){
        next = start;
        i = 0;
        do {
            i++;
            next = (next % h) * w + next / h;
        } while (next > start);
        if (next < start || i == 1) continue;

        tmp = m[next = start];
        do {
            i = (next % h) * w + next / h;
            m[next] = (i == start) ? tmp : m[i];
            next = i;
        } while (next > start);
    }
}

void transpose_inplace_double_complex(double complex *m, int w, int h)
{
    int start, next, i;
    double complex tmp;

    for (start = 0; start <= w * h - 1; start++){
        next = start;
        i = 0;
        do {
            i++;
            next = (next % h) * w + next / h;
        } while (next > start);
        if (next < start || i == 1) continue;

        tmp = m[next = start];
        do {
            i = (next % h) * w + next / h;
            m[next] = (i == start) ? tmp : m[i];
            next = i;
        } while (next > start);
    }
}
