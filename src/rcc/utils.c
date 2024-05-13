/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
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

#include <complex.h>
#include <ctype.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <utils.h>


/**
 * time[sec] -> number of days, hours, minutes, seconds, milliseconds
 * t - total wall time [sec]
 * Returns:
 * d - days, h - hours, m - minutes, s - seconds, ms - milliseconds
 */
void decompose_time(double t, int *d, int *h, int *m, int *s, int *ms)
{
    // d, h, m in seconds
    const double d2s = 24 * 3600.0, h2s = 3600.0, m2s = 60.0;

    *d = floor(t / d2s);
    t -= *d * d2s;
    *h = floor(t / h2s);
    t -= *h * h2s;
    *m = floor(t / m2s);
    t -= *m * m2s;
    // now: units of remainder are [seconds]
    *s = floor(t);
    t -= *s;
    *ms = floor(t * 1000);
}


/**
 * Signum math function.
 */
int sgn(double x)
{
    if (x > 0) { return 1; }
    if (x < 0) { return -1; }
    return 0;
}


/**
 * checks if string 's' contains only whitespaces
 */
int is_whitespace_string(const char *s)
{
    while (*s != '\0') {
        if (!isspace((unsigned char) *s)) {
            return 0;
        }
        s++;
    }
    return 1;
}


/**
 * imax
 *
 * Finds maximum value in an array of integer numbers
 */
int imax(size_t count, int *values)
{
    int themax = values[0];

    for (size_t i = 1; i < count; ++i) {
        themax = values[i] > themax ? values[i] : themax;
    }

    return themax;
}


int max(int x, int y)
{
    return (x >= y) ? x : y;
}


int min(int x, int y)
{
    return (x <= y) ? x : y;
}


/**
 * Finds maximum value in an array of numbers of type 'size_t'
 */
size_t size_t_max(size_t count, size_t *values)
{
    size_t themax = values[0];

    for (size_t i = 1; i < count; ++i) {
        themax = values[i] > themax ? values[i] : themax;
    }

    return themax;
}


/**
 * compares two integer arrays.
 * returns 0 if equal, > 0 if a > b; else < 0 (the same as the strcmp).
 * a[n], b[n]
 */
int intcmp(size_t n, int *a, int *b)
{
    return memcmp(a, b, sizeof(int) * n);
}


/**
 * power for integer numbers
 */
int64_t int_pow(int base, int exp)
{
    int64_t result = 1;
    while (exp) {
        if (exp & 1) {
            result *= base;
        }
        exp /= 2;
        base *= base;
    }
    return result;
}


/**
 * returns maximum of three integers
 */
int imax3(int a, int b, int c)
{
    int m = a;
    (m < b) && (m = b); //these are not conditional statements.
    (m < c) && (m = c); //these are just boolean expressions.
    return m;
}


/**
 * Returns unique identifier of the program run
 */
size_t run_id()
{
    static size_t id = 0;

    if (id == 0) {
        id = (size_t) time(NULL);
    }

    return id;
}


/**
 * Search integer in the array of integers.
 * Using the simplest O(N) algorithm.
 * Returns index of the element 'x' or -1 if not found.
 */
int int_search(size_t count, int *arr, int x)
{
    for (size_t i = 0; i < count; i++) {
        if (arr[i] == x) {
            return i;
        }
    }
    return -1;
}


/**
 * Checks if given number is in range: x in [a,b]
 */
int in_range(double x, double a, double b)
{
    if (a <= x && x <= b) {
        return 1;
    }
    else {
        return 0;
    }
}


/**
 * Copy integer arrays.
 */
void intcpy(int *dst, int *src, size_t n)
{
    memcpy(dst, src, sizeof(int) * n);
}


/**
 * replaces all occurences of symbol 'x' with symbol 'y'
 */
void str_replace(char *s, char x, char y)
{
    char *p = s;
    while (*p++) {
        if (*p == x) {
            *p = y;
        }
    }
}


int cmpfunc_int(const void *a, const void *b)
{
    return (*(int *) a - *(int *) b);
}


int cmpfunc_double(const void *a, const void *b)
{
    if (*(double *) a > *(double *) b) {
        return 1;
    }
    else if (*(double *) a < *(double *) b) {
        return -1;
    }
    else {
        return 0;
    }
}


double sum_doubles(size_t n, double *v)
{
    double sum = 0.0;

    for (size_t i = 0; i < n; i++) {
        sum += v[i];
    }

    return sum;
}


void get_real_parts(size_t n, double complex *src, double *dst)
{
    for (size_t i = 0; i < n; i++) {
        dst[i] = creal(src[i]);
    }
}


void get_imag_parts(size_t n, double complex *src, double *dst)
{
    for (size_t i = 0; i < n; i++) {
        dst[i] = cimag(src[i]);
    }
}


/**
 * sort array in acsending order
 * TODO: optimize for n=2,3 (without qsort)
 */
void isort(int n, int *src, int *dst)
{
    intcpy(dst, src, n);

    qsort(dst, n, sizeof(int), cmpfunc_int);
}


/**
 * counts occurrence of the given character in the string
 */
size_t strcount(char *s, char c)
{
    size_t count = 0;

    while (*s) {
        if (*s++ == c) {
            count++;
        }
    }

    return count;
}


/**
 * prints line of hyphens, only for decoration
 */
void print_hyphens(int offs, int len)
{
    for (int i = 0; i < offs; i++) {
        printf(" ");
    }
    for (int i = 0; i < len; i++) {
        printf("-");
    }
    printf("\n");
}


/**
 * prints string with current date and time
 */
void print_asctime()
{
    time_t tm = time(0);
    printf(" %s", asctime(localtime(&tm)));
}


/**
 * gets real parts of the array of complex numbers.
 * arrays should be pre-allocated.
 *
 * @param n
 * @param array_complex
 * @param array_real
 */
void array_get_real_part(size_t n, double complex *array_complex, double *array_real)
{
    for (size_t i = 0; i < n; i++) {
        array_real[i] = creal(array_complex[i]);
    }
}


/**
 * An alternative implementation of asprintf
 */
int cc_asprintf(char **str, char *fmt, ...)
{
    va_list argp;

    // get length of the string
    va_start(argp, fmt);
    char one_char[1];
    int len = vsnprintf(one_char, 1, fmt, argp);
    if (len < 1) {
        fprintf(stderr, "An encoding error occurred. Setting the input pointer to NULL.\n");
        *str = NULL;
        return len;
    }
    va_end(argp);

    // allocate memory for the string
    *str = malloc(len + 1);
    if (!str) {
        fprintf(stderr, "Couldn't allocate %i bytes.\n", len + 1);
        return -1;
    }

    // save data to the string
    va_start(argp, fmt);
    vsnprintf(*str, len + 1, fmt, argp);
    va_end(argp);

    return len;
}


/*
 * parses string of integer numbers, converts it into array of integers.
 * output array of integers must be pre-allocated.
 * returns 0 on success, -1 otherwise
 *
 * example: "1234" -> [1, 2, 3, 4]
 */
int str_to_int_array(char *str, int *array)
{
    size_t len = strlen(str);

    for (int i = 0; i < len; i++) {
        char c = str[i];
        if (isdigit(c)) {
            array[i] = c - '0';
        }
        else {
            return -1;
        }
    }

    return 0;
}


/*
 * converts array of positive integer numbers in range [0,9] into a string.
 * memory for the string must be pre-allocated.
 * returns 0 on success, -1 otherwise
 *
 * example: [1, 2, 3, 4] -> "1234\0"
 * (string length = 5 symbols, including '\0')
 */
int int_array_to_str(int n, const int *array, char *str)
{
    for (int i = 0; i < n; i++) {
        if (array[i] < 0 || array[i] > 9) {
            return -1;
        }
        str[i] = (char) ('0' + array[i]);
    }
    str[n] = '\0';

    return 0;
}
