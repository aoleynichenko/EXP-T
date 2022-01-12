/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2022 The EXP-T developers.
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

#include "cubic_spline.h"

/*
 * Cubic spline interpolation.
 *
 * The algorithm is adapted from:
 * W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery,
 * Numerical Recipes in C, Second Edition.
 * Cambridge University Press, 1992
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "errquit.h"

void sort_pairs_ascending(int n, double *x, double *y);

void spline_find_piece(cubic_spline_t *spline, double x, int *low, int *high);


/**
 * Constructs and returns the cubic spline object.
 * Natural spline is assumed (boundary conditions: y''[0] = y''[n-1] = 0.0)
 */
cubic_spline_t *construct_cubic_spline(int n, double *x, double *y)
{
    cubic_spline_t *spline = (cubic_spline_t *) calloc(1, sizeof(cubic_spline_t));

    // new spline structure
    spline->n = n;
    spline->x = (double *) calloc(n, sizeof(double));
    spline->y = (double *) calloc(n, sizeof(double));
    spline->y2 = (double *) calloc(n, sizeof(double));

    memcpy(spline->x, x, sizeof(double) * n);
    memcpy(spline->y, y, sizeof(double) * n);

    // node points of the spline must be sorted in the order: x[0] < x[1] < ... < x[n-1]
    sort_pairs_ascending(n, spline->x, spline->y);

    // temporary array
    double *u = (double *) calloc(n, sizeof(double));

    // we will construct natural spline
    spline->y2[0] = 0.0;
    u[0] = 0.0;

    // decomposition loop of the tridiagonal algorithm
    // 'y2' and 'u' are used for temporary storage of the decomposed factors

    for (int i = 1; i < n - 1; i++) {
        double sig = (spline->x[i] - spline->x[i - 1]) / (spline->x[i + 1] - spline->x[i - 1]);
        double p = sig * spline->y2[i - 1] + 2.0;
        spline->y2[i] = (sig - 1.0) / p;
        u[i] = (spline->y[i + 1] - spline->y[i]) / (spline->x[i + 1] - spline->x[i]) -
               (spline->y[i] - spline->y[i - 1]) / (spline->x[i] - spline->x[i - 1]);
        u[i] = (6.0 * u[i] / (spline->x[i + 1] - spline->x[i - 1]) - sig * u[i - 1]) / p;
    }

    // the "natural" upper boundary condition
    double qn = 0.0;
    double un = 0.0;

    // backsubstitution loop of the tridiaonal algorithm

    spline->y2[n - 1] = (un - qn * u[n - 2]) / (qn * spline->y2[n - 2] + 1.0);
    for (int k = n - 2; k >= 0; k--) {
        spline->y2[k] = spline->y2[k] * spline->y2[k + 1] + u[k];
    }

    free(u);

    return spline;
}


/**
 * Evaluates the value of the spline at the given point 'x'
 */
double evaluate_spline(cubic_spline_t *spline, double x)
{
    int i1, i2;
    spline_find_piece(spline, x, &i1, &i2);

    double h = spline->x[i2] - spline->x[i1];
    double a = (spline->x[i2] - x) / h;
    double b = (x - spline->x[i1]) / h;

    // evaluate cubic spline polynomial
    double y = a * spline->y[i1] + b * spline->y[i2] +
               ((a * a * a - a) * spline->y2[i1] + (b * b * b - b) * spline->y2[i2]) * h * h / 6.0;

    return y;
}


/**
 * Deallocation of the cubic spline object
 */
void delete_spline(cubic_spline_t *spline)
{
    if (spline != NULL) {
        free(spline->x);
        free(spline->y);
        free(spline->y2);
        free(spline);
    }
}


/**
 * Calculates the first derivative at the point 'x'
 */
double spline_first_derivative(cubic_spline_t *spline, double x)
{
    double x_min = spline->x[0];
    double x_max = spline->x[spline->n - 1];
    if (x < x_min || x > x_max) {
        errquit("spline_first_derivative(): argument of the spline (%g) is out of range [%g,%g]", x, x_min, x_max);
    }

    int i1, i2;
    spline_find_piece(spline, x, &i1, &i2);

    double h = spline->x[i2] - spline->x[i1];
    double a = (spline->x[i2] - x) / h;
    double b = (x - spline->x[i1]) / h;

    return (spline->y[i2] - spline->y[i1]) / h +
           (h / 6.0) * (spline->y2[i1] * (1 - 3.0*a*a) + spline->y2[i2] * (3.0*b*b - 1));
}


/**
 * Calculates the second derivative at the point 'x'
 */
double spline_second_derivative(cubic_spline_t *spline, double x)
{
    double x_min = spline->x[0];
    double x_max = spline->x[spline->n - 1];
    if (x < x_min || x > x_max) {
        errquit("spline_second_derivative(): argument of the spline (%g) is out of range [%g,%g]", x, x_min, x_max);
    }

    int i1, i2;
    spline_find_piece(spline, x, &i1, &i2);

    double h = spline->x[i2] - spline->x[i1];
    double a = (spline->x[i2] - x) / h;
    double b = (x - spline->x[i1]) / h;

    return spline->y2[i1] * a + spline->y2[i2] * b;
}


/**
 * Calculates the third derivative at the point 'x'
 */
double spline_third_derivative(cubic_spline_t *spline, double x)
{
    double x_min = spline->x[0];
    double x_max = spline->x[spline->n - 1];
    if (x < x_min || x > x_max) {
        errquit("spline_second_derivative(): argument of the spline (%g) is out of range [%g,%g]", x, x_min, x_max);
    }

    int i1, i2;
    spline_find_piece(spline, x, &i1, &i2);

    double h = spline->x[i2] - spline->x[i1];

    return (spline->y2[i2] - spline->y2[i1]) / h;
}


/**
 * find right place in the table by means of bisection
 */
void spline_find_piece(cubic_spline_t *spline, double x, int *index_low, int *index_high)
{
    int low = 0;
    int high = spline->n - 1;
    while (high - low > 1) {
        int k = (high + low) >> 1;
        if (spline->x[k] > x) {
            high = k;
        } else {
            low = k;
        }
    }

    *index_low = low;
    *index_high = high;
}


/**
 * Returns the minimum point (xmin,ymin) for the given spline-interpolated function
 */
void find_spline_minimum(cubic_spline_t *spline, double *xmin, double *ymin, int print_level)
{
    double const x_thresh = 1e-5;
    double const y_thresh = 1e-8;
    double const max_grad = 1e-4;
    int maxiter = 50;
    double alpha = 0.5;

    if (print_level > 0) {
        printf(" > minimization of the spline function\n");
        printf("\n");
        printf("   optimization parameters:\n");
        printf("   x threshold      = %g\n", x_thresh);
        printf("   y threshold      = %g\n", y_thresh);
        printf("   max gradient     = %g\n", max_grad);
        printf("   max iterations   = %d\n", maxiter);
        printf("   speed of descent = %g\n", alpha);
    }

    // find initial guess for the iterative minimization
    *xmin = spline->x[0];
    *ymin = spline->y[0];
    for (int i = 1; i < spline->n; i++) {
        if (spline->y[i] < *ymin) {
            *xmin = spline->x[i];
            *ymin = spline->y[i];
        }
    }

    if (print_level > 0) {
        printf("\n");
        printf("   initial guess:\n");
        printf("   x0 = %20.12f\n", *xmin);
        printf("   y0 = %20.12f\n", *ymin);
        printf("\n");

        printf("   it.           xmin                ymin                 grad\n");
    }

    // the gradient descent algorithm is used to locate the minimum
    int it;
    for (it = 0; it < maxiter; it++) {

        double grad = spline_first_derivative(spline, *xmin);
        double xnew = *xmin - alpha * grad;
        double ynew = evaluate_spline(spline, xnew);

        if (print_level > 0) {
            printf("   %3d%20.12f%20.12f%20.12f\n", it, xnew, ynew, grad);
        }

        if (fabs(*xmin - xnew) < x_thresh &&
            fabs(*ymin - ynew) < y_thresh &&
            fabs(grad) < max_grad) {
            *xmin = xnew;
            *ymin = ynew;
            break;
        }

        *xmin = xnew;
        *ymin = ynew;
    }

    if (print_level > 0) {
        if (it == maxiter) {
            printf("   not converged\n");
        } else {
            printf("   converged\n");
        }
    }

    if (print_level > 0) {
        printf("\n");
        printf("   xmin = %20.12f\n", *xmin);
        printf("   ymin = %20.12f\n", *ymin);
        printf("\n");
    }
}


/**
 * Sorts pairs of points (x[i],y[i]) by the x[i] values (in an ascending order)
 */
void sort_pairs_ascending(int n, double *x, double *y)
{
    double xtmp, ytmp;
    size_t i, j;

    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {

            if (x[i] <= x[j]) {
                continue;
            }

            xtmp = x[i];
            x[i] = x[j];
            x[j] = xtmp;

            ytmp = y[i];
            y[i] = y[j];
            y[j] = ytmp;
        }
    }
}















