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
 * Diagram contractions.
 * Are based on low-level operations (see diagram.h)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cuda_code.h"
#include "engine.h"
#include "platform.h"
#include "datamodel.h"
#include "error.h"
#include "linalg.h"
#include "options.h"
#include "symmetry.h"
#include "timer.h"
#include "utils.h"

#if defined BLAS_MKL
  #include "mkl.h"
#endif

#include "omp.h"

/*
 * type of multiplication:
 * D -- diagram is stored on Disk
 * M -- diagram is stored in Memory
 */
enum {
    MULT_M_MM,
    MULT_M_DM,
    MULT_M_MD,
    MULT_M_DD,
    MULT_D_MM,
    MULT_D_MD,
    MULT_D_DM,
    MULT_D_DD
};

static void mult_check_creation_annihilation(diagram_t *dg1, diagram_t *dg2, int ncontr);

static void mult_check_quasiparticles(diagram_t *dg1, diagram_t *dg2, int ncontr);

static void mult_check_valence_t3space(diagram_t *dg1, diagram_t *dg2, int ncontr);

static diagram_t *mult_product_template(diagram_t *dg1, diagram_t *dg2, int ncontr, int perm_unique);

static int mult_type(diagram_t *op1, diagram_t *op2, diagram_t *prod);

void mult_algorithm_m_mm_openmp_external(diagram_t *op1, diagram_t *op2, diagram_t *tgt, int ncontr);

void mult_algorithm_m_mm_openmp_internal(diagram_t *op1, diagram_t *op2, diagram_t *tgt, int ncontr);

void mult_algorithm_m_dm(diagram_t *op1, diagram_t *op2, diagram_t *tgt, int ncontr);

void target_order(int *ord1, int rk1, int *ord2, int rk2, int *ord3, int rk3);

void supmat_dims(block_t *b1, block_t *b2, int ncontr, int *m, int *n, int *k);

void mulblocks_lapack(double complex *A, int ma, int na,
                      double complex *B, int mb, int nb,
                      double beta, double complex *C);

void mulblocks(block_t *op1, block_t *op2, block_t *prod, int ncontr, int nthreads);


/**
 * Performs contraction of two diagrams:
 * target = fact1 * name2
 * 'target', 'name1', 'name2' are allowed to coincide.
 *
 * Arguments:
 *   name1    symbolic name of the diagram -- first term
 *   name2    symbolic name of the diagram -- second term
 *   target   symbolic name of the resulting diagram
 *   mmm      three numbers:
 *            how many lines we should save in the first multiplier (m1), in the
 *            second (m2) and how many lines will be contracted (m3).
 *            mmm = [m1][m2][m3].
 *            Example: m1 = 1, m2 = 3, m3 = 1  =>  mmm = 131
 */
void mult(char *name1, char *name2, char *target, int ncontr)
{
    diagram_t *dg1, *dg2, *dg_prod;

    dg1 = diagram_stack_find(name1);
    if (dg1 == NULL) {
        errquit("mult(): diagram '%s' not found", name1);
    }

    dg2 = diagram_stack_find(name2);
    if (dg2 == NULL) {
        errquit("mult(): diagram '%s' not found", name2);
    }

    dg_prod = diagram_mult(dg1, dg2, ncontr, (target[0] == '$') ? 1 : 0);
    strcpy(dg_prod->name, target);

    // save new diagram to stack
    // (replace the old diagram named 'target' if needed)
    if (diagram_stack_find(target) != NULL) {
        diagram_stack_replace(target, dg_prod);
    }
    else {
        diagram_stack_push(dg_prod);
    }
}


/**
 * Evaluates contraction of two diagrams: name2 * name3 -> namet.
 *
 * Arguments:
 *   dg1, dg2 -- source diagrams to be contracted
 *   mmm   -- how many lines we should save in the first multiplier (m1), in the
 *            second (m2) and how many lines will be contracted (m3).
 *            mmm = [m1][m2][m3].
 *            Example: m1 = 1, m2 = 3, m3 = 1  =>  mmm = 131
 * Returns:
 *   target diagram
 */
diagram_t *diagram_mult(diagram_t *dg1, diagram_t *dg2, int ncontr, int perm_unique)
{
    timer_new_entry("mult", "Diagram contraction (mult)");
    timer_start("mult");
    timer_new_entry("mult_mmm", "mult M <- M x M");
    timer_new_entry("mult_mdm", "mult M <- D x M");

    mult_check_quasiparticles(dg1, dg2, ncontr);
    mult_check_valence_t3space(dg1, dg2, ncontr);
    mult_check_creation_annihilation(dg1, dg2, ncontr);

    diagram_t *tgt = mult_product_template(dg1, dg2, ncontr, perm_unique);

    // choose the appropriate algorithm with the least number of I/O operations
    int algo = mult_type(dg1, dg2, tgt);
    switch (algo) {
        case MULT_M_MM:
        case MULT_D_MM:
        case MULT_D_DM:
            timer_start("mult_mmm");
            if (cc_opts->openmp_algorithm == CC_OPENMP_ALGORITHM_EXTERNAL) {
                mult_algorithm_m_mm_openmp_external(dg1, dg2, tgt, ncontr);
            }
            else {
                mult_algorithm_m_mm_openmp_internal(dg1, dg2, tgt, ncontr);
            }
            timer_stop("mult_mmm");
            break;
        case MULT_M_DM:
            timer_start("mult_mdm");
            mult_algorithm_m_dm(dg1, dg2, tgt, ncontr);
            timer_stop("mult_mdm");
            break;
        default:
            diagram_summary(dg1);
            diagram_summary(dg2);
            printf(" ncontr = %d\n", ncontr);
            errquit("unimplemented contraction algorithm number %d\n", algo);
            break;
    }

    timer_stop("mult");

    return tgt;
}


void restore_unique_blocks(diagram_t *dg)
{
    for (size_t ib = 0; ib < dg->n_blocks; ib++) {
        block_t *block = dg->blocks[ib];

        if (block->is_unique == 0) {
            restore_block(dg, block);
        }
    }
}


void destroy_unique_blocks(diagram_t *dg)
{
    for (size_t ib = 0; ib < dg->n_blocks; ib++) {
        block_t *block = dg->blocks[ib];

        if (block->is_unique == 0) {
            destroy_block(block);
        }
    }
}


/*
 * sequence of loops:
 * for block C in product:
 *      for block A in operand-1:
 *          for block B in operand-2:
 *              C += A * B
 */
void mult_algorithm_m_mm_openmp_external(diagram_t *op1, diagram_t *op2, diagram_t *tgt, int ncontr)
{
    int rk1 = op1->rank;
    int rk2 = op2->rank;

    restore_unique_blocks(op1);
    restore_unique_blocks(op2);

    #pragma omp parallel for schedule(dynamic) num_threads(cc_opts->nthreads)
    for (size_t ib3 = 0; ib3 < tgt->n_blocks; ib3++) {
        block_t *b3 = tgt->blocks[ib3];
        block_load(b3);
        if (b3->is_unique == 0) {
            continue;
        }

        for (size_t ib1 = 0; ib1 < op1->n_blocks; ib1++) {
            block_t *b1 = op1->blocks[ib1];
            if (intcmp(rk1 - ncontr, b1->spinor_blocks, b3->spinor_blocks) != 0) {
                continue;
            }
            block_load(b1);

            for (size_t ib2 = 0; ib2 < op2->n_blocks; ib2++) {
                block_t *b2 = op2->blocks[ib2];
                if (intcmp(rk2 - ncontr, b2->spinor_blocks, b3->spinor_blocks + rk1 - ncontr) != 0) {
                    continue;
                }
                if (intcmp(ncontr, b1->spinor_blocks + rk1 - ncontr, b2->spinor_blocks + rk2 - ncontr) != 0) {
                    continue;
                }

                block_load(b2);
                mulblocks(b1, b2, b3, ncontr, 1);
                block_unload(b2);
            }
            block_unload(b1);
        }
        block_store(b3);
    }

    destroy_unique_blocks(op1);
    destroy_unique_blocks(op2);
}


/*
 * sequence of loops:
 * for block C in product:
 *      for block A in operand-1:
 *          for block B in operand-2:
 *              C += A * B
 */
void mult_algorithm_m_mm_openmp_internal(diagram_t *op1, diagram_t *op2, diagram_t *tgt, int ncontr)
{
    int rk1 = op1->rank;
    int rk2 = op2->rank;

    for (size_t ib3 = 0; ib3 < tgt->n_blocks; ib3++) {
        block_t *b3 = tgt->blocks[ib3];
        block_load(b3);
        if (b3->is_unique == 0) {
            continue;
        }

        for (size_t ib1 = 0; ib1 < op1->n_blocks; ib1++) {
            block_t *b1 = op1->blocks[ib1];
            if (intcmp(rk1 - ncontr, b1->spinor_blocks, b3->spinor_blocks) != 0) {
                continue;
            }
            if (b1->is_unique == 0) {
                restore_block(op1, b1);
            }
            block_load(b1);

            for (size_t ib2 = 0; ib2 < op2->n_blocks; ib2++) {
                block_t *b2 = op2->blocks[ib2];
                if (intcmp(rk2 - ncontr, b2->spinor_blocks, b3->spinor_blocks + rk1 - ncontr) != 0) {
                    continue;
                }
                if (intcmp(ncontr, b1->spinor_blocks + rk1 - ncontr, b2->spinor_blocks + rk2 - ncontr) != 0) {
                    continue;
                }
                if (b2->is_unique == 0) {
                    restore_block(op2, b2);
                }

                block_load(b2);
                mulblocks(b1, b2, b3, ncontr, cc_opts->nthreads);
                block_unload(b2);

                if (b2->is_unique == 0) {
                    destroy_block(b2);
                }
            }
            block_unload(b1);
            if (b1->is_unique == 0) {
                destroy_block(b1);
            }
        }
        block_store(b3);
    }
}



/*
 * sequence of loops:
 * for block A in operand-1:
 *     for block B in operand-2:
 *         for block C in product:
 *             C += A * B
 */
void mult_algorithm_m_dm(diagram_t *op1, diagram_t *op2, diagram_t *tgt, int ncontr)
{
    int rk1 = op1->rank;
    int rk2 = op2->rank;

    for (size_t ib1 = 0; ib1 < op1->n_blocks; ib1++) {
        block_t *b1 = op1->blocks[ib1];
        block_load(b1);
        if (b1->is_unique == 0) {
            restore_block(op1, b1);
        }

        for (size_t ib2 = 0; ib2 < op2->n_blocks; ib2++) {
            block_t *b2 = op2->blocks[ib2];
            if (intcmp(ncontr, b1->spinor_blocks + rk1 - ncontr, b2->spinor_blocks + rk2 - ncontr) != 0) {
                continue;
            }
            block_load(b2);
            if (b2->is_unique == 0) {
                restore_block(op2, b2);
            }

            for (size_t ib3 = 0; ib3 < tgt->n_blocks; ib3++) {
                block_t *b3 = tgt->blocks[ib3];

                if (intcmp(rk1 - ncontr, b1->spinor_blocks, b3->spinor_blocks) != 0) {
                    continue;
                }
                if (intcmp(rk2 - ncontr, b2->spinor_blocks, b3->spinor_blocks + rk1 - ncontr) != 0) {
                    continue;
                }
                if (b3->is_unique == 0) {
                    continue;
                }

                block_load(b3);
                mulblocks(b1, b2, b3, ncontr, cc_opts->nthreads);
                block_unload(b3);
            }
            if (b2->is_unique == 0) {
                destroy_block(b2);
            }
            block_unload(b2);
        }
        if (b1->is_unique == 0) {
            destroy_block(b1);
        }
        block_store(b1);
    }
}


/**
 * general interface for contraction of two blocks.
 * encapsulates invocations of library-specific threading operations
 * (in case of OpenMP) or CUDA BLAS library (in case of CUDA)
 *
 * @note all locks must be preloaded!
 */
void mulblocks(block_t *op1, block_t *op2, block_t *prod, int ncontr, int nthreads)
{
    int M, N, K;

    supmat_dims(op1, op2, ncontr, &M, &N, &K);

    double complex *A = op1->buf;
    double complex *B = op2->buf;
    double complex *C = prod->buf;

    omp_set_num_threads(1);

    if (!cc_opts->cuda_enabled) {

        // enable internal threading
#if defined BLAS_MKL
        mkl_set_num_threads_local(nthreads);
#elif defined BLAS_OPENBLAS
        openblas_set_num_threads(nthreads);
#else
        // printf("else!\n");
#endif

        mulblocks_lapack(A, N, K, B, M, K, 1.0, C);

#if defined BLAS_MKL
        mkl_set_num_threads_local(1);
#elif defined BLAS_OPENBLAS
        openblas_set_num_threads(1);
#endif
    }
    else {
#ifdef CUDA_FOUND
        // parallel CUBLAS on GPU
            int err = mulblocks_cuda(carith, A, N, K, B, M, K, 1.0, C);
            if (err > 0) {
                errquit("Error in mult (CUDA error)");
            }
#endif /* CUDA_FOUND */
    }

    omp_set_num_threads(nthreads);
}


/**
 *   Wrapper for LAPACK's subroutine zgemm().
 *   It is much more convenient to use the CBLAS interface than invoke dgemm()
 *   directly from the Fortran code since all CC diagrams are stored in the
 *   row-major order.
 *   C = A * B + beta * C
 */
void mulblocks_lapack(double complex *A, int ma, int na,
                      double complex *B, int mb, int nb,
                      double beta, double complex *C)
{
    int m = ma;
    int n = mb;
    int k = na;
    double complex zalpha = 1.0 + 0.0 * I;
    double complex zbeta = beta + 0.0 * I;

    xgemm(WORKING_TYPE, "N", "T", m, n, k, &zalpha, A, k, B, k, &zbeta, C, n);
}


/**
 * determines dimensions of supermatrices:
 *   block b1 -> m x k
 *   block b2 -> n x k
 *
 * @param ncontr number of lines to be contracted
 * @note blocks are considered as "complementary" to each other
 * (dimensions to be contracted should have the same ranges of indices)
 */
void supmat_dims(block_t *b1, block_t *b2, int ncontr, int *m, int *n, int *k)
{
    int rk1 = b1->rank;
    int rk2 = b2->rank;
    int **ix1 = b1->indices;
    int **ix2 = b2->indices;
    int N = 1, M = 1, K = 1;

    for (int i = 0; i < rk1 - ncontr; i++) {
        N *= ix1[i][0];
    }
    for (int i = 0; i < rk2 - ncontr; i++) {
        M *= ix2[i][0];
    }
    for (int i = 0; i < ncontr; i++) {
        K *= ix1[rk1 - ncontr + i][0];
    }

    *m = M, *n = N, *k = K;
}


/**
 * compare (check) quasiparticle labels of two diagrams to be contracted.
 * Allowed pairs of quasiparticles are P-P and H-H
 *
 * @param dg1 pointer to the first diagram
 * @param dg2 pointer to the second diagram
 * @param ncontr number of lines to be contracted
 */
static void mult_check_quasiparticles(diagram_t *dg1, diagram_t *dg2, int ncontr)
{
    int rk1 = dg1->rank;
    int rk2 = dg2->rank;

    for (int i = 0; i < ncontr; i++) {
        if (dg1->qparts[rk1 - 1 - i] != dg2->qparts[rk2 - 1 - i]) {
            printf("dimensions to be contracted are not consistent:\n");
            diagram_summary(dg1);
            diagram_summary(dg2);
            errquit("error in mult: inconsistent diagrams");
        }
    }
}


/**
 * compare (check) valence labels of two diagrams to be contracted.
 * Allowed pairs of valence labels are G-G and A-A (G - general, A - active)
 *
 * @param dg1 pointer to the first diagram
 * @param dg2 pointer to the second diagram
 * @param ncontr number of lines to be contracted
 */
static void mult_check_valence_t3space(diagram_t *dg1, diagram_t *dg2, int ncontr)
{
    int rk1 = dg1->rank;
    int rk2 = dg2->rank;

    for (int i = 0; i < ncontr; i++) {
        if (dg1->valence[rk1 - 1 - i] != dg2->valence[rk2 - 1 - i]) {
            printf("dimensions to be contracted are not consistent (by valence):\n");
            diagram_summary(dg1);
            diagram_summary(dg2);
            errquit("error in mult: inconsistent diagrams");
        }
        if (dg1->t3space[rk1 - 1 - i] != dg2->t3space[rk2 - 1 - i]) {
            printf("dimensions to be contracted are not consistent (by t3space):\n");
            diagram_summary(dg1);
            diagram_summary(dg2);
            errquit("error in mult: inconsistent diagrams");
        }
    }
}


/**
 * check diagrams to be contracted for proper pairs creation-annihilation.
 * only pairs Creation-Annihilation are allowed.
 *
 * @param dg1 pointer to the first diagram
 * @param dg2 pointer to the second diagram
 * @param ncontr number of lines to be contracted
 */
static void mult_check_creation_annihilation(diagram_t *dg1, diagram_t *dg2, int ncontr)
{
    char creat_annih_1[CC_DIAGRAM_MAX_RANK];
    char creat_annih_2[CC_DIAGRAM_MAX_RANK];
    int rk1 = dg1->rank;
    int rk2 = dg2->rank;

    // A -- annihilation, C -- creation
    for (int i = 0; i < rk1; i++) {
        if (dg1->order[i] <= rk1 / 2) {
            creat_annih_1[i] = 'A';
        }
        else {
            creat_annih_1[i] = 'C';
        }
    }
    for (int i = 0; i < rk2; i++) {
        if (dg2->order[i] <= rk2 / 2) {
            creat_annih_2[i] = 'A';
        }
        else {
            creat_annih_2[i] = 'C';
        }
    }
    for (int i = 0; i < ncontr; i++) {
        if (creat_annih_1[rk1 - 1 - i] == creat_annih_2[rk2 - 1 - i]) {
            printf("dimensions to be contracted are not consistent:\n");
            diagram_summary(dg1);
            diagram_summary(dg2);
            printf("creation operators (C) must be paired with annihilation operators (A)\n");
            printf("diagram-1 (%s): ", dg1->name);
            for (int j = 0; j < rk1; j++) {
                printf("%c", creat_annih_1[j]);
            }
            printf("\n");
            printf("diagram-2 (%s): ", dg2->name);
            for (int j = 0; j < rk2; j++) {
                printf("%c", creat_annih_2[j]);
            }
            printf("\n");
            errquit("error in mult: inconsistent diagrams");
        }
    }
}


static diagram_t *mult_product_template(diagram_t *dg1, diagram_t *dg2, int ncontr, int perm_unique)
{
    char qparts[CC_DIAGRAM_MAX_RANK + 1];
    char valence[CC_DIAGRAM_MAX_RANK + 1];
    char t3space[CC_DIAGRAM_MAX_RANK + 1];
    char order[CC_DIAGRAM_MAX_RANK + 1];
    int tgt_order[CC_DIAGRAM_MAX_RANK];
    int rk1, rk2, rk3;

    rk1 = dg1->rank;
    rk2 = dg2->rank;
    rk3 = rk1 + rk2 - 2 * ncontr;

    // construct qparts, valence and order for the resulting diagram
    for (int i = 0; i < rk1 - ncontr; i++) {
        qparts[i] = dg1->qparts[i];
        valence[i] = dg1->valence[i] + '0';
        t3space[i] = dg1->t3space[i] + '0';
    }
    for (int i = 0; i < rk2 - ncontr; i++) {
        qparts[rk1 - ncontr + i] = dg2->qparts[i];
        valence[rk1 - ncontr + i] = dg2->valence[i] + '0';
        t3space[rk1 - ncontr + i] = dg2->t3space[i] + '0';
    }
    target_order(dg1->order, rk1, dg2->order, rk2, tgt_order, rk3);
    for (int i = 0; i < rk3; i++) {
        order[i] = tgt_order[i] + '0';
    }
    qparts[rk3] = '\0';
    valence[rk3] = '\0';
    t3space[rk3] = '\0';
    order[rk3] = '\0';

    // symmetry of the product operator
    int irrep_1 = dg1->symmetry;
    int irrep_2 = dg2->symmetry;
    int irrep_prod = mulrep2_abelian(irrep_1, irrep_2);

    // construct empty diagram with proper structure
    return diagram_new("mult-intermediate", qparts, valence, t3space, order, perm_unique, irrep_prod);
}


/**
 * Returns the type of the algorithm to be used for tensor contraction.
 * The order of three nested loops completely determines performance.
 * Recall: M -- Memory, D -- Disk
 *
 * @param op1 operand 1
 * @param op2 operand 2
 * @param prod resulting diagram
 * @return type of multiplication algorithm (MULT_M_MM, MULT_M_DM, ...)
 */
static int mult_type(diagram_t *op1, diagram_t *op2, diagram_t *prod)
{
    int type_prod = diagram_get_storage_type(prod) == CC_DIAGRAM_ON_DISK ? 1 : 0;
    int type_op1 = diagram_get_storage_type(op1) == CC_DIAGRAM_ON_DISK ? 1 : 0;
    int type_op2 = diagram_get_storage_type(op2) == CC_DIAGRAM_ON_DISK ? 1 : 0;

    int algorithms[2][2][2];
    algorithms[0][0][0] = MULT_M_MM;
    algorithms[0][1][0] = MULT_M_DM;
    algorithms[0][0][1] = MULT_M_MD;
    algorithms[0][1][1] = MULT_M_DD;
    algorithms[1][0][0] = MULT_D_MM;
    algorithms[1][1][0] = MULT_D_DM;
    algorithms[1][0][1] = MULT_D_MD;
    algorithms[1][1][1] = MULT_D_DD;

    return algorithms[type_prod][type_op1][type_op2];
}


/***************** auxiliary routines for tensor contractions *****************/

typedef struct {
    char yp;   // annihilation (y - Y(U)nichtozhenie) or creation (p - P(R)ozhdenie)
    int8_t ie;   // number of associated electron (particle)
} idxdata_t;


// comparator for struct idxdata_t (required by qsort)
int idxdata_t_less(const void *op1, const void *op2)
{
    idxdata_t *dat1 = (idxdata_t *) op1;
    idxdata_t *dat2 = (idxdata_t *) op2;
    // как числа, 'y' > 'p' и поэтому пойдет справа, мы эту ситуацию инвертируем знаком '-'
    if (dat1->yp > dat2->yp) {  // Y > P (уничтожения старше, они идут слева)
        return -1;
    }
    else if (dat1->yp < dat2->yp) {
        return 1;
    }
    else {
        // сравниваем номера частиц (электронов)
        return dat1->ie - dat2->ie;
    }
}


void target_order(int *ord1, int rk1, int *ord2, int rk2, int *ord3, int rk3)
{
    // number of dimensions ot be contracted
    int ncontr = (rk1 + rk2 - rk3) / 2;

    idxdata_t s1[CC_DIAGRAM_MAX_RANK];
    idxdata_t s2[CC_DIAGRAM_MAX_RANK];
    idxdata_t s3[CC_DIAGRAM_MAX_RANK];
    idxdata_t s3sorted[CC_DIAGRAM_MAX_RANK];
    int mapping[CC_DIAGRAM_MAX_RANK * 2];
    int i, j;


    // construct metainfo about dimensions using the order array
    // operand 1
    for (i = 0; i < rk1; i++) {
        if (ord1[i] <= rk1 / 2) {
            s1[i].yp = 'y';
            s1[i].ie = ord1[i];
        }
        else {
            s1[i].yp = 'p';
            s1[i].ie = ord1[i] - rk1 / 2;
        }
    }
    // operand 2
    for (i = 0; i < rk2; i++) {
        if (ord2[i] <= rk2 / 2) {
            s2[i].yp = 'y';
            s2[i].ie = ord2[i] + rk1 / 2; // чтобы гарантировать уникальные номера частиц
        }
        else {
            s2[i].yp = 'p';
            s2[i].ie = ord2[i] - rk2 / 2 + rk1 / 2;
        }
    }

    // соответствие номеров электронов для сворачиваемых измерений
    for (i = 0; i < 2 * CC_DIAGRAM_MAX_RANK; i++) {
        mapping[i] = 0;
    }
    for (i = 0; i < ncontr; i++) {
        mapping[s1[rk1 - ncontr + i].ie] = s2[rk2 - ncontr + i].ie;
        mapping[s2[rk2 - ncontr + i].ie] = s1[rk1 - ncontr + i].ie;
    }

    // строим строку из несвернутых индексов -- результат
    for (i = 0; i < rk1 - ncontr; i++) {
        s3[i] = s1[i];
    }
    for (i = 0; i < rk2 - ncontr; i++) {
        s3[rk1 - ncontr + i] = s2[i];
        // подменяем электроны, пришедшие из второго операнда,
        // с помощью построенного ранее соответствия номеров частиц
        if (mapping[s2[i].ie] != 0) {
            s3[rk1 - ncontr + i].ie = mapping[s2[i].ie];
        }
    }

    // упорядочиваем строку к "нормальной форме" сортировкой
    for (i = 0; i < rk3; i++) {
        s3sorted[i] = s3[i];
    }
    qsort(s3sorted, rk3, sizeof(idxdata_t), idxdata_t_less);

    // найдем перестановку s3sorted -> s3
    // просто смотрим в какую позицию j уехал индекс с позиции i
    for (i = 0; i < rk3; i++) {
        for (j = 0; j < rk3; j++) {
            if (s3sorted[i].yp == s3[j].yp && s3sorted[i].ie == s3[j].ie) {
                ord3[j] = i + 1;  // +1 потому что строки перестановок содержат только 1, 2, ...
                break;
            }
        }
    }
}
