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

/*******************************************************************************
 * dgstack.c
 *
 * Operations with the "diagram stack".
 *
 * 2018-2020 Alexander Oleynichenko
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <options.h>

#include "datamodel.h"
#include "error.h"
//#include "spinors.h"
#include "utils.h"

#define CC_MAX_STACK_DEPTH 1024

static diagram_t *dg_stack[CC_MAX_STACK_DEPTH];
static int dg_stack_top = 0;   // next free position

/**
 * pushes diagram 'dg' on the top of the diagram stack
 * @return pointer to the diagram pushed to the stack
 */
diagram_t *diagram_stack_push(diagram_t *dg)
{
    if (dg_stack_top == CC_MAX_STACK_DEPTH) {
        errquit("diagram_stack_push: stack overflow "
                "(please, increase the CC_MAX_STACK_DEPTH parameter, currently == %d)",
                CC_MAX_STACK_DEPTH);
    }

    dg_stack[dg_stack_top++] = dg;
    return dg;
}


/**
 * replaces diagram named 'name' with another one (dg)
 */
diagram_t *diagram_stack_replace(char *name, diagram_t *dg)
{
    diagram_t *d_old;
    int i;

    i = diagram_stack_find_index(name);
    if (i == -1) {
        return NULL;
    }

    d_old = dg_stack[i];
    dg_stack[i] = dg;

    // delete dg
    diagram_delete(d_old);

    return dg;
}


/**
 * @return top of the diagram stack
 */
dg_stack_pos_t get_stack_pos()
{
    return dg_stack_top;
}


#include "memory.h"

void print_mem_usage()
{
    double curr_usage = cc_get_current_memory_usage() / (1024.0 * 1024.0 * 1024.0);
    double peak_usage = cc_get_peak_memory_usage() / (1024.0 * 1024.0 * 1024.0);
    printf(" $ current memory usage = %.3f Gb; peak memory usage = %.3f Gb\n", curr_usage, peak_usage);
}

/**
 * all diagrams with index >= pos are deleted
 */
void restore_stack_pos(dg_stack_pos_t pos)
{
    int i;

    for (i = pos; i < dg_stack_top; i++){
        diagram_delete(dg_stack[i]);
    }
    dg_stack_top = pos;
    // pos = next free position
}


/**
 * find diagram in the stack by name
 * @return index of the diagram
 */
int diagram_stack_find_index(char *name)
{
    int i;

    for (i = 0; i < dg_stack_top; i++){
        if (strcmp(name, dg_stack[i]->name) == 0) {
            return i;
        }
    }

    return -1;
}


/**
 * removes diagram from the stack and deallocate all resources
 * associated with this diagram.
 * NOTE: this function breaks positions in the diagram stack!
 */
void diagram_stack_erase(char *name)
{
    int i;
    diagram_t *diag = NULL;

    for (i = 0; i < dg_stack_top; i++){
        if (strcmp(name, dg_stack[i]->name) == 0) {
            diag = dg_stack[i];

            // shift diagrams in the stack by 1
            for (int j = i; j < dg_stack_top-1; j++) {
                dg_stack[j] = dg_stack[j+1];
            }
            dg_stack_top--;

            break;
        }
    }

    if (diag != NULL) {
        diagram_delete(diag);
    }
}


/**
 * find diagram in the stack by name
 * @return pointer to the diagram found
 */
diagram_t *diagram_stack_find(char *name)
{
    int i;

    for (i = 0; i < dg_stack_top; i++){
        if (strcmp(name, dg_stack[i]->name) == 0) {
            return dg_stack[i];
        }
    }

    return NULL;
}


/**
 * print stack of diagrams with short info about them
 */
void diagram_stack_print()
{
    size_t isb;
    size_t nsb_mem, nsb_disk;
    size_t disk_used, ram_used;
    diagram_t *dg;
    block_t *sb;
    char s_qparts[CC_DIAGRAM_MAX_RANK + 1];
    char s_valence[CC_DIAGRAM_MAX_RANK + 1];
    char s_order[CC_DIAGRAM_MAX_RANK + 1];

    if (dg_stack_top == 0) {
        printf("diagram_stack_print: stack is empty\n");
        return;
    }

    printf("\ndiagram stack:\n");
    printf("     <name>      iii     iiv     iiu       #sb mem   #sb disk   #sb tot  size, GB\n");

    for (int i = 0; i < dg_stack_top; i++){
        dg = dg_stack[i];
        diagram_get_quasiparticles(dg, s_qparts);
        diagram_get_valence(dg, s_valence);
        diagram_get_order(dg, s_order);

        // collect info about blocks
        nsb_mem = 0;
        nsb_disk = 0;
        int nsb_uniq = 0;
        for (isb = 0; isb < dg->n_blocks; isb++){
            sb = dg->blocks[isb];
            if (sb->storage_type == CC_DIAGRAM_ON_DISK) {
                nsb_disk++;
            }
            else{
                nsb_mem++;
            }
            if (sb->is_unique) {
                nsb_uniq++;
            }
        }

        diagram_get_memory_used(dg_stack[i], &ram_used, &disk_used);
        printf("[%2d] %-12s%-8s%-8s%-8s%10ld%10ld%10ld%10.3f%4d/%d\n", i, dg_stack[i]->name, s_qparts, s_valence, s_order,
               nsb_mem, nsb_disk, dg->n_blocks, (ram_used + disk_used) / (1024.0 * 1024.0 * 1024.0), nsb_uniq, dg_stack[i]->n_blocks);
    }
    printf("end of diagram stack\n");
}
