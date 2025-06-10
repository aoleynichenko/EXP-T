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

/*
 * Operations with the "diagram stack".
 */

#include <stdio.h>
#include <string.h>

#include "dgstack.h"
#include "error.h"

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
                "(please, increase the CC_MAX_STACK_DEPTH parameter in file %s, line %d; currently == %d)",
                __FILE__, __LINE__, CC_MAX_STACK_DEPTH);
    }

    dg_stack[dg_stack_top++] = dg;
    return dg;
}


/**
 * replaces diagram named 'name' with another one (dg)
 */
diagram_t *diagram_stack_replace(char *name, diagram_t *dg)
{
    int i = diagram_stack_find_index(name);
    if (i == -1) {
        return NULL;
    }

    diagram_t *d_old = dg_stack[i];
    dg_stack[i] = dg;

    // destroy unused diagram
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


/**
 * all diagrams with index >= pos are deleted
 */
void restore_stack_pos(dg_stack_pos_t pos)
{
    for (int i = pos; i < dg_stack_top; i++) {
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
    for (int i = 0; i < dg_stack_top; i++) {
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
    diagram_t *diag = NULL;

    for (int i = 0; i < dg_stack_top; i++) {
        if (strcmp(name, dg_stack[i]->name) == 0) {
            diag = dg_stack[i];

            // shift diagrams in the stack by 1
            for (int j = i; j < dg_stack_top - 1; j++) {
                dg_stack[j] = dg_stack[j + 1];
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
    for (int i = 0; i < dg_stack_top; i++) {
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
    printf("\n diagram stack:\n");
    printf(" ----------------------------------------------------------------------------------------------\n");
    printf("       <name>      iii     iiv     iiu       #sb mem   #sb disk   #sb tot  size, GB     #unique\n");
    printf(" ----------------------------------------------------------------------------------------------\n");

    for (int idg = 0; idg < dg_stack_top; idg++) {
        diagram_t *dg = dg_stack[idg];

        /*
         * list of quasiparticles, valence lines and order of indices
         */
        char s_qparts[CC_DIAGRAM_MAX_RANK + 1];
        char s_valence[CC_DIAGRAM_MAX_RANK + 1];
        char s_order[CC_DIAGRAM_MAX_RANK + 1];
        diagram_get_quasiparticles(dg, s_qparts);
        diagram_get_valence(dg, s_valence);
        diagram_get_order(dg, s_order);

        /*
         * collect info about blocks
         */
        size_t count_blocks_in_mem = 0;
        size_t count_blocks_on_disk = 0;
        size_t count_unique_blocks = 0;

        for (int isb = 0; isb < dg->n_blocks; isb++) {
            block_t *block = dg->blocks[isb];
            if (block->storage_type == CC_DIAGRAM_ON_DISK) {
                count_blocks_on_disk++;
            }
            else {
                count_blocks_in_mem++;
            }
            if (block->is_unique) {
                count_unique_blocks++;
            }
        }

        /*
         * collect info about memory usage
         */
        size_t ram_used = 0;
        size_t disk_used = 0;
        diagram_get_memory_used(dg_stack[idg], &ram_used, &disk_used);
        double mem_used_gb = (double) (ram_used + disk_used) / (1024.0 * 1024.0 * 1024.0);

        printf(" [%3d] %-12s%-8s%-8s%-8s%10ld%10ld%10ld%10.3f%8ld/%ld\n",
               idg, dg_stack[idg]->name, s_qparts, s_valence, s_order,
               count_blocks_in_mem, count_blocks_on_disk, dg->n_blocks, mem_used_gb, count_unique_blocks,
               dg_stack[idg]->n_blocks);
    }

    printf(" ----------------------------------------------------------------------------------------------\n\n");
}


/**
 * Checks is a diagram exists.
 */
void assert_diagram_exists(char *name)
{
    diagram_t *dg = diagram_stack_find(name);
    if (dg == NULL) {
        errquit("diagram '%s' doesn't exist", name);
    }
}


/**
 * Renames diagram in the stack.
 */
void rename_diagram(char *old_name, char *new_name)
{
    assert_diagram_exists(old_name);
    diagram_t *dg = diagram_stack_find(old_name);
    strcpy(dg->name, new_name);
}

