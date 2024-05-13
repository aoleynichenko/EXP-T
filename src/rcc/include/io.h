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
 * Platform-independent I/O Unix-like interface.
 */


#ifndef CC_IO_H_INCLUDED
#define CC_IO_H_INCLUDED

#include <stddef.h>
#include <stdint.h>

typedef struct io_stat {
    size_t n_written;
    size_t n_read;
    size_t n_open;
    size_t n_closed;
    size_t n_created;
    size_t n_removed;
} io_stat_t;

int io_open(char *path, char *mode);

int io_close(int fd);

int64_t io_read(int fd, void *buf, size_t count);

int64_t io_write(int fd, const void *buf, size_t count);

int io_file_exists(char *filename);

int io_directory_exists(char *dirname);

int io_mkdir(char *path);

int io_chdir(char *path);

char *io_getcwd(char *path, size_t maxsz);

int io_remove(char *pathname);

size_t io_fsize(char *pathname);

size_t io_write_compressed(int fd, const void *buf, size_t count);

size_t io_read_compressed(int fd, void *buf, size_t count);

void print_compression_stats();

void io_statistics(io_stat_t *st);

#endif /* CC_IO_H_INCLUDED */
