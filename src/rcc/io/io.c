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

/*******************************************************************************
 *                                I/O interface
 *
 * Routines for binary (and I hope fast) input/output.
 * This interface is system-independent.
 *
 * For Unix-like systems the following routines are simply wrappers
 * for Unix system calls (open, close, read, write).
 *
 ******************************************************************************/

#include "io.h"

#include <errno.h>
#include <string.h>
#include <stdio.h>

#include "error.h"
#include "memory.h"

#define CHUNK_SIZE (1024*1024*1024)

// Linux version
#if defined(__linux__) || defined(__APPLE__)

#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

static struct io_stat IO_STAT = {0};


// mode = "w", "r"
int io_open(char *path, char *mode)
{
    mode_t md;
    int ret;
    int creat = 0;

    if (strcmp(mode, "w") == 0) {
        md = O_CREAT /*| O_TRUNC*/ | O_WRONLY;
        creat = 1;
        ret = open(path, md, S_IRUSR | S_IWUSR);
    }
    else if (strcmp(mode, "a") == 0) {
        md = O_CREAT | O_APPEND | O_WRONLY;
        creat = 1;
        ret = open(path, md, S_IRUSR | S_IWUSR);
    }
    else if (strcmp(mode, "r") == 0) {
        md = O_RDONLY;
        ret = open(path, md);
    }
    else {
        printf("io_open(): wrong mode '%s'", mode);
        return -1;
    }

    if (ret == -1) {
        //printf("path = '%s'\n", path);
        //printf("io_open(): %s", strerror(errno));
    }

    if (creat) {
        IO_STAT.n_created++;
    }
    IO_STAT.n_open++;
    return ret;
}


int io_close(int fd)
{
    int status = close(fd);
    if (status == -1) {
        errquit("io_close(): %s", strerror(errno));
    }

    IO_STAT.n_closed++;
    return status;
}


int64_t io_read(int fd, void *buf, size_t count)
{
    /*ssize_t status = read(fd, buf, count);
    if (status != count) {
        errquit("io_read(): %s", strerror(errno));
    }

    IO_STAT.n_read += status;
    return status;*/
    size_t count_total = count;

    while (count > 0) {
        size_t nr = (count > CHUNK_SIZE) ? CHUNK_SIZE : count;
        count -= nr;

        ssize_t status = read(fd, buf, nr);
        if (status != nr) {
            errquit("io_read(): %s", strerror(errno));
        }

        buf += nr;
        IO_STAT.n_read += status;
    }

    return count_total;
}


int64_t io_write(int fd, const void *buf, size_t count)
{
    size_t count_total = count;

    while (count > 0) {
        size_t nw = (count > CHUNK_SIZE) ? CHUNK_SIZE : count;
        count -= nw;

        int status = write(fd, buf, nw);
        if (status != nw) {
            errquit("io_write(): %s", strerror(errno));
        }

        buf += nw;
        IO_STAT.n_written += status;
    }

    return count_total;
}


int io_file_exists(char *filename)
{
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
}


// returns 1 if exists, 0 otherwise
int io_directory_exists(char *dirname)
{
    struct stat st = {0};

    if (stat(dirname, &st) == -1) {
        return 0;
    }
    else {
        return 1;
    }
}


// returns zero on success, or -1 if an error occurred (in which case, errno is set appropriately).
int io_mkdir(char *path)
{
    int s = mkdir(path, 0777);
    if (s == -1) {
        printf("path = '%s'\n", path);
        errquit("io_mkdir(): %s", strerror(errno));
    }
    return s;
}


// On success, zero is returned.  On error, -1 is returned, and errno is set appropriately.
int io_chdir(char *path)
{
    int s = chdir(path);
    if (s == -1) {
        printf("path = '%s'\n", path);
        errquit("io_chdir(): %s", strerror(errno));
    }
    return s;
}


char *io_getcwd(char *buf, size_t maxsz)
{
    char *s = getcwd(buf, maxsz);

    if (s == NULL) {
        errquit("io_getcwd(): (errno = %d, strerror = %s)", errno, strerror(errno));
    }

    return s;
}


/**
 * remove a file or directory
 */
int io_remove(char *pathname)
{
    int s;

    s = remove(pathname);
    if (s == -1) {
        errquit("io_remove(): (errno = %d, strerror = %s)", errno, strerror(errno));
    }

    IO_STAT.n_removed++;
    return s;
}


/**
 * file size (in bytes)
 */
size_t io_fsize(char *pathname)
{
    struct stat st = {0};

    if (stat(pathname, &st) == -1) {
        return 0;
    }
    else {
        return st.st_size;
    }
}


void io_statistics(io_stat_t *st)
{
    *st = IO_STAT;
}


#else
#error "only Linux and Mac OS X platforms are currently supported"
#endif

/*******************************************************************************
 *
 * Data compression -- interface to data compression libraries.
 * Implemented bindings: LZ4
 * All subroutines perform compression only if data compression is required by
 * the 'compress' option (see options.h)
 * References:
 * LZ4   https://github.com/lz4/lz4 (Copyright (C) 2011-present, Yann Collet)
 *
 * 2019 Alexander Oleynichenko
 ******************************************************************************/

#include <math.h>
#include <stdio.h>

#include "lz4.h"
#include "options.h"
#include "utils.h"

static void grow_zbuf(size_t count);

// buffer for compressed data
static char *zbuf;
static size_t zbuf_len = 0;

// for collecting statistics
static int compress_module_is_initialized = 0;
static int histogram[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // size = 20
static double min_compression_ratio = 100.0;
static double max_compression_ratio = 0.0;
static double mean_compression_ratio = 0.0;
static double mean_compression_ratio_10 = 0.0; // for values of ratios <= 10
static int n_compressions = 0;  // for evaluation of mean compression ratio
static int n_compressions_10 = 0;  // for values <= 10


/*******************************************************************************
 * io_write_compressed
 *
 * Writes data to disk. Performs data compression if required (see options.h,
 * cc_options_t->compress field).
 * NOTE: compression is not performed if buffer length <= 16
 * Arguments:
 *   fd     file descriptor
 *   buf    data buffer
 *   count  number of bytes in the buffer
 * Returns:
 *   number of bytes written (<= count, < if data were compressed)
 ******************************************************************************/
size_t io_write_compressed(int fd, const void *buf, size_t count)
{
    int src_size, dst_size;

    // if no compression is required
    if (cc_opts->compress == 0 || count <= 16) {
        return io_write(fd, buf, count);
    }
    // else compress data ...

    src_size = count;
    if (src_size >= LZ4_MAX_INPUT_SIZE) {
        errquit("io_write_compressed(): src_size (%ld) > LZ4_MAX_INPUT_SIZE (%ld). "
                "Please, decrease the tize size or contact A. Oleynichenko for "
                "increasing this limit (alexvoleynichenko@gmail.com)", src_size, LZ4_MAX_INPUT_SIZE);
    }

    dst_size = LZ4_compressBound(src_size);

    if (dst_size > zbuf_len) {
        grow_zbuf(dst_size);
    }

    dst_size = LZ4_compress_default(buf, zbuf, src_size, dst_size);

    // save statistics
    double comp_ratio = ((double) src_size) / dst_size;
    if (comp_ratio < min_compression_ratio) {
        min_compression_ratio = comp_ratio;
    }
    if (comp_ratio > max_compression_ratio) {
        max_compression_ratio = comp_ratio;
    }
    if (0.0 <= comp_ratio && comp_ratio <= 10.0) {
        histogram[(int) ceil(comp_ratio * 2)] += 1;  // step for histogram == 0.5
        n_compressions_10 += 1;
        mean_compression_ratio_10 += comp_ratio;
    }
    mean_compression_ratio += comp_ratio;
    n_compressions += 1;

    // write: [actual length of compressed data] [data]
    io_write(fd, &dst_size, sizeof(int));
    return io_write(fd, zbuf, dst_size);
}


/*******************************************************************************
 * io_read_compressed
 *
 * Reads data from disk. Performs data decompression if required (see options.h,
 * cc_options_t->compress field).
 * NOTE: compression is not performed if buffer length <= 16
 * Arguments:
 *   fd     file descriptor
 *   buf    data buffer (output buffer)
 *   count  number of bytes in the buffer
 * Returns:
 *   number of bytes obtained after decompression (must be == count)
 ******************************************************************************/
size_t io_read_compressed(int fd, void *buf, size_t count)
{
    int compressed_size;
    int dst_capacity;

    // if no decompression is needed
    if (cc_opts->compress == 0 || count <= 16) {
        return io_read(fd, buf, count);
    }
    // else decompress data ...

    dst_capacity = count;

    size_t s = io_read(fd, &compressed_size, sizeof(int));

    if (compressed_size >= LZ4_MAX_INPUT_SIZE) {
        errquit("io_read_compressed(): compressed_size (%ld) > LZ4_MAX_INPUT_SIZE (%ld). "
                "Please, decrease the tize size or contact A. Oleynichenko for "
                "increasing this limit (alexvoleynichenko@gmail.com)", compressed_size, LZ4_MAX_INPUT_SIZE);
    }

    // increase zbuf ?
    if (compressed_size > zbuf_len) {
        grow_zbuf(compressed_size);
    }

    io_read(fd, zbuf, compressed_size);

    return LZ4_decompress_safe(zbuf, buf, compressed_size, dst_capacity);
}


/*******************************************************************************
 * print_compression_stats
 *
 * Prints statistical data about compression.
 ******************************************************************************/
void print_compression_stats()
{
    if (cc_opts->compress == CC_COMPRESS_NONE) {
        return;
    }

    printf(" Compression statistics:\n");
    printf("   min compression ratio = %.3f\n", min_compression_ratio);
    printf("   max compression ratio = %.3f\n", max_compression_ratio);
    printf("   mean compression ratio = %.3f\n", mean_compression_ratio / n_compressions);
    printf("   mean compression ratio (<=10) = %.3f\n", mean_compression_ratio_10 / n_compressions_10);
    printf("   distribution of compression ratios:\n");

    // build histogram
    int maxval = imax(20, histogram);
    for (int i = 0; i < 20; i++) {
        int nbars = ((double) histogram[i] / maxval) * 50;
        printf("   [%4.1f -%4.1f] |", 0.5 * i, 0.5 * i + 0.5);
        for (int j = 0; j < nbars; j++) {
            printf("=");
        }
        printf("\n");
    }
}


/*******************************************************************************
 * grow_zbuf
 *
 * Reallocates memory for the temporary buffer used for (de)compression.
 ******************************************************************************/
void grow_zbuf(size_t count)
{
    if (zbuf_len == 0) {
        zbuf = (char *) cc_malloc(count);
        zbuf_len = count;
    }
    else {  // realloc
        char *zbuf2 = (char *) cc_malloc(count);
        memcpy(zbuf2, zbuf, zbuf_len);
        cc_free(zbuf);
        zbuf = zbuf2;
        zbuf_len = count;
    }
}
