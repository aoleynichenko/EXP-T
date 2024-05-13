/**
 * LIBUNF - tools for accessing Fortran binary unformatted files from projects
 * written in the C programming language
 *
 * 2024 Alexander Oleynichenko
 * alexvoleynichenko@gmail.com
 */

#ifndef LIBUNF_H_INCLUDED
#define LIBUNF_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

typedef enum {
    UNF_ACCESS_SEQUENTIAL,
    UNF_ACCESS_DIRECT,
    UNF_ACCESS_STREAM
} unf_access_t;

typedef enum {
    UNF_POS_BEGIN,
    UNF_POS_CURRENT,
    UNF_POS_END
} unf_position_t;

enum {
    UNF_ERROR = -1,
    UNF_SUCCESS = 0,
};

typedef struct {
    FILE *file_ptr;
    int access;
    int record_len; // is used only for direct-access files
    int error_flag;
} unf_file_t;

unf_file_t *unf_open(const char *path, const char *mode, unf_access_t access, ...);

int unf_close(unf_file_t *file);

int unf_write(unf_file_t *file, char *fmt, ...);

int unf_write_rec(unf_file_t *file, int rec, char *fmt, ...);

int unf_read(unf_file_t *file, char *fmt, ...);

int unf_read_rec(unf_file_t *file, int rec, char *fmt, ...);

int unf_next_rec_size(unf_file_t *file);

int unf_seek(unf_file_t *file, unf_position_t pos, int offset);

int unf_rewind(unf_file_t *file);

int unf_backspace(unf_file_t *file);

int unf_skip(unf_file_t *file);

int unf_eof(unf_file_t *file);

int unf_error(unf_file_t *file);

#ifdef __cplusplus
}
#endif

#endif // LIBUNF_H_INCLUDED
