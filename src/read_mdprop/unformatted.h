//
// Created by Alexander Oleynichenko on 29.12.2023.
//

#ifndef UNFORMATTED_H
#define UNFORMATTED_H

enum {
    TYPE_LEN_8 = 8,
    TYPE_LEN_16 = 16,
    TYPE_LEN_32 = 32,
    TYPE_LEN_64 = 64
};

enum {
    TYPE_CHAR,
    TYPE_INTEGER_4,
    TYPE_INTEGER_8,
    TYPE_REAL_4,
    TYPE_REAL_8,
    TYPE_COMPLEX_4,
    TYPE_COMPLEX_8,
};

typedef struct {
    int file_descriptor;
} unformatted_file_t;

unformatted_file_t *unformatted_open(const char *path, const char *mode);

int unformatted_close(unformatted_file_t *file);

int unformatted_write(unformatted_file_t *file, char *fmt, ...);

int unformatted_read(unformatted_file_t *file, const char *fmt, ...);

int unformatted_next_record_size(unformatted_file_t *file);

void unformatted_rewind_last_record(unformatted_file_t *file);

#endif //UNFORMATTED_H
