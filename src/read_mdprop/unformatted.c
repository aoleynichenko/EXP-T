//
// Created by Alexander Oleynichenko on 29.12.2023.
//

#include <ctype.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#include "unformatted.h"


/**
 * Opens an unformatted file indicated by filename and returns a file stream
 * associated with that file. mode is used to determine the file access mode:
 * "r" - open a file for reading
 * "w" - create a file for writing
 * "a" - append to a file
 * If successful, returns a pointer to the object that controls the opened file stream.
 * On error, returns a null pointer (and errno is set in this case).
 */
unformatted_file_t *unformatted_open(const char *path, const char *mode)
{
    int flags;
    mode_t rights;
    int file_descr;

    if (strcmp(mode, "w") == 0) {
        flags = O_CREAT /*| O_TRUNC*/ | O_WRONLY;
        rights = S_IRUSR | S_IWUSR;
        file_descr = open(path, flags, rights);
    }
    else if (strcmp(mode, "a") == 0) {
        flags = O_CREAT | O_APPEND | O_WRONLY;
        rights = S_IRUSR | S_IWUSR;
        file_descr = open(path, flags, rights);
    }
    else if (strcmp(mode, "r") == 0) {
        flags = O_RDONLY;
        file_descr = open(path, flags);
    }
    else {
        // set errno: wrong argument
        return NULL;
    }

    if (file_descr == -1) {
        // set errno
        return NULL;
    }

    unformatted_file_t *file = (unformatted_file_t *) calloc(1, sizeof(unformatted_file_t));
    if (file == NULL) {
        // set errno
        return NULL;
    }

    file->file_descriptor = file_descr;

    return file;
}


/**
 * Closes the given unformatted file.
 * Whether or not the operation succeeds, the stream is no longer associated with a file.
 * The behavior is undefined if the value of the pointer stream is used after
 * unformatted_close() returns.
 * Return value: 0 on success, EOF otherwise
 */
int unformatted_close(unformatted_file_t *file)
{
    if (file == NULL) {
        // set errno
        return EOF;
    }

    int status = close(file->file_descriptor);
    if (status == -1) {
        // set errno
        return EOF;
    }

    free(file);

    return 0; // success
}


int unformatted_write(unformatted_file_t *file, char *fmt, ...)
{
    if (file == NULL) {
        return 0;
    }
}

// c, i4, i8, r4, r8, z4, z8
// c[], i4[], i8[], r4[], r8[], z4[], z8[]
// example:
// c[],i4,r8,i8[],r4[],i4
int unformatted_read(unformatted_file_t *file, const char *fmt, ...)
{
    if (file == NULL) {
        return 0;
    }

    /*
     * read size of the record in bytes
     */
    int record_size_1 = 0;
    int err = read(file->file_descriptor, &record_size_1, sizeof(int32_t));
    if (err != sizeof(int32_t)) {
        return 0;
    }

    //printf("record size: %d\n", record_size_1);

    int n_arguments_read = 0;
    int n_bytes_read = 0;
    va_list ap;
    va_start(ap, fmt);

    //
    //printf("fmt = %s\n", fmt);

    for (const char *p = fmt; *p;) {
        int type;
        int type_len;

        //printf("remainder = %s\n", p);

        if (*p == ',') {
            // nothing to do, just skip
            p++;
            continue;
        }

        int num_repeat = 1;
        if (isdigit(*p)) {
            num_repeat = 0;
            while (isdigit(*p)) {
                num_repeat = 10 * num_repeat + (*p - '0');
                p++;
            }
            //printf("num_repeat = %d\n", num_repeat);
        }

        if (*p == 'c') {
            //printf("type char, *p == %c\n", *p);
            // 'c'
            type = TYPE_CHAR;
            type_len = sizeof(char); // default: only one character
            p++;

            // if needed - Fortran string
            if (isdigit(*p)) {
                type_len = 0;
                while (isdigit(*p)) {
                    type_len = 10 * type_len + (*p - '0');
                    p++;
                }
            }
            //printf("Fortran string len = %d\n", type_len);
        }
        else if (*p == 'i' && *(p + 1) == '4') {
            // 'i4'
            type = TYPE_INTEGER_4;
            type_len = sizeof(int32_t);
        }
        else if (*p == 'i' && *(p + 1) == '8') {
            // 'i8'
            type = TYPE_INTEGER_8;
            type_len = sizeof(int64_t);
        }
        else if (*p == 'r' && *(p + 1) == '4') {
            // 'r4'
            type = TYPE_REAL_4;
            type_len = sizeof(float);
        }
        else if (*p == 'r' && *(p + 1) == '8') {
            // 'r8'
            type = TYPE_REAL_8;
            type_len = sizeof(double);
        }
        else if (*p == 'z' && *(p + 1) == '4') {
            // 'z4'
            type = TYPE_COMPLEX_4;
            type_len = sizeof(float _Complex);
        }
        else if (*p == 'z' && *(p + 1) == '8') {
            // 'z8'
            type = TYPE_COMPLEX_8;
            type_len = sizeof(double _Complex);
        }
        else {
            errno = EINVAL;
            break;
        }

        if (type != TYPE_CHAR) {
            p += 2;
            //printf("remainder = %s\n", p);
        }

        /*
         * array or not?
         * array dimension can be specified by either integer-4 '[i4]' or integer-8 '[i8]' number
         */
        int is_array = 0;
        int array_dim_type = TYPE_INTEGER_4;
        if (*p == '[' && *(p + 1) == 'i' && *(p + 2) == '4' && *(p + 3) == ']') {
            is_array = 1;
            array_dim_type = TYPE_INTEGER_4;
            p += 4;
        }
        else if (*p == '[' && *(p + 1) == 'i' && *(p + 2) == '8' && *(p + 3) == ']') {
            is_array = 1;
            array_dim_type = TYPE_INTEGER_8;
            p += 4;
        }

        for (int i_repeat = 0; i_repeat < num_repeat; i_repeat++) {
            /*
             * get pointer to data to be written
             */
            void *data_ptr = va_arg(ap, void *);

            /*
             * if needed, get array length
             */
            size_t array_dim = 1;
            if (is_array && array_dim_type == TYPE_INTEGER_4) {
                int32_t *array_dim_ptr = va_arg(ap, int32_t *);
                if (array_dim_ptr == NULL) {
                    errno = EINVAL; // wrong argument
                    break;
                }
                array_dim = *array_dim_ptr;
            }
            else if (is_array && array_dim_type == TYPE_INTEGER_8) {
                int64_t *array_dim_ptr = va_arg(ap, int64_t *);
                if (array_dim_ptr == NULL) {
                    errno = EINVAL; // wrong argument
                    break;
                }
                array_dim = *array_dim_ptr;
            }

            /*
             * read entry from the unformatted file
             */
            size_t n_bytes = array_dim * type_len;
            err = read(file->file_descriptor, data_ptr, n_bytes);
            if (err != n_bytes) {
                break;
            }

            n_bytes_read += n_bytes;
            n_arguments_read++;
        }
    }

    va_end(ap);

    //printf("arguments read: %d\n", n_arguments_read);
    //printf("num bytes read: %d\n", n_bytes_read);

    /*
     * rewind to the end of the record if needed
     */
    if (n_bytes_read != record_size_1) {
        off_t offset = record_size_1 - n_bytes_read;
        //printf("do lseek %d bytes\n", offset);
        lseek(file->file_descriptor, offset, 1);
    }

    /*
     * read size of the record in bytes
     */
    int record_size_2 = 0;
    err = read(file->file_descriptor, &record_size_2, sizeof(int32_t));
    if (err != sizeof(int32_t)) {
        return 0;
    }

    /*
     * two sizes of a record must coincide
     */
    if (record_size_1 != record_size_2) {
        return 0;
    }

    //printf("%d %d\n", record_size_1, record_size_2);

    return n_arguments_read;
}


int unformatted_next_record_size(unformatted_file_t *file)
{
    if (file == NULL) {
        return 0;
    }

    /*
     * read size of the record in bytes
     */
    int record_size = 0;
    int err = read(file->file_descriptor, &record_size, sizeof(int32_t));
    if (err != sizeof(int32_t)) {
        return 0;
    }

    off_t offset = -4;
    lseek(file->file_descriptor, offset, 1);

    return record_size;
}


void unformatted_rewind_last_record(unformatted_file_t *file)
{
    off_t offset = -4;
    int err = lseek(file->file_descriptor, offset, 1);
    if (err == -1) {
        return;
    }

    int32_t record_size = 0;
    read(file->file_descriptor, &record_size, sizeof(int32_t));

    offset = -(2 * sizeof(int32_t) + record_size);
    lseek(file->file_descriptor, offset, 1);
}
