/**
 * LIBUNF - tools for accessing Fortran binary unformatted files from projects
 * written in the C programming language
 *
 * 2024 Alexander Oleynichenko
 * alexvoleynichenko@gmail.com
 */

#include <assert.h>
#include <complex.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>

#include "libunf.h"

enum {
    TYPE_CHAR,
    TYPE_INTEGER_1,
    TYPE_INTEGER_2,
    TYPE_INTEGER_4,
    TYPE_INTEGER_8,
    TYPE_REAL_4,
    TYPE_REAL_8,
    TYPE_COMPLEX_4,
    TYPE_COMPLEX_8,
};

static int try_read_bytes(unf_file_t *file, char *fmt, size_t *n_bytes_read, int *n_args_read, va_list ap);

static int try_write_bytes(unf_file_t *file, char *fmt, size_t *n_bytes_written, int *n_args_written, va_list ap);

static int fmt_get_type_size(char **fmt, int *data_type, int *type_size, int *num_repeats);

static int seek_backward(unf_file_t *file);

static int seek_forward(unf_file_t *file);


/**
 * Opens an unformatted file indicated by filename and returns a file stream
 * associated with that file.
 *
 * mode is used to determine the file access mode:
 * "r" - open a file for reading
 * "w" - create a file for writing
 * "a" - append to a file
 *
 * access:
 * UNF_ACCESS_SEQUENTIAL
 * UNF_ACCESS_DIRECT
 * UNF_ACCESS_STREAM
 *
 * If successful, returns a pointer to the object that controls the opened file stream.
 * On error, returns a null pointer (and errno is set in this case).
 */
unf_file_t *unf_open(const char *path, const char *mode, unf_access_t access, ...)
{
    FILE *file = NULL;
    int record_len = 0;

    if (!(access == UNF_ACCESS_SEQUENTIAL ||
          access == UNF_ACCESS_DIRECT ||
          access == UNF_ACCESS_STREAM)) {
        errno = EINVAL;
        return NULL;
    }

    // get record length in the special case of the direct unformatted I/O
    if (access == UNF_ACCESS_DIRECT) {
        va_list ap;
        va_start(ap, access);
        record_len = va_arg(ap, int);
        va_end(ap);

        if (record_len <= 0) {
            errno = EINVAL;
            return NULL;
        }
    }

    // parse mode and open file if success
    if (strcmp(mode, "w") == 0) {
        file = fopen(path, "wb");
    }
    else if (strcmp(mode, "a") == 0) {
        file = fopen(path, "ab");
    }
    else if (strcmp(mode, "r") == 0) {
        file = fopen(path, "rb");
    }
    else {
        errno = EINVAL;
        return NULL;
    }

    if (file == NULL) {
        return NULL;
    }

    unf_file_t *unf_file = (unf_file_t *) calloc(1, sizeof(unf_file_t));
    if (unf_file == NULL) {
        fclose(file);
        return NULL;
    }

    unf_file->file_ptr = file;
    unf_file->access = access;
    unf_file->record_len = record_len;
    unf_file->error_flag = 0;

    return unf_file;
}


/**
 * Closes the given unformatted file.
 * Whether or not the operation succeeds, the stream is no longer associated with a file.
 * The behavior is undefined if the value of the pointer stream is used after
 * unf_close() returns.
 * Returns UNF_SUCCESS upon success, UNF_ERROR otherwise
 */
int unf_close(unf_file_t *file)
{
    if (file == NULL) {
        errno = EINVAL;
        return UNF_ERROR;
    }

    int status = fclose(file->file_ptr);
    if (status == EOF) {
        return UNF_ERROR;
    }

    free(file);

    return UNF_SUCCESS; // success
}


/**
 * Writes data (the next record) to sequential and stream access files.
 * Arrays must be passed to the function by pointer, scalars and array dimensions
 * must be passed by value.
 * Format string is used to specify argument data types.
 * Format string consists of type specifiers separated by commas.
 *
 * Available type specifiers are:
 * "c" - character
 * "c<N>" - Fortran string (not null-terminated), N stands for the lenght in bytes
 * "i1", "i2", "i4", "i8" - integer numbers from integer(1) = int8_t to integer(8) = int64_t
 * "r4", "r8" - real numbers, real(4) = float and real(8) = double, respectively
 * "z4", "z8" - complex numbers, complex(4) = float _Complex and complex(8) = double _Complex, respectively
 *
 * Arrays are specified by addition of the "[i4]" or "[i8]" formats,
 * where "i4" (or "i8") denotes the type of the variable used to specify array length.
 *
 * Returns the number of receiving arguments successfully written
 * (which may be zero in case a matching failure occurred before
 * the first receiving argument was assigned).
 *
 * The unf_eof() and unf_error() functions and the errno variable are to be used to handle errors.
 */
int unf_write(unf_file_t *file, char *fmt, ...)
{
    if (file == NULL ||
        file->access == UNF_ACCESS_DIRECT) {
        errno = EINVAL;
        return 0;
    }

    // only for sequential files: size of the record in bytes
    // here: just template, to be overwritten in future
    if (file->access == UNF_ACCESS_SEQUENTIAL) {
        int32_t record_size = 0;
        size_t n_written = fwrite(&record_size, 1, sizeof(int32_t), file->file_ptr);
        if (n_written != sizeof(int32_t)) {
            file->error_flag = 1;
            return 0;
        }
    }

    // write target bytes.
    // the same code for both sequential and stream access mode.
    size_t n_bytes_written = 0;
    int n_arguments_written = 0;

    va_list ap;
    va_start(ap, fmt);
    int err = try_write_bytes(file, fmt, &n_bytes_written, &n_arguments_written, ap);
    va_end(ap);

    if (err == UNF_ERROR) {
        file->error_flag = 1;
        return n_arguments_written;
    }

    // only for sequential files: size of the record in bytes
    if (file->access == UNF_ACCESS_SEQUENTIAL) {
        // begin entry
        int err = fseek(file->file_ptr, -sizeof(int32_t) - n_bytes_written, SEEK_CUR);
        if (err != 0) {
            file->error_flag = 1;
            return n_arguments_written;
        }

        size_t n_written = fwrite(&n_bytes_written, 1, sizeof(int32_t), file->file_ptr);
        if (n_written != sizeof(int32_t)) {
            file->error_flag = 1;
            return n_arguments_written;
        }

        // end entry
        err = fseek(file->file_ptr, n_bytes_written, SEEK_CUR);
        if (err != 0) {
            file->error_flag = 1;
            return n_arguments_written;
        }

        n_written = fwrite(&n_bytes_written, 1, sizeof(int32_t), file->file_ptr);
        if (n_written != sizeof(int32_t)) {
            file->error_flag = 1;
            return n_arguments_written;
        }
    }

    return n_arguments_written;
}


/**
 * Writes data to the record number 'rec' of a direct access file.
 * Arrays must be passed to the function by pointer, scalars and array dimensions
 * must be passed by value.
 * Format string is used to specify argument data types; see the description of the
 * unf_write() function for details.
 *
 * Returns the number of receiving arguments successfully written
 * (which may be zero in case a matching failure occurred before
 * the first receiving argument was assigned).
 *
 * The unf_eof() and unf_error() functions and the errno variable are to be used to handle errors.
 */
int unf_write_rec(unf_file_t *file, int rec, char *fmt, ...)
{
    if (file == NULL ||
        file->access != UNF_ACCESS_DIRECT ||
        file->record_len <= 0 ||
        rec < 1) {
        errno = EINVAL;
        return 0;
    }

    // find the required entry
    off_t offset = (rec - 1) * file->record_len;
    int err = fseek(file->file_ptr, offset, SEEK_SET);
    if (err != 0) {
        file->error_flag = 1;
        return 0;
    }

    // write target bytes.
    // the same code for both sequential and stream access mode.
    size_t n_bytes_written = 0;
    int n_arguments_written = 0;

    va_list ap;
    va_start(ap, fmt);
    err = try_write_bytes(file, fmt, &n_bytes_written, &n_arguments_written, ap);
    va_end(ap);

    if (n_bytes_written > file->record_len) {
        file->error_flag = 1;
    }
    else if (n_bytes_written < file->record_len) {
        // if needed: write zeros to preserve correct alignment of records
        for (int i = 0; i < file->record_len - n_bytes_written; i++) {
            fputc(0, file->file_ptr);
        }
    }

    if (err == UNF_ERROR) {
        file->error_flag = 1;
    }

    return n_arguments_written;
}


/**
 * Reads data (the next record) from sequential and stream access files.
 * All arguments, both arrays and scalars, must be passed to the function by pointer.
 * Format string is used to specify argument data types.
 * Format string consists of type specifiers separated by commas.
 *
 * Available type specifiers are:
 * "c" - character
 * "c<N>" - Fortran string (not null-terminated), N stands for the lenght in bytes
 * "i1", "i2", "i4", "i8" - integer numbers from integer(1) = int8_t to integer(8) = int64_t
 * "r4", "r8" - real numbers, real(4) = float and real(8) = double, respectively
 * "z4", "z8" - complex numbers, complex(4) = float _Complex and complex(8) = double _Complex, respectively
 *
 * Arrays are specified by addition of the "[i4]" or "[i8]" formats,
 * where "i4" (or "i8") denotes the type of the variable used to specify array length.
 *
 * Returns the number of receiving arguments successfully assigned
 * (which may be zero in case a matching failure occurred before
 * the first receiving argument was assigned).
 *
 * The unf_eof() and unf_error() functions and the errno variable are to be used to handle errors.
 */
int unf_read(unf_file_t *file, char *fmt, ...)
{
    if (file == NULL ||
        file->access == UNF_ACCESS_DIRECT) {
        errno = EINVAL;
        return 0;
    }

    // only for sequential files: size of the record in bytes
    int32_t record_size = 0;
    if (file->access == UNF_ACCESS_SEQUENTIAL) {
        size_t n_read = fread(&record_size, 1, sizeof(int32_t), file->file_ptr);
        if (n_read != sizeof(int32_t)) {
            file->error_flag = 1;
            return 0;
        }
    }

    // read target bytes.
    // the same code for both sequential and stream access mode.
    size_t n_bytes_read = 0;
    int n_arguments_read = 0;

    va_list ap;
    va_start(ap, fmt);
    int err = try_read_bytes(file, fmt, &n_bytes_read, &n_arguments_read, ap);
    va_end(ap);

    if (err == UNF_ERROR) {
        file->error_flag = 1;
        return n_arguments_read;
    }

    // only for sequential files: check the record size and rewind to the end of the entry

    if (file->access == UNF_ACCESS_SEQUENTIAL) {
        // rewind to the end of the record if needed
        if (n_bytes_read != record_size) {
            off_t offset = record_size - n_bytes_read;
            fseek(file->file_ptr, offset, SEEK_CUR);
        }

        // read size of the record in bytes
        int32_t record_size_2 = 0;
        size_t n_read = fread(&record_size_2, 1, sizeof(int32_t), file->file_ptr);
        if (n_read != sizeof(int32_t)) {
            file->error_flag = 1;
            return n_arguments_read;
        }

        // two sizes of a record must coincide
        if (record_size != record_size_2) {
            file->error_flag = 1;
            return n_arguments_read;
        }
    }

    return n_arguments_read;
}


/**
 * Reads data from record number 'rec'.
 * For direct access files only.
 * All arguments, both arrays and scalars, must be passed to the function by pointer.
 * Format string is used to specify argument data types; see the description of the
 * unf_read() function for details.
 *
 * Returns the number of receiving arguments successfully assigned
 * (which may be zero in case a matching failure occurred before
 * the first receiving argument was assigned).
 *
 * The unf_eof() and unf_error() functions and the errno variable are to be used to handle errors.
 */
int unf_read_rec(unf_file_t *file, int rec, char *fmt, ...)
{
    if (file == NULL ||
        file->access != UNF_ACCESS_DIRECT ||
        file->record_len <= 0 ||
        rec < 1) {
        errno = EINVAL;
        return 0;
    }

    // find the required entry
    off_t offset = (rec - 1) * file->record_len;
    int err = fseek(file->file_ptr, offset, SEEK_SET);
    if (err != 0) {
        file->error_flag = 1;
        return 0;
    }

    // read target bytes
    size_t n_bytes_read = 0;
    int n_arguments_read = 0;

    va_list ap;
    va_start(ap, fmt);
    err = try_read_bytes(file, fmt, &n_bytes_read, &n_arguments_read, ap);
    va_end(ap);

    if (n_bytes_read > file->record_len) {
        file->error_flag = 1;
    }

    if (err == UNF_ERROR) {
        file->error_flag = 1;
    }

    return n_arguments_read;
}


/**
 * Returns size of the next record (in bytes).
 * For sequential access files only.
 */
int unf_next_rec_size(unf_file_t *file)
{
    if (file == NULL ||
        file->access != UNF_ACCESS_SEQUENTIAL) {
        errno = EINVAL;
        return 0;
    }

    /*
     * read size of the record in bytes
     */
    int record_size = 0;
    size_t err = fread(&record_size, 1, sizeof(int32_t), file->file_ptr);
    if (err != sizeof(int32_t)) {
        return 0;
    }

    off_t offset = -4;
    fseek(file->file_ptr, offset, SEEK_CUR);

    return record_size;
}


/**
 * Sets the record position indicator for the sequential unformatted file
 * to the value pointed to by offset.
 * This function is similar to fseek(), but positions of records are used
 * rather than position of bytes.
 *
 * Position 'pos' to which offset is added can have one of the following values:
 * UNF_POS_BEGIN     from beginning of the file
 * UNF_POS_CURRENT   the current file position
 * UNF_POS_END       end of the file
 *
 * Offset is a number of records to shift the position relative to origin.
 *
 * Returns UNF_SUCCESS upon success, UNF_ERROR otherwise.
 */
int unf_seek(unf_file_t *file, unf_position_t pos, int offset)
{
    if (file == NULL) {
        return UNF_ERROR;
    }

    if (file->access != UNF_ACCESS_SEQUENTIAL) {
        return UNF_ERROR;
    }

    // set the position for further seeking
    if (pos == UNF_POS_BEGIN) {
        int err = fseek(file->file_ptr, 0, SEEK_SET);
        if (err != 0) {
            return UNF_ERROR;
        }
    }
    else if (pos == UNF_POS_END) {
        int err = fseek(file->file_ptr, 0, SEEK_END);
        if (err != 0) {
            return UNF_ERROR;
        }
    }
    else {
        // nothing to do
    }

    // shift position step-by-step
    if (offset == 0) {
        // nothing to do
    }
    else if (offset > 0) {
        while (offset) {
            int status = seek_forward(file);
            if (unf_eof(file) || status == UNF_ERROR) {
                return UNF_ERROR;
            }
            offset--;
        }
    }
    else { // offset < 0
        while (offset) {
            int status = seek_backward(file);
            if (status == UNF_ERROR) {
                return UNF_ERROR;
            }
            offset++;
        }
    }

    return UNF_SUCCESS;
}


/**
 * Positions the unformatted file to its initial point.
 * The next unf_read() call will read the first record.
 * For sequential files only.
 *
 * Returns UNF_SUCCESS upon success, UNF_ERROR otherwise.
 */
int unf_rewind(unf_file_t *file)
{
    return unf_seek(file, UNF_POS_BEGIN, 0);
}


/**
 * Positions the specified file to just before the preceding record.
 * For sequential files only.
 *
 * Returns UNF_SUCCESS upon success, UNF_ERROR otherwise.
 */
int unf_backspace(unf_file_t *file)
{
    return unf_seek(file, UNF_POS_CURRENT, -1);
}


/**
 * Skips the next record.
 * For sequential files only.
 *
 * Returns UNF_SUCCESS upon success, UNF_ERROR otherwise.
 */
int unf_skip(unf_file_t *file)
{
    return unf_seek(file, UNF_POS_CURRENT, 1);
}


/**
 * Checks if the end of the given unformatted file has been reached.
 *
 * Returns nonzero value if the end of the file has been reached, otherwise 0.
 */
int unf_eof(unf_file_t *file)
{
    return feof(file->file_ptr);
}


/**
 * Checks the given unformatted file for errors.
 *
 * Returns nonzero value if the file stream has errors occurred, 0 otherwise.
 */
int unf_error(unf_file_t *file)
{
    if (file->error_flag) {
        return UNF_ERROR;
    }

    return ferror(file->file_ptr) ? UNF_ERROR : UNF_SUCCESS;
}


/*
 *
 * auxiliary functions
 *
 */


static int try_read_bytes(unf_file_t *file, char *fmt, size_t *n_bytes_read, int *n_args_read, va_list ap)
{
    assert(file != NULL);

    *n_args_read = 0;
    *n_bytes_read = 0;

    char *p = fmt;
    while (*p) {
        // nothing to do, just skip
        if (*p == ',') {
            p++;
            continue;
        }

        int num_repeats = 1;
        int type_size;
        int data_type;
        if (fmt_get_type_size(&p, &data_type, &type_size, &num_repeats) == UNF_ERROR) {
            return UNF_ERROR;
        }

        /*
         * array or not?
         * array dimension can be specified by either integer-4 '[i4]' or integer-8 '[i8]' number
         */
        int is_array = 0;
        int array_dim_type = TYPE_INTEGER_4;
        if (strncmp(p, "[i4]", 4) == 0) {
            is_array = 1;
            array_dim_type = TYPE_INTEGER_4;
            p += 4;
        }
        else if (strncmp(p, "[i8]", 4) == 0) {
            is_array = 1;
            array_dim_type = TYPE_INTEGER_8;
            p += 4;
        }

        for (int i_repeat = 0; i_repeat < num_repeats; i_repeat++) {
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
                    return UNF_ERROR;
                }
                array_dim = *array_dim_ptr;
            }
            else if (is_array && array_dim_type == TYPE_INTEGER_8) {
                int64_t *array_dim_ptr = va_arg(ap, int64_t *);
                if (array_dim_ptr == NULL) {
                    errno = EINVAL; // wrong argument
                    return UNF_ERROR;
                }
                array_dim = *array_dim_ptr;
            }

            /*
             * read entry from the unformatted file
             * or skip it, if the data pointer is NULL
             */
            size_t n_bytes = array_dim * type_size;
            if (data_ptr != NULL) {
                size_t err = fread(data_ptr, 1, n_bytes, file->file_ptr);
                if (err != n_bytes) {
                    return UNF_ERROR;
                }
            }
            else {
                int err = fseek(file->file_ptr, (off_t) n_bytes, SEEK_CUR);
                if (err != 0) {
                    return UNF_ERROR;
                }
            }

            *n_bytes_read += n_bytes;
            *n_args_read += 1;
        }
    }

    return UNF_SUCCESS;
}


static int try_write_bytes(unf_file_t *file, char *fmt, size_t *n_bytes_written, int *n_args_written, va_list ap)
{
    assert(file != NULL);

    *n_args_written = 0;
    *n_bytes_written = 0;

    char *p = fmt;
    while (*p) {
        // nothing to do, just skip
        if (*p == ',') {
            p++;
            continue;
        }

        int num_repeats = 1;
        int type_size;
        int data_type;

        if (fmt_get_type_size(&p, &data_type, &type_size, &num_repeats) == UNF_ERROR) {
            return UNF_ERROR;
        }

        /*
         * array or not?
         * array dimension can be specified by either integer-4 '[i4]' or integer-8 '[i8]' number
         */
        int is_array = 0;
        int array_dim_type = TYPE_INTEGER_4;
        if (strncmp(p, "[i4]", 4) == 0) {
            is_array = 1;
            array_dim_type = TYPE_INTEGER_4;
            p += 4;
        }
        else if (strncmp(p, "[i8]", 4) == 0) {
            is_array = 1;
            array_dim_type = TYPE_INTEGER_8;
            p += 4;
        }

        for (int i_repeat = 0; i_repeat < num_repeats; i_repeat++) {
            /*
             * get pointer to data to be written
             */
            char data_char = 0;
            int8_t data_int1 = 0;
            int16_t data_int2 = 0;
            int32_t data_int4 = 0;
            int64_t data_int8 = 0;
            float data_real4 = 0;
            double data_real8 = 0;
            float _Complex data_complex4 = 0.0 + 0.0 * _Complex_I;
            double _Complex data_complex8 = 0.0 + 0.0 * _Complex_I;
            void *data_ptr = NULL;

            if (is_array || (data_type == TYPE_CHAR && type_size > 1)) {
                data_ptr = va_arg(ap, void *);
            }
            else {
                if (data_type == TYPE_CHAR) {
                    data_char = va_arg(ap, char);
                    data_ptr = &data_char;
                }
                else if (data_type == TYPE_INTEGER_1) {
                    data_int1 = va_arg(ap, int8_t);
                    data_ptr = &data_int1;
                }
                else if (data_type == TYPE_INTEGER_2) {
                    data_int2 = va_arg(ap, int16_t);
                    data_ptr = &data_int2;
                }
                else if (data_type == TYPE_INTEGER_4) {
                    data_int4 = va_arg(ap, int32_t);
                    data_ptr = &data_int4;
                }
                else if (data_type == TYPE_INTEGER_8) {
                    data_int8 = va_arg(ap, int64_t);
                    data_ptr = &data_int8;
                }
                else if (data_type == TYPE_REAL_4) {
                    data_real4 = va_arg(ap, float);
                    data_ptr = &data_real4;
                }
                else if (data_type == TYPE_REAL_8) {
                    data_real8 = va_arg(ap, double);
                    data_ptr = &data_real8;
                }
                else if (data_type == TYPE_COMPLEX_4) {
                    data_complex4 = va_arg(ap, float _Complex);
                    data_ptr = &data_complex4;
                }
                else if (data_type == TYPE_COMPLEX_8) {
                    data_complex8 = va_arg(ap, double _Complex);
                    data_ptr = &data_complex8;
                }
            }

            /*
             * if needed, get array length
             */
            size_t array_dim = 1;
            if (is_array && array_dim_type == TYPE_INTEGER_4) {
                array_dim = va_arg(ap, int32_t);
            }
            else if (is_array && array_dim_type == TYPE_INTEGER_8) {
                array_dim = va_arg(ap, int64_t);
            }

            /*
             * write entry to the unformatted file
             * or write zeros, if the data pointer is NULL
             */
            size_t n_bytes = array_dim * type_size;
            if (data_ptr != NULL) {
                size_t err = fwrite(data_ptr, 1, n_bytes, file->file_ptr);
                if (err != n_bytes) {
                    return UNF_ERROR;
                }
            }
            else {
                for (int i = 0; i < n_bytes; i++) {
                    fputc(0, file->file_ptr);
                }
                if (ferror(file->file_ptr)) {
                    return UNF_ERROR;
                }
            }

            *n_bytes_written += n_bytes;
            *n_args_written += 1;
        }
    }

    return UNF_SUCCESS;
}


// returns 0 if success, -1 if error
static int fmt_get_type_size(char **fmt, int *data_type, int *type_size, int *num_repeats)
{
    char *p = *fmt;

    // try to get number of repeats
    long parsed_int = strtol(p, &p, 10);
    if (parsed_int < 0 || errno == ERANGE) {
        return -1;
    }
    else if (parsed_int > 0) {
        *num_repeats = parsed_int;
    }
    else {
        *num_repeats = 1;
    }


    if (*p == 'c') {
        // 'c'
        *type_size = sizeof(char); // default: only one character
        *data_type = TYPE_CHAR;
        p++;

        // try to read length - is it a Fortran string?
        parsed_int = strtol(p, &p, 10);
        if (parsed_int < 0 || errno == ERANGE) {
            return -1;
        }
        else if (parsed_int > 0) {
            *type_size = parsed_int;
        }
    }
    else {
        // one of numeric types
        if (strncmp(p, "i1", 2) == 0) {
            // 'i1'
            *type_size = sizeof(int8_t);
            *data_type = TYPE_INTEGER_1;
        }
        else if (strncmp(p, "i2", 2) == 0) {
            // 'i2'
            *type_size = sizeof(int16_t);
            *data_type = TYPE_INTEGER_2;
        }
        else if (strncmp(p, "i4", 2) == 0) {
            // 'i4'
            *type_size = sizeof(int32_t);
            *data_type = TYPE_INTEGER_4;
        }
        else if (strncmp(p, "i8", 2) == 0) {
            // 'i8'
            *type_size = sizeof(int64_t);
            *data_type = TYPE_INTEGER_8;
        }
        else if (strncmp(p, "r4", 2) == 0) {
            // 'r4'
            *type_size = sizeof(float);
            *data_type = TYPE_REAL_4;
        }
        else if (strncmp(p, "r8", 2) == 0) {
            // 'r8'
            *type_size = sizeof(double);
            *data_type = TYPE_REAL_8;
        }
        else if (strncmp(p, "z4", 2) == 0) {
            // 'z4'
            *type_size = sizeof(float _Complex);
            *data_type = TYPE_COMPLEX_4;
        }
        else if (strncmp(p, "z8", 2) == 0) {
            // 'z8'
            *type_size = sizeof(double _Complex);
            *data_type = TYPE_COMPLEX_8;
        }
        else {
            errno = EINVAL;
            return UNF_ERROR;
        }

        p += 2;
    }

    *fmt = p;

    return UNF_SUCCESS;
}


/**
 * Seeks 1 record forward (auxiliary function).
 * For sequential files only.
 * Returns UNF_SUCCESS upon success, UNF_ERROR otherwise.
 */
static int seek_forward(unf_file_t *file)
{
    // read size of the record in bytes
    int32_t record_size = 0;
    size_t err = fread(&record_size, 1, sizeof(int32_t), file->file_ptr);
    if (err != sizeof(int32_t)) {
        return UNF_ERROR;
    }

    int status = fseek(file->file_ptr, record_size, SEEK_CUR);
    if (status != 0) {
        return UNF_ERROR;
    }

    int32_t record_size_2 = 0;
    err = fread(&record_size_2, 1, sizeof(int32_t), file->file_ptr);
    if (err != sizeof(int32_t)) {
        return UNF_ERROR;
    }

    if (record_size != record_size_2) {
        return UNF_ERROR;
    }

    return UNF_SUCCESS;
}


/**
 * Seeks 1 record backward (auxiliary function).
 * For sequential files only.
 * Returns UNF_SUCCESS upon success, UNF_ERROR otherwise.
 */
static int seek_backward(unf_file_t *file)
{
    // if we are at the beginning of the file,
    // the following fseek() call will simply return error
    // (without setting error flag however)

    // access the last record size
    off_t offset = -4;
    int err = fseek(file->file_ptr, offset, SEEK_CUR);
    if (err != 0) {
        return UNF_ERROR;
    }

    int32_t record_size = 0;
    size_t n_read = fread(&record_size, 1, sizeof(int32_t), file->file_ptr);
    if (n_read != sizeof(int32_t)) {
        return UNF_ERROR;
    }

    offset = -(2 * sizeof(int32_t) + record_size);
    err = fseek(file->file_ptr, offset, SEEK_CUR);
    if (err != 0) {
        return UNF_ERROR;
    }

    return UNF_SUCCESS;
}
