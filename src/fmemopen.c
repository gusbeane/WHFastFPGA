/**
 * @file    fmemopen.c
 * @brief   Implementation of fmemopen for old MacOSX and Windows.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section LICENSE
 * Copyright (c) 2023 Hanno Rein
 * 
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Most of the code below is taken from glibc. 
 * See https://github.com/lattera/glibc/blob/master/libio/fmemopen.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include "fmemopen.h"

#ifdef __MACH__

// This code here is a workaround for old MacOS versions which do 
// not come with fmemopen. (Note that as a result of fmemopen 
// missing, conda-forge builds fail).
// The following code has been slightly modified from
// https://raw.githubusercontent.com/Cibiv/PDA/master/fmemopen.c

struct fmem {
    size_t pos;
    size_t size;
    char *buffer;
};

static int readfn(void *handler, char *buf, int size) {
    struct fmem *mem = handler;
    size_t available = mem->size - mem->pos;

    if (size > available) {
        size = (int)available;
    }
    memcpy(buf, mem->buffer + mem->pos, sizeof(char) * size);
    mem->pos += size;

    return size;
}

static int writefn(void *handler, const char *buf, int size) {
    struct fmem *mem = handler;
    size_t available = mem->size - mem->pos;

    if (size > available) {
        size = (int)available;
    }
    memcpy(mem->buffer + mem->pos, buf, sizeof(char) * size);
    mem->pos += size;

    return size;
}

static fpos_t seekfn(void *handler, fpos_t offset, int whence) {
    size_t pos;
    struct fmem *mem = handler;

    switch (whence) {
        case SEEK_SET: pos = offset; break;
        case SEEK_CUR: pos = mem->pos + offset; break;
        case SEEK_END: pos = mem->size + offset; break;
        default: return -1;
    }

    if (pos > mem->size) {
        return -1;
    }

    mem->pos = pos;
    return (fpos_t)pos;
}

static int closefn(void *handler) {
    free(handler);
    return 0;
}

FILE *reb_fmemopen(void *buf, size_t size, const char *mode) {
    // This data is released on fclose.
    struct fmem* mem = (struct fmem *) malloc(sizeof(struct fmem));

    // Zero-out the structure.
    memset(mem, 0, sizeof(struct fmem));

    mem->size = size;
    mem->buffer = buf;

    // funopen's man page: https://developer.apple.com/library/mac/#documentation/Darwin/Reference/ManPages/man3/funopen.3.html
    return funopen(mem, readfn, writefn, seekfn, closefn);
}

#else
// Just an alias for everyone else
FILE *reb_fmemopen(void *buf, size_t len, const char *type){
    return fmemopen(buf, len, type);
}

#endif
