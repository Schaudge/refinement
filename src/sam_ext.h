/* the extension of sam file was copied or modified from samtools software.
 * padding.c -- depad subcommand.

    Copyright (C) 2011, 2012 Broad Institute.
    Copyright (C) 2014-2016, 2019-2020 Genome Research Ltd.
    Portions copyright (C) 2012, 2013 Peter Cock, The James Hutton Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <htslib/kstring.h>
#include <htslib/sam.h>

static int replace_cigar(bam1_t *b, uint32_t n, uint32_t *cigar) {
    size_t diff = 0;
    if (n != b->core.n_cigar) {
        size_t o = b->core.l_qname + b->core.n_cigar * 4;
        if (n > b->core.n_cigar) {
            diff = (n - b->core.n_cigar) * 4;
            if ((INT_MAX - b->l_data)/4 < (n - b->core.n_cigar)) {
                fprintf(stderr, "[depad] ERROR: BAM record too big\n");
                return -1;
            }
            if (b->l_data + diff > b->m_data) {
                b->m_data = b->l_data + diff;
                kroundup32(b->m_data);
                auto *tmp = (uint8_t*)realloc(b->data, b->m_data);
                if (!tmp) {
                    fprintf(stderr, "[depad] ERROR: Memory allocation failure.\n");
                    return -1;
                }
                b->data = tmp;
            }
        } else {
            diff = -(int)((b->core.n_cigar - n) * 4);
        }
        memmove(b->data + b->core.l_qname + n * 4, b->data + o, b->l_data - o);
        b->core.n_cigar = n;
    }

    memcpy(b->data + b->core.l_qname, cigar, n * 4);
    b->l_data += diff;

    return 0;
}

#define write_cigar(_c, _n, _m, _v) do { \
        if (_n == _m) { \
            _m = _m? _m<<1 : 4; \
            _c = (uint32_t*)realloc(_c, _m * 4); \
            if (!(_c)) { \
                fprintf(stderr, "[depad] ERROR: Memory allocation failure.\n"); \
                return -1; \
            } \
        } \
        _c[_n++] = (_v); \
    } while (0)

// Functions to swap partial serial sequences and qualities in a bam record, written by Schaudge King!
static inline int swap_front_serial_seq(uint8_t * seq_ptr, uint32_t break_pos, uint32_t len) {
    auto block_size = (len + 1) >> 1;
    auto * swap_ptr = (uint8_t*) malloc(block_size);
    if (!swap_ptr) {
        fprintf(stderr, "Memory allocation for sequences swap space failure..\n");
        return -1;
    }
    memmove(swap_ptr, seq_ptr, block_size);
    auto offset = (break_pos + 1) >> 1;
    if (len % 2 == 0) {
        memmove(seq_ptr, seq_ptr + block_size, offset);
        if (break_pos % 2 == 0)
            memmove(seq_ptr + offset, swap_ptr, block_size);
        else {
            seq_ptr[offset - 1] = (seq_ptr[offset - 1] & 0xf0) | (swap_ptr[0] >> 4);
            for (auto i = 0; i < block_size - 1; i++)
                seq_ptr[offset + i] = (swap_ptr[i] << 4) | (swap_ptr[i + 1] >> 4);
            auto last_pos = offset + block_size - 1;
            seq_ptr[last_pos] = (swap_ptr[block_size - 1] << 4) | (seq_ptr[last_pos] & 0x0f);
        }
    } else {
        for (auto i = 0; i < offset; i++)
            seq_ptr[i] = (seq_ptr[block_size + i - 1] << 4) | (seq_ptr[block_size + i] >> 4);
        if (break_pos % 2 == 0) {
            memmove(seq_ptr + offset, swap_ptr, block_size - 1);
            auto last_pos = offset + block_size - 1;
            seq_ptr[last_pos] = (swap_ptr[block_size - 1] & 0xf0) | (seq_ptr[last_pos] & 0x0f);
        } else {
            seq_ptr[offset - 1] = (seq_ptr[offset - 1] & 0xf0) | (swap_ptr[0] >> 4);
            for (auto i = 0; i < block_size - 1; i++)
                seq_ptr[offset + i] = (swap_ptr[i] << 4) | (swap_ptr[i + 1] >> 4);
        }
    }
    free(swap_ptr);
    return 0;
}
