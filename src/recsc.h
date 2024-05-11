#ifndef RECSC_H
#define RECSC_H

#include <stdio.h>
#include <stdlib.h>
#include <tuple>
#include <vector>
#include <string>
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "options.h"

using namespace std;

class Recsc {

    static const int kmer_size{5}; // k-mer size for prefix/postfix match to dual reference sequence
    static const int max_sc_mm{2}; // maximum soft-clip mismatch for correction

public:
    Recsc(Options *opt);
    ~Recsc();

    void correct2();

private:
    Options * options;
    samFile * in;
    hts_idx_t * bam_idx;
    bam_hdr_t * bam_header;
    hts_itr_t * bam_iter;
    samFile * output;
    htsThreadPool p = {NULL, 0};
    std::vector<std::tuple<std::string, std::string, std::string>> regions;
    std::vector<std::tuple<int32_t, hts_pos_t, hts_pos_t>> numerical_regions;
    const char tagName[2] = {'M', 'D'};
};

#endif
