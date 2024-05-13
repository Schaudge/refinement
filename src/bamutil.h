#ifndef BAM_UTIL_H
#define BAM_UTIL_H

#include <cstdio>
#include <cstdlib>
#include <map>
#include "util.h"
#include "htslib/sam.h"

using namespace std;

class BamUtil {
public:
    BamUtil();
    ~BamUtil();

public:
    static string getQName(const bam1_t *b);
    static string getSeq(const bam1_t *b);
    static string getQual(const bam1_t *b);
    static string getCigar(const bam1_t *b);
    static char fourbits2base(uint8_t val);
    static uint8_t base2fourbits(char base);
    static void dump(bam1_t *b);
    static bool isPartOf(bam1_t *part, bam1_t *whole, bool isLeft);
    static void dumpHeader(bam_hdr_t* hdr);
    static bool isPrimary(bam1_t *b);
    static bool isProperPair(bam1_t *b);
    static int getRightRefPos(bam1_t *b);

    static inline unsigned int charToOp(const char cigar_char) {
        return cigar_char == 'M' ? BAM_CMATCH : (cigar_char == 'S' ? BAM_CSOFT_CLIP : (cigar_char == 'I' ? BAM_CINS : BAM_CDEL));
    }

    static inline tuple<bool, bool> readRealignQualify(bam1_t *b) {
        if (b->core.n_cigar < 2) {
            return make_tuple(false, false);
        } else {
            auto *cigar_data = (uint32_t *) bam_get_cigar(b);
            auto last_idx = b->core.n_cigar - 1;
            for (auto idx = 1; idx < last_idx; idx++) {
                if (bam_cigar_opchr(cigar_data[idx]) != 'M' && bam_cigar_opchr(cigar_data[last_idx]) != 'I' && bam_cigar_opchr(cigar_data[last_idx]) != 'D')
                    return make_tuple(false, false);
            }
            if (bam_cigar_opchr(cigar_data[0]) == 'M' && bam_cigar_opchr(cigar_data[last_idx]) == 'S') {
                return make_tuple(true, false);
            } else if (bam_cigar_opchr(cigar_data[0]) == 'S') {
                if (bam_cigar_opchr(cigar_data[last_idx]) == 'M') {
                    return make_tuple(false, true);
                } else if (bam_cigar_opchr(cigar_data[last_idx]) == 'S') {
                    if (bam_cigar_oplen(cigar_data[0]) + 2 < bam_cigar_oplen(cigar_data[last_idx]))
                        return make_tuple(true, false);
                    else if (bam_cigar_oplen(cigar_data[0]) > bam_cigar_oplen(cigar_data[last_idx]) + 2)
                        return make_tuple(false, true);
                }
            }
            return make_tuple(false, false);
        }
    }

    static inline tuple<string, hts_pos_t, string, string> parseSATag(uint8_t * sa_data_ptr) {
        string sa_tag{bam_aux2Z(sa_data_ptr)};
        size_t sa_idx_end = sa_tag.find_first_of(',');
        auto sa_ref = sa_tag.substr(0, sa_idx_end);
        auto sa_idx_start = sa_idx_end + 1;
        sa_idx_end = sa_tag.find_first_of(',', sa_idx_start);
        auto sa_pos = stol(sa_tag.substr(sa_idx_start, sa_idx_end - sa_idx_start)) - 1;
        sa_idx_start = sa_idx_end + 1;
        sa_idx_end = sa_tag.find_first_of(',', sa_idx_start);
        auto sa_strand = sa_tag.substr(sa_idx_start, sa_idx_end - sa_idx_start);
        sa_idx_start = sa_idx_end + 1;
        sa_idx_end = sa_tag.find_first_of(',', sa_idx_start);
        auto sa_cigar = sa_tag.substr(sa_idx_start, sa_idx_end - sa_idx_start);
        return make_tuple(sa_ref, sa_pos, sa_strand, sa_cigar);
    }

    static inline int outputAndEraseReads(std::map<size_t, std::vector<bam1_t*>> & position_reads_dict,
                                          samFile * & bam_fp, bam_hdr_t * & bam_header_t, size_t previous_pos) {
        bool written_reads_error = false;
        if (!position_reads_dict.empty()) {
            for (auto erase_iter = position_reads_dict.begin(); erase_iter != position_reads_dict.end();) {
                if (previous_pos < 1 || erase_iter->first <= previous_pos) {
                    for (auto &read: erase_iter->second) {
                        if (!written_reads_error)
                            if (sam_write1(bam_fp, bam_header_t, read) < 0)
                                written_reads_error = true;
                        bam_destroy1(read);
                    }
                    erase_iter = position_reads_dict.erase(erase_iter);
                } else break;
            }
        }
        return written_reads_error ? -1 : 0;
    }

    static string getZTag(bam1_t* b, const char tagName[2]) {
        auto tagData = (uint8_t*) bam_aux_get(b, tagName);
        return bam_aux2Z(tagData);
    }

    static void updateITag(bam1_t* b, const char tagName[2], const int newTagVal) {
        // update edit distance (NM) tag info
        // const char tagName[2] = {'N', 'M'};
        auto tagData = (uint8_t*) bam_aux_get(b, tagName);
        // int tagValue = bam_aux2i(tagData);
        char tagType = * tagData;
        if (tagType == 'C' && newTagVal >= 0 && newTagVal <= 255) {
            tagData[1] = newTagVal;
        }
    }

};

#endif
