#include <algorithm>
#include <sstream>
#include <ranges>
#include <cmath>
#include <regex>
#include <unordered_map>
#include "recsc.h"
#include "sam_ext.h"
#include "bamutil.h"
#include "bioio.hpp"

Recsc::Recsc(Options *opt) {
    options = opt;
    in = nullptr;
    bam_idx = nullptr;
    bam_header = nullptr;
    bam_iter = nullptr;
    output = nullptr;
}

Recsc::~Recsc(){
    if (in != nullptr) {
       if (sam_close(in) < 0) {
           cerr << "ERROR: failed to close " << options->input << endl;
           exit(-1);
       }
    }
    if (bam_idx != nullptr) {
        hts_idx_destroy(bam_idx);
    }
    if (bam_header != nullptr) {
        bam_hdr_destroy(bam_header);
    }
    if (bam_iter != nullptr) {
        hts_itr_destroy(bam_iter);
    }
    if (output != nullptr) {
        if (sam_close(output) < 0) {
            cerr << "ERROR: failed to close " << options->output << endl;
            exit(-1);
        }
    }
    if (p.pool)
        hts_tpool_destroy(p.pool);
}

void Recsc::correct2() {
    in = sam_open(options->input.c_str(), "r");
    if (!in) {
        cerr << "ERROR: failed to open bam file " << options->input << endl;
        exit(-1);
    }
    bam_idx = sam_index_load(in, options->input.c_str());
    if (!bam_idx) {
        cerr << "ERROR: failed to open bam index file" << endl;
        exit(-1);
    }
    if (ends_with(options->output, "bam"))
        output = sam_open(options->output.c_str(), "wb");
    else
        output = sam_open(options->output.c_str(), "w");
    if (!output) {
        cerr << "ERROR: failed to open output " << options->output << endl;
        exit(-1);
    }
    bam_header = sam_hdr_read(in);
    if (bam_header == nullptr || bam_header->n_targets == 0) {
        cerr << "ERROR: this SAM file has no header " << options->input << endl;
        exit(-1);
    }
    if (sam_hdr_write(output, bam_header) < 0) {
        cerr << "failed to write header" << endl;
        exit(-1);
    }
    if (ends_with(options->output, "bam") && sam_idx_init(output, bam_header, 0, options->output_idx.c_str()) < 0) {
        cerr << "failed to auto create the bam index" << endl;
        exit(-1);
    }
    if (options->nthreads > 1) {
        if (!(p.pool = hts_tpool_init(options->nthreads))) {
            cerr << "Error for creating thread pool" << endl;
            exit(-1);
        }
        hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
        hts_set_opt(output, HTS_OPT_THREAD_POOL, &p);
    }

    ifstream fa_index {options->refFile + ".fai"};
    if (!fa_index) {
        cerr << "ERROR: failed to open reference index file (.fai) " << options->output << endl;
        exit(-1);
    }
    const bioio::FastaIndex ref_index = bioio::read_fasta_index(fa_index);

    bool error_exit_flag{false};
    // read bed/interval file into merged region vector for loop
    string line, chr, region_start, region_end;
    ifstream region_input(options->regionFile);
    while (getline(region_input, line)) {
        if (line[0] == '@' || line[0] == '#') continue;
        std::istringstream isl(line);
        isl >> chr >> region_start >> region_end;
        isl.ignore(16);
        if (!numerical_regions.empty()) {
            int32_t current_tid = sam_hdr_name2tid(bam_header, chr.c_str());
            size_t last_region_idx = numerical_regions.size() - 1;
            int32_t previous_tid = get<0>(numerical_regions[last_region_idx]);
            int64_t previous_end = get<2>(numerical_regions[last_region_idx]);
            if (current_tid == previous_tid && stoi(region_start) - previous_end < 10000) {
                hts_pos_t pre_start = get<1>(numerical_regions[numerical_regions.size() - 1]);
                numerical_regions[numerical_regions.size() - 1] = make_tuple(current_tid, pre_start, stol(region_end));
            } else
                numerical_regions.emplace_back(current_tid, stol(region_start), stol(region_end));
        } else
            numerical_regions.emplace_back(sam_hdr_name2tid(bam_header, chr.c_str()), stol(region_start), stol(region_end));
    }
    region_input.close();

    long start{0}, previous{0}, riter_idx{0}, pre_iter_idx{-1};
    size_t ref_genome_start{0}, genome_ref_len{300};
    bool font_matched, back_matched, capture_region;
    unsigned int range = options->refRange * 2;
    string genome_sequence, realigned_ref_seq, contig_name;
    bam1_t *b = bam_init1();
    auto region_size = numerical_regions.size();
    std::unordered_map<string, uint16_t> has_realigned_reads;
    std::map<size_t, std::vector<bam1_t*>> back_moved_reads;
    while (sam_read1(in, bam_header, b) >= 0) {
        while (riter_idx < region_size && b->core.tid > get<0>(numerical_regions[riter_idx])) {
            if (BamUtil::outputAndEraseReads(back_moved_reads, output, bam_header, 0) != 0) {
                error_exit_flag = true;
                goto bam_destroy_for_free;
            }
            ++riter_idx;
        }
        while (riter_idx < region_size && b->core.tid == get<0>(numerical_regions[riter_idx]) &&
                b->core.pos > get<2>(numerical_regions[riter_idx])) ++riter_idx;
        if (riter_idx > pre_iter_idx && riter_idx < region_size) {
            // generate the reverse complement sequence from genome (.fa/.fasta) file for each region only once!
            auto riter_region_start = get<1>(numerical_regions[riter_idx]);
            ref_genome_start = riter_region_start > 150 ? riter_region_start - 150 : 1;
            genome_ref_len = get<2>(numerical_regions[riter_idx]) - ref_genome_start + 150;
            contig_name = sam_hdr_tid2name(bam_header, get<0>(numerical_regions[riter_idx]));
            genome_sequence = bioio::read_fasta_contig(options->refFile, ref_index.at(contig_name), ref_genome_start, genome_ref_len);
            reverse(genome_sequence.begin(), genome_sequence.end());
            transform(genome_sequence.begin(), genome_sequence.end(), genome_sequence.begin(), complement);
            pre_iter_idx = riter_idx;
        }
        capture_region = riter_idx < region_size && b->core.tid == get<0>(numerical_regions[riter_idx]) &&
                         b->core.pos > get<1>(numerical_regions[riter_idx]) &&
                         b->core.pos < get<2>(numerical_regions[riter_idx]);
        // previous genome realignment reads output first
        if (BamUtil::outputAndEraseReads(back_moved_reads, output, bam_header, b->core.pos) != 0) {
            error_exit_flag = true;
            goto bam_destroy_for_free;
        }
        unsigned int read_len = b->core.l_qseq;
        auto *tagData = (uint8_t *) bam_aux_get(b, tagName);
        if (capture_region && (read_len > options->endPosition) && (b->core.mtid >= 0) && tagData != nullptr) {
            // determine this read whether a local chimeric alignment that maybe a potential middle insert/deletion?
            string read_name = bam_get_qname(b);
            if (has_realigned_reads.find(read_name) != has_realigned_reads.end()) {
                auto previous_flag = has_realigned_reads[read_name];
                if (previous_flag & b->core.flag) continue;
            }
            // case 1: font/back end soft-clip read
            auto *cigar_data = (uint32_t *) bam_get_cigar(b);
            char first_operator = bam_cigar_opchr(cigar_data[0]);
            uint32_t first_op_len = bam_cigar_oplen(cigar_data[0]);
            char last_operator = bam_cigar_opchr(cigar_data[b->core.n_cigar - 1]);
            uint32_t last_op_len = bam_cigar_oplen(cigar_data[b->core.n_cigar - 1]);
            auto [is_match_softclip_reads, is_softclip_match_reads] = BamUtil::readRealignQualify(b);
            auto* sa_tag_ptr = (uint8_t*) bam_aux_get(b, "SA");
            if (sa_tag_ptr != nullptr && (is_match_softclip_reads || is_softclip_match_reads)) {
                auto read_strand = (b->core.flag & BAM_FREVERSE) == 0 ? "+" : "-";
                auto [sa_ref, sa_pos, sa_strand, sa_cigar_string] = BamUtil::parseSATag(sa_tag_ptr);   // Structured binding
                if (sa_ref == contig_name && sa_strand == read_strand) {
                    auto [is_clip_match_sa, is_match_clip_sa, padding_len, sa_cigar_count] = clipMatchConvert(sa_cigar_string);
                    if (is_match_softclip_reads && is_clip_match_sa) {
                        auto read_end_pos = bam_endpos(b);
                        auto front_read_len = read_len - padding_len;
                        if (likely(b->core.pos <= sa_pos && sa_pos < (read_end_pos + options->maxdel) && (read_end_pos + 17 < sa_pos || sa_pos < read_end_pos + 10))) {
                            auto new_cigar_len = b->core.n_cigar + sa_cigar_count;
                            auto* new_cigar = (uint32_t*) malloc(new_cigar_len * 4);
                            if (unlikely(!new_cigar)) {
                                error_exit_flag = true;
                                goto bam_destroy_for_free;
                            }
                            auto latest_genome_pos {b->core.pos};
                            size_t preceding_read_len{0}, new_cigar_idx{0};
                            bool indel_has_padded{false};
                            for (; new_cigar_idx < b->core.n_cigar - 1; ++new_cigar_idx) {    // skip the last 'S' operator
                                auto current_op_len = bam_cigar_oplen(cigar_data[new_cigar_idx]);
                                auto current_op_chr = bam_cigar_opchr(cigar_data[new_cigar_idx]);
                                if (latest_genome_pos + current_op_len >= sa_pos || preceding_read_len + current_op_len >= front_read_len) {
                                    auto under_padding_seq_len = front_read_len - preceding_read_len;
                                    auto remained_genome_distance  = sa_pos - latest_genome_pos;
                                    auto padded_len = under_padding_seq_len > remained_genome_distance ? remained_genome_distance : under_padding_seq_len;
                                    if (latest_genome_pos + padded_len >= sa_pos) {
                                        if (padded_len > 0) {
                                            new_cigar[new_cigar_idx] = (padded_len << BAM_CIGAR_SHIFT) | (cigar_data[new_cigar_idx] & BAM_CIGAR_MASK);
                                            ++new_cigar_idx;
                                        }
                                        new_cigar[new_cigar_idx] = bam_cigar_gen(under_padding_seq_len - padded_len, BAM_CINS);
                                    } else {
                                        new_cigar[new_cigar_idx] = (padded_len << BAM_CIGAR_SHIFT) | (cigar_data[new_cigar_idx] & BAM_CIGAR_MASK);
                                        new_cigar[++new_cigar_idx] = bam_cigar_gen(remained_genome_distance - padded_len, BAM_CDEL);
                                    }
                                    indel_has_padded = true;
                                    break;
                                } else
                                    new_cigar[new_cigar_idx] = cigar_data[new_cigar_idx];
                                if (current_op_chr == 'M') {
                                    latest_genome_pos += current_op_len;
                                    preceding_read_len += current_op_len;
                                } else if (current_op_chr == 'S' || current_op_chr == 'I')
                                    preceding_read_len += current_op_len;
                                else if (current_op_chr == 'D')
                                    latest_genome_pos += current_op_len;
                            }
                            size_t mismatch_seq_len{0};
                            if (!indel_has_padded) {   // complex indels case, insertions + deletions
                                auto padded_len = front_read_len - preceding_read_len;
                                if (sa_pos < latest_genome_pos + 11 && padded_len > 11) {         // long insertion
                                    mismatch_seq_len = sa_pos - latest_genome_pos;
                                    new_cigar[new_cigar_idx] = bam_cigar_gen(padded_len - mismatch_seq_len, BAM_CINS);
                                } else if (sa_pos > latest_genome_pos + 17 && padded_len < 17) {  // long deletion
                                    mismatch_seq_len = padded_len;
                                    new_cigar[new_cigar_idx] = bam_cigar_gen(sa_pos - latest_genome_pos - mismatch_seq_len, BAM_CDEL);
                                } else {                                                          // should never occur
                                    BamUtil::dump(b);
                                    free(new_cigar);
                                    goto bam_write;
                                }
                            }
                            size_t back_padding_start{0};
                            bool first_clip_skipped{false};
                            for (auto pi = 1; pi < sa_cigar_string.size(); ++pi) {
                                auto sa_op_char = sa_cigar_string[pi];
                                if (sa_op_char == 'S' && !first_clip_skipped) {
                                    first_clip_skipped = true;
                                    back_padding_start = pi + 1;
                                } else if (unlikely(first_clip_skipped && mismatch_seq_len > 0 && sa_op_char == 'M')) {
                                    auto op_len = stol(sa_cigar_string.substr(back_padding_start, pi - back_padding_start));
                                    back_padding_start = pi + 1;
                                    new_cigar[++new_cigar_idx] = bam_cigar_gen(op_len + mismatch_seq_len, BAM_CMATCH);
                                    mismatch_seq_len = 0;
                                } else if (first_clip_skipped && isalpha(sa_op_char)) {
                                    auto op_len = stol(sa_cigar_string.substr(back_padding_start, pi - back_padding_start));
                                    back_padding_start = pi + 1;
                                    new_cigar[++new_cigar_idx] = bam_cigar_gen(op_len, BamUtil::charToOp(sa_op_char));
                                }
                            }
                            bam_aux_del(b, sa_tag_ptr);
                            if (replace_cigar(b, new_cigar_idx + 1, new_cigar) != 0) {
                                free(new_cigar);
                                error_exit_flag = true;
                                goto bam_destroy_for_free;
                            }
                            if (has_realigned_reads.find(read_name) == has_realigned_reads.end()) {
                                if (b->core.flag & BAM_FSECONDARY)
                                    b->core.flag &= (~ BAM_FSECONDARY);
                                has_realigned_reads[read_name] = b->core.flag & (BAM_FREAD1 | BAM_FREAD2);
                            } else
                                has_realigned_reads[read_name] = BAM_FREAD1 | BAM_FREAD2;
                            free(new_cigar);
                            goto bam_write;
                        }
                    } else if (is_softclip_match_reads && is_match_clip_sa) {       // special case for duplication (insertion)
                        vector<uint32_t> reverse_cigar;
                        size_t sa_idx_start{0}, consumed_read_len{0}, ni{0};
                        auto new_start_pos = sa_pos;
                        auto read_end_pos = bam_endpos(b);
                        if (likely(b->core.pos <= sa_pos && sa_pos + padding_len < read_end_pos)) {
                            auto new_cigar_len = b->core.n_cigar + sa_cigar_count;
                            auto* new_cigar = (uint32_t*) malloc(new_cigar_len * 4);
                            if (unlikely(!new_cigar)) {
                                error_exit_flag = true;
                                goto bam_destroy_for_free;
                            }
                            for (auto ci = 1; ci < sa_cigar_string.size() - 2; ++ci) {   // skip the last 'S' operator
                                char op_char = sa_cigar_string[ci];
                                if (isalpha(op_char)) {
                                    auto op_len = stol(sa_cigar_string.substr(sa_idx_start, ci - sa_idx_start));
                                    new_cigar[ni++] = bam_cigar_gen(op_len, BamUtil::charToOp(op_char));
                                    if (op_char != 'D')
                                        consumed_read_len += op_len;
                                    if (op_char == 'M' || op_char == 'D')
                                        sa_pos += op_len;
                                    sa_idx_start = ci + 1;
                                }
                            }
                            for (auto i = b->core.n_cigar - 1; i > 0; --i) {    // skip the first 'S' operator for raw CIGAR
                                auto reverse_op_char = bam_cigar_opchr(cigar_data[i]);
                                auto reverse_op_len = bam_cigar_oplen(cigar_data[i]);
                                if (reverse_op_char != 'S' && read_end_pos - reverse_op_len <= sa_pos) {
                                    b->core.pos = new_start_pos;
                                    bam_aux_del(b, sa_tag_ptr);
                                    auto reverse_extend_len = read_end_pos - sa_pos;
                                    auto insert_len = reverse_op_char == 'D' ? read_len - consumed_read_len : read_len - consumed_read_len - reverse_extend_len;
                                    new_cigar[ni++] = bam_cigar_gen(insert_len, BAM_CINS);
                                    new_cigar[ni++] = (reverse_extend_len << BAM_CIGAR_SHIFT) | (cigar_data[i] & BAM_CIGAR_MASK);
                                    if (!reverse_cigar.empty())
                                        for (auto c : std::views::reverse(reverse_cigar))
                                            new_cigar[ni++] = c;
                                    if (replace_cigar(b, ni, new_cigar) != 0) {
                                        free(new_cigar);
                                        error_exit_flag = true;
                                        goto bam_destroy_for_free;
                                    }
                                    if (has_realigned_reads.find(read_name) == has_realigned_reads.end()) {
                                        if (b->core.flag & BAM_FSECONDARY)
                                            b->core.flag &= (~ BAM_FSECONDARY);
                                        has_realigned_reads[read_name] = b->core.flag & (BAM_FREAD1 | BAM_FREAD2);
                                    } else
                                        has_realigned_reads[read_name] = BAM_FREAD1 | BAM_FREAD2;
                                    back_moved_reads[new_start_pos].emplace_back(bam_dup1(b));
                                    break;
                                }
                                reverse_cigar.emplace_back(cigar_data[i]);
                                if (reverse_op_char != 'D')
                                    consumed_read_len += reverse_op_len;
                                if (reverse_op_char == 'M' || reverse_op_char == 'D')
                                    read_end_pos -= reverse_op_len;
                            }
                            free(new_cigar);
                            continue;
                        } else
                            BamUtil::dump(b);
                    }
                }
            }
            // case 2: read end sequence with mismatched base
            string mismatch_del_string{bam_aux2Z(tagData)};
            size_t last_char_index = mismatch_del_string.size() - 1;
            if (last_char_index < 3) {   // < 10 bp end mismatch !
                font_matched = back_matched = false;
            } else if (isdigit(mismatch_del_string[last_char_index]) && isdigit(mismatch_del_string[last_char_index-1]) &&
                       (mismatch_del_string[last_char_index-2] == 'A' || mismatch_del_string[last_char_index-2] == 'C' ||
                        mismatch_del_string[last_char_index-2] == 'G' || mismatch_del_string[last_char_index-2] == 'T')) {
                font_matched = false;
                back_matched = stoi(mismatch_del_string.substr(last_char_index - 1)) <= options->endPosition;
            } else if (isdigit(mismatch_del_string[last_char_index]) &&
                       (mismatch_del_string[last_char_index-1] == 'A' || mismatch_del_string[last_char_index-1] == 'C' ||
                        mismatch_del_string[last_char_index-1] == 'G' || mismatch_del_string[last_char_index-1] == 'T')) {
                font_matched = false; back_matched = true;
            } else if (isdigit(mismatch_del_string[0]) && isdigit(mismatch_del_string[1]) &&
                       (mismatch_del_string[2] == 'A' || mismatch_del_string[2] == 'C' ||
                        mismatch_del_string[2] == 'G' || mismatch_del_string[2] == 'T')) {
                font_matched = stoi(mismatch_del_string.substr(0, 2)) <= options->endPosition;
                back_matched = false;
            } else if (isdigit(mismatch_del_string[0]) && (mismatch_del_string[1] == 'A' ||
                mismatch_del_string[1] == 'C' || mismatch_del_string[1] == 'G' || mismatch_del_string[1] == 'T')) {
                font_matched = true; back_matched = false;
            } else {
                font_matched = false; back_matched = false;
            }
            // TODO: The LL and RR orientation paired reads to be corrected
            if (first_operator == 'S' || last_operator == 'S' || back_matched || font_matched) {
                // we prefer the soft-clip alignment first, then for the back/font end mismatch!
                start = b->core.pos > options->refRange ? b->core.pos - options->refRange : 0;
                bool back_end_crt = last_operator == 'S' || (first_operator != 'S' && back_matched);
                if (back_end_crt) {
                    for (int i = 0; i < sizeof(cigar_data); ++i) {
                        char cigar_operator = bam_cigar_opchr(cigar_data[i]);
                        if (cigar_operator == 'M')
                            start += (bam_cigar_oplen(cigar_data[i]));
                    }
                }
                // the calculation of the start position for reverse realigned reference sequence
                if (abs(start - previous) > 3) {
                    if (start + range <= ref_genome_start || start >= ref_genome_start + genome_ref_len) {
                        realigned_ref_seq = "";
                    } else if (start <= ref_genome_start) {
                        size_t reverse_pos = start + range >= ref_genome_start + genome_ref_len ? 0 :
                                             ref_genome_start + genome_ref_len - start - range;
                        realigned_ref_seq = genome_sequence.substr(reverse_pos);
                    } else if (start + range > ref_genome_start + genome_ref_len) {
                        realigned_ref_seq = genome_sequence.substr(0, ref_genome_start + genome_ref_len - start);
                    } else {
                        realigned_ref_seq = genome_sequence.substr(ref_genome_start + genome_ref_len - start - range, range);
                    }
                    previous = start;
                }
                // cout << realigned_ref_seq << endl;

                string read_seq = BamUtil::getSeq(b);
                unsigned int read_start_pos = back_end_crt ? read_len - kmer_size : 1;
                string kmer_seq = read_seq.substr(read_start_pos, kmer_size);
                // adjust the dynamic k-mer size for non-unique matched reference sequences
                size_t match_count = substr_count(realigned_ref_seq, kmer_seq);
                if (match_count > 1 && back_end_crt) {
                    while (match_count > 1 && read_start_pos > kmer_size) {
                        kmer_seq = read_seq.substr(--read_start_pos);
                        match_count = substr_count(realigned_ref_seq, kmer_seq);
                    }
                } else if (match_count > 1) {
                    size_t extended_kmer_size{kmer_size};
                    while (match_count > 1 && extended_kmer_size * 2 < read_len) {
                        kmer_seq = read_seq.substr(1, ++extended_kmer_size);
                        match_count = substr_count(realigned_ref_seq, kmer_seq);
                    }
                }
                // we use a mask quality (to 0x1) method to hide the enzyme cutting padding sequence!
                if (match_count > 0) {
                    size_t matched_ref_pos = realigned_ref_seq.find(kmer_seq);
                    if (last_operator == 'S') {
                        int sc_mismatch{0};
                        while (sc_mismatch <= max_sc_mm && 0 < read_start_pos && read_start_pos < read_len
                               && matched_ref_pos < range) {
                            if (realigned_ref_seq[matched_ref_pos--] != read_seq[read_start_pos--])
                                sc_mismatch += (read_start_pos + last_op_len > read_len ? 1 : 3);
                        }
                    } else if (first_operator == 'S') {
                        int sc_mismatch{0};
                        while (sc_mismatch <= max_sc_mm && 0 < read_start_pos && read_start_pos < read_len
                               && matched_ref_pos < range) {
                            if (realigned_ref_seq[matched_ref_pos++] != read_seq[read_start_pos++])
                                sc_mismatch += (read_start_pos < first_op_len ? 1 : 3);
                        }
                    } else if (back_matched) {
                        while (0 < read_start_pos && read_start_pos < read_len && matched_ref_pos < range
                               && realigned_ref_seq[matched_ref_pos] == read_seq[read_start_pos]) {
                            --matched_ref_pos;
                            --read_start_pos;
                        }
                    } else {
                        while (0 < read_start_pos && read_start_pos < read_len && matched_ref_pos < range
                               && realigned_ref_seq[matched_ref_pos] == read_seq[read_start_pos]) {
                            ++matched_ref_pos;
                            ++read_start_pos;
                        }
                    }
                    uint8_t * out_quality = bam_get_qual(b);
                    if (back_end_crt)
                        for (auto i = read_start_pos; i < read_len; ++i) out_quality[i] = 0x1;
                    else
                        for (auto i = 0; i < read_start_pos; ++i) out_quality[i] = 0x1;
                }
            }
        }
bam_write:
        if (sam_write1(output, bam_header, b) < 0) {
            error_exit_flag = true;
            goto bam_destroy_for_free;
        }
    }
bam_destroy_for_free:
    bam_destroy1(b);
    if (BamUtil::outputAndEraseReads(back_moved_reads, output, bam_header, 0) != 0)
        error_exit("bam/sam writing failed, free and exiting ...");
    if (ends_with(options->output, "bam") && sam_idx_save(output) < 0)
        error_exit("bam index writing failed, exiting ...");
    if (error_exit_flag)
        error_exit("some error occurred in memory allocation or writing bam/sam, exiting ...");
}
