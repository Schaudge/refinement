#include "bamutil.h"
#include <sstream>

BamUtil::BamUtil(){
}

BamUtil::~BamUtil(){
}

void BamUtil::dump(bam1_t *b) {
    cerr << b->core.tid << ":" << b->core.pos << ", M:" << b->core.mtid << ":" << b->core.mpos << " TLEN:" << b->core.isize << " ID:" << b->id << endl;
    cerr << getQName(b) << " " << getCigar(b) << endl;
    cerr << getSeq(b) << endl;
    cerr << getQual(b) << endl;
}

string BamUtil::getQName(const bam1_t *b) {
    return bam_get_qname(b);
}

string BamUtil::getQual(const bam1_t *b) {
    uint8_t *data = bam_get_qual(b);
    int len = b->core.l_qseq;
    string s(len, '\0');
    for (int i=0; i < len; i++) {
        s[i] = (char)(data[i] + 33);
    }
    return s;
}

string BamUtil::getSeq(const bam1_t *b) {
    uint8_t *data = bam_get_seq(b);
    int len = b->core.l_qseq;
    string s(len, '\0');
    for(int i=0; i<len; i++) {
        char base;
        if(i%2 == 1)
            base = fourbits2base(data[i/2] & 0xF);
        else
            base = fourbits2base((data[i/2]>>4) & 0xF);
        s[i] = base;
    }
    return s;
}

//Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G, 8 for T and 15 for N.
char BamUtil::fourbits2base(uint8_t val) {
    switch(val) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            cerr << "ERROR: Wrong base with value "<< (int)val << endl ;
            return 'N';
    }
}

uint8_t BamUtil::base2fourbits(char base) {
    switch(base) {
        case 'A':
            return 1;
        case 'C':
            return 2;
        case 'G':
            return 4;
        case 'T':
            return 8;
        case 'N':
            return 15;
        default:
            cerr << "ERROR: Wrong base "<< base << endl ;
            return 15;
    }
}

/*
@discussion In the CIGAR array, each element is a 32-bit integer. The
 lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
 length of a CIGAR.
*/

string BamUtil::getCigar(const bam1_t *b) {
    auto *data = (uint32_t *) bam_get_cigar(b);
    int cigarNum = b->core.n_cigar;
    stringstream ss;
    for (auto i = 0; i < cigarNum; i++) {
        uint32_t val = data[i];
        char op = bam_cigar_opchr(val);
        uint32_t len = bam_cigar_oplen(val);
        ss << len << op;
    }
    return ss.str();
}

bool BamUtil::isPartOf(bam1_t *part, bam1_t *whole, bool isLeft) {
    auto *cigarPart = (uint32_t *) bam_get_cigar(part);
    int cigarNumPart = part->core.n_cigar;
    auto *cigarWhole = (uint32_t *) bam_get_cigar(whole);
    int cigarNumWhole = whole->core.n_cigar;

    if (cigarNumWhole < cigarNumPart)
        return false;

    for (int i=0; i < cigarNumPart; i++) {
        uint32_t valPart = cigarPart[i];
        // if it's right aligned, backward
        if (!isLeft)
            valPart = cigarPart[cigarNumPart - i - 1];
        char opPart = bam_cigar_op(valPart);
        uint32_t lenPart = bam_cigar_oplen(valPart);

        uint32_t valWhole = cigarWhole[i];
        // if it's right aligned, backward
        if(!isLeft)
            valWhole = cigarWhole[cigarNumWhole - i - 1];
        char opWhole = bam_cigar_op(valWhole);
        uint32_t lenWhole = bam_cigar_oplen(valWhole);

        if (opPart != opWhole)
            return false;

        if (lenPart > lenWhole)
            return false;

        if (lenPart < lenWhole) {
            // length mismatch is only allowed in the last bases
            if (i != cigarNumPart-1) {
                // we only allow one CLIP in the end
                if (i != cigarNumPart-2)
                    return false;
                // check for next, is it CLIP ?
                int next = i+1;
                uint32_t valPartNext = cigarPart[next];
                // if it's right aligned, backward
                if (!isLeft)
                    valPartNext = cigarPart[cigarNumPart - next - 1];
                char opPartNext = bam_cigar_op(valPartNext);
                uint32_t lenPartNext = bam_cigar_oplen(valPartNext);
                if (opPartNext != BAM_CHARD_CLIP)
                    return false;
            }
        }
    }

    return true;
}

void BamUtil::dumpHeader(bam_hdr_t* hdr) {
    cerr << hdr->n_targets << " contigs in the bam file:" << endl;
    int dumped = 0;
    while (dumped < hdr->n_targets) {
        char *targetName = hdr->target_name[dumped];
        int targetLen = hdr->target_len[dumped];
        string name(targetName);
        cerr << targetName << ": " << targetLen << " bp" << endl;
        dumped++;
    }
    cerr << endl;
}

/* bam_cigar_type returns a bit flag with:
 *   bit 1 set if the cigar operation consumes the query
 *   bit 2 set if the cigar operation consumes the reference
 *
 * For reference, the unobfuscated truth table for this function is:
 * BAM_CIGAR_TYPE  QUERY  REFERENCE
 * --------------------------------
 * BAM_CMATCH      1      1
 * BAM_CINS        1      0
 * BAM_CDEL        0      1
 * BAM_CREF_SKIP   0      1
 * BAM_CSOFT_CLIP  1      0
 * BAM_CHARD_CLIP  0      0
 * BAM_CPAD        0      0
 * BAM_CEQUAL      1      1
 * BAM_CDIFF       1      1
 * BAM_CBACK       0      0
 * --------------------------------
 */

const int     QUERY_CONSUM[16] = {1, 1, 0, 0, 1, 0, 0, 1, 1, 0};
const int REFERENCE_CONSUM[16] = {1, 0, 1, 1, 0, 0, 0, 1, 1, 0};

bool BamUtil::isPrimary(bam1_t *b) {
    if (b->core.flag & BAM_FSECONDARY || b->core.flag & BAM_FSUPPLEMENTARY)
        return false;
    else
        return true;
}

bool BamUtil::isProperPair(bam1_t *b) {
    return b->core.flag & BAM_FPROPER_PAIR;
}

int BamUtil::getRightRefPos(bam1_t *b) {
    if (b->core.pos < 0)
        return -1;
    return b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
}
