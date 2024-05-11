#include "options.h"
#include "util.h"

Options::Options() {
    input = "";
    output = "";
    output_idx = "";
    refFile = "";
    regionFile = "";

    nthreads = 1;
    endPosition = 12;
    refRange = 300;

    debug = false;
}

bool Options::validate() {
    if (input.empty()) {
        error_exit("input should be specified by --in1");
    } else {
        check_file_valid(input);
    }

    if (output.size() > 4 && output.substr(output.size()-3) == "bam")
        output_idx = output.substr(0, output.size() - 4) + ".bai";

    if (refFile.empty()) {
        cerr << "reference fasta file should be specified by -r" << endl;
        exit(-1);
    }

    if (regionFile.empty()) {
        cerr << "region bed/interval file should be specified by -g" << endl;
        exit(-1);
    }

    if (endPosition > 15) {
        error_exit("ratio_threshold cannot be greater than 15.");
    } else if (endPosition < 3) {
        error_exit("ratio_threshold cannot be less than 3.");
    }

    if (refRange > 1000) {
        error_exit("ratio_threshold cannot be greater than 1000.");
    } else if(refRange < 30) {
        error_exit("ratio_threshold cannot be less than 30.");
    }

    return true;
}
