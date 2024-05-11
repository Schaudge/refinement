/*
 * correcting dual artificial sequences introduced by enzyme digestion
 * Written By Schaudge King
 * 2021-06-09
 * MIT LICENSE
 */
#include <stdio.h>
#include <time.h>
#include <sstream>
#include "cmdline.h"
#include "common.h"
#include "recsc.h"
#include "options.h"

using namespace std;

string command;

int main(int argc, char* argv[]){

    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cerr << "recsc " << VERSION_NUMBER << endl;
        return 0;
    }

    cmdline::parser cmd;
    // input/output
    cmd.add<string>("in", 'i', "input sorted bam/sam file. STDIN will be read from if it's not specified", true, "-");
    cmd.add<string>("out", 'o', "output bam/sam file. STDOUT will be written to if it's not specified", false, "-");
    cmd.add<string>("ref", 'r', "reference fasta file name (should be an uncompressed .fa/.fasta file)", true, "");
    cmd.add<string>("region", 'g', "(sorted) bed/interval file to specify the capturing region, none by default", true, "");
    cmd.add<int>("thread", 'p', "the number of threads for processing bam file", false, 1);

    // thresholds
    cmd.add<int>("max_end_distance", 'm', "The maximum length for the reads position with mismatch base, we trigger the correct strategy to run. The value should less than 15.", false, 12);
    cmd.add<int>("near_ref_range", 'n', "The range size for the realigned reference sequence around the checked genome position. The value should less than 1000.", false, 300);

    // debugging
    cmd.add("debug", 0, "output some debug information to STDERR.");

    cmd.parse_check(argc, argv);

    Options opt;
    opt.input = cmd.get<string>("in");
    opt.output = cmd.get<string>("out");
    opt.refFile = cmd.get<string>("ref");
    opt.regionFile = cmd.get<string>("region");
    opt.nthreads = cmd.get<int>("thread");
    opt.endPosition = cmd.get<int>("max_end_distance");
    opt.refRange = cmd.get<int>("near_ref_range");
    opt.debug = cmd.exist("debug");

    opt.validate();
    
    time_t t1 = time(NULL);

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    Recsc recsc(&opt);
    //recsc.correct();
    recsc.correct2();

    time_t t2 = time(NULL);
    cerr << command << "\n";
    cerr << "recsc v" << VERSION_NUMBER << ", time used: " << (t2)-t1 << " seconds." << endl;

    return 0;
}
