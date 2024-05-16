#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include "htslib/sam.h"

using namespace std;


class Options {

public:
    Options();
    bool validate();

public:
    string input;
    string output;
    string output_idx;
    string refFile;
    string regionFile;
    bool debug;

    int nthreads;
    // thresholds
    int endPosition;
    int refRange;
    int maxdel;
};

#endif