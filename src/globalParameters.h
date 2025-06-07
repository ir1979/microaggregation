#ifndef _GLOBAL_PARAMETERS_H_
#define _GLOBAL_PARAMETERS_H_

class globalParameters {
public:
    static int debugLevel;
    static bool ignoreDLD;
    static bool standard01;
    static bool sdid;
    static double sdidRatio;
    static const int INC_QUANTIL = 10;
};

// Initialize static members
int globalParameters::debugLevel = 0;
bool globalParameters::ignoreDLD = false;
bool globalParameters::standard01 = true;
bool globalParameters::sdid = false;
double globalParameters::sdidRatio = 0.0;

#endif
