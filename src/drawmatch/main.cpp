#include "opts.h"
#include <iostream>
#include <fstream>
#include "dataf/dataf.h"
#include "mrcimg/mrc2img.h"
#include "modelmatch/match_core.h"

// #define TESTDEVELOP

using namespace std;
// using namespace ann_1_1_char;

int main(int argc, char **argv) {

    EX_TRACE("\nMARKERAUTO --version 1.5.3\n")

    struct options opts;
    opts.diameter = -1;
    opts.verbose = 0;
    opts.rotation_angle = 0;

    if(GetOpts(argc, argv, &opts) <= 0) {
        EX_TRACE("***WRONG INPUT.\n");
        return -1;
    }

    util::MrcStack mrcs;
    mrcs.Open(opts.input);

    util::FiducialStack fidstack;

    util::ImgMatchVector imvector;
    HSetVector hset;
	
	imvector.ReadVectorByFolder("matches");
    hset.ReadVectorByFolder("transmx");
	
	IplImage* tmplt;
    util::SeriesReadFromFile(&tmplt, "avgtmplt");
	fidstack.ReadFidsByFile("fids.txt");

    ModelMatch::Test(mrcs, fidstack, tmplt->width, imvector, 0.5f, "matches_ill");

    mrcs.Close();
}

