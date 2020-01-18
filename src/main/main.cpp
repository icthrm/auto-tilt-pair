#include "opts.h"
#include <iostream>
#include <fstream>
#include "dataf/dataf.h"
#include "modelmatch/match_core.h"

using namespace std;
// using namespace ann_1_1_char;

int main(int argc, char **argv) {
    struct options opts;
    opts.diameter = -1;
    opts.verbose = 0;
    opts.rotation_angle = 0;
	
	GetOpts(argc, argv, &opts);

    util::FiducialStack fidstack;

    util::ImgMatchVector imvector;
    HSetVector hset;
	
	vector<float> angles;
    if(!util::ReadAnglesByName(opts.inputangle, &angles)) {
        std::cout<<"Can't open tilt angle file."<<endl;
        return -1;
    }
    
    IplImage* tmplt;
    util::SeriesReadFromFile(&tmplt, "avgtmplt");
	fidstack.ReadFidsByFile("fids.txt");
    EX_TIME_BEGIN("\n%sDo MatchMain", _WAVE)
    ModelMatch::MatchMain(fidstack, angles, &imvector, &hset, 0.85*tmplt->width, true, false);
    imvector.WriteVectorByFolder("matches");
    hset.WriteVectorByFolder("transmx");

    cvReleaseImage(&tmplt);

    EX_TIME_END("Do MatchMain")
}
