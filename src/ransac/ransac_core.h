#ifndef RANSAC_CORE_H__
#define RANSAC_CORE_H__

#include "xform.h"
#include "dataf/dataf.h"
#include <vector>

class Ransac{
public:
    static void RansacMain(const util::ImgMatchVector& imvector, util::ImgMatchVector* inlier_vec, util::ImgMatchVector* outlier_vec = NULL)
    {
        for(int i = 0; i < imvector.Size(); i++){
            util::img_match tmp;
            tmp.idx1 = imvector[i].idx1;
            tmp.idx2 = imvector[i].idx2;
            inlier_vec->PushBack(tmp);
            outlier_vec->PushBack(tmp);

            CvMat* M = RansacXForm(imvector[i].pairs, lsq_homog, 4, 0.01, homog_xfer_err, 2.8/*3.0*/,
                                   &((*inlier_vec)[i].pairs), &((*outlier_vec)[i].pairs));
            std::cout<<"["<<imvector[i].idx1<<","<<imvector[i].idx2<<"]: "<<imvector[i].pairs.size()
                     <<" & "<<(*inlier_vec)[i].pairs.size()<<" & "<<(*outlier_vec)[i].pairs.size()<<std::endl;
            cvReleaseMat(&M);
        }
    }
};

#endif