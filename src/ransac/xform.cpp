/*
  This file contains definitions for functions to compute transforms from
  image feature correspondences

  Copyright (C) 2006-2010  Rob Hess <hess@eecs.oregonstate.edu>

  @version 1.1.2-20100521
*/

#include "xform.h"
#include "util/exception.h"
#include "dataf/keypoint.h"
#include <cxcore.h>
#include <stdlib.h>
#include <ctime>

/************************* Local Function Prototypes *************************/

static int CalcMinInliers(int n, int m, double p_badsupp, double p_badxform);
static inline double log_factorial(int);
static inline double dist_sq_2D(CvPoint2D32f p1, CvPoint2D32f p2);
static int* DrawRansacSampleIndex(bool* flags, int ori_size, int sam_size);
static void ExtractCorresPts(const std::vector<util::pair>& pairs, int* idxs, int sam_size, CvPoint2D32f** pts, CvPoint2D32f** mpts);
static int FindConsensus(const std::vector<util::pair>& pairs, CvMat* M, RansacErrFUNC err_fn, double err_tol, int** consensus);
static inline void ReleaseMem(CvPoint2D32f* pts1, CvPoint2D32f* pts2, int* sam_idxs);

#define __utilpt2cv2d(P) cvPoint2D32f((P).x, (P).y)

/********************** Functions prototyped in model.h **********************/

/**
  Calculates a best-fit image transform from image feature correspondences
  using RANSAC.

  @param features an array of features; only features with a non-NULL match
    of type mtype are used in homography computation
  @param n number of features in feat
  @param mtype determines which of each feature's match fields to use
    for model computation; should be one of FEATURE_FWD_MATCH,
    FEATURE_BCK_MATCH, or FEATURE_MDL_MATCH; if this is FEATURE_MDL_MATCH,
    correspondences are assumed to be between a feature's img_pt field
    and its match's mdl_pt field, otherwise correspondences are assumed to
    be between the the feature's img_pt field and its match's img_pt field
  @param xform_fn pointer to the function used to compute the desired
    transformation from feature correspondences
  @param m minimum number of correspondences necessary to instantiate the
    model computed by xform_fn
  @param p_badxform desired probability that the final transformation
    returned by RANSAC is corrupted by outliers (i.e. the probability that
    no samples of all inliers were drawn)
  @param err_fn pointer to the function used to compute a measure of error
    between putative correspondences and a computed model
  @param err_tol correspondences within this distance of a computed model are
    considered as inliers
  @param inliers if not NULL, output as an array of pointers to the final
    set of inliers
  @param n_in if not NULL and \a inliers is not NULL, output as the final
    number of inliers

  @return Returns a transformation matrix computed using RANSAC or NULL
    on error or if an acceptable transform could not be computed.
*/
CvMat* RansacXForm(const std::vector<util::pair>& pairs, RansacXFormFUNC xform_fn,
                   int min_req, double p_badxform, RansacErrFUNC err_fn,
                   double err_tol, std::vector<util::pair>* inliers, std::vector<util::pair>* outliers)
{
    int* sample, * consensus, * consensus_max = NULL;
    int consen_size, consen_min, consen_max = 0;

    CvPoint2D32f* pts, * mpts;
    int ori_size = pairs.size();
    bool* flags = new bool[ori_size];

    CvMat* M = NULL;
    double p, in_frac = RANSAC_INLIER_FRAC_EST;
    int k = 0;

    if(ori_size < min_req) {
        EX_TRACE("Warning: not enough matches to compute xform(%d < %d)", ori_size, min_req)
        goto end;
    }

    srandom(time(NULL));

    consen_min = CalcMinInliers(ori_size, min_req, RANSAC_PROB_BAD_SUPP, p_badxform);
    p = pow(1.0 - pow(in_frac, min_req), k);

    while(p > p_badxform) {
        sample = DrawRansacSampleIndex(flags, ori_size, min_req);

        ExtractCorresPts(pairs, sample, min_req, &pts, &mpts);
        M = xform_fn(pts, mpts, min_req);

        if(!M) {
            goto iteration_end;
        }

        consen_size = FindConsensus(pairs, M, err_fn, err_tol, &consensus);

        if(consen_size > consen_max) {
            if(consensus_max) {
                delete [] consensus_max;
            }
            consensus_max = consensus;
            consen_max = consen_size;
            in_frac = (double)consen_max/ori_size;
        }
        else {
            delete [] consensus;
        }
        cvReleaseMat(&M);

iteration_end:
        ReleaseMem(pts, mpts, sample);
        p = pow(1.0 - pow(in_frac, min_req), ++k);
    }

    /* calculate final transform based on best consensus set */
    if(consen_max >= consen_min) {
        ExtractCorresPts(pairs, consensus_max, consen_max, &pts, &mpts);
        M = xform_fn(pts, mpts, consen_max);
        consen_size = FindConsensus(pairs, M, err_fn, err_tol, &consensus);
        cvReleaseMat(&M);
        ReleaseMem(pts, mpts, consensus_max);
        ExtractCorresPts(pairs, consensus, consen_size, &pts, &mpts);
        M = xform_fn(pts, mpts, consen_size);

        if(inliers) {
            for(int i = 0; i < consen_size; i++) {
                inliers->push_back(pairs[consensus[i]]);
            }
        }

        if(outliers) {
            int* ptr = consensus;
            int* end = consensus+consen_size;

            for(int i = 0; i < ori_size; i++) {
                if(ptr != end) {
                    if(i < *ptr) {
                        outliers->push_back(pairs[i]);
                    }
                    else if(i == *ptr) {
                        ptr++;
                    }
                }
                else {
                    outliers->push_back(pairs[i]);
                }
            }
        }

        ReleaseMem(pts, mpts, consensus);
    }
    else if(consensus_max) {
        delete [] consensus_max;
    }

end:
    delete [] flags;
    return M;
}



/*
  Calculates a planar homography from point correspondeces using the direct
  linear transform.  Intended for use as a ransac_xform_fn.

  @param pts array of points
  @param mpts array of corresponding points; each pts[i], i=0..n-1,
    corresponds to mpts[i]
  @param n number of points in both pts and mpts; must be at least 4

  @return Returns the 3x3 planar homography matrix that transforms points
    in pts to their corresponding points in mpts or NULL if fewer than 4
    correspondences were provided
*/
CvMat* dlt_homog(CvPoint2D32f* pts, CvPoint2D32f* mpts, int n)
{
    CvMat* H, * A, * VT, * D, h, v9;
    double _h[9];
    int i;

    if(n < 4)
        return NULL;

    /* set up matrices so we can unstack homography into h; Ah = 0 */
    A = cvCreateMat(2*n, 9, CV_64FC1);
    cvZero(A);
    for(i = 0; i < n; i++)
    {
        cvmSet(A, 2*i, 3, -pts[i].x);
        cvmSet(A, 2*i, 4, -pts[i].y);
        cvmSet(A, 2*i, 5, -1.0 );
        cvmSet(A, 2*i, 6, mpts[i].y * pts[i].x);
        cvmSet(A, 2*i, 7, mpts[i].y * pts[i].y);
        cvmSet(A, 2*i, 8, mpts[i].y);
        cvmSet(A, 2*i+1, 0, pts[i].x);
        cvmSet(A, 2*i+1, 1, pts[i].y);
        cvmSet(A, 2*i+1, 2, 1.0 );
        cvmSet(A, 2*i+1, 6, -mpts[i].x * pts[i].x);
        cvmSet(A, 2*i+1, 7, -mpts[i].x * pts[i].y);
        cvmSet(A, 2*i+1, 8, -mpts[i].x);
    }
    D = cvCreateMat(9, 9, CV_64FC1);
    VT = cvCreateMat(9, 9, CV_64FC1);
    cvSVD(A, D, NULL, VT, CV_SVD_MODIFY_A + CV_SVD_V_T);
    v9 = cvMat(1, 9, CV_64FC1, NULL);
    cvGetRow(VT, &v9, 8);
    h = cvMat(1, 9, CV_64FC1, _h);
    cvCopy(&v9, &h, NULL);
    h = cvMat(3, 3, CV_64FC1, _h);
    H = cvCreateMat(3, 3, CV_64FC1);
    cvConvert(&h, H);

    cvReleaseMat(&A);
    cvReleaseMat(&D);
    cvReleaseMat(&VT);
    return H;
}



/**
  Calculates a least-squares planar homography from point correspondeces.

  @param pts array of points
  @param mpts array of corresponding points; each pts[i], i=1..n, corresponds
    to mpts[i]
  @param n number of points in both pts and mpts; must be at least 4

  @return Returns the 3 x 3 least-squares planar homography matrix that
    transforms points in pts to their corresponding points in mpts or NULL if
    fewer than 4 correspondences were provided
*/
CvMat* lsq_homog(CvPoint2D32f* pts, CvPoint2D32f* mpts, int n)
{
    CvMat* H, * A, * B, X;
    double x[9];
    int i;

    if(n < 4)
    {
        fprintf(stderr, "Warning: too few points in lsq_homog(), %s line %d\n",
                __FILE__, __LINE__);
        return NULL;
    }

    /* set up matrices so we can unstack homography into X; AX = B */
    A = cvCreateMat(2*n, 8, CV_64FC1);
    B = cvCreateMat(2*n, 1, CV_64FC1);
    X = cvMat(8, 1, CV_64FC1, x);
    H = cvCreateMat(3, 3, CV_64FC1);
    cvZero(A);
    for(i = 0; i < n; i++)
    {
        cvmSet(A, i, 0, pts[i].x);
        cvmSet(A, i+n, 3, pts[i].x);
        cvmSet(A, i, 1, pts[i].y);
        cvmSet(A, i+n, 4, pts[i].y);
        cvmSet(A, i, 2, 1.0);
        cvmSet(A, i+n, 5, 1.0);
        cvmSet(A, i, 6, -pts[i].x * mpts[i].x);
        cvmSet(A, i, 7, -pts[i].y * mpts[i].x);
        cvmSet(A, i+n, 6, -pts[i].x * mpts[i].y);
        cvmSet(A, i+n, 7, -pts[i].y * mpts[i].y);
        cvmSet(B, i, 0, mpts[i].x);
        cvmSet(B, i+n, 0, mpts[i].y);
    }
    cvSolve(A, B, &X, CV_SVD);
    x[8] = 1.0;
    X = cvMat(3, 3, CV_64FC1, x);
    cvConvert(&X, H);

    cvReleaseMat(&A);
    cvReleaseMat(&B);
    return H;
}



/*
  Calculates the transfer error between a point and its correspondence for
  a given homography, i.e. for a point x, it's correspondence x', and
  homography H, computes d(x', Hx)^2.

  @param pt a point
  @param mpt pt's correspondence
  @param H a homography matrix

  @return Returns the transfer error between pt and mpt given H
*/
double homog_xfer_err(CvPoint2D32f pt, CvPoint2D32f mpt, CvMat* H)
{
    CvPoint2D32f xpt = persp_xform_pt(pt, H);

    return sqrt(dist_sq_2D(xpt, mpt));
}

double dist_sq_2D(CvPoint2D32f p1, CvPoint2D32f p2)
{
    double x_diff = p1.x - p2.x;
    double y_diff = p1.y - p2.y;

    return x_diff * x_diff + y_diff * y_diff;
}

/*
  Performs a perspective transformation on a single point.  That is, for a
  point (x, y) and a 3 x 3 matrix T this function returns the point
  (u, v), where

  [x' y' w']^T = T * [x y 1]^T,

  and

  (u, v) = (x'/w', y'/w').

  Note that affine transforms are a subset of perspective transforms.

  @param pt a 2D point
  @param T a perspective transformation matrix

  @return Returns the point (u, v) as above.
*/
CvPoint2D32f persp_xform_pt(CvPoint2D32f pt, CvMat* T)
{
    CvMat XY, UV;
    double xy[3] = { pt.x, pt.y, 1.0 }, uv[3] = { 0 };
    CvPoint2D32f rslt;

    cvInitMatHeader(&XY, 3, 1, CV_64FC1, xy, CV_AUTOSTEP);
    cvInitMatHeader(&UV, 3, 1, CV_64FC1, uv, CV_AUTOSTEP);
    cvMatMul(T, &XY, &UV);
    rslt = cvPoint2D32f(uv[0] / uv[2], uv[1] / uv[2]);

    return rslt;
}

/**
  Calculates the minimum number of inliers as a function of the number of
  putative correspondences.  Based on equation (7) in

  Chum, O. and Matas, J.  Matching with PROSAC -- Progressive Sample Consensus.
  In <EM>Conference on Computer Vision and Pattern Recognition (CVPR)</EM>,
  (2005), pp. 220--226.

  @param n number of putative correspondences
  @param m min number of correspondences to compute the model in question
  @param p_badsupp prob. that a bad model is supported by a correspondence
  @param p_badxform desired prob. that the final transformation returned is bad

  @return Returns the minimum number of inliers required to guarantee, based
    on p_badsupp, that the probability that the final transformation returned
    by RANSAC is less than p_badxform
*/
int CalcMinInliers(int n, int m, double p_badsupp, double p_badplane)
{
    double pi, sum;
    int i, j;

    for(j = m+1; j <= n; j++) {
        sum = 0;
        for(i = j; i <= n; i++) {
            pi = (i-m) * log(p_badsupp) + (n-i+m) * log(1.0 - p_badsupp) +
                 log_factorial(n - m) - log_factorial(i - m) -
                 log_factorial(n - i);
            /*
             * Last three terms above are equivalent to log(n-m choose i-m)
             */
            sum += exp(pi);
            if(sum < p_badplane) {
                break;
            }
        }
        if(sum < p_badplane) {
            break;
        }
    }
    return j;
}

/** Calculates the natural log of the factorial of a number @param n number @return Returns log(n!) */
inline double log_factorial(int n)
{
    double f = 0;
    int i;

    for(i = 1; i <= n; i++)
        f += log(i);

    return f;
}

int* DrawRansacSampleIndex(bool* flags, int ori_size, int sam_size)
{
    memset(flags, 0, sizeof(bool)*ori_size);

    int* sam_idxs = new int[sam_size];
    int ranx;

    for(int i = 0; i < sam_size; i++) {
        do {
            ranx = random()%ori_size;
        } while(flags[ranx]);
        sam_idxs[i] = ranx;
        flags[ranx] = true;
    }

    return sam_idxs;
}

void ExtractCorresPts(const std::vector< util::pair >& pairs, int* idxs, int sam_size, CvPoint2D32f** pts, CvPoint2D32f** mpts)
{
    CvPoint2D32f* _pts, * _mpts;
    _pts = new CvPoint2D32f[sam_size];
    _mpts = new CvPoint2D32f[sam_size];

    for(int i = 0; i < sam_size; i++) {
        _pts[i] = __utilpt2cv2d(pairs[idxs[i]].first);
        _mpts[i] = __utilpt2cv2d(pairs[idxs[i]].second);
    }

    *pts = _pts;
    *mpts = _mpts;
}

/**For a given model and error function, finds a consensus from a set of feature correspondences.
  @param err_fn error function used to measure distance from M
  @param err_tol correspondences within this distance of M are added to the consensus set
  @param consensus output as an array of index in the pairs set
  @return Returns the number of points in the consensus set
*/
int FindConsensus(const std::vector<util::pair>& pairs, CvMat* M, RansacErrFUNC err_fn, double err_tol, int** consensus)
{
    CvPoint2D32f pt, mpt;
    double err;

    int* consensus_idxs = new int[pairs.size()];
    int index = 0;

    for(int i = 0; i < pairs.size(); i++) {
        pt = __utilpt2cv2d(pairs[i].first);
        mpt = __utilpt2cv2d(pairs[i].second);
        err = err_fn(pt, mpt, M);
        if(err <= err_tol) {
            consensus_idxs[index++] = i;
        }
    }

    *consensus = consensus_idxs;
    return index;
}

/** Releases memory and reduces code size above
  @param pts1 an array of points
  @param pts2 an array of points
  @param features an array of pointers to features; can be NULL
*/
void ReleaseMem(CvPoint2D32f* pts1, CvPoint2D32f* pts2, int* sam_idxs)
{
    delete [] pts1;
    delete [] pts2;
    if(sam_idxs) {
        delete [] sam_idxs;
    }
}
