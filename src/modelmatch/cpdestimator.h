#ifndef CPDESTIMATOR_H__
#define CPDESTIMATOR_H__

#include "dataf/dataf.h"
#include "util/exception.h"
#include "mrcimg/mrc2img.h"
#include "mrcimg/img_util.h"
#include "ifgt.h"
#include "levmar.h"
#include "clapack.h"
#include "cminpack/minpack.h"
#include "searchtree.h"
#include "util/matrix.h"

/*
#define DEFAULT_LAMBDA		3
#define DEFAULT_BETA		3*/

class CPDEstimator{
public:
	static const double COVERAGE_REFRESH_RATIO = 0.7;
private:
	static const int MAX_CPD_ITERATION = 200;
	static const int DEFAULT_NORMALIZE = 1;
	static const double OUTLIERS_WEIGHT = 0.1;//0.1
	static const double ERROR_TOL = 1e-6;
	static const int DEFAULT_SIGMA2 = -1;
	static const int MIN_PT_NUM = 6;
	static const int KDTREE_NK = 5;
	
	static const double DEFAULT_BETA = 3.0;		// Default value for beta (non-rigid). 
	static const double DEFAULT_LAMBDA = 3.0;			// Default value for lambda (non-rigid).
	
class Normalizer{
private:
	static void CalculateCenter(const std::vector<util::point2d>& pts, util::point2d& center);
	static void MovePointSetToCenter(const util::point2d& center, std::vector<util::point2d>& pts);
	static void CalculateScale(const std::vector<util::point2d>& pts, float* scale);
public:
	static void ScalePointSet(const float scale, std::vector<util::point2d>& pts);
	static void GlobalNormalize(std::vector<util::point2d>& fixed, std::vector<util::point2d>& moving, 
								util::point2d& fixed_reverse_center, util::point2d& moving_reverse_center, float* scale);
	static void ReverseNormalization(const util::point2d& reverse_center, float scale, std::vector<util::point2d>& points);
};

// Probability matrices produced by comparing two data sets with a `GaussTransform`.
struct Probabilities{
    // The probability matrix, multiplied by the identity matrix.
    std::vector<double> p1;
    // The probability matrix, transposes, multiplied by the identity matrix.
    std::vector<double> pt1;
    // The probability matrix multiplied by the fixed points.
    std::vector<util::point2d> px;
    // The total error.
    double l;
	
public:
	Probabilities(){}
};

public:
/**
 *     |p0 p1|       |p4|
 * R = |p2 p3|   t = |p5|
 * 
 * if it is rigid transform: x'=s*R*x+t, R is constranted by rigid;
 * if it is affine transform: x'=R*x+t, R is not constranted, scale is not used;
 * 
 */
struct Params{
    double p[6];
    double scale;
	double sigma2;
};

struct NonrigidParams{
	std::vector<std::vector<float> > m_g;
	std::vector<std::vector<float> > m_w;
	double sigma2;
	double m_lambda;
	double m_beta;
};

private:
	double default_sigma2;
	double cpd_err_tol;		//error tolerance
	double outliers_weight;
	int max_iter;
	mx::Matrix<2, 3, double> initial_trans_mat;
	
	const std::vector<util::point2d>& ori_fixed;
	const std::vector<util::point2d>& ori_moving;
	SearchTree ftree;					// served as the search tree for fixed point set 
	std::vector<util::point2d> fixed;	// will be normalized in the construct method
	std::vector<util::point2d> moving;	// will be normalized in the construct method
	util::point2d fixed_center;			// save the old center of fixed point set
	util::point2d moving_center;		// save the old center of moving point set
	float scale;						// save the scale 
	std::vector<util::point2d> points;
	std::vector<int> correspondence;
	Params fin_params;						// for fine align
	NonrigidParams nr_params;
	double pts_dist_tol;
	
	static double DefaultSigma2(const std::vector<util::point2d>& fixed, const std::vector<util::point2d>& moving);
	static double GetMaxDiscernibleSigma(const std::vector<util::point2d>& pts1, const std::vector<util::point2d>& pts2);
	
	static void ComputeDirectGaussTransform(const std::vector<util::point2d>& fixed, 
									   const std::vector<util::point2d>& moving, double sigma2, double outliers, Probabilities& probabilities);
	static void ComputeFastGaussTransform(const std::vector<util::point2d>& fixed, 
												   const std::vector<util::point2d>& moving, double sigma2, double outliers, Probabilities& probabilities);
	static void CalculateCorrespondence(const std::vector< util::point2d >& fixed, 
										const std::vector< util::point2d >& moving, double sigma2, double outliers, std::vector<int>& correspondence);
	
	static void AffineParamsProduct(double h, double alpi1, double alpi2, double t0, double t1, Params& tparams);
	
	void InitializeNonrigidParams();
	void ModifyProbabilities(const NonrigidParams& nr_params, Probabilities& probabilities);
	
	double ComputeDirectGaussCorrelation(const std::vector<util::point2d>& fixed, const std::vector<util::point2d>& moving, double sigma2);
	double ComputeTreeGaussCorrelation(const std::vector<util::point2d>& fixed, const std::vector<util::point2d>& moving, double sigma2);
	
	void MaximumStepAsAffineTransform(const std::vector<util::point2d>& fixed, const std::vector<util::point2d>& moving, const Probabilities& probabilities, double sigma2, Params& params);
	
	void MaximumStepAsNonrigidTransform(const std::vector<util::point2d>& fixed, const std::vector<util::point2d>& moving, 
										const Probabilities& probabilities, double sigma2, NonrigidParams& result);
	
	static void UpdatePointSet(const std::vector<util::point2d>& ori, const Params& params, std::vector<util::point2d>& result);
	static void UpdatePointSet(const std::vector<util::point2d>& ori, const mx::Matrix<2, 3, double>& A, std::vector<util::point2d>& result);
	static void UpdatePointSetDrift(const std::vector<util::point2d>& ori, const NonrigidParams& params, std::vector<util::point2d>& result);
	
	bool GaussianMixedModelMatch(double dist_err_tol, std::vector<std::pair<int, int> >& matches, bool eigenlimit = true);
	/* should not be used solo */
	void MiniOptimizeByAffineGMM(const std::vector<util::point2d>& fixed, std::vector<util::point2d>& moving);
	
	void ComputeTransformByLeastSquares(const std::vector<util::point2d>& fixed, const std::vector<util::point2d>& moving, 
										const std::vector<std::pair<int, int> >& matches, Params& params);
	
	/** a wrap with old interface */
	void CalculateMatch(const std::vector<util::point2d>& fids1, const std::vector<util::point2d>& fids2, 
						   std::vector<std::pair<int, int> >& rawmatches, std::vector<double>& dists, double dist_err_tol);
	void FromCorrespondenceToMatch(const std::vector< util::point2d >& fids1, const std::vector< util::point2d >& fids2, std::vector<int>& correspondence,
									   std::vector< std::pair<int, int> >& rawmatches, double dist_err_tol) const;
	
public:
	CPDEstimator(const std::vector<util::point2d>& __fids1, const std::vector<util::point2d>& __fids2);
	
	void PointDriftEstimation(double dist_err_tol, std::vector<std::pair<util::point2d, util::point2d> >& matchvec, 
							  std::vector<mx::Matrix<2, 3, double> >& hset, bool test = false, bool eigenlimit = true);
	void DistributionCorrection(double dist_err_tol, std::vector<std::pair<util::point2d, util::point2d> >& matchvec);
};

#endif
