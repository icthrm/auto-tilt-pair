#ifndef R4PESTIMATOR_H__
#define R4PESTIMATOR_H__

#include "dataf/dataf.h"
#include <ANN/ANN.h>
#include "4pcs.h"
#include "util/exception.h"
#include "mrcimg/mrc2img.h"
#include "mrcimg/img_util.h"
#include "levmar.h"
#include "clapack.h"
#include "cminpack/minpack.h"
#include "util/matrix.h"

extern "C" {
#include "f2c.h"
#include "clapack.h"
}


// #define MIN_PT_NUM		6

typedef std::pair<fpcs::Quadrilateral, fpcs::Quadrilateral> QuadPair;

class Ran4PEstimator{
public:
	static const double COVERAGE_REFRESH_RATIO = 0.95;
private:
	static const int MIN_PT_NUM = 6;
	static const double GRAPH_BASE_SAFE_TOL = 0.3;//0.3
	static const double ERROR_TOL = 0.05;
	static const double ERROR_TOL_MIN_DIST = 1.5/0.01;
	static const double RAN_COVERAGE = 0.25;
	static const double AREA_DIFF_TOL = 0.15;
	
class HIterator{	
private:
	static Ran4PEstimator* r4per;
	static mx::Matrix<2, 3, double> HH;
	static std::vector<std::pair<int, int> >* matches;
private:
	static void lmderfun(double* p, double* x, int m, int n, void *data);
	static void lmderjac(double* p, double *jac, int m, int n, void *data);
	static bool SplitToSecondRoutionData(double* p, std::vector<util::point2d>& fids1, std::vector<int>& idx1, std::vector<util::point2d>& fids2, 
		std::vector<int>& idx2, std::vector<std::pair<int, int> >& matches, std::vector<std::pair<int, int> >& submm, mx::Matrix<2, 3, double>& HH, float dist_tol);
	
public:
	static void ReEstimateTransform(std::vector<std::vector<std::pair<int, int> > >& mmset, std::vector<mx::Matrix<2, 3, double> >& hset);
	static void SetInitialInput(Ran4PEstimator* __r4per, mx::Matrix<2, 3, double>& __HH, std::vector<std::pair<int, int> >* __matches);
};
	
struct PointPair{
public:
	const util::point2d* p1;
	const util::point2d* p2;

public:	
	int idx1;
	int idx2;

	util::point2d p1top2;		//the vector of (p2-p1)
	float squared_length;
	
public:
	PointPair(int index1, int index2, const util::point2d& point1, const util::point2d& point2):idx1(index1), idx2(index2), p1(&point1), p2(&point2){
		p1top2.x = p2->x - p1->x;
		p1top2.y = p2->y - p1->y;
		squared_length = util::L2(*(p1), *(p2));
	}
	
	const util::point2d& P1() const{return *p1;}
	
	const util::point2d& P2() const{return *p2;}
	
	void FractionPoint(float fraction, util::point2d& fract);	//a=p1+fraction*(p2-p1)
};

private:
	std::vector<util::point2d>& fids1;			//point set P; select a base Quadrilateral from P
	std::vector<util::point2d>& fids2;			//point set Q; generate a reference tree form Q
	
	std::vector<util::point2d> unmatched_fids1;
	std::vector<util::point2d> unmatched_fids2;
	
	std::vector<PointPair> ppairs1;
	std::vector<PointPair> ppairs2;
	
	double graph_baseline_safe_tol;
	double graph_diameter;
	double min_baseline;
	double pts_dist_tol;
	
	double cos_beta1, cos_beta2;					// cosine of the rotation angle

private:
	void PickRandomQuadrilateral(const std::vector<PointPair>& ppairs, fpcs::Quadrilateral& quad, float* invariant1, float* invariant2);
	
	void GeneratePointPairs(const std::vector<util::point2d>& fids, std::vector<PointPair>& pointpairs);
	void FilterPointPairs(float max_dist, float min_dist, std::vector<PointPair>& pointpairs);
	void PairDistanceMeanStd(const std::vector<PointPair>& pointpairs, float* dist_mean, float* dist_std);
	void PrepareLocalLimitedPointPairs(const std::vector<util::point2d>& fids, std::vector<PointPair>& pointpairs);
	static void ApplyTransformToAll(const std::vector<util::point2d>& fids1, std::vector<util::point2d>& fidtrans, const mx::Matrix<2, 3, double>& H);
	
	void SearchQuadrilateralPairs(const fpcs::Quadrilateral& quad, float invariant1, float invariant2, float area_diff_tol, float dist_err_tol, std::vector<QuadPair>& pairs);
	
	/**affine transform **/
	void EstimateTransform(const std::vector<util::point2d>& fids1, 
								 const std::vector<util::point2d>& fids2, const QuadPair& qpair, mx::Matrix<2, 3, double>& H);
	void EstimateTransform(const std::vector<util::point2d>& fids1, const std::vector<util::point2d>& fids2, 
						   const std::vector<std::pair<int, int> >& matches, mx::Matrix<2, 3, double>& H);
	
	void CalculateMatch(const std::vector<util::point2d>& fids1, const std::vector<util::point2d>& fids2, 
						   std::vector<std::pair<int, int> >& rawmatches, std::vector<float>& dists, double dist_err_tol);
	void EigenValue2X2(const mx::Matrix<2, 3, double>& H, float* eval1, float* eval2);
	void RansacMatch(double dist_err_tol, std::vector<std::pair<int, int> >& matches, std::vector<std::pair<int, int> >& ambiguity, mx::Matrix<2, 3, double>& HH, bool eigenlimit = true);
	void FindUnmatchedAmbiguity(const std::vector<std::vector<std::pair<int, int> > >& mmset, const std::vector<mx::Matrix<2, 3, double> >& hset, double search_radius);
	
public:
	Ran4PEstimator(std::vector<util::point2d>& __fids1, std::vector<util::point2d>& __fids2);
	~Ran4PEstimator();
	
	void SetSupplementData(double beta1, double beta2);	
	
	std::vector<util::point2d>& UmFids1();
	std::vector<util::point2d>& UmFids2();
	
	/** The estimated transform hset should be applied to fids1 **/
	void AffineTransformEstimation(double dist_err_tol, std::vector<std::pair<util::point2d, util::point2d> >& matchvec, 
								   std::vector<mx::Matrix<2, 3, double> >& hset, bool test = false, bool eigenlimit = true);
	void DistributionCorrection(double dist_err_tol, std::vector<std::pair<util::point2d, util::point2d> >& matchvec);
};


#endif

