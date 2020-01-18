#ifndef IFGT_H__
#define IFGT_H__

#include <vector>
#include <dataf/dataf.h>

namespace ifgt{

#define TRUNCATION_NUMBER_UL  200
	
// The results from k-means clustering.
class Clustering{
public:
    double max_radius;						// The maximum cluster radius.
    std::vector<int> indices;				// The cluster membership ids for each points.
    std::vector<util::point2d> centers;		// The centers of each cluster.
    std::vector<int> npoints;				// The number of points in each cluster.
    std::vector<double> radii;				// The radius of each cluster.
public:
	//k-means clustering, specifying the starting cluster centers.
	Clustering(){};
	Clustering(const std::vector<util::point2d>& points, double epsilon, const std::vector<util::point2d>& icenters);
	Clustering(const std::vector<util::point2d>& points, int nclusters, double epsilon);
};

class IFGaussianTransform{
private:
	const std::vector<util::point2d>& source;
	double bandwidth;
	double epsilon;
    int truncation_number;
    int p_max_total;
    std::vector<double> constant_series;
    std::vector<double> ry_square;
	
	// The number of clusters that should be used for the IFGT.
    int nclusters;
    // The cutoff radius.
    double cutoff_radius;
	Clustering clustering;
	
private:
    void CalculateMonomials(const double* d, std::vector<double>& monomials) const;
    void CalculateConstantSeries(std::vector<double>& constant_series);
	
public:
	IFGaussianTransform(const std::vector<util::point2d>& source, double bandwidth, double eps);
// 	~IFGaussianTransform();
	
	void CalculateTransformWith(const std::vector<util::point2d>& target, std::vector<double>& result);
	void CalculateTransformWith(const std::vector<util::point2d>& target, const std::vector<double>& weights, std::vector<double>& result);
};

}

#endif
