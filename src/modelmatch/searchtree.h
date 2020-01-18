#ifndef SEARCHTREE_H__
#define SEARCHTREE_H__

#include <ANN/ANN.h>
#include <vector>
#include "dataf/dataf.h"

// #define ERROR_TOL  0.05

class SearchTree{
private:
	int nk;
	int psdim;
	double pseps;
	int psNPts; // actual number of data points for kdtree
    ANNpoint psQryPt; // query point
    ANNpointArray psDataPts; // data points
    ANNidxArray psNNIdx; // near neighbor indices
    ANNdistArray psDists; // near neighbor distances
    ANNkd_tree* psKdTree; // search structure
public:
	SearchTree():psKdTree(NULL), psQryPt(NULL), psDataPts(NULL){}
	SearchTree(const std::vector<util::point2d>& fixed, bool alloc_buffer = true);
	~SearchTree();
	void GetKNearestPoints(const util::point2d& quary, int k, int* indices, double* squared_distances);
	void GetRadiusNearestPoints(const util::point2d& quary, double squared_radius, int* k, int* indices, double* squared_distances);
	int* IndicesArray();
	double* DistanceArray();
};

#endif
