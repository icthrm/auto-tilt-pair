#include "searchtree.h"


SearchTree::SearchTree(const std::vector< util::point2d >& fixed, bool alloc_buffer):nk(fixed.size()), psdim(2), pseps(0.01), psNPts(fixed.size())
{
    psQryPt = annAllocPt(psdim); // query point
    psDataPts = annAllocPts(psNPts, psdim); // data points
	if(alloc_buffer){
		psNNIdx = new ANNidx[nk]; // near neighbor indices
		psDists = new ANNdist[nk]; // near neighbor distances
	}
	else{
		psNNIdx = NULL;
		psDists = NULL;
	}

    for(int i = fixed.size(); i--;){
		psDataPts[i][0] = fixed[i].x;
		psDataPts[i][1] = fixed[i].y;
	}
	
	psKdTree = new ANNkd_tree(psDataPts, psNPts, psdim);
}

double* SearchTree::DistanceArray()
{
	return psDists;
}

int* SearchTree::IndicesArray()
{
	return psNNIdx;
}

void SearchTree::GetKNearestPoints(const util::point2d& quary, int k, int* indices, double* squared_distances)
{
	psQryPt[0] = quary.x;
	psQryPt[1] = quary.y;
	psKdTree->annkSearch(psQryPt, k, indices, squared_distances, pseps);
}

void SearchTree::GetRadiusNearestPoints(const util::point2d& quary, double radius, int* k, int* indices, double* squared_distances)
{
	psQryPt[0] = quary.x;
	psQryPt[1] = quary.y;
	*k = psKdTree->annkFRSearch(psQryPt, radius, psNPts, indices, squared_distances, pseps);
}

SearchTree::~SearchTree()
{
	if(psNNIdx) delete [] psNNIdx;
    if(psDists) delete [] psDists;
    if(psKdTree) delete psKdTree;
    if(psDataPts) annDeallocPts(psDataPts);
    if(psQryPt) annDeallocPt(psQryPt);
}

