#ifndef FPCS_CORE_H__
#define FPCS_CORE_H__

#include "dataf/dataf.h"
#include <vector>
#include <ctime>
#include <cstdlib>
#include "ANN/ANN.h"

namespace fpcs{

/**
 * v[4]; a, b, c, d;
 * a   c
 * 	\ /
 * 	 e
 *  / \
 * d   b           invarant1=|a-e|/|a-b|; invarant2=|c-e|/|c-d|
 */
struct Quadrilateral{
    int v[4];
	double S; 		//value of area
	Quadrilateral(){}
    Quadrilateral(int vertex0, int vertex1, int vertex2, int vertex3){
        v[0] = vertex0;
        v[1] = vertex1;
        v[2] = vertex2;
        v[3] = vertex3;
    }
    int& operator[](int idex){
		return v[idex];
    }
    
    const int& operator[](int idex) const{
		return v[idex];
    }
    
    const Quadrilateral& operator=(const Quadrilateral& quad){
		memcpy(v, quad.v, sizeof(int)*4);
		return *this;
    }
};

void GenerateRandomQuadIndex(int size, Quadrilateral* quad);

// void ReorderAsAnticlockwise(const util::point2d& pi, util::point2d& p1, util::point2d& p2, util::point2d& p3);

/**
 *  QuadrilateralInvariant: compute the invariants of Quadrilateral
 *  p1, p2, q1, q2
 */
bool QuadrilateralInvariant(const util::point2d& p1, const util::point2d& p2, const util::point2d& q1, 
						   const util::point2d& q2, int* index, float* invariant1, float* invariant2, float* baseline1, float* baseline2);

}

std::ostream& operator<<(std::ostream& o, const fpcs::Quadrilateral& quad);

#endif