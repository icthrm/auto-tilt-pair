#include "4pcs.h"
#include <cmath>
#include <numeric>
#include <matrix/matrix.h>

std::ostream& operator<<(std::ostream& o, const fpcs::Quadrilateral& quad)
{
	o<<"["<<quad[0]<<"]["<<quad[1]<<"]["<<quad[2]<<"]["<<quad[3]<<"]";
	return o;
}

void fpcs::GenerateRandomQuadIndex(int size, Quadrilateral* quad){
	(*quad)[0] = rand()%size;
	int idx = 1;
	while(true){
generatebeginrqervfe124334321:
		int tmp = rand()%size;
		for(int i = 0; i < idx; i++){
			if(tmp == (*quad)[i]){
				goto generatebeginrqervfe124334321;
			}
		}
		(*quad)[idx] = tmp;
		if(++idx >= 4){
			break;
		}
	}
}

struct ptvec
{
	util::point2d p;
	float angle;
	bool is_start;
	
	bool operator > (const ptvec& pv2) const{
        return this->angle > pv2.angle;
    }
    
    bool operator < (const ptvec& pv2) const{
        return this->angle < pv2.angle;
    }
};

// void ReorderAsAnticlockwise(const util::point2d& pi, util::point2d& p1, util::point2d& p2, util::point2d& p3)
// {
// 	ptvec pv[4];
// 	memset(pv, 0, sizeof(ptvec)*4);
// 	
// 	if(((p1.y-pi.y)*(p2.x-pi.x)-(p1.x-pi.x)*(p2.y-pi.y))*((p1.y-pi.y)*(p3.x-pi.x)-(p1.x-pi.x)*(p3.y-pi.y)) < 0){
// 		pv[1].p = p1;
// 		pv[0].p = p2; pv[2].p = p3;
// 	}
// 	else if(((p2.y-pi.y)*(p1.x-pi.x)-(p2.x-pi.x)*(p1.y-pi.y))*((p2.y-pi.y)*(p3.x-pi.x)-(p2.x-pi.x)*(p3.y-pi.y)) < 0){
// 		pv[1].p = p2;
// 		pv[0].p = p1; pv[2].p = p3;
// 	}
// 	else{
// 		pv[1].p = p3;
// 		pv[0].p = p1; pv[2].p = p2;
// 	}
// 	
// 	util::point2d center;
// 	INTERSECT_P(pi, pv[1].p, pv[0].p, pv[2].p, center);
// // 	GetIntersectionP(pi, pv[1].p, pv[0].p, pv[2].p, center);
// 	pv[3].p = pi;
// 	pv[3].is_start = true;
// 	
// 	for(int i = 0; i < 4; i++){
// 		pv[i].angle = atan2((pv[i].p.y-center.y), (pv[i].p.x-center.x));
// 		if(pv[i].angle < 0){
// 			pv[i].angle += 6.28318530717959;
// 		}
// 	}
// 	
// 	std::sort(pv, pv+4);
// 	
// 	if(pv[0].is_start){ p1 = pv[1].p; p2 = pv[2].p; p3 = pv[3].p;}
// 	else if(pv[1].is_start){ p1 = pv[2].p; p2 = pv[3].p; p3 = pv[0].p;}
// 	else if(pv[2].is_start){ p1 = pv[3].p; p2 = pv[0].p; p3 = pv[1].p;}
// 	else{ p1 = pv[0].p; p2 = pv[1].p; p3 = pv[2].p;}
// }

bool fpcs::QuadrilateralInvariant(const util::point2d& p1, const util::point2d& p2, const util::point2d& q1, 
								  const util::point2d& q2, int* index, float* invariant1, float* invariant2, float* baseline1, float* baseline2){
    if(invariant1 == NULL || invariant2 == NULL || index == NULL){
        return false;
	}
// 	std::cout<<"("<<p1.x<<","<<p1.y<<")"<<"("<<p2.x<<","<<p2.y<<")"<<"("<<q1.x<<","<<q1.y<<")"<<"("<<q2.x<<","<<q2.y<<")"<<std::endl;
	
	float p1p2x = (p1.x-p2.x), p1p2y = (p1.y-p2.y);
	float p1q1x = (p1.x-q1.x), p1q1y = (p1.y-q1.y);
	float p1q2x = (p1.x-q2.x), p1q2y = (p1.y-q2.y);
	float p2q1x = (p2.x-q1.x), p2q1y = (p2.y-q1.y);  //
	float p2q2x = (p2.x-q2.x), p2q2y = (p2.y-q2.y);  
	float q1q2x = (q1.x-q2.x), q1q2y = (q1.y-q2.y);

	//LINEp1p2(q1)
	float lp1p2q1 = (-p1p2y*p1q1x+p1q1y*p1p2x);
	//LINEp1p2(q2)
	float lp1p2q2 = (-p1p2y*p1q2x+p1q2y*p1p2x);
	
	//LINEq1q2(p1)
	float lq1q2p1 = (q1q2y*p1q1x-p1q1y*q1q2x);
	//LINEq1q2(p2)
	float lq1q2p2 = (q1q2y*p2q1x-p2q1y*q1q2x);
	
	
	//LINEp1q1(p2)
	float lp1q1p2 = -lp1p2q1;
	//LINEp1q1(q2)
	float lp1q1q2 = (-p1q1y*p1q2x+p1q2y*p1q1x);
	
	//LINEp2q2(p1)
	float lp2q2p1 = (p2q2y*p1p2x-p1p2y*p2q2x);
	//LINEp2q2(q1)
	float lp2q2q1 = (-p2q2y*p2q1x+p2q1y*p2q2x);
		
	//LINEp1q2(p2)
	float lp1q2p2 = -lp1p2q2;
	//LINEp1q2(q1)
	float lp1q2q1 = -lp1q1q2;
	
	//LINEp2q1(p1)
	float lp2q1p1 = (p2q1y*p1p2x-p1p2y*p2q1x);
	//LINEp2q1(q2)
	float lp2q1q2 = -lp2q2q1;

// #define LINEab(c,c)		\
// 	(-aby*acx+acy*abx)
// #define LINEab(d,d)		\
// 	(-aby*adx+ady*abx)
	
	
#define INTERSECT_P(a, b, c, d, e)		\
	{e.y = ((c.y-d.y)*(a.y-b.y)*(c.x-a.x)+(c.y-d.y)*(a.x-b.x)*a.y-(a.y-b.y)*(c.x-d.x)*c.y)/((c.y-d.y)*(a.x-b.x)-(a.y-b.y)*(c.x-d.x)); \
	e.x = (c.x-d.x)*(e.y-c.y)/(c.y-d.y)+c.x;}
#define NORMAMINSB(a, b, tmp)	\
	{float tmk;	\
	tmk = (a.x-b.x); tmk = tmk*tmk;		\
	tmp = (a.y-b.y); tmp = tmp*tmp;	tmp += tmk;	\
	}
#define INVARIANT(a, b, c, d, e, inva1, inva2, basel1, basel2)		\
	{float s1, s2, t1, t2; float tmp1, tmp2;	\
	NORMAMINSB(a, e, tmp1); NORMAMINSB(a, b, tmp2); inva1 = sqrt(tmp1/tmp2); basel1 = sqrt(tmp2);	\
	NORMAMINSB(c, e, tmp1); NORMAMINSB(c, d, tmp2); inva2 = sqrt(tmp1/tmp2); basel2 = sqrt(tmp2);	\
	}
#define REASSIGNINDEX(index, i1, i2, i3, i4)	\
	index[0] = i1; index[1] = i2; index[2] = i3; index[3] = i4;

	
	util::point2d e;
	float s1, s2, t1, t2;
	
	if(lp1p2q1*lp1p2q2 < 0 && lq1q2p1*lq1q2p2 < 0){		//p1-p2;q1-q2
		INTERSECT_P(p1, p2, q1, q2, e);
		INVARIANT(p1, p2, q1, q2, e, *invariant1, *invariant2, *baseline1, *baseline2);
		REASSIGNINDEX(index, 0, 1, 2, 3)
	}
	else if(lp1q1p2*lp1q1q2 < 0 && lp2q2p1*lp2q2q1 < 0){		//p1-q1;p2-q2
		INTERSECT_P(p1, q1, p2, q2, e);
		INVARIANT(p1, q1, p2, q2, e, *invariant1, *invariant2, *baseline1, *baseline2);
		REASSIGNINDEX(index, 0, 2, 1, 3)
	}
	else if(lp1q2p2*lp1q2q1 < 0 && lp2q1p1*lp2q1q2 < 0){		//p1-q2;p2-q1
		INTERSECT_P(q2, p1, q1, p2, e);
		INVARIANT(q2, p1, q1, p2, e, *invariant1, *invariant2,  *baseline1, *baseline2);
		REASSIGNINDEX(index, 0, 3, 1, 2)
	}
	else{
		*invariant1 = -1;
		*invariant2 = -1;
		return false;
	}

#undef INTERSECT_P
#undef NORMAMINSB
#undef INVARIANT
#undef REASSIGNINDEX

    return true;
}