#ifndef KEYPOINT_H__
#define KEYPOINT_H__

#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include "util/exception.h"

namespace util{

/** max feature descriptor length */
#define FEATURE_MAX_D 128
#define Round(X)	int((X)+0.5)

struct RECT{
	int left;
	int top;
	int right;
	int bottom;
};

struct _point{
	union{
	float v[2];
	struct{float x;
    float y;};
	};
    _point(float _x, float _y):x(_x), y(_y){}
    _point(){}
    bool operator != (const _point& pt_){return !(x == pt_.x && y == pt_.y);}
};

typedef struct{
    float x, y; 		///< x,y location in a picture(original)
    float sigma; 		///< scale sigma; to original picture
    float orient;		///< orientation of this key point
}keypt_t;

typedef struct{
    keypt_t kpt;				///< key point
    int d;                      		///< descriptor length
    unsigned char descr[FEATURE_MAX_D];  	///< descriptor 
}feature;

inline static float SquareDistance(const util::feature& f1, const util::feature& f2)
{
    assert(f1.d == f2.d);
    float tmp = 0;
    
    for(int i = 0; i < f1.d; i++){
	float diff = f1.descr[i]-f2.descr[i];
	tmp += diff*diff;
    }
    return tmp;
}

inline static void PrintKeyPoint(const keypt_t& kpt, std::ostream& o)
{
    o<<"["<<kpt.x<<"]["<<kpt.y<<"]:sigma("<<kpt.sigma<<") orient("<<kpt.orient<<")"<<std::endl;
}

inline static void InputKeyPoint(keypt_t* kpt, std::istream& in)
{
    char ch;
    in>>ch>>kpt->x>>ch>>ch>>kpt->y>>ch>>ch>>ch>>ch>>ch>>ch>>ch
    >>ch>>kpt->sigma>>ch>>ch>>ch>>ch>>ch>>ch>>ch>>ch>>kpt->orient>>ch;
}

inline static void PrintFeature(const feature& feat, std::ostream& o)
{
    PrintKeyPoint(feat.kpt, o);
    o<<"DIM"<<feat.d<<" ";
    for(int i = 0; i < feat.d; i++){
	o<<"["<<(float)feat.descr[i]<<"]";
    }
    o<<std::endl;
}

inline static void InputFeature(feature* feat, std::istream& in)
{
    char ch;
    InputKeyPoint(&(feat->kpt), in);
    in>>ch>>ch>>ch>>feat->d;
    for(int i = 0; i < feat->d; i++){
	in>>ch>>*((int*)&(feat->descr[i]))>>ch;
    }
}

inline static int WriteFeaturesByName(const std::vector<feature>& feats, const char* filename) throw(ex::Exception)
{
    std::ofstream out(filename); 
    if(!out.good()) {
        ex::EX_THROW("Can't Create File");
    }
    for(int i = 0; i < feats.size(); i++){
	PrintFeature(feats[i], out);
    }
    out<<"#";
    return feats.size();
}

inline static int ReadFeaturesByName(std::vector<feature>* feats, const char* filename) throw(ex::Exception)
{
    std::ifstream in(filename); 
    if(!in.good()) {
        ex::EX_THROW("Can't Open File");
    }

    while(in.peek()!='#' && in.good()){
	feature f;
	InputFeature(&f, in);
	feats->push_back(f);
	in.ignore(4096, '\n');
    }
   
    return feats->size();
}

}
#endif