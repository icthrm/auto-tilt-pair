#ifndef DATAFORMAT_H__
#define DATAFORMAT_H__

#include <vector>
#include <cstdio>
#include <cassert>
#include <algorithm>
#include <string>
#include <climits>
#include <cmath>
#include <cxcore.hpp>
#include <cv.hpp>
#include "keypoint.h"
#include "util/exception.h"

namespace util {

typedef _point point2d;

inline float L2(const point2d& p1, const point2d& p2)
{
	float difx = p1.x-p2.x;
	float dify = p1.y-p2.y;
	
	return difx*difx+dify*dify;
}

inline float CalTriangleArea(const util::point2d& p1, const util::point2d& p2, const util::point2d& p3)
{
    float s = ((p3.x-p1.x)*(p2.y-p1.y) - (p2.x-p1.x)*(p3.y-p1.y))*.5;

    return s>0 ? s : -s;
}


inline bool IsPolygonConvex(const util::point2d& p1, const util::point2d& p2, const util::point2d& p3, const util::point2d& p4)
{
	float z1, z2, z3, z4;
 
	z1 = ((p2.x - p1.x) * (p4.y - p1.y) - (p4.x - p1.x) * (p2.y - p1.y));
	z2 = ((p4.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p4.y - p1.y));
	z3 = ((p4.x - p2.x) * (p3.y - p2.y) - (p3.x - p2.x) * (p4.y - p2.y));
	z4 = ((p3.x - p2.x) * (p1.y - p2.y) - (p1.x - p2.x) * (p3.y - p2.y));
 
	return (z1 * z2 > 0) && (z3 * z4 > 0);
}

inline float CalQuadrangleArea(const util::point2d& p1, const util::point2d& p2, const util::point2d& p3, const util::point2d& p4)
{
    return (CalTriangleArea(p1, p2, p3)+CalTriangleArea(p2, p3, p4)+CalTriangleArea(p1, p2, p4)+CalTriangleArea(p1, p3, p4))*.5;
}

/** @brief Stack of fiducial markers*/
class FiducialStack {
private:
    int size;
    std::vector<util::point2d>* vfidp;
    int width, height;
    float ratio;

public:
    FiducialStack():size(0), vfidp(NULL), width(-1), height(-1), ratio(1) {}

    FiducialStack(int _size):size(_size) {
        vfidp = new std::vector<util::point2d>[size];
    }

    ~FiducialStack() {
        if(vfidp) {
            delete [] vfidp;
        }
    }

    void ReSize(int  _size) {
        if(size != _size) {
            delete [] vfidp;
            size = _size;
            vfidp = new std::vector<util::point2d>[size];
        }
    }
    void SetRatio(float r) {
        ratio = r;
    }

    void Clear() {
        delete [] vfidp;
        vfidp = new std::vector<util::point2d>[size];
    }

    void SetWxH(int _width, int _height) {
        width = _width;
        height = _height;
    }

    void Release() {
        if(vfidp) {
            delete [] vfidp;
        }
        vfidp = NULL;
    }

    std::vector<util::point2d>& V(int index) {
        assert(index>=0&&index<size);
        return vfidp[index];
    }
    const std::vector<util::point2d>& V(int index) const {
        assert(index>=0&&index<size);
        return vfidp[index];
    }

    int Size() const {
        return size;
    }
    int Width() const {
        return width;
    }
    int Height() const {
        return height;
    }
    float Ratio() const {
        return ratio;
    }

    void WriteFidsByFile(const char* filename) const;
    bool ReadFidsByFile(const char* filename);
};

typedef std::pair<util::_point, util::_point> pair;

struct img_match {
    int idx1, idx2;		//the index of first image and sencond image; which is matched
    std::vector<pair> pairs;
    size_t size() const {
        return pairs.size();
    }
};

/** @brief MatchPairStack; the output of keymatch process*/
class ImgMatchVector {
private:
    std::vector<img_match>* match_vector;
public:
    ImgMatchVector()
    {
        match_vector = new std::vector<img_match>();
    }

    ~ImgMatchVector()
    {
        if(match_vector) {
            delete match_vector;
        }
    }

    void Clear()
    {
        match_vector->clear();
    }

    void Release()
    {
        if(match_vector) {
            delete match_vector;
        }
        match_vector = NULL;
    }

    int Size() const {
        return match_vector->size();
    }

    img_match& operator[](int index)
    {
        assert(index>=0&&index<match_vector->size());
        return (*match_vector)[index];
    }

    const img_match& operator[](int index) const
    {
        assert(index>=0&&index<match_vector->size());
        return (*match_vector)[index];
    }

    img_match& V(int index) {
        assert(index>=0&&index<match_vector->size());
        return (*match_vector)[index];
    }
    const img_match& V(int index) const {
        assert(index>=0&&index<match_vector->size());
        return (*match_vector)[index];
    }

    img_match& MallocNewMatch()
    {
        img_match data;
        match_vector->push_back(data);
        return (*match_vector)[match_vector->size()-1];
    }

    void PushBack(const img_match& data)
    {
        match_vector->push_back(data);
    }
    
    img_match* GetMatchSetWithIdx(int idx1, int idx2, bool& noex){							//!!!!!!!!!!!!!!!!!!!!!!!!!!!
		for(int i = 0; i < match_vector->size(); i++){
			if((*match_vector)[i].idx1 == idx1 && (*match_vector)[i].idx2 == idx2){
				noex = true;
				return &((*match_vector)[i]);
			}
			if((*match_vector)[i].idx1 == idx2 && (*match_vector)[i].idx2 == idx1){
				noex = false;
				return &((*match_vector)[i]);
			}
		}
		return NULL;
	}

    static void ReadPairs(img_match* pairs, std::istream& in);
    void CoordinateTransform(int width, int height);
    void PreRotate(float angle);		//angle in degree
    void PrintPairs(int index, std::ostream& o) const;
    void WriteVectorByFolder(const char* folderpath) const;
    void ReadVectorByFolder(const char* folderpath);
};

#define ISEQUAL(x, y)		fabs((x)-(y)) < 0.001

static bool ReadAnglesByName(const char* name, std::vector<float>* angles)
{
    std::ifstream in(name);
    if(!in.good()) {
        return false;
    }

    while(in.good()) {
        float angle;
        in>>angle;
        if(in.fail()) {
            break;
        }
        angles->push_back(angle);
    }
    in.close();
    return true;
}

}
#endif
