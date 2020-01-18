#ifndef MODEL_MATCH_CORE_H__
#define MODEL_MATCH_CORE_H__

#include "dataf/dataf.h"
#include <vector>
#include <ANN/ANN.h>
#include "ransac/xform.h"
#include "4pcs.h"
#include "r4pestimator.h"
#include "cpdestimator.h"
#include "util/exception.h"
#include "mrcimg/mrc2img.h"
#include "mrcimg/img_util.h"

#define MIN_PT_NUM		6
#define DIS_ERR_TOL		5
#define DIS_ERR_TOL4096		10
#define SUCCESS_TOL		0.999999

struct h_set{
    int idx1, idx2;		//the index of first image and sencond image; which is matched
    std::vector<mx::Matrix<2, 3, double> > h;
	h_set():idx1(-1),idx2(-1){}
	
	h_set(const h_set& __h):idx1(__h.idx1), idx2(__h.idx2){
		h = __h.h;
	}
	
    ~h_set(){
		h.clear();
    }
    size_t size() const {
        return h.size();
    }
};

/** @brief Transform matrix vector; aid for the the output of keymatch process*/
class HSetVector {
private:
    std::vector<h_set>* hset_vector;
public:
    HSetVector()
    {
        hset_vector = new std::vector<h_set>();
    }

    ~HSetVector()
    {
        if(hset_vector) {
            delete hset_vector;
        }
    }

    void Clear()
    {
        hset_vector->clear();
    }

    void Release()
    {
        if(hset_vector) {
            delete hset_vector;
        }
        hset_vector = NULL;
    }
	
	h_set* GetHSetWithIdx(int idx1, int idx2){
		for(int i = 0; i < hset_vector->size(); i++){
			if((*hset_vector)[i].idx1 == idx1 && (*hset_vector)[i].idx2 == idx2){
				return &((*hset_vector)[i]);
			}
		}
		return NULL;
	}
	
    int Size() const {
        return hset_vector->size();
    }

    h_set& operator[](int index)
    {
        assert(index>=0&&index<hset_vector->size());
        return (*hset_vector)[index];
    }

    const h_set& operator[](int index) const
    {
        assert(index>=0&&index<hset_vector->size());
        return (*hset_vector)[index];
    }

    h_set& V(int index) {
        assert(index>=0&&index<hset_vector->size());
        return (*hset_vector)[index];
    }
    
    const h_set& V(int index) const {
        assert(index>=0&&index<hset_vector->size());
        return (*hset_vector)[index];
    }

    h_set& MallocNewHSet()
    {
        h_set data;
        hset_vector->push_back(data);
        return (*hset_vector)[hset_vector->size()-1];
    }

    void PushBack(const h_set& data)
    {
        hset_vector->push_back(data);
    }

    static void ReadTransforms(h_set* hset, std::istream& in);
	static void WriteTransforms(const h_set& hset, std::ostream& out);
    void WriteVectorByFolder(const char* folderpath) const;
    void ReadVectorByFolder(const char* folderpath);
};

class ModelMatch{
private:
	static void PairMatch(std::vector<util::point2d>& __fids1, std::vector<util::point2d>& __fids2, double beta1, double beta2, double dist_err_tol, 
						  std::vector<std::pair<util::point2d, util::point2d> >& matchvec, std::vector<mx::Matrix<2, 3, double> >& hset);
public:	
	static void MatchMain(util::FiducialStack& fstack, const std::vector<float>& angles, util::ImgMatchVector* imvector, HSetVector* hset, float dist_err_tol, bool eigenlimit = true, bool do_test = false); //dist_err_tol = half radius
	static void Test(util::MrcStack& mrcr, const util::FiducialStack& fidsk, float diameter, const util::ImgMatchVector& imvector, float ratio, const char* folder = "matches_ill");
};

#endif
