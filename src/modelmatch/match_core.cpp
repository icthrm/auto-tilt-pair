#include "match_core.h"
#include <opencv/cv.h>
#include <iomanip>
#include <sys/stat.h>
#include <unistd.h>
#include <levmar.h>
#include "util/matrix.h"

void PrintMat(CvMat * A)
{
    for(int i = 0; i < A->rows; i++){
        for(int j = 0; j < A->cols; j++){
            //if(CV_MAT_DEPTH(A->type)==CV_64F)
			std::cout<<cvGetReal2D(A,i,j)<<"\t";
		}
		std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void ModelMatch::PairMatch(std::vector< util::point2d >& __fids1, std::vector< util::point2d >& __fids2, double beta1, double beta2, double dist_err_tol, 
						   std::vector<std::pair<util::point2d, util::point2d> >& matchvec, std::vector<mx::Matrix<2, 3, double> >& hset)
{
	Ran4PEstimator ran4PEstimator(__fids1, __fids2);
	ran4PEstimator.SetSupplementData(beta1, beta2);
	ran4PEstimator.AffineTransformEstimation(dist_err_tol, matchvec, hset, false, true);	//ran4p use half diameter
	
	if(!ran4PEstimator.UmFids1().size()|!ran4PEstimator.UmFids2().size()){		//there are no obvious warp
	}
	else{
		CPDEstimator cpdEstimator(ran4PEstimator.UmFids2(), ran4PEstimator.UmFids1());		//the first is fixed and the second is moving; therefore, it is inverse with 4p method
		cpdEstimator.PointDriftEstimation(dist_err_tol, matchvec, hset, false, false);		//CPD use full diameter
	}
}

void ModelMatch::MatchMain(util::FiducialStack& fstack, const std::vector<float>& angles, util::ImgMatchVector* imvector, HSetVector* hsetvec, float dist_err_tol, bool eigenlimit, bool do_test)
{
	EX_TRACE("Used distance threshold = %.2f\n", dist_err_tol)
	
	#define STEP_ARRAY_SIZE			2
    assert(STEP_ARRAY_SIZE < fstack.Size());

    imvector->Clear();
	hsetvec->Clear();
	
    int step_length[STEP_ARRAY_SIZE];
    for(int i = 0; i < STEP_ARRAY_SIZE; i++){
        step_length[i] = i+1;
    }
	
    int turn = 0;
	
	int sampling = fstack.Width()*fstack.Height()/(dist_err_tol*dist_err_tol*4)*.5;
	
    while(turn < STEP_ARRAY_SIZE){
        for(int i = 0; i+step_length[turn] < fstack.Size(); i++){
            int idx1 = i;//14;//6;//60;//
            int idx2 = i+step_length[turn];//16;//105;//
            
			if(fstack.V(idx1).size() > fstack.V(idx2).size()){		//this strategy is for Ran4PEstimator
				int tmp = idx1;
				idx1 = idx2;
				idx2 = tmp;
			}
			
            EX_TIME_BEGIN("#\nMatching Point Set (MRC[%d] & MRC[%d])", idx1, idx2)
            util::img_match& imatch = imvector->MallocNewMatch();
            imatch.idx1 = idx1;
            imatch.idx2 = idx2;
			
			h_set& hset = hsetvec->MallocNewHSet();
			hset.idx1 = idx1;
			hset.idx2 = idx2;
			
			PairMatch(fstack.V(idx1), fstack.V(idx2), angles[idx1], angles[idx2], dist_err_tol, imatch.pairs, hset.h);
			
			//if(imatch.size() < 5 || (imatch.size() < fstack.V(idx1).size()*CPDEstimator::COVERAGE_REFRESH_RATIO*.5 &&  
				//imatch.size() < fstack.V(idx2).size()*CPDEstimator::COVERAGE_REFRESH_RATIO*.5) ){
			if(imatch.size() < 5 || (imatch.size() < 9 && imatch.size() < fstack.V(idx1).size()*CPDEstimator::COVERAGE_REFRESH_RATIO*.25 &&  
				imatch.size() < fstack.V(idx2).size()*CPDEstimator::COVERAGE_REFRESH_RATIO*.25) ){
				imatch.pairs.clear();
				hset.h.clear();
			}

            EX_TIME_END("MRC[%d] & MRC[%d]: %ld/(%ld,%ld) Pairs found", idx1, idx2, imatch.size(), fstack.V(idx1).size(), fstack.V(idx2).size())
// 			return;
        }
        turn++;
    }
}

static void ImagePairs2Canvas(const IplImage* img1, const IplImage* img2, IplImage** canvas)
{
	*canvas = cvCreateImage(cvSize(img1->width+img2->width, img1->height), IPL_DEPTH_8U, 3);
	cvSetZero(*canvas);
	cvSetImageROI(*canvas, cvRect(0, 0, img1->width, img1->height));
    cvAdd(img1, *canvas, *canvas, NULL);
    cvSetImageROI(*canvas, cvRect(img1->width, 0, img1->width+img2->width, img2->height));
    cvAdd(img2, *canvas, *canvas, NULL);
    cvResetImageROI(*canvas);
}

void ModelMatch::Test(util::MrcStack& mrcr, const util::FiducialStack& fidsk, float diameter, const util::ImgMatchVector& imvector, float ratio, const char* folder)
{
    EX_TIME_BEGIN("%sMatch Testing", _DASH)

	diameter = diameter*ratio;
	
	int r = int(diameter+.5), t = int(.1*r+.5);if(t <= 0){t = 1;}
	CvScalar red = CV_RGB(255, 0, 0), green = CV_RGB(0, 255, 0);
	
    for(int i = 0; i < imvector.Size(); i++){
		IplImage* p[2], *tmp;
        tmp = mrcr.GetIplImage(imvector[i].idx1);
        p[0] = cvCreateImage(cvSize(tmp->width*ratio, tmp->height*ratio), IPL_DEPTH_32F, 1);
        cvResize(tmp, p[0], CV_INTER_CUBIC);
        cvReleaseImage(&tmp);
		util::ConvertGrayMRC2RGB(p[0], &tmp);
		cvReleaseImage(&p[0]); p[0] = tmp; 

		const std::vector<util::point2d>& fids1 = fidsk.V(imvector[i].idx1);
		for(int k = 0; k < fids1.size(); k++){
			cvCircle(p[0], cvPoint(int(fids1[k].x*ratio+.5), int(fids1[k].y*ratio+.5)), r, red, t);
		}
		
        tmp = mrcr.GetIplImage(imvector[i].idx2);
        p[1] = cvCreateImage(cvSize(tmp->width*ratio, tmp->height*ratio), IPL_DEPTH_32F, 1);
        cvResize(tmp, p[1], CV_INTER_CUBIC);
        cvReleaseImage(&tmp);
		util::ConvertGrayMRC2RGB(p[1], &tmp);
		cvReleaseImage(&p[1]); p[1] = tmp;
		
		const std::vector<util::point2d>& fids2 = fidsk.V(imvector[i].idx2);
		for(int k = 0; k < fids2.size(); k++){
			cvCircle(p[1], cvPoint(int(fids2[k].x*ratio+.5), int(fids2[k].y*ratio+.5)), r, red, t);
		}
		
		IplImage* canvas;
		ImagePairs2Canvas(p[0], p[1], &canvas);

		std::vector<std::pair<util::point2d, util::point2d> > vpair;
		for(int j = 0; j < imvector[i].pairs.size(); j++){
			util::point2d pair1, pair2;
			pair1.x = imvector[i].pairs[j].first.x*ratio;
			pair1.y = imvector[i].pairs[j].first.y*ratio;
			pair2.x = imvector[i].pairs[j].second.x*ratio;
			pair2.y = imvector[i].pairs[j].second.y*ratio;
			vpair.push_back(std::make_pair(pair1, pair2));
		}
		
        for(int i = 0; i < vpair.size(); i++){
			util::DrawCircle(canvas, vpair[i].first.x, vpair[i].first.y, r, green, t+1);
			util::DrawCircle(canvas, p[1]->width+vpair[i].second.x, vpair[i].second.y, r, green, t+1);
			util::DrawLine(canvas, vpair[i].first, util::point2d(p[1]->width+vpair[i].second.x, vpair[i].second.y), green, t+1);
		}
		
        if(access(folder,0) == -1) {		//create file folder
            mkdir(folder,0777);
        }
        std::ostringstream oss;
        oss <<folder<<"/"<<"("<<imvector[i].idx1<<")&("<<imvector[i].idx2<<").png";
        try {
//             util::SaveImage(canvas, oss.str().c_str());
			std::cout<<"Saving "<<oss.str()<<std::endl;
			cvSaveImage(oss.str().c_str(), canvas);
        } catch(ex::Exception& e) {
            EX_TRACE("%s\n", e.Msg())
        }
        cvReleaseImage(&p[0]);
        cvReleaseImage(&p[1]);
		cvReleaseImage(&canvas);
    }
    EX_TIME_END("Match Testing")
}

void HSetVector::WriteTransforms(const h_set& hset, std::ostream& out)
{
	out<<hset.idx1<<"\t"<<hset.idx2<<"\t"<<hset.size()<<std::endl;
	for(int i = 0; i < hset.size(); i++){
		for(int m = 0; m < 2; m++){
			for(int n = 0; n < 3; n++){
				out<<hset.h[i].V(m,n)<<"\t";
			}
			out<<std::endl;
		}
		out<<std::endl<<std::endl;
	}
}

void HSetVector::ReadTransforms(h_set* hset, std::istream& in)
{
	int hsize;
	in>>hset->idx1>>hset->idx2>>hsize;
	for(int i = 0; i < hsize; i++){
		mx::Matrix<2, 3, double> h;
		for(int m = 0; m < 2; m++){
			for(int n = 0; n < 3; n++){
				in>>h.V(m,n);
			}
		}
		hset->h.push_back(h);
	}
}

void HSetVector::WriteVectorByFolder(const char* folderpath) const
{
    if(access(folderpath,0) == -1) {		//create file folder
        mkdir(folderpath,0777);
    }
    std::ostringstream ooo;
    ooo <<folderpath<<"/attributes";
    std::ofstream out(ooo.str().c_str());
    out<<"Z:"<<hset_vector->size()<<"\n";
    out.close();
    for(int i = 0; i < hset_vector->size(); i++) {
        std::ostringstream oss;
        oss <<folderpath<<"/"<<i<<".txt";
        std::ofstream o(oss.str().c_str());
        try {
            WriteTransforms((*hset_vector)[i], o);
        } catch(ex::Exception& e) {
            EX_TRACE("%s\n", e.Msg())
        }
    }
}

void HSetVector::ReadVectorByFolder(const char* folderpath)
{
    Clear();
    std::cout <<std::setprecision(8)<<std::endl;
    std::ostringstream iii;
    iii <<folderpath<<"/attributes";
    std::ifstream in(iii.str().c_str());
    char ch;
    int _size;
    in>>ch>>ch>>_size;
    in.close();

    for(int i = 0; i < _size; i++) {
        std::ostringstream oss;
        oss <<folderpath<<"/"<<i<<".txt";
        std::ifstream in(oss.str().c_str());
        if(!in.good()) {
            ex::EX_THROW("Can't Open File");
        }
        h_set& hset = MallocNewHSet();
        ReadTransforms(&hset, in);
        in.close();
    }
}
