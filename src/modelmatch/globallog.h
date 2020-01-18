#ifndef GLOBALLOG__
#define GLOBALLOG__

// #define __DEBUG__

#ifdef __DEBUG__
#define TRACELOG
#endif

#ifdef TRACELOG
#include <dataf/dataf.h>
#include <iostream>

class GLog
{
private:
	std::vector<util::point2d> unmatched_fids1;
	std::vector<util::point2d> unmatched_fids2;
	
	std::vector<std::pair<util::point2d, util::point2d> > tmatches;		//transformed matches
public:
	void PushBackTransformedPairs(const std::pair<util::point2d, util::point2d>& match)
	{
		tmatches.push_back(match);
	}
	
	void PushBackTransformedPairs(const util::point2d& pt1, const util::point2d& pt2)
	{
		tmatches.push_back(std::make_pair(pt1, pt2));
	}
	
	void PushBackUnmatchedFids1(const util::point2d& pt){
		unmatched_fids1.push_back(pt);
	}
	
	void PushBackUnmatchedFids2(const util::point2d& pt){
		unmatched_fids2.push_back(pt);
	}
	
	void PushBackUnmatchedFids1(const std::vector<util::point2d>& pts){
		for(int i = 0; i < pts.size(); i++){
			unmatched_fids1.push_back(pts[i]);
		}
	}
	
	void PushBackUnmatchedFids2(const std::vector<util::point2d>& pts){
		for(int i = 0; i < pts.size(); i++){
			unmatched_fids2.push_back(pts[i]);
		}
	}
	
	size_t MSize(){
		return tmatches.size();
	}
	
	void PrintTMPairs(std::ostream& o=std::cout) const
	{
		for(int i = 0; i < tmatches.size(); i++) {
			o<<tmatches[i].first.x<<" "<<tmatches[i].first.y<<"\t\t"<<tmatches[i].second.x<<" "<<tmatches[i].second.y<<"\n";
		}
		o<<std::endl;
	}
	
	void PrintFids1(std::ostream& o=std::cout) const
	{
		for(int i = 0; i < tmatches.size(); i++) {
			o<<tmatches[i].first.x<<" "<<tmatches[i].first.y<<"\n";
		}
		for(int i = 0; i < unmatched_fids1.size(); i++){
			o<<unmatched_fids1[i].x<<" "<<unmatched_fids1[i].y<<"\n";
		}
		o<<std::endl;
	}
	
	void PrintFids2(std::ostream& o=std::cout) const
	{
		for(int i = 0; i < tmatches.size(); i++) {
			o<<tmatches[i].second.x<<" "<<tmatches[i].second.y<<"\n";
		}
		for(int i = 0; i < unmatched_fids2.size(); i++){
			o<<unmatched_fids2[i].x<<" "<<unmatched_fids2[i].y<<"\n";
		}
		o<<std::endl;
	}
};

extern GLog glog;

#endif

#endif