#include "dataf.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <iostream>
#include "keypoint.h"
#include <sys/stat.h>
#include <unistd.h>
#include <iomanip>

void util::FiducialStack::WriteFidsByFile(const char* filename) const {
    std::ofstream out(filename);
    if(!out.good()) {
        ex::EX_THROW("Can't Create File");
    }
    std::cout <<std::setprecision(8)<<std::endl;
    out<<width<<"\t"<<height<<"\t"<<size<<"\t"<<ratio<<std::endl;
    for(int i = 0 ; i < size ; i++) {
        out<<"frame "<<i<<std::endl;
        out<<"num = "<<vfidp[i].size() <<std::endl;
        for(int j = 0 ; j < vfidp[i].size() ; j++) {
            float x = vfidp[i][j].x;
            float y = vfidp[i][j].y;
            out<<x<<"\t"<<y<<std::endl;
        }
    }
    return;
}

bool util::FiducialStack::ReadFidsByFile(const char* filename) {
    std::string s;
    std::stringstream ss;
    std::string first_str;
    char ch;

    std::ifstream fin(filename);

    if(!fin.good()) {
        std::cout<<"Unable to open "<<filename<<std::endl;
        return false;
    }

    getline(fin , s);

    ss.str("");
    ss<<s;
    ss>>width>>height>>size>>ratio;

    if(vfidp) {
        delete [] vfidp;
    }
    vfidp = new std::vector<util::point2d>[size];

    while(getline(fin ,s)) {
        ss.clear();
        ss<<s;
        ss>>first_str;
        int feats_num;
        if(first_str == "frame") {
            int frame;
            ss>>frame;
            getline(fin ,s);
            ss.clear();
            ss<<s;
            ss>>first_str;

            if(first_str == "num") {
                ss>>ch>>feats_num;
                for(int i=0; i< feats_num ; i++) {
                    util::point2d pt;
                    fin>>pt.x>>pt.y;
                    vfidp[frame].push_back(pt);
                }
            }
        }
    }
    return true;
}

void util::ImgMatchVector::ReadPairs(img_match* imatch, std::istream& in)
{
    char ch;
    in>>imatch->idx1>>imatch->idx2;
    while(in.peek()!='#' && in.good()) {
        pair p;
        in>>ch>>p.first.x>>ch>>p.first.y>>ch>>ch>>ch>>p.second.x>>ch>>p.second.y>>ch;
        imatch->pairs.push_back(p);
        in.ignore(4096, '\n');
    }
}


void util::ImgMatchVector::PrintPairs(int index, std::ostream& o) const
{
    o<<(*match_vector)[index].idx1<<" "<<(*match_vector)[index].idx2<<std::endl;
    for(int i = 0; i < (*match_vector)[index].pairs.size(); i++) {
        o<<"("<<((*match_vector)[index].pairs)[i].first.x<<","<<((*match_vector)[index].pairs)[i].first.y<<")&"
         <<"("<<((*match_vector)[index].pairs)[i].second.x<<","<<((*match_vector)[index].pairs)[i].second.y<<")\n";
    }
    o<<"#";
}

void util::ImgMatchVector::CoordinateTransform(int width, int height)
{
    for(int index = 0; index < (*match_vector).size(); index++) {
        for(int i = 0; i < (*match_vector)[index].pairs.size(); i++) {
            ((*match_vector)[index].pairs)[i].first.x -= width/2;
            ((*match_vector)[index].pairs)[i].first.y -= height/2;
            ((*match_vector)[index].pairs)[i].second.x -= width/2;
            ((*match_vector)[index].pairs)[i].second.y -= height/2;
        }
    }
}

void util::ImgMatchVector::PreRotate(float angle)
{
    for(int index = 0; index < (*match_vector).size(); index++) {
        for(int i = 0; i < (*match_vector)[index].pairs.size(); i++) {
            ((*match_vector)[index].pairs)[i].first.x = cos(angle)*((*match_vector)[index].pairs)[i].first.x-sin(angle)*((*match_vector)[index].pairs)[i].first.y;
            ((*match_vector)[index].pairs)[i].first.y = sin(angle)*((*match_vector)[index].pairs)[i].first.x+cos(angle)*((*match_vector)[index].pairs)[i].first.y;
            ((*match_vector)[index].pairs)[i].second.x = cos(angle)*((*match_vector)[index].pairs)[i].second.x-sin(angle)*((*match_vector)[index].pairs)[i].second.y;
            ((*match_vector)[index].pairs)[i].second.y = sin(angle)*((*match_vector)[index].pairs)[i].second.x+cos(angle)*((*match_vector)[index].pairs)[i].second.y;
        }
    }
}

void util::ImgMatchVector::WriteVectorByFolder(const char* folderpath) const
{
    if(access(folderpath,0) == -1) {		//create file folder
        mkdir(folderpath,0777);
    }
    std::ostringstream ooo;
    ooo <<folderpath<<"/attributes";
    std::ofstream out(ooo.str().c_str());
    out<<"Z:"<<match_vector->size()<<"\n";
    out.close();
    for(int i = 0; i < match_vector->size(); i++) {
        std::ostringstream oss;
        oss <<folderpath<<"/"<<i<<".txt";
        std::ofstream o(oss.str().c_str());
        try {
            PrintPairs(i, o);
        } catch(ex::Exception& e) {
            EX_TRACE("%s\n", e.Msg())
        }
    }
}

void util::ImgMatchVector::ReadVectorByFolder(const char* folderpath)
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
        img_match& imatch = MallocNewMatch();
        ReadPairs(&imatch, in);
        in.close();
    }
}
