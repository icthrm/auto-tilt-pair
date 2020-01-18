#include "mrc2img.h"

bool util::MrcStack::Open(const char* filename) {
    in.open(filename);
    if(!in.good()) {
        return false;
    }
    header = new MRCheader();

    in.read((char*) header, sizeof(MrcHeader));
//     if ( ! ( header->cmap[0]=='M'&&header->cmap[1]=='A'&&header->cmap[2]=='P' ) ) {
//         in.close();
//         delete header;
//         header = NULL;
//         return false;
//     }

    name = std::string(filename);
    return true;
}

IplImage* util::MrcStack::ReadIplImage(size_t index)
{
    assert(header != NULL && index >= 0 && index < header->nz);
    assert(in.good());

    const int& nx = header->nx;			//the coordinate is different(have not been tested)
    const int& ny = header->ny;

    CvSize size;
    size.width = nx;
    size.height = ny;
    IplImage *img = cvCreateImage(size, IPL_DEPTH_32F, 1);

    int bufsize = nx*ny;

    switch(header->mode) {
    case MRC_MODE_BYTE: {
        unsigned char* tmpbuf = new unsigned char[bufsize];
        for(size_t y = 0; y < ny; y++){
            in.seekg(sizeof(MrcHeader) +header->next+ (index*nx*ny+y*nx) *sizeof(unsigned char), std::ios::beg);
            in.read((char*) & (tmpbuf[y*nx]), nx*sizeof(unsigned char));
        }
        
        unsigned char* src = tmpbuf;
        for(size_t y = 0; y < img->height; y++){
			float* start = (float*)(img->imageData+y*img->widthStep);
			for(int x = 0; x < img->width; x++){
				*start++ = *src++/255.0f;
			}
		}
        
        delete [] tmpbuf;
        break;
    }
    case MRC_MODE_SHORT: {
        short* tmpbuf = new short[bufsize];
        for(size_t y = 0; y < ny; y++) {
            in.seekg(sizeof(MrcHeader) +header->next+ (index*nx*ny+y*nx) *sizeof(short), std::ios::beg);
            in.read((char*) & (tmpbuf[y*nx]),nx*sizeof(short));
        }
        
        short* src = tmpbuf;
        for(size_t y = 0; y < img->height; y++){
			float* start = (float*)(img->imageData+y*img->widthStep);
			for(int x = 0; x < img->width; x++){
				*start++ = *src++;
			}
		}
        
        delete [] tmpbuf;
        break;
    }
    case MRC_MODE_FLOAT: {
        for(size_t y = 0; y < ny; y++) {
            in.seekg(sizeof(MrcHeader)+header->next+(index*nx*ny+y*nx)*sizeof(float), std::ios::beg);
            in.read((char*)(img->imageData+y*img->widthStep),nx*sizeof(float));
        }

        break;
    }
    default:
        cvReleaseImage(&img);
        return NULL;
    }

    return img;
}

void util::MrcStack::DoCaching()
{
	if(!slices){
		slices = new IplImage*[header->nz];
	}
	for(int i = 0; i < header->nz; i++){
		slices[i] = ReadIplImage(i);
	}
	
	is_cashed = true;
}

void util::MrcStack::FreeCache()
{
	if(!slices){
		return;
	}
	
	for(int i = 0; i < header->nz; i++){
		cvReleaseImage(&(slices[i]));
	}
	delete [] slices;
	slices = NULL;
	is_cashed = false;
}

IplImage* util::MrcStack::GetIplImage(int index) {
	if(!is_cashed){
		return ReadIplImage(index);
	}
	else{
		return cvCloneImage(slices[index]);
	}
}

IplImage* const& util::MrcStack::operator[](int idx)
{
	return slices[idx];
}

IplImage* const& util::MrcStack::operator[](int idx) const
{
	return slices[idx];
}

void util::MrcStack::PrintHeader(std::ostream& o) const {
    o<<"\n nx: "<<header->nx;         /*  # of Columns                  */
    o<<"\n ny: "<<header->ny;         /*  # of Rows                     */
    o<<"\n nz: "<<header->nz;         /*  # of Sections.                */
    o<<"\n mode: "<<header->mode;       /*  given by #define MRC_MODE...  */

    o<<"\n nxstart: "<<header->nxstart;    /*  Starting point of sub image.  */
    o<<"\n nystart: "<<header->nystart;
    o<<"\n nzstart: "<<header->nzstart;

    o<<"\n mx: "<<header->mx;         /* Grid size in x, y, and z       */
    o<<"\n my: "<<header->my;
    o<<"\n mz: "<<header->mz;

    o<<"\n xlen: "<<header->xlen;       /* length of x element in um.     */
    o<<"\n ylen: "<<header->ylen;       /* get scale = xlen/nx ...        */
    o<<"\n zlen: "<<header->zlen;

    o<<"\n alpha: "<<header->alpha;      /* cell angles, ignore */
    o<<"\n beta: "<<header->beta;
    o<<"\n gamma: "<<header->gamma;

    o<<"\n mapc: "<<header->mapc;       /* map coloumn 1=x,2=y,3=z.       */
    o<<"\n mapr: "<<header->mapr;       /* map row     1=x,2=y,3=z.       */
    o<<"\n maps: "<<header->maps;       /* map section 1=x,2=y,3=z.       */

    o<<"\n amin: "<<header->amin;
    o<<"\n amax: "<<header->amax;
    o<<"\n amean: "<<header->amean;

    o<<"\n ispg: "<<header->ispg;       /* image type */
    o<<"\n nsymbt: "<<header->nsymbt;     /* space group number */


    /* 64 bytes */

    o<<"\n next: "<<header->next;
    o<<"\n sizeof header: "<<sizeof(MRCheader);
    o<<"\n creatid: "<<header->creatid;  /* Creator id, hvem = 1000, DeltaVision = -16224 */


    o<<"\n blank: "<<header->blank;

    o<<"\n nint: "<<header->nint;
    o<<"\n nreal: "<<header->nreal;
    o<<"\n sub: "<<header->sub;
    o<<"\n zfac: "<<header->zfac;

    o<<"\n min2: "<<header->min2;
    o<<"\n max2: "<<header->max2;
    o<<"\n min3: "<<header->min3;
    o<<"\n max3: "<<header->max3;
    o<<"\n min4: "<<header->min4;
    o<<"\n max4: "<<header->max4;


    o<<"\n idtype: "<<header->idtype;
    o<<"\n lens: "<<header->lens;
    o<<"\n nd1: "<<header->nd1;     /* Devide by 100 to get o<<header-> value. */
    o<<"\n nd2: "<<header->nd2;
    o<<"\n vd1: "<<header->vd1;
    o<<"\n vd2: "<<header->vd2;
    o<<"\n tiltangles: "<<header->tiltangles[0]<<" "<<header->tiltangles[1]<<" "
     <<header->tiltangles[2]<<" "<<header->tiltangles[3]<<" "<<header->tiltangles[4]<<" "
     <<header->tiltangles[5]<<" ";  /* 0,1,2 = original:  3,4,5 = current */


    o<<"\n xorg: "<<header->xorg;
    o<<"\n yorg: "<<header->yorg;
    o<<"\n zorg: "<<header->zorg;
    o<<"\n cmap: "<<header->cmap[0]<<header->cmap[1]<<header->cmap[2]<<header->cmap[3];
    o<<"\n stamp: "<<header->stamp[0]<<header->stamp[1]<<header->stamp[2]<<header->stamp[3];
    o<<"\n rms: "<<header->rms;

    for(int i = 0; i < header->nlabl; i++) {
        o<<"\n labels["<<i<<"]: "<<header->labels[i];
    }
    o.flush();
}


