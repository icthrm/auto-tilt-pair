#include <cxcore.hpp>
#include <cv.hpp>

union pp{
	struct {	
double x;
double y;
double z;
	};
	double p[3];
};

// 
// double* independent_variables;
// double* oberved_values;
// 
int main(int argc, char* argv[])
{
	pp t;
	t.x = 1;
	t.y = 2;
	t.z = 3;
	std::cout<<t.p[0]<<"\t"<<t.p[1]<<"\t"<<t.p[2]<<"\t"<<std::endl;;
//     int m=10; //no. of observation variables
//     int n=3; //no. of paremeters 
//     independent_variables = new double[m]; //this is the vector of independent variables 
//     oberved_values = new double[m]; //this is the vector of observed variables to be fitted
//     for (int i=0;i<m;i++)
//     {
//       independent_variables[i]=i*2+5;
//       double xi=independent_variables[i];
//       double yi=16*xi*xi*xi + 12*xi*xi + 15;
//       oberved_values[i]=yi;
//     }
//     double* x=new double[n]; //initial estimate of parameters vector
//     x[0]=0.1;
//     x[1]=0.1;
//     x[2]=0.1;
// 
//     double* fvec=new double[m]; //no need to populate 
//     double ftol=1e-08; //tolerance
//     double xtol=1e-08; //tolerance
//     double gtol=1e-08; //tolerance
//     int maxfev=400; //maximum function evaluations
//     double epsfcn=1e-08; //tolerance
//     double* diag=new double[n]; //some internal thing
//     int mode=1; //some internal thing
//     double factor=1; // a default recommended value
//     int nprint=0; //don't know what it does
//     int info=0; //output variable
//     int nfev=0; //output variable will store no. of function evals
//     double* fjac=new double[m*n]; //output array of jacobian
//     int ldfjac=m; //recommended setting
//     int* ipvt=new int[n]; //for internal use
//     double* qtf=new double[n]; //for internal use
//     double* wa1=new double[n]; //for internal use
//     double* wa2=new double[n]; //for internal use
//     double* wa3=new double[n]; //for internal use
//     double* wa4=new double[m]; //for internal use
// 
//     lmdif_( m, n, x, fvec, ftol,
//            xtol, gtol, maxfev, epsfcn, diag,  mode,  factor,
//            nprint,  &info, &nfev, fjac,ldfjac, ipvt, qtf,
//            wa1, wa2, wa3, wa4);
// 
//     //the below is output result. compare it with the values used to generate 
//     //observed variables
//     double a=x[0];
//     double b=x[1];
//     double c=x[2];
// 
	return 0;

}