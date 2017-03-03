//file class_integral_Gauss_3d_complex.h
#include <iostream>
#include <math.h>
#include <cmath>
#include <complex>
#include <limits>

using namespace std;
typedef complex< double > dcomplex;
const dcomplex II = dcomplex( 0.0, 1.0 );

#define PI M_PI

class integral_Gauss_3d_complex//members by default are private
{
      typedef dcomplex (*pfn)(double,double,double,double *); 	//define a function pointer that is function of three ind. variables 
								//plus an additional array containing every other parameters of the function
      double h;//integral bound, the integration is over a cube i.e -h < x < h , -h < y < h , -h < z < h
      double tol;                           //tolerance for the adaptive method
      double *varargin;                     //Contains all others parameters to be passed to the integrand
      pfn f;                                //integrand function
      double x;
      double y;
      double z;
      int min_recursion;
      int recursion;
      dcomplex init_value;

public:                                     //public members
       integral_Gauss_3d_complex(double bound,double tolerance,double *param, int size ,pfn integrand,double xc,double yc,double zc)//constructor. size is the size of varargin
       {
            h=bound; tol = numeric_limits<double>::max(); f=integrand; x = xc; y = yc; z = zc;
            varargin = new double[size];
            min_recursion = 3;
            recursion = 1;
            for (int r = 0; r < size; r++ )
            {
                varargin[r]=param[r];
            }
            init_value = adaptative_integration_Gauss_3d();
            tol = tolerance*abs(init_value);
       }
       ~integral_Gauss_3d_complex() {delete[] varargin;} //destructor

       dcomplex adaptative_integration_Gauss_3d();
private:
       dcomplex Gauss_3d_non_adapt(pfn ,double,double ,double ,double ,double *);
       dcomplex adapt_int_Gauss(pfn,double ,double,double,double,double,double *,int,int);
};

dcomplex integral_Gauss_3d_complex::Gauss_3d_non_adapt(pfn f,double h,double x0,double y0,double z0,double *varargin)
{
 //the integration is over a cube i.e -h < x < h , -h < y < h , -h < z < h
 //and is normalized by h^3, the point (x0,y0,z0) is the center of the cube


       double w1 = 40.0/361.0;
       double w2 = 121.0/2888.0;

       double a = sqrt(19.0/30.0)*h;
       double b = sqrt(19.0/33.0)*h;

       dcomplex temp_w1 = f(x0+a,y0,z0,varargin) + f(x0-a,y0,z0,varargin) + f(x0,y0+a,z0,varargin) + f(x0,y0-a,z0,varargin) + f(x0,y0,z0+a,varargin) + f(x0,y0,z0-a,varargin);

       dcomplex temp_w2 = f(x0+b,y0+b,z0+b,varargin) + f(x0+b,y0+b,z0-b,varargin) + f(x0+b,y0-b,z0+b,varargin) + f(x0+b,y0-b,z0-b,varargin) + f(x0-b,y0+b,z0+b,varargin) + f(x0-b,y0+b,z0-b,varargin) + f(x0-b,y0-b,z0+b,varargin) + f(x0-b,y0-b,z0-b,varargin);

       return (w1*temp_w1 + w2*temp_w2);

}


dcomplex integral_Gauss_3d_complex::adaptative_integration_Gauss_3d()
{
      return adapt_int_Gauss(f,h,x,y,z,tol,varargin,recursion,min_recursion);
}

dcomplex integral_Gauss_3d_complex::adapt_int_Gauss(pfn f,double h1,double x0,double y0,double z0,double tol1,double *varargin,int recursion,int min_recursion)
{
       dcomplex one_trapezoid_area = Gauss_3d_non_adapt(f,h1,x0,y0,z0,varargin);

       dcomplex two_trapezoid_area = (Gauss_3d_non_adapt(f,h1/2,x0+h1/2,y0+h1/2,z0+h1/2,varargin) + Gauss_3d_non_adapt(f,h1/2,x0+h1/2,y0+h1/2,z0-h1/2,varargin) +
Gauss_3d_non_adapt(f,h1/2,x0+h1/2,y0-h1/2,z0+h1/2,varargin)  + Gauss_3d_non_adapt(f,h1/2,x0+h1/2,y0-h1/2,z0-h1/2,varargin)+
Gauss_3d_non_adapt(f,h1/2,x0-h1/2,y0+h1/2,z0+h1/2,varargin) + Gauss_3d_non_adapt(f,h1/2,x0-h1/2,y0+h1/2,z0-h1/2,varargin) +
Gauss_3d_non_adapt(f,h1/2,x0-h1/2,y0-h1/2,z0+h1/2,varargin)  + Gauss_3d_non_adapt(f,h1/2,x0-h1/2,y0-h1/2,z0-h1/2,varargin))/8.0;

        if (abs(two_trapezoid_area - one_trapezoid_area) < tol1 && recursion >= min_recursion)//( fabs(real(two_trapezoid_area - one_trapezoid_area)) < tol1 && fabs(imag(two_trapezoid_area - one_trapezoid_area)) < tol1  ) 
        {

           return two_trapezoid_area;
        }
        else
        {
            tol1 = 8.0*tol1;

            return (adapt_int_Gauss(f,h1/2,x0-h1/2,y0+h1/2,z0+h1/2,tol1/sqrt(8),varargin,recursion+1,min_recursion) + adapt_int_Gauss(f,h1/2,x0-h1/2,y0+h1/2,z0-h1/2,tol1/sqrt(8),varargin,recursion+1,min_recursion) + 
            adapt_int_Gauss(f,h1/2,x0-h1/2,y0-h1/2,z0+h1/2,tol1/sqrt(8),varargin,recursion+1,min_recursion)  + adapt_int_Gauss(f,h1/2,x0-h1/2,y0-h1/2,z0-h1/2,tol1/sqrt(8),varargin,recursion+1,min_recursion)+ 
            adapt_int_Gauss(f,h1/2,x0+h1/2,y0+h1/2,z0+h1/2,tol1/sqrt(8),varargin,recursion+1,min_recursion) + adapt_int_Gauss(f,h1/2,x0+h1/2,y0+h1/2,z0-h1/2,tol1/sqrt(8),varargin,recursion+1,min_recursion) + 
            adapt_int_Gauss(f,h1/2,x0+h1/2,y0-h1/2,z0+h1/2,tol1/sqrt(8),varargin,recursion+1,min_recursion)  + adapt_int_Gauss(f,h1/2,x0+h1/2,y0-h1/2,z0-h1/2,tol1/sqrt(8),varargin,recursion+1,min_recursion))/8.0;
        }
}
