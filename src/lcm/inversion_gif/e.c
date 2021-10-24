/*


https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Parameter_plane#Plane_types




family of one parameter functions ( complex quadratic polynomial)


parameter plane for each type ( parameter)

parameter ( plane) types
* plain plane = c plane ( basic , reference plane)
* inverted c plane 
* c_Myrberg_type : "The point 1.40115 is called the "Myreberg point" of the Mandelbrot set. The sequence of circles attached to the right of the main cardioid get smaller and smaller and approach this point. That point is not the end of the Mandelbrot set since there's a path leading off to the right. Inverting on that point makes all these circles larger and larger instead of smaller and smaller. Exploring this inverted plane can be quite interesting. The original cardioid is turned around and distorted a bit. It appears near the center of this image. The big circle to its left is the inversion of the small circle to the right of the original cardioid. The little bit of a line moving off to the right of the image is the end of the path in the µ-plane ending at z = 2. "




Alternate parameter planes : 
The collection of quadratic polynomials can be parameterized in different ways which lead to different shapes for the Mandelbrot sets.






  Adam Majewski
  adammaj1 aaattt o2 dot pl  // o like oxygen not 0 like zero 
  
  
  console program in c programing language 
===============================================================





  
  ==============================================
  
  
  Structure of a program or how to analyze the program 
  
  
  
  Creating graphic:
  * memory array
  * save it to the disk as a pgm file
  * convert pgm file to png usnigng Image Magic convert
  
  
  creating image
  * rectangle from complex plane: p= plane
  * map it to the c plane: for each pixel of plane compute c or lambda using  map_parameter
  
  
  
  

   
  ==========================================

  
  ---------------------------------
  indent d.c 
  default is gnu style 
  -------------------



  c console progam 
  
	export  OMP_DISPLAY_ENV="TRUE"	
  	gcc d.c -lm -Wall -march=native -fopenmp
  	time ./a.out > b.txt


  gcc e.c -lm -Wall -march=native -fopenmp


  time ./a.out

  time ./a.out >a.txt
  
  ./g.sh
  
  ============================
  
  gcc e.c -lm -Wall -march=native -fopenmp -pg
  gprof ./a.out > p.txt
  
   
   

  ----------------------
  
 real	0m19,809s
user	2m26,763s
sys	0m0,161s


  

*/

#include <stdio.h>
#include <stdlib.h>		// malloc
#include <string.h>		// strcat
#include <math.h>		// M_PI; needs -lm also
#include <complex.h> 		// complex numbers : https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c
#include <omp.h>		// OpenMP




// https://sourceforge.net/p/predef/wiki/Standards/

#if defined(__STDC__)
#define PREDEF_STANDARD_C_1989
#if defined(__STDC_VERSION__)
#if (__STDC_VERSION__ >= 199409L)
#define PREDEF_STANDARD_C_1994
#endif
#if (__STDC_VERSION__ >= 199901L)
#define PREDEF_STANDARD_C_1999
#endif
#endif
#endif







/* --------------------------------- global variables and consts ------------------------------------------------------------ */


// each typedef should have different range !!!



/* Representation FunctionType 
	https://mrob.com/pub/muency/representationfunction.html
	function defining relation between data and the image
*/
typedef enum  {
			LSM =100, 
			DEM = 101, 
			Unknown = 102, 
			BD = 103, 
			MBD = 104, 
			SAC, 
			DLD, 
			ND, 
			NP, 
			POT, 
			Blend
		
		} RepresentationFunctionTypeT; 


#define FMAX 2 // number of Family Types; 
typedef enum  { 
			c_type = 0, 
			lambda_type = 1
		
		} FamilyTypeT; 
// 


#define PMAX 4 // number of Transformation types; see ProjectionTypeT !!!!!!!!!!!!!!!!!!!!!!!!!!
/*

*/
typedef enum  { //
			identity = 10, // https://en.wikipedia.org/wiki/Identity_function
			inversion = 11, 
			exponentiation = 12, 
			moebius = 13
		} ProjectionTypeT; 



// virtual 2D array and integer ( screen) coordinate
// Indexes of array starts from 0 not 1 
//unsigned int ix, iy; // var
static unsigned int ixMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int ixMax;	//
static unsigned int iWidth;	// horizontal dimension of array

static unsigned int iyMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int iyMax;	//

static unsigned int iHeight = 1000;	//  
// The size of array has to be a positive constant integer 
static unsigned int iSize;	// = iWidth*iHeight; 


// ----------memmory 1D arrays ==================
// unsigned char = for 1 byte ( 8 bit) colors 
unsigned char *data;
unsigned char *edge;

 

// unsigned int i; // var = index of 1D array
//static unsigned int iMin = 0; // Indexes of array starts from 0 not 1
static unsigned int iMax;	// = i2Dsize-1  = 
// The size of array has to be a positive constant integer 
// unsigned int i1Dsize ; // = i2Dsize  = (iMax -iMin + 1) =  ;  1D array with the same size as 2D array

/*


Image Width = 3.0000000000000000 in world coordinate
PixelWidth = 0.0015007503751876 
plane radius = 1.5000000000000000 
plane center = -0.7500000000000000 +0.0000000000000000 
xMin  = -2.2500000000000000 	 xMax = 0.7500000000000000 
yMin  = -1.5000000000000000 	 yMax = 1.5000000000000000 
File LSM_c_2000_-0.750000_1.500000.pgm saved . Comment = one parameter family of complex quadratic polynomial, parameter plane ;  LSM_c 


Image Width = 5.4000000000000004 in world coordinate
PixelWidth = 0.0027013506753377 
plane radius = 2.7000000000000002 
plane center = 1.3300000000000001 +0.0000000000000000 
xMin  = -1.3700000000000001 	 xMax = 4.0300000000000002 
yMin  = -2.7000000000000002 	 yMax = 2.7000000000000002 
File LSM_c_inverted_2000_1.330000_2.700000.pgm saved . Comment = one parameter family of complex quadratic polynomial, parameter plane ;  LSM_c_inverted 



Image Width = 10.0000000000000000 in world coordinate
PixelWidth = 0.0050025012506253 
plane radius = 5.0000000000000000 
plane center = 4.0000000000000000 +0.0000000000000000 
xMin  = -1.0000000000000000 	 xMax = 9.0000000000000000 
yMin  = -5.0000000000000000 	 yMax = 5.0000000000000000 
File LSM_c_parabola_2000_4.000000_5.000000.pgm saved . Comment = one parameter family of complex quadratic polynomial, parameter plane ;  LSM_c_parabola 


Image Width = 801.3999999999999773 in world coordinate
PixelWidth = 0.4009004502251126 
plane radius = 400.6999999999999886 
plane center = 1.3300000000000001 +0.0000000000000000 
xMin  = -399.3700000000000045 	 xMax = 402.0299999999999727 
yMin  = -400.6999999999999886 	 yMax = 400.6999999999999886 
ratio of image  = 1.000000 ; it should be 1.000 ...
Maximal number of iterations = iterMax_LSM = 2000 




Mandelbrot Set (in the 1/(mu+.25) plane) : x in [-9.08763557,3.41706117]; y in [-5.62095252,5.74695361].
Mandelbrot Set (in the 1/(mu-1.40115) plane) : x in [-6.89980824,6.49615956]; y in [-6.71624278,6.67972502].
Mandelbrot Set (in the 1/(mu-2) plane) : x in [-2.32859532,-0.17140468]; y in [-0.93790897,0.93790897].
Mandelbrot Set (in the 1/lambda plane) : x in [-1.15856298,1.05217571]; y in [-1.10369313,1.10704555].
Mandelbrot Set (in the 1/(lambda-1) plane) : x in [-0.64246706,0.66013131]; y in [-1.08440356,1.08659373].


const double CxMin= -1.8;
const double CxMax=  3.8;
const double CyMin= -1.55;
const double CyMax=  1.55;

*/
// see set_plane				c 	lambda


const double plane_radii[FMAX] = 		{1.5, 	1.5};

const complex double plane_centers[PMAX] = 	{-0.5, 0.0};


						//c	lambda
const complex double critical_points[FMAX] = 	{0.0, 	0.5};

const double  DisplayAspectRatio  = 1.0; // https://en.wikipedia.org/wiki/Aspect_ratio_(image)




const complex double cf = - 1.401155; //the Feigenbaum point -1.401155



// parameter plane 
double xMin ;	//-0.05;
double xMax ;	//0.75;
double yMin ;	//-0.1;
double yMax ;	//0.7;


double PixelWidth;	// =(CxMax-CxMin)/ixMax;
double PixelHeight;	// =(CyMax-CyMin)/iyMax;
double ratio;






const int iterMax_LSM = 1000;
const int iterMax_DEM = 2500;

// EscapeRadius for bailout test 
double ER = 2000.0;	






// dem
double MinBoundaryWidth = 0.03; // fixed value. To do computing it for every pixel ??       
double BoundaryWidth = 2.0; // % of image width  



/* colors = shades of gray from 0 to 255 */
unsigned char iColorOfExterior = 250;
unsigned char iColorOfInterior = 200;
unsigned char iColorOfBoundary = 0;
unsigned char iColorOfUnknown = 30;





/* ------------------------------------------ functions -------------------------------------------------------------*/



//------------------complex numbers -----------------------------------------------------




// from screen to world coordinate ; linear mapping
// uses global cons
static inline double Give_x (const int ix)
{
  return (xMin + ix * PixelWidth);
}

// uses global cons
static inline double Give_y (const int iy) {
  return (yMax - iy * PixelHeight);
}				// reverse y axis


static inline complex double Give_p (const int ix, const int iy)
{
  double x = Give_x (ix);
  double y = Give_y (iy);

  return x + y * I;




}


complex double fm( const double complex z , const complex double m ){

	return m*z*(1.0-z);
}


complex double fc( const double complex z , const complex double c ){

	return z*z +c;
}



/* complex function. 


*/
complex double f(const FamilyTypeT FamilyType, const double complex z0 , const complex double p ) {

  	complex double z = z0;
  
  	switch(FamilyType){
	
		case c_type :		{z = z*z + p;  break;}	 // complex quadratic polynomial, p is changed in give_parameter function
		
		case lambda_type: 	{z = p*z*(1.0-z);  break;} // p is changed in give_parameter function
	
	
		default: {z = z*z + p; }
	
	
	}
  
  return  z;
}
	
	

// projection from p to c or lambda 
complex double map_parameter(const ProjectionTypeT ProjectionType, const complex double parameter, const complex double translation){

	
	complex double p; 
	// plane transformation 
	switch(ProjectionType){
	
		case identity :{p = translation + parameter;  break;} // first translation and then identity 
		
		case inversion :{p = translation + 1.0/parameter; break;} // first translation then inverion, 2 transformations
		
		case exponentiation :{p = translation + cexp(parameter) ; break;} // here one can change cf to get different image 
		
		
		default: {p = parameter;}
	
	
	}
	



  
	
  return p;


}




complex double give_parameter(const ProjectionTypeT ProjectionType,  const complex double translation, const int ix, const int iy){

	// initial value of parameter
	complex double parameter= Give_p(ix,iy);
  	

	parameter = map_parameter(ProjectionType, parameter, translation);

  
	
  	return parameter;
	
	


}





// uses global var
int set_plane(const FamilyTypeT FamilyType, const ProjectionTypeT ProjectionType){

	complex double center = plane_centers[FamilyType];
	double radius = plane_radii[FamilyType];
	
	if (ProjectionType != exponentiation)
		{	

			

  			xMin = creal(center) - radius*DisplayAspectRatio;	
  			xMax = creal(center) + radius*DisplayAspectRatio;	//0.75;
  			yMin = cimag(center) - radius;	// inv
  			yMax = cimag(center) + radius;	//0.7;
  		}
  		
  		else {
  		
  				
  			xMax = 0.7; //  gives  0.5089024742041425 after transformation
  			xMin = xMax - 2.0*radius*DisplayAspectRatio; // 
  			yMin = cimag(center) - radius;	// inv
  			yMax = cimag(center) + radius;	//0.7;
  		}
  	
  	return 0;

}



void print_local_info(const RepresentationFunctionTypeT RepresentationFunctionType, const FamilyTypeT FamilyType, const ProjectionTypeT ProjectionType, const double translation){

	// view rectangle 
	printf ("Image Width = %.16f in world coordinate\n", xMax - xMin);
  	printf ("PixelWidth = %.16f \n", PixelWidth);
  	
  	printf ("plane radius = %.16f \n", plane_radii[FamilyType]);
  	complex double c = plane_centers[FamilyType];
  	printf ("plane center = %.16f %+.16f \n", creal(c),cimag(c) );
  	printf ("\tplane before transformation = p-plane\n" );
  	printf ("xMin  = %.16f \t xMax = %.16f \n", xMin, xMax );
  	printf ("yMin  = %.16f \t yMax = %.16f \n", yMin, yMax );
  	printf ("\tplane after transformation ( projection = modified c-plane  \n" );
  	printf ("xMin  = %.16f \t xMax = %.16f \n", creal(map_parameter(ProjectionType,xMin, translation)) , creal(map_parameter(ProjectionType,xMax, translation)) );
  	printf ("yMin  = %.16f \t yMax = %.16f \n", cimag(map_parameter(ProjectionType,yMin*I, translation)), cimag(map_parameter(ProjectionType,yMax*I, translation)) );
  	
  	
  	
    	printf ("ratio of image  = %f ; it should be 1.000 ...\n", ratio);
  
  	
  	 // map_parameter(const ProjectionTypeT ProjectionType, const complex double parameter)
  	// image corners in world coordinate
  	// center and radius
  	// center and zoom
  	// GradientRepetition
  	
  	
  	printf ("Maximal number of iterations = iterMax_LSM = %d \n", iterMax_LSM);
  	printf ("Maximal number of iterations = iterMax_DEM = %d \n", iterMax_DEM);
  	printf (" BoundaryWidth*iWidth/2000.0 = %f \n", BoundaryWidth*iWidth/2000.0);
 	printf ("MinimalBoundaryWidth = %.16f = %f pixels = %f * image width\n", MinBoundaryWidth, MinBoundaryWidth/PixelWidth, MinBoundaryWidth/(xMax - xMin) );
  	printf("\n\n");
  
  	
  //

}




int local_setup(const RepresentationFunctionTypeT RepresentationFunctionType, const FamilyTypeT FamilyType, const ProjectionTypeT ProjectionType, const double translation){

	set_plane(FamilyType, ProjectionType);	
	
	
  	/* Pixel sizes of the initial plane, before transformation !!!! */
  	PixelWidth = (xMax - xMin) / ixMax;	//  ixMax = (iWidth-1)  step between pixels in world coordinate 
	PixelHeight = (yMax - yMin) / iyMax;
  	ratio = ((xMax - xMin) / (yMax - yMin)) / ((double) iWidth / (double) iHeight);	// it should be 1.000 ...
  	
  	MinBoundaryWidth = 0.0300000000000000; // PixelWidth*BoundaryWidth*iWidth/2000.0; //0.01*cabs(parameter);
  	  	
	print_local_info( RepresentationFunctionType, FamilyType, ProjectionType, translation);
	return 0;
};






/* -----------  array functions = drawing -------------- */

/* gives position of 2D point (ix,iy) in 1D array  ; uses also global variable iWidth */
static inline unsigned int Give_i (const int ix, const int iy)
{
  return ix + iy * iWidth;
}


// ***********************************************************************************************
// ********************** edge detection usung Sobel filter ***************************************
// ***************************************************************************************************

// from Source to Destination
int ComputeBoundaries(const unsigned char S[], unsigned char D[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
  /* sobel filter */
  unsigned char G, Gh, Gv; 
  // boundaries are in D  array ( global var )
 
  // clear D array
  memset(D, iColorOfExterior, iSize*sizeof(*D)); // 
 
  // printf(" find boundaries in S array using  Sobel filter\n");   
#pragma omp parallel for schedule(dynamic) private(i,iY,iX,Gv,Gh,G) shared(iyMax,ixMax)
  for(iY=1;iY<iyMax-1;++iY){ 
    for(iX=1;iX<ixMax-1;++iX){ 
      Gv= S[Give_i(iX-1,iY+1)] + 2*S[Give_i(iX,iY+1)] + S[Give_i(iX-1,iY+1)] - S[Give_i(iX-1,iY-1)] - 2*S[Give_i(iX-1,iY)] - S[Give_i(iX+1,iY-1)];
      Gh= S[Give_i(iX+1,iY+1)] + 2*S[Give_i(iX+1,iY)] + S[Give_i(iX-1,iY-1)] - S[Give_i(iX+1,iY-1)] - 2*S[Give_i(iX-1,iY)] - S[Give_i(iX-1,iY-1)];
      G = sqrt(Gh*Gh + Gv*Gv);
      i= Give_i(iX,iY); /* compute index of 1D array from indices of 2D array */
      if (G==0) {D[i]=255;} /* background */
      else {D[i]=0;}  /* boundary */
    }
  }
 
   
 
  return 0;
}



// copy from Source to Destination
int CopyBoundaries(const unsigned char S[],  unsigned char D[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
 
 
  fprintf(stderr, "copy boundaries from S array to D array \n");
  for(iY=1;iY<iyMax-1;++iY)
    for(iX=1;iX<ixMax-1;++iX)
      {i= Give_i(iX,iY); if (S[i]==0) D[i]=0;}
 
 
 
  return 0;
}




// ***************************************************************************************************************************
// ************************** LSM*****************************************
// ****************************************************************************************************************************

unsigned char ComputeColorOfLSM(const FamilyTypeT FamilyType, complex double p){

	int nMax = iterMax_LSM;
  	
  	unsigned char iColor;
	
  	int n;
  	
  	complex double z = critical_points[FamilyType];

  	for (n=0; n < nMax; n++){ //forward iteration
	
    		if (cabs(z) > ER) break; // esacping
    	
  		
  		z = fm(z,p); //  for speed only one family here without switch 	
   		//z = f(FamilyType, z,p); /* forward iteration : complex quadratic polynomial */ 
  	}
  
  	if (n ==nMax)
  		{iColor = 0;} // interior = non escaping set
  		else iColor = 255 - 255.0 * ((double) n)/60; // nMax or lower walues in denominator ; exterior = escaping set
  
  
  	return iColor;


}


// ***************************************************************************************************************************
// ************************** DEM = exterior DE Method where DE = Distance Estimation  only for z^+c family !!!! ************
// ****************************************************************************************************************************


double Give_DE_c(double complex C )
{
  int i=0; // iteration 
   
   
  double complex Z= 0.0; // initial value for iteration Z0
  double R; // =radius = cabs(Z)
  double D; 
  double complex dC = 0; // derivative
  double de; // = 2 * z * log(cabs(z)) / dc;
  int iMax = iterMax_DEM;
  
   
    
  // iteration = computing the orbit
  for(i=0;i<iMax;i++)
    { 
    	// only for c family 
      dC = 2 * Z * dC + 1.0; 
      Z= Z*Z+C; // https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/qpolynomials
      
            
      R = cabs(Z);
      if(R > ER) break; // exterior of M set
   
      
    } // for(i=0
   
   
  if (i == iMax) D = -1.0; // interior 
    else { // exterior
      de = 2.0 * R * log(R) / cabs(dC) ; // 
    
      if (de < MinBoundaryWidth) D = de; //FP_ZERO; //  boundary
         else  D = 1.0; // exterior
    }
    
  return D; 
}
 
 
 


double Give_DE_m(double complex M )
{
  int i=0; // iteration 
   
   
  double complex Z = 0.5; // initial value for iteration Z0
  double R; // =radius = cabs(Z)
  double D; 
  double complex dM = 1.0;; // derivative
  double de; // = 2 * z * log(cabs(z)) / dc;
  int iMax = iterMax_DEM;
  
	
	
    
  // iteration = computing the orbit
  for(i=0;i<iMax;i++)
    { 
    	// only for m family 
      dM = (1 - 2*Z)*dM*M + Z -Z*Z; // (1-2*An)*Dn*m-An^2+An
      Z  = M*Z*(1-Z); // https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/qpolynomials
      
            
      R = cabs(Z);
      if ( R > ER) break; // exterior of M set
   
      
    } // for(i=0
   
   
  if (i == iMax) 
  	{D = -1.0;} // interior 
  	else { // exterior
  		de = 2.0 * R * log(R) / cabs(dM) ; // 
    	
  		if (de < MinBoundaryWidth) 
  			{D = FP_ZERO;} //  boundary
         		else  D = 1.0; // exterior
    }
    
  return D; 
}
  
 
 
 

unsigned char ComputeColorOfDE( complex double p){

	//int nMax = iterMax_DEM;
  	
  	unsigned char iColor;
	
  	//int n;
  	
  	
  	
  	double de = Give_DE_m( p );
  	// 
  	if (de < 0.0)
     		{iColor = iColorOfInterior;}/*  interior of Mandelbrot set = inside_color =  */
                               
    		else // exterior and boundary
     		{      
       			if (de == FP_ZERO) 
       				{iColor = iColorOfBoundary; }// boundary     
          			else iColor = iColorOfExterior;// exterior
     		};

  
  	return iColor;


}


 
 
/* ==================================================================================================
 ============================= Draw functions ===============================================================
 =====================================================================================================
*/ 
unsigned char ComputeColor(const RepresentationFunctionTypeT RepresentationFunctionType, const FamilyTypeT FamilyType, const complex double parameter){

	unsigned char iColor= 0;
	
	
	
	switch(RepresentationFunctionType){
	
		case LSM :{iColor = ComputeColorOfLSM(FamilyType, parameter); break;}
		
		case DEM : {iColor = ComputeColorOfDE( parameter); break; } // only lambda family 
		
		
	
		default: {}
	
	
	}
	
	return iColor;



}




complex double GiveParameterAndComputePixelWidth(const FamilyTypeT FamilyType, const ProjectionTypeT ProjectionType, const double translation, const int ix, const int iy){


	complex double parameter = give_parameter(ProjectionType, translation, ix, iy);
	
	
	
	
	return parameter; 

}


// plots  raster point (ix,iy) = computes it's color and save it to the array A
int DrawPoint (const RepresentationFunctionTypeT RepresentationFunctionType, const FamilyTypeT FamilyType, const ProjectionTypeT ProjectionType, const double translation, const int ix, const int iy, unsigned char A[])
{
	complex double parameter = GiveParameterAndComputePixelWidth( FamilyType,  ProjectionType, translation,  ix, iy);
	
	unsigned char iColor = ComputeColorOfLSM(FamilyType, parameter);
	//ComputeColorOfDE( parameter); // for speed only one family here without switch 
	// ComputeColor(RepresentationFunctionType, FamilyType, parameter);
  
	unsigned int i = Give_i (ix, iy);	/* compute index of 1D array from indices of 2D array */
  
  	A[i] = iColor ; // 
  
  	return 0;
}




// fill array 
// uses global var :  ...
// scanning complex plane 
int DrawImage (const RepresentationFunctionTypeT RepresentationFunctionType, const FamilyTypeT FamilyType, const ProjectionTypeT ProjectionType, const double translation, unsigned char A[])
{
  	unsigned int ix, iy;		// pixel coordinate 
  	
  	local_setup(RepresentationFunctionType, FamilyType, ProjectionType, translation);

  	fprintf(stderr, "compute image RepresentationFunctionType = %d ProjectionType = %d translation = %f\n", RepresentationFunctionType, ProjectionType, translation);
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(ix,iy) shared(A, ixMax , iyMax)
  	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix)
      			{DrawPoint(RepresentationFunctionType, FamilyType, ProjectionType, translation, ix, iy, A);}	//  
  		}

  return 0;
}













 
// *******************************************************************************************
// ********************************** save A array to pgm file ****************************
// *********************************************************************************************

int SaveImage(const unsigned char A[], const char *shortName )
{

  FILE *fp;
  const unsigned int MaxColorComponentValue = 255;	/* color component is coded from 0 to 255 ;  it is 8 bit color file */
  
  
  
  // https://programmerfish.com/create-output-file-names-using-a-variable-in-c-c/
  char fileName[512];
  const char* fileType = ".pgm";
  sprintf(fileName,"%s%s", shortName, fileType); // 
  
  
  
  char long_comment[200];
  sprintf (long_comment, "one parameter family of complex quadratic polynomial, parameter plane ");





  // save image array to the pgm file 
  fp = fopen (fileName, "wb");	// create new file,give it a name and open it in binary mode 
  fprintf (fp, "P5\n # %s\n %u %u\n %u\n", long_comment, iWidth, iHeight, MaxColorComponentValue);	// write header to the file
  size_t rSize = fwrite (A, sizeof(A[0]), iSize, fp);	// write whole array with image data bytes to the file in one step 
  fclose (fp);

  // info 
  if ( rSize == iSize) 
  	{
  		printf ("File %s saved ", fileName);
  		if (long_comment == NULL || strlen (long_comment) == 0)
    		printf ("\n");
  			else { printf (". Comment = %s \n", long_comment); }
  	}
  	else {printf("wrote %zu elements out of %u requested\n", rSize,  iSize);}
  	
  	
   
  
  return 0;
}




const char* GiveName(const double translation)
{

	static char Name[512];
	  	
    	sprintf(Name,"%f",  translation);
    	
    	
    
    	return Name;
}




int MakeImages( const FamilyTypeT FamilyType, const ProjectionTypeT ProjectionType, const double translation){

	const char *Name;

	DrawImage(LSM, FamilyType, ProjectionType, translation, data);
	Name = GiveName(translation);
	//SaveImage(data, Name); 
	
	ComputeBoundaries(data,edge);
	//Name = GiveName(translation);
	SaveImage(edge, Name); 
	
	//CopyBoundaries(edge, data);
	//shortName = GiveName("LSCM",  ProjectionType);
	//SaveImage(data, shortName); 
	
	printf("==========================================================================================================================\n\n\n\n");

	return 0;

}




/*

********************************************* info 

*/




int PrintInfoAboutProgam()
{
	
  
  // 
  printf (" \n");
  

  printf("gcc version: %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__); // https://stackoverflow.com/questions/20389193/how-do-i-check-my-gcc-c-compiler-version-for-my-eclipse
  // OpenMP version is displayed in the console 
  return 0;
}






// *****************************************************************************
//;;;;;;;;;;;;;;;;;;;;;;  program setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// **************************************************************************************

int setup ()
{

  fprintf (stderr, "setup start\n");
  
  
  
  
  
	
  /* 2D array ranges */
  
  iWidth = iHeight* DisplayAspectRatio;
  iSize = iWidth * iHeight;	// size = number of points in array 
  // iy
  iyMax = iHeight - 1;		// Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].
  //ix

  ixMax = iWidth - 1;

  /* 1D array ranges */
  // i1Dsize = i2Dsize; // 1D array with the same size as 2D array
  iMax = iSize - 1;		// Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].
  
  
  
	
  
  
   	
  /* create dynamic 1D arrays for colors ( shades of gray ) */
  data = malloc (iSize * sizeof (unsigned char));
  edge = malloc (iSize * sizeof (unsigned char));
   
  //
 
  	
  if (data == NULL || edge == NULL ){
    fprintf (stderr, " Setup error : Could not allocate memory");
    return 1;
  }

  
  
  
  fprintf (stderr," end of setup \n");
	
  return 0;

} // ;;;;;;;;;;;;;;;;;;;;;;;;; end of the setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




int end(){


  fprintf (stderr," allways free memory (deallocate )  to avoid memory leaks \n"); // https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
  free (data);
  free(edge);
  
 
  PrintInfoAboutProgam();
  return 0;

}









// ********************************************************************************************************************
/* -----------------------------------------  main   -------------------------------------------------------------*/
// ********************************************************************************************************************

int main () {
  
  
  
	setup ();
	
	double t0 = -4.0 ; // translation
	int nMax = 100; // nimber of images 
	double dt = 8.0/ nMax;
	// double m; 
	double t;
	
	
	FamilyTypeT family = 1;
	ProjectionTypeT projection = 11;
	
	
	// translation 	for lamba from 0 to 2
	// tramslation for c from ? to ? 
	for (int n = 0; n <= nMax; ++n){
		
		t =  t0 + n*dt;
		MakeImages(family, projection,t);
		}
		//MakeImages(family, projection, 0.0);
	
	
	
	
	
	
	
  	end();
  	
  	//printf(" dt = %f = %f * pixelWidth = %f * MinBoundaryWidth\n", dt, dt/PixelWidth, dt / MinBoundaryWidth);

  	return 0;
}
