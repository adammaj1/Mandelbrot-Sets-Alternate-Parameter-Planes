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
  
  
  ============== Image X ========================
  
  DrawImageOfX -> DrawPointOfX -> ComputeColorOfX 
  
  first 2 functions are identical for every X
  check only last function =  ComputeColorOfX
  which computes color of one pixel !
  
  

   
  ==========================================

  
  ---------------------------------
  indent d.c 
  default is gnu style 
  -------------------



  c console progam 
  
	export  OMP_DISPLAY_ENV="TRUE"	
  	gcc d.c -lm -Wall -march=native -fopenmp
  	time ./a.out > b.txt


  gcc d.c -lm -Wall -march=native -fopenmp


  time ./a.out

  time ./a.out >a.txt
  
  
  ./g.sh

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

#define NoP 8 // number of parameter types; see ParameterTypeT



/* Representation FunctionType 
	https://mrob.com/pub/muency/representationfunction.html
	function defining relation between data and the image
*/
typedef enum  {LSM , DEM, Unknown, BD, MBD , SAC, DLD, ND, NP, POT, Blend
		
		} RepresentationFunctionTypeT; 
// FunctionTypeT FunctionType;


/*

*/
typedef enum  {c_type, c_inverted_type, c_parabola_type, c_Myrberg_type, c_inverted_2_type, lambda_type , lambda_inverted_type, lambda_inverted_1_type
		
		} ParameterTypeT; 



// virtual 2D array and integer ( screen) coordinate
// Indexes of array starts from 0 not 1 
//unsigned int ix, iy; // var
static unsigned int ixMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int ixMax;	//
static unsigned int iWidth;	// horizontal dimension of array

static unsigned int iyMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int iyMax;	//

static unsigned int iHeight = 5000;	//  
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




*/
// see set_plane				c 	c_inv 	c_par 	c_Myr 	c_2	lambda	lambd_inv


const double plane_radii[NoP] = 		{1.5, 	2.7, 	5.0, 	400.7, 	5.0, 	3.2,	1.12,	1.1};

const complex double plane_centers[NoP] = 	{-0.75, 1.33, 	4.0, 	1.33, 	2.0,	1.0,	0.0,	0.0};


const complex double critical_points[NoP] = 	{0.0, 	0.0, 	0.0, 	0.0, 	0.0, 	0.5, 	0.5,	0.5};

const double  DisplayAspectRatio  = 1.0; // https://en.wikipedia.org/wiki/Aspect_ratio_(image)



// parameter plane 
double xMin ;	//-0.05;
double xMax ;	//0.75;
double yMin ;	//-0.1;
double yMax ;	//0.7;


double PixelWidth;	// =(CxMax-CxMin)/ixMax;
double PixelHeight;	// =(CyMax-CyMin)/iyMax;
double ratio;













const int iterMax_LSM = 2000;


// EscapeRadius for bailout test 
double ER = 200.0;	


double BoundaryWidth = 3.0; // % of image width  
double distanceMax; //distanceMax = BoundaryWidth*PixelWidth;


/* colors = shades of gray from 0 to 255 */
unsigned char iColorOfExterior = 250;
unsigned char iColorOfInterior = 200;
unsigned char iColorOfBoundary = 0;
unsigned char iColorOfUnknown = 30;





/* ------------------------------------------ functions -------------------------------------------------------------*/

// complex function
complex double f(const ParameterTypeT ParameterType, const double complex zn , const complex double p ) {

  	complex double z = zn;
  
  	switch(ParameterType){
	
		case c_type :		 // complex quadratic polynomial
		
		case c_inverted_type :	 //  p is inverted in give_parameter function
		
		case c_parabola_type :	 // p is changed in give_parameter function
		
		case c_Myrberg_type :	 // p is changed in give_parameter function
			
		case c_inverted_2_type :	{z = z*z + p;  break;} // p is changed in give_parameter function
		
		case lambda_type :		 // p is changed in give_parameter function
		
		case lambda_inverted_type :		// p is changed in give_parameter function
		
		case lambda_inverted_1_type :		{z = p*z*(1.0-z);  break;} // p is changed in give_parameter function
	
	
		default: {}
	
	
	}
  
  return  z;
}
	

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


/*
https://mathcs.clarku.edu/~djoyce/julia/MandelbrotColorPlane.java

public static Complex convert (Complex z, int oldPlane, int newPlane) {
    if (oldPlane == newPlane)
       return new Complex(z);
    // first convert from the old plane to either the lambda or the mu plane
    if (oldPlane != MU && oldPlane != LAMBDA) {
      z = z.reciprocal();
      if (oldPlane == RECIPMUPLUSFOURTH)
        z = z.minus(0.25);
      else if (oldPlane == RECIPLAMBDAMINUSONE)
        z = z.plus(1.0);
      else if (oldPlane == RECIPMUMINUSMYER)
        z = z.plus(MYERBERG);
    } // if
    // next, convert to mu or lambda as necessary
    if (oldPlane==LAMBDA || oldPlane==RECIPLAMBDA || oldPlane==RECIPLAMBDAMINUSONE) {
      if (newPlane!=LAMBDA && newPlane!=RECIPLAMBDA && newPlane!=RECIPLAMBDAMINUSONE) {
        // convert lambda to mu.  mu = (lambda/2)^2 - (lambda/2)
        z = z.over(2.0);
        z = z.times(z).minus(z);
      } // if
    } else {
      if (newPlane==LAMBDA || newPlane==RECIPLAMBDA || newPlane==RECIPLAMBDAMINUSONE) {
        // convert mu to lambda.  lambda = 1 + sqrt(1+4mu)
        z = z.times(4.0).plus(1.0);
        z = z.sqrt().plus(1.0);
      } // if
    } // if/else
    // finally, convert to the new plane
    switch (newPlane) {
      case MU: return z;
      case LAMBDA: return z;
      case RECIPMU: return z.reciprocal();
      case RECIPMUPLUSFOURTH: return z.plus(0.25).reciprocal();
      case RECIPLAMBDA: return z.reciprocal();
      case RECIPLAMBDAMINUSONE: return z.minus(1.0).reciprocal();
      case RECIPMUMINUSMYER: return z.minus(MYERBERG).reciprocal();
    } // switch
    return z;  // never used, but makes compiler happy
  } // convert







c parameter is a reference plane 
so all planes are compared to c 
for simplicity

*/
complex double give_parameter(const ParameterTypeT ParameterType, const int ix, const int iy){


	complex double parameter= Give_p(ix,iy);
  	


	switch(ParameterType){
	
		case c_type :{  break;}
		
		case c_inverted_type :{parameter = 1.0/parameter; break;}
		
		case c_parabola_type :{parameter = 0.25+ 1.0/parameter; break;}
		
		case c_Myrberg_type :{parameter = -1.401155 - 1.0/parameter; break;}
		
		case c_inverted_2_type :{parameter = -2.0 + 1.0/parameter; break;}
		
		case lambda_type :{  break;}
		
		case lambda_inverted_type :{parameter = 1.0/parameter; break;}
		
		case lambda_inverted_1_type :{parameter =1.0+ 1.0/parameter; break;}
		
		default: {}
	
	
	}
	



  
	
  return parameter;
	
	


}





// uses global var
int set_plane(const ParameterTypeT ParameterType){

	complex double center = plane_centers[ParameterType];
	double radius = plane_radii[ParameterType];

  	xMin = creal(center) - radius*DisplayAspectRatio;	
  	xMax = creal(center) + radius*DisplayAspectRatio;	//0.75;
  	yMin = cimag(center) - radius;	// inv
  	yMax = cimag(center) + radius;	//0.7;
  	
  	return 0;

}



void print_local_info(const RepresentationFunctionTypeT RepresentationFunctionType, const ParameterTypeT ParameterType){

	// view rectangle 
	printf ("Image Width = %.16f in world coordinate\n", xMax - xMin);
  	printf ("PixelWidth = %.16f \n", PixelWidth);
  	printf ("plane radius = %.16f \n", plane_radii[ParameterType]);
  	complex double c = plane_centers[ParameterType];
  	printf ("plane center = %.16f %+.16f \n", creal(c),cimag(c) );
  	printf ("xMin  = %.16f \t xMax = %.16f \n", xMin, xMax );
  	printf ("yMin  = %.16f \t yMax = %.16f \n", yMin, yMax );
    	printf ("ratio of image  = %f ; it should be 1.000 ...\n", ratio);
  
  	
  	 
  	// image corners in world coordinate
  	// center and radius
  	// center and zoom
  	// GradientRepetition
  	
  	
  	printf ("Maximal number of iterations = iterMax_LSM = %d \n", iterMax_LSM);
 
  	printf("\n");
  
  	
  //

}




int local_setup(const RepresentationFunctionTypeT RepresentationFunctionType, const ParameterTypeT ParameterType){

	set_plane(ParameterType);	
	
	
  	/* Pixel sizes */
  	PixelWidth = (xMax - xMin) / ixMax;	//  ixMax = (iWidth-1)  step between pixels in world coordinate 
	PixelHeight = (yMax - yMin) / iyMax;
  	ratio = ((xMax - xMin) / (yMax - yMin)) / ((double) iWidth / (double) iHeight);	// it should be 1.000 ...
  	  	
	print_local_info( RepresentationFunctionType, ParameterType);
	return 0;
};


// ****************** DYNAMICS = trap tests ( target sets) ****************************







// bailout test
// z escapes when 
// abs(z)> ER or cabs2(z)> ER2 
// https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Julia_set#Boolean_Escape_time
// this function is not used !!!! dead code 

int Escapes(const complex double z){
 // here target set (trap) is the exterior  circle with radsius = ER ( EscapeRadius) 
  // with ceter = origin z= 0
  // on the Riemann sphere it is a circle with point at infinity as a center  
   
  if (cabs(z)>ER) return 1;
  return 0;
}








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

unsigned char ComputeColorOfLSM(const ParameterTypeT ParameterType, complex double p){

	int nMax = iterMax_LSM;
  	
  	unsigned char iColor;
	
  	int n;
  	
  	complex double z = critical_points[ParameterType];

  	for (n=0; n < nMax; n++){ //forward iteration
	
    		if (cabs(z) > ER) break; // esacping
    	
  			
   		z = f(ParameterType, z,p); /* forward iteration : complex quadratic polynomial */ 
  	}
  
  	if (n ==nMax)
  		{iColor = 0;} // interior = non escaping set
  		else iColor = 255 - 255.0 * ((double) n)/60; // nMax or lower walues in denominator ; exterior = escaping set
  
  
  	return iColor;


}

 
 
/* ==================================================================================================
 ============================= Draw functions ===============================================================
 =====================================================================================================
*/ 
unsigned char ComputeColor(const RepresentationFunctionTypeT RepresentationFunctionType, const ParameterTypeT ParameterType, const complex double parameter){

	unsigned char iColor= 0;
	
	
	
	switch(RepresentationFunctionType){
	
		case LSM :{iColor = ComputeColorOfLSM(ParameterType, parameter); break;}
		
		
		
		
	
		default: {}
	
	
	}
	
	return iColor;



}


// plots  raster point (ix,iy) = computes it's color and save it to the array A
int DrawPoint (const RepresentationFunctionTypeT RepresentationFunctionType, const ParameterTypeT ParameterType, const int ix, const int iy, unsigned char A[])
{
	complex double parameter = give_parameter(ParameterType, ix, iy);
  
	unsigned char iColor = ComputeColor(RepresentationFunctionType, ParameterType, parameter);
  
	unsigned int i = Give_i (ix, iy);	/* compute index of 1D array from indices of 2D array */
  
  	A[i] = iColor ; // 
  
  	return 0;
}




// fill array 
// uses global var :  ...
// scanning complex plane 
int DrawImage (const RepresentationFunctionTypeT RepresentationFunctionType, const ParameterTypeT ParameterType, unsigned char A[])
{
  	unsigned int ix, iy;		// pixel coordinate 
  	
  	local_setup(RepresentationFunctionType, ParameterType);

  	fprintf(stderr, "compute image RepresentationFunctionType = %d ParameterType = %d\n", RepresentationFunctionType, ParameterType);
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(ix,iy) shared(A, ixMax , iyMax)
  	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix)
      			DrawPoint(RepresentationFunctionType, ParameterType, ix, iy, A);	//  
  }

  return 0;
}
















 
// *******************************************************************************************
// ********************************** save A array to pgm file ****************************
// *********************************************************************************************

int SaveImage(const unsigned char A[], const ParameterTypeT ParameterType, const char *shortName , const char *comment)
{

  FILE *fp;
  const unsigned int MaxColorComponentValue = 255;	/* color component is coded from 0 to 255 ;  it is 8 bit color file */
  
  int n = iHeight;
  double center = plane_centers[ParameterType];
  double r = plane_radii[ParameterType];
  
  // https://programmerfish.com/create-output-file-names-using-a-variable-in-c-c/
  char fileName[512];
  const char* fileType = ".pgm";
  sprintf(fileName,"%s_%d_%f_%f%s", shortName, n, center,r, fileType); // 
  
  
  
  char long_comment[200];
  sprintf (long_comment, "one parameter family of complex quadratic polynomial, parameter plane ;  %s", comment);





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
	
	//
	DrawImage(LSM, c_type, data);
	SaveImage(data, c_type, "LSM_c","LSM_c"); 
	
	ComputeBoundaries(data,edge);
	SaveImage(edge, c_type, "LCM_c","LCM_c"); 
	
	CopyBoundaries(edge, data);
	SaveImage(data, c_type, "LSCM_c","LSCM_c"); 
  
 	// 
  
  	DrawImage(LSM, c_inverted_type, data);
	SaveImage(data, c_inverted_type, "LSM_c_inverted","LSM_c_inverted"); 
	
	
	ComputeBoundaries(data,edge);
	SaveImage(edge, c_inverted_type, "LCM_c_inverted","LCM_c_inverted"); 
	
	CopyBoundaries(edge, data);
	SaveImage(data, c_inverted_type, "LSCM_c_inverted","LSCM_c_inverted"); 
	
	//
	
	
	
	DrawImage(LSM, c_parabola_type, data);
	SaveImage(data, c_parabola_type, "LSM_c_parabola","LSM_c_parabola"); 
	
	ComputeBoundaries(data,edge);
	SaveImage(edge, c_parabola_type, "LCM_c_parabola","LCM_c_parabola"); 
	
	CopyBoundaries(edge, data);
	SaveImage(data, c_parabola_type, "LSCM_c_parabola","LSCM_c_parabola"); 
	
	//
	
	
	
	DrawImage(LSM, c_Myrberg_type, data);
	SaveImage(data, c_Myrberg_type, "LSM_c_Myrberg","LSM_c_Myrberg"); 
	
	ComputeBoundaries(data,edge);
	SaveImage(edge, c_Myrberg_type, "LCM_c_Myrberg","LCM_c_Myrberg"); 
	
	CopyBoundaries(edge, data);
	SaveImage(data, c_Myrberg_type, "LSCM_c_Myrberg","LSCM_c_Myrberg"); 
	
	//
	
  
  	DrawImage(LSM, c_inverted_2_type, data);
	SaveImage(data, c_inverted_2_type, "LSM_c_inverted_2","LSM_c_inverted_2"); 
	
	ComputeBoundaries(data,edge);
	SaveImage(edge, c_inverted_2_type, "LCM_c_inverted_2","LCM_c_inverted_2"); 
	
	CopyBoundaries(edge, data);
	SaveImage(data, c_inverted_2_type, "LSCM_c_inverted_2","LSCM_c_inverted_2"); 
	
	//
	
	
	
	DrawImage(LSM, lambda_type, data);
	SaveImage(data, lambda_type, "LSM_lambda","LSM_lambda"); 
	
	ComputeBoundaries(data,edge);
	SaveImage(edge, lambda_type, "LCM_lambda","LCM_lambda"); 
	
	CopyBoundaries(edge, data);
	SaveImage(data, lambda_type, "LSCM_lambda","LSCM_lambda"); 
  
 	// 
	
	
	DrawImage(LSM, lambda_inverted_type, data);
	SaveImage(data, lambda_inverted_type, "LSM_lambda_inverted","LSM_lambda_inverted"); 
	
	ComputeBoundaries(data,edge);
	SaveImage(edge, lambda_inverted_type, "LCM_lambda_inverted","LCM_lambda_inverted"); 
	
	CopyBoundaries(edge, data);
	SaveImage(data, lambda_inverted_type, "LSCM_lambda_inverted","LSCM_lambda_inverted"); 
  
 	// 
	
	
	DrawImage(LSM, lambda_inverted_1_type, data);
	SaveImage(data, lambda_inverted_1_type, "LSM_lambda_inverted_1","LSM_lambda_inverted_1"); 
	
	ComputeBoundaries(data,edge);
	SaveImage(edge, lambda_inverted_1_type, "LCM_lambda_inverted_1","LCM_lambda_inverted_1"); 
	
	CopyBoundaries(edge, data);
	SaveImage(data, lambda_inverted_1_type, "LSCM_lambda_inverted_1","LSCM_lambda_inverted_1"); 
 
  	end();

  	return 0;
}
