# Mandelbrot-Sets-Alternate-Parameter-Planes


Original page : [Julia and Mandelbrot Sets. Alternate Parameter Planes by David E. Joyce © 1994.](https://mathcs.clarku.edu/~djoyce/julia/altplane.html)


Images of [complex quadratic polynomials](https://en.wikipedia.org/wiki/Complex_quadratic_polynomial) under
* [Geometric transformation](https://en.wikipedia.org/wiki/Geometric_transformation)
* [complex mapping](https://www.falstad.com/complexviewer/): [A complex function f is the transformation function which maps the complex plane z to w](http://home.iitk.ac.in/~saiwal/engineering/complex-mappings/)








## z^2+c family

Here initial parameter p is transformed:   


    c = t(p)
    

    
  
    
    
So the whole plane (more precisely  rectangle viewport) is transformed. For example for c_exponnetial 

```
	plane before transformation 
xMin  = -8.9000000000000004 	 xMax = 0.7000000000000000 
yMin  = -2.3999999999999999 	 yMax = 2.3999999999999999 
	plane after transformation 
xMin  = -1.4010186110735179 	 xMax = 0.6125977074704767 
yMin  = -0.6754631805511510 	 yMax = 0.6754631805511510 
```
    
    
###  z^2 + c

First example is the standard plane, simple: c = p

  
![](./img/LCM_c_5000_-0.750000_1.500000_600.png)

Description by  [David E. Joyce](https://mathcs.clarku.edu/~djoyce/julia/altplane.html)  
>  the most notable feature of the set is a cardioid. 


### z^2 + (k +1/p)

Here 

c = k + 1/p

so this transformation consist of 2 transformations: 
* translation by k 
* the inversion with respect to the unit circle centered at zero


Translation is performed before inversion. It moves the set with respect to unit circle and thus changes the shape after the whole transformation


See [video by Max Million](https://www.youtube.com/watch?v=JiAIELoh4YU)

####  z^2 + 1/p

Second example is [simple inversion](https://mathworld.wolfram.com/Inversion.html) about a unit circle centered at origin  

Here c = 1/p ( and k = 0 so no translatio)  

  
![](./img/LCM_c_inverted_5000_1.330000_2.700000_600.png)



Description by  [Nikola Ubavić](https://ubavic.rs/work/julia_and_mandelbrot_set/?lang=sr)  
> composition of the inversion with respect to the unit circle centered at zero, and the conjugation (axial symmetry with respect to the real line). 
> The cardioid from the boundary of the Mandelbrot set in the "standard" parameterization corresponds to the tear-shaped curve in the 1/c parameterization.

Description by  [David E. Joyce](https://mathcs.clarku.edu/~djoyce/julia/altplane.html)  
> The inverse of the cardioid is the exterior of a teardrop shape. The circles on the outside of the cardioid are inverted to circles on the inside of the teardrop. The cusp of the cardioid becomes the cusp of the teardrop.


#### z^2 -2.0+1.0/p

Here : c =  -2.0+1.0/p ( and k = -2)

![](./img/LCM_c_inverted_2_5000_2.000000_5.000000_600.png)


####  z^2 + 0.25+1.0/p

Here: c =  0.25 + 1.0/p ( and k = 0.25)


  
  
![](./png/LCM_c_parabola_5000_4.000000_5.000000_600.png)

Description by  [Nikola Ubavić](https://ubavic.rs/work/julia_and_mandelbrot_set/?lang=sr)  
> If a 1/4 translation is performed before the inversion, then the cardioid is imaged in a parabola.


> The cusp of the cardioid is then moved to the origin, and that inverts to infinity. 
> the cardioid is sent to the outside of a parabola. 
> The circles on the outside of the cardioid are inverted to circles on the inside of the parabola.


Albert Chern: 
> Since a cardioid is a Möbius inversion of a parabola, you find that the Mandelbrot set is actually a fractal growing along a parabola. 


#### z^2 -1.401155 - 1.0/p

Here: c =  -1.401155 - 1.0/p

  

![](./img/LCM_c_Myrberg_5000_1.330000_400.700000_600.png)



### Exponential map

Here : c = cf + e^p

where 
* cf is the Feigenbaum point -1.401155 = Accumulation point of period-doubling





[Period doubling cascade](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/1over2_family): 
* one can see 7 hyperbolic components  for periods 2^0, 2^1, 2^2, 2^3, 2^4, 2^5, 2^6
* note that each componnet has tha same size 
   
![](./img/LCM_c_log_5000_-8.400000_2.400000_600.png)



## m*z*(z-1) = logistic family


  m*z*(1.0-z)
  
so here m = p  
  
  
![](./img/LCM_lambda_5000_1.000000_3.200000_600.png)


  z*(1.0-z)/p
  
here m = 1/p

![](./img/LCM_lambda_inverted_5000_0.000000_1.120000_600.png)


  (1 + 1/p)*z*(1.0-z)
  
here m = 1 + 1/p


![](./img/LCM_lambda_inverted_1_5000_0.000000_1.100000_600.png)  


Animated gif ( click the image to see the animation)  

![](./img/inversion.gif)  

The code for creating gif file is in [src/lcm/inversion_gif/ directory](src/lcm/inversion_gif/)





Compare
* [changing the inversion offset for the lambda plane mandelbrot by Max Million](https://www.youtube.com/watch?v=Xob9xZ2dBIs). Here m = k + 1/p where k is changing from 0 to 2 




# Code


## Families ( forms)  of the complex quadratic polynomial 
* [wikipedia](https://en.wikipedia.org/wiki/Complex_quadratic_polynomial#Forms)


```c
complex double f(const FamilyTypeT FamilyType, const double complex z0 , const complex double p ) {

  	complex double z = z0;
  
  	switch(FamilyType){
	
		case c_type :		{z = z*z + p;  break;}	 // complex quadratic polynomial, p is changed in give_parameter function
		
		case lambda_type: 	{z = p*z*(1.0-z);  break;} // p is changed in give_parameter function
	
	
		default: {z = z*z + p; }
	
	
	}
  
  return  z;
}
```

## Algorithms
* LCM = Level Curves Method = first compute level sets of escape time, then find boundaries of level sets  ( =  level curves ). Here Sobel filter is used  
* DEM = Distance Estimation Method 
  * [FF](https://fractalforums.org/programming/11/exterior-distance-estimation-for-logistic-map/4458)



## Plane transformations
    
```c
//    transformation = projection = map  from p    to c or lambda parameter
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
```
 


# To do
* [Mercator Mandelbrot Maps by: Anders Sandberg](https://www.flickr.com/photos/arenamontanus/albums/72157615740829949) : a logarithmic projection around a point c0: z-> (log(|z-c0|), arg(z-c0)). The idea for this kind of logarithmic map of the set is from [David Madore](http://www.madore.org/~david/math/mandelbrot.html). While it is a bit of a stretch to call it a Mercator projection, it sounds better than just a log scale map. If the Mandelbrot set is assumed to be about meter-sized these zooms reach the size of an atomic nucleus.  
* transverse Mercator
  * [The North Transept Transfigurational Dreamwarp by Robert Lee](https://www.flickr.com/photos/robertlee99/50289710431/) : 
* Riemann sphere
  * [mathematica.stackexchange question: how-to-draw-fractal-images-of-iteration-functions-on-the-riemann-sphere](https://mathematica.stackexchange.com/questions/15047/how-to-draw-fractal-images-of-iteration-functions-on-the-riemann-sphere)
* [Conformal_Warping_Around_The_Cardioid_In_The_Mandelbrot_Set](https://commons.wikimedia.org/wiki/File:Conformal_Warping_Around_The_Cardioid_In_The_Mandelbrot_Set.gif)
* [distorted Mandelbrot by Ben's Fractals ](https://www.youtube.com/watch?v=xp-teLc9bxk)
* Moebius by Claude
  * [interpolating_moebius_transformations](https://mathr.co.uk/blog/2015-02-06_interpolating_moebius_transformations.html)
  * [Moebius by  Claude Heiland-Allen](https://archive.org/details/MandelbrotMoebiusExperiments/mandelbrot_moebius_experiment_1.mkv)
  * [moebius_transformations_and_the_riemann_sphere](https://mathr.co.uk/blog/2017-12-01_moebius_transformations_and_the_riemann_sphere.html)
  * [Arbitrary Moebius Transformation of the Poincare Disk by Vladimir Bulatov](http://bulatov.org/math/1001/moebius-animation1.html)
* [Artistic Models of the Hyperbolic Geometry by Vladimir Bulatov](http://bulatov.org/math/1003/index.html)
* [Zⁿ+ A∙Sin(Cᵐ) by Jens-Peter Christensen](https://www.youtube.com/watch?v=po2t_yMWABs&list=RDGMEMJQXQAmqrnmK1SEjY_rKBGAVMpo2t_yMWABs&start_radio=1)

# Git

create a new repository on the command line

```git
echo "# Mandelbrot-Sets-Alternate-Parameter-Planes" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:adammaj1/Mandelbrot-Sets-Alternate-Parameter-Planes.git
git push -u origin main
```


## Local repo
```
~/Dokumenty/mandelbrot_planes$ 
```




## Subdirectory

```git
mkdir png
git add *.png
git mv  *.png ./png
git commit -m "move"
git push -u origin main
```
then link the images:

```txt
![](./png/n.png "description") 

```

```git
gitm mv -f 
```

[Remove a file/directory from a Git repository without deleting it from the local filesystem](https://stackoverflow.com/questions/1143796/remove-a-file-from-a-git-repository-without-deleting-it-from-the-local-filesyste)

```git
git rm --cached myfile.log
```

Single directory and all it's content

```git
git rm --cached -r ./png
git commit -m "png"
git push -u origin main

```


Rename directory

```git
git add ./lcm/
git mv ./lcm ./png 
git commit -m "png"
git push -u origin main

```

## Github
* [GitHub Flavored Markdown Spec](https://github.github.com/gfm/)
* [md cheat sheet](http://mdcheatsheet.com/)
* [CommonMark Spec](https://spec.commonmark.org)
* [Markdown parser ](https://markdown-it.github.io/)
