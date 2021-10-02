# Mandelbrot-Sets-Alternate-Parameter-Planes


[Julia and Mandelbrot Sets. Alternate Parameter Planes by David E. Joyce © 1994.](https://mathcs.clarku.edu/~djoyce/julia/altplane.html)


Images of [complex quadratic polynomials](https://en.wikipedia.org/wiki/Complex_quadratic_polynomial)


## z^2+p family

Here 


    p = t(c)
    
###  z^2 + c

First example is simple " p = c

   z^2 + c


![](./png/LCM_c_5000_-0.750000_1.500000_600.png)

Description by  [David E. Joyce](https://mathcs.clarku.edu/~djoyce/julia/altplane.html)  
>  the most notable feature of the set is a cardioid 



### z^2 + 1/c

Second example iss simple inversion. Here p = 1/c

  
![](./png/LCM_c_inverted_5000_1.330000_2.700000_600.png)



Description by  [Nikola Ubavić](https://ubavic.rs/work/julia_and_mandelbrot_set/?lang=sr)  
> composition of the inversion with respect to the unit circle centered at zero, and the conjugation (axial symmetry with respect to the real line). 
> The cardioid from the boundary of the Mandelbrot set in the "standard" parameterization corresponds to the tear-shaped curve in the 1/c parameterization.


> The inverse of the cardioid is the exterior of a teardrop shape. The circles on the outside of the cardioid are inverted to circles on the inside of the teardrop. The cusp of the cardioid becomes the cusp of the teardrop.


### z^2 -2.0+1.0/c

Third example : p =  -2.0+1.0/c

![](./png/LCM_c_inverted_2_5000_2.000000_5.000000_600.png)


### z^2 + 0.25+1.0/c

Fourth example : p =  0.25 + 1.0/c


  
  
![](./png/LCM_c_parabola_5000_4.000000_5.000000_600.png)

Description by  [Nikola Ubavić](https://ubavic.rs/work/julia_and_mandelbrot_set/?lang=sr)  
> If a 1/4 translation is performed before the inversion, then the cardioid is imaged in a parabola.


> The cusp of the cardioid is then moved to the origin, and that inverts to infinity. 
> the cardioid is sent to the outside of a parabola. 
> The circles on the outside of the cardioid are inverted to circles on the inside of the parabola.


### z^2 -1.401155 - 1.0/c

Fifth example: p =  -1.401155 - 1.0/c

  

![](./png/LCM_c_Myrberg_5000_1.330000_400.700000_600.png)



### Exponential map

Where
* c = c0 + cexp(x + i y) 
* x + i y = clog(c - c0)
* with scaling by 2pi/width







## p*z*(z-1) = logistic family


  m*z*(1.0-z)
  
![](./png/LCM_lambda_5000_1.000000_3.200000_600.png)


  z*(1.0-z)/m

![](./png/LCM_lambda_inverted_5000_0.000000_1.120000_600.png)


  (1 + 1/m)*z*(1.0-z)


![](./png/LCM_lambda_inverted_1_5000_0.000000_1.100000_600.png)  







# To do
* [Mercator Mandelbrot Maps by: Anders Sandberg](https://www.flickr.com/photos/arenamontanus/albums/72157615740829949) : a logarithmic projection around a point c0: z-> (log(|z-c0|), arg(z-c0)). The idea for this kind of logarithmic map of the set is from [David Madore](http://www.madore.org/~david/math/mandelbrot.html). While it is a bit of a stretch to call it a Mercator projection, it sounds better than just a log scale map. If the Mandelbrot set is assumed to be about meter-sized these zooms reach the size of an atomic nucleus.  
* transverse Mercator
  * [The North Transept Transfigurational Dreamwarp by Robert Lee](https://www.flickr.com/photos/robertlee99/50289710431/) : 


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