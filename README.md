# Mandelbrot-Sets-Alternate-Parameter-Planes


[Julia and Mandelbrot Sets. Alternate Parameter Planes by David E. Joyce © 1994.](https://mathcs.clarku.edu/~djoyce/julia/altplane.html)


Images of [complex quadratic polynomials](https://en.wikipedia.org/wiki/Complex_quadratic_polynomial)


## z^2+p family

   z^2 + c


![](./png/LCM_c_5000_-0.750000_1.500000_600.png)


  z^2 + 1/c
  
![](./png/LCM_c_inverted_5000_1.330000_2.700000_600.png)


  z^2 -2.0+1.0/c


![](./png/LCM_c_inverted_2_5000_2.000000_5.000000_600.png)


  z^2 + 0.25+1.0/c
  
![](./png/LCM_c_parabola_5000_4.000000_5.000000_600.png)


  z^2 -1.401155 - 1.0/c

![](./png/LCM_c_Myrberg_5000_1.330000_400.700000_600.png)

## p*z*(z-1) = logistic family


  m*z*(1.0-z)
  
![](./png/LCM_lambda_5000_1.000000_3.200000_600.png)


  z*(1.0-z)/m

![](./png/LCM_lambda_inverted_5000_0.000000_1.120000_600.png)


  (1 + 1/m)*z*(1.0-z)


![](./png/LCM_lambda_inverted_1_5000_0.000000_1.100000_600.png)  


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

