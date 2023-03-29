# Team12: The Micromechanics Project 



## Students

- Andreas Devogel: andreas.devogel@student.kuleuven.be
- Zirui Yan: zirui.yan@student.kuleuven.be 

## Description

The aim of this project is to compute the strain and the stress of a 3-dimensional composite material cuboid based on the periodic Lippman-Schiwinger equation. 

The algorithm we use is mainly based on [1] A fast numerical method for computing the linear and nonlinear mechanical properties of composites.

## Used libraries

FFTW: FFTW (Fastest Fourier Transform in the West) is a library for computing the discrete Fourier transform (DFT) in one or more dimensions, of both real and complex data.

Eigen: Eigen is a C++ template library for linear algebra.

GTest: Google Test (also known as GTest) is a C++ testing framework developed by Google.

Plotly: Plotly is an open-source data visualization python library that allows users to create interactive charts, graphs, and other visualizations in a web-based environment.

## How to setup and use this code

Clone this repository by following command:

```
git clone git@gitlab.kuleuven.be:math-eng/h0t46a/2023/team12.git
```
or
```
git clone https://gitlab.kuleuven.be/math-eng/h0t46a/2023/team12.git
```
Download cmake first if you do not have it on your computer. Since we've made cmakelists.txt for you, you can use cmake to build all the code by following command:
```
cmake -B Build
cd Build
cmake --build .
``` 
Then there are three executables in the Build folder. They are `main.exe`, `IOTest.exe`, `algorithmTest.exe`.

### main.exe

`main.exe` reads in the properties of a 3-dimensional composite material cuboid and the average strain of this cuboid and outputs either the stress or the strain after computation.
`main.exe` takes 7 command line arguments: in_mat, in_strain, out, e_or_s, ij, tol, maxit.

- in_mat: material-input file 
- in_strain: average-strain-input file
- out: output file
- e_or_s: indicate you want to ouput the strain or the stress
- ij: which component of strain tensor or stress tensor should be outut. ij could be one of xx, yy, zz, xy, yx, xz, zx, yz, zy
- tol: tolerance of error
- maxit: max iteration number

If you are in the team12 folder, you use command like this:
```
./main "../Input/test_mat.in" "../Input/test_e.in" "../Output/demo.out" e xx 1e-16 30
```

### IOTest.exe

`IOTest.exe` tests `IO.hpp`which contains the function for reading input files and writing output files.

If you are at the folder team12, you can simply execute  algorithmTest.exe by following command:

```
cd Build
./IOTest

```


### algorithmTest.exe

`algorithmTest.exe` tests each member function of `micromechanics.hpp` which contains a class named micromechanics to do the algorithm given in [1].

If you are at the folder team12, you can simply execute  algorithmTest.exe by following command:

```
cd Build
./algorithmTest

```



## Visualization

We use a python library named `plotly` for visualization.
Before using it, you have to install it and two other necessary libraries:
```
pip install plotly==5.13.1
pip install kaleido
pip install packaging

```

`visualization.py` takes commend line arguments: input file, output file, [options]. 

options: `-disp`, `-slice`

### without options

Commend line arguments: input file, output file

- input file
- output file: If output file != None, it will generate picture in object directory.

If you are in the Src folder, you use command like this:

```
python visualization.py "../Output/test_vis.out" "../Output/test_vis_noopt.png"
```

### -disp
Commend line arguments: input file, output file, -disp

- input file
- output file: If output file != None, it will both generate picture in object directory and in browser; if output file == None, it will only generate picture in browser.

Displays interactive visualization in browser. This option can be combined with the other options, however, if specified, it should always be the first. 

If you are in the Src folder, you use command like this::

```
python visualization.py "../Output/test_vis.out" "../Output/test_vis_noopt.png" -disp
```

### -slice
Commend line arguments: input file, output file, [-disp],-slice, dimension, slices, opacity

- dimension: x, y or z
- slices: a variable number of real numbers specifying the position of the slice in the specified dimension should be in the range of the read in material lengths
- opacity: the opacity of the visualization real number between 0 and 1, with 1 being non-opaque and 0 being completely seetrough

If you are in the Src folder, you use command like one of these:

```
python visualization.py "../Output/test_vis.out" "../Output/test_vis_slicex.png" -slice x 0 1 2 0.7
python visualization.py "../Output/test_vis.out" "../Output/test_vis_slicey.png" -slice y 0 1 2 0.7
python visualization.py "../Output/test_vis.out" "../Output/test_vis_slicez.png" -slice z 0 1.5 3 0.7
```




## Reference

- [1] H. Moulinec, Pierre Suquet. A fast numerical method for computing the linear and nonlinear mechanical properties of composites. Comptes Rendus de l’Académie des sciences. Série II. Mécanique, physique, chimie, astronomie, Gauthier-Villars, 1994. hal-03019226
- [2] H. Moulinec, Pierre Suquet. A numerical method for computing the overall response of nonlinear
composites with complex microstructure. Computer Methods in Applied Mechanics and Engineering,
1998, 157 (1-2), 10.1016/S0045-7825(97)00218-1. hal-01282728



