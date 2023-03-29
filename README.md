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
Doweload cmake first if you do not have it on your computer. Since we've made cmakelists.txt for you, you can use cmake to build all the code by following command:
```
cmake -B Build
cd Build
cmake --build .
``` 
Then there are three executables in the Build folder. They are `main.exe`, `IOTest.exe`, `algorithmTest.exe`.

### main.exe

`main.exe` reads in the properties of a 3-dimensional composite material cuboid and the average strain of this cuboid and outputs either the stress or the strain after computation.
`main.exe` takes 

### IOTest.exe
`IOTest.exe` tests `IO.hpp`. 

### algorithmTest.exe
`algorithmTest.exe` tests `micromechanics.hpp` which contains a class named micromechanics to do the algorithm given in [1].

## Visualization

We use a python library named `plotly` for visualization.
Before using it, you have to install it and two other necessary libraries:
```
pip install plotly==5.13.1
pip install kaleido
pip install packaging

```
`visualization.py` takes 2 commend line arguments. The first one is the data 
``` 
cd Src
python visualization.py "../Output/test.out" "../Output/test.png"
```


## Reference

- [1] H. Moulinec, Pierre Suquet. A fast numerical method for computing the linear and nonlinear mechanical properties of composites. Comptes Rendus de l’Académie des sciences. Série II. Mécanique, physique, chimie, astronomie, Gauthier-Villars, 1994. hal-03019226
- [2] H. Moulinec, Pierre Suquet. A numerical method for computing the overall response of nonlinear
composites with complex microstructure. Computer Methods in Applied Mechanics and Engineering,
1998, 157 (1-2), 10.1016/S0045-7825(97)00218-1. hal-01282728



