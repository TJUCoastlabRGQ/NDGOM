# NDG-FEM

## Introduction

This software uses Nodal Discontinuous Galerkin Finite Element Methods (NDG-FEM) to solve advection, advection-diffusion and various hydraulic problems in one, two and three dimension. 
The software is written in Matlab and C languages and uses OpenMP for parallelization.
The simulation results can be exported in the nc format and vtk format, while the visualization is also supported by Matlab.

## Install

Part of the software codes are written in C Mex functions. 
Before the users start their simulation, these C functions need to be compiled first.

### 1. Windows

* Install TDM-GCC software, including OpenMP package
* configure TDM-GCC as the default compiler for Matlab
* run script `NdgSetup.m`

### 2. OS X

* Install Intel C compiler
* configure `icc` as the default compiler for Matlab
* run script `NdgSetup.m`

## Usages

Use following commands to run a two-dimensional advection probelm

```
adv = ConstAdvUniformMesh2d( 1, 60, enumStdCellType.Quad );
adv.matSolve;
pos = makeNdgPostProcessFromNdgPhys( adv );
pos.drawResult( 1, 1, 0 )
```

![](http://ww1.sinaimg.cn/large/7a1c18a8ly1frtnec7hwfj21st1cmq9i.jpg)

## Acknowledgment

This software is inspired by _Nodal Discontinuous Galerkin Methods Algorithms, Analysis, and Applications Texts in Applied Mathematics_ (Hesthaven and Tim Warburton, 2008), please refer to [NodalDG](http://www.caam.rice.edu/~timwar/Book/NodalDG.html) for more details.
