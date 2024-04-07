# 2021 UPDATE

A new version (DISCCOFAN) that handles more advanced operations and can process more complex data sets is available at https://github.com/sgazagnes/disccofan. This version is not updated anymore


# Project Title

Distributed Component Forests in 2-D: Hierarchical Image Representations Suitable for Tera-Scale Images.

From: International Journal of Pattern Recognition and Artificial Intelligence (2019): 
https://www.worldscientific.com/doi/pdf/10.1142/S0218001419400123


## Purpose

This program filters 2D data set using distributed-memory distributed component tree forests. 


# Getting Started


## Prerequisites

Libraries needed: OpenMPI, FreeImage, CFITSIO, HDF5

## Installing

The compilation is done with 

```
make
```

You might want to adapt the Makefile if some libraries are installed in specific locations.


## Running the program

To run a simple area filtering:

```
mpirun -np [NUM_PROC] ./area -g [X],[Y] --inprefix [imgprefix] --intype [imgtype] --lambda [L] 
```

To run a DAP segmentation:
```
mpirun -np [NUM_PROC] ./dap -g [X],[Y] --inprefix [imgprefix] --intype [imgtype] 
```


To run a 1D Pattern spectra
```
mpirun -np [NUM_PROC] ./pattern -g [X],[Y] --inprefix [imgprefix] --intype [imgtype] 
```


The grid option (**-g**) is REQUIRED and defines the grid division of the data set. If the dataset is divided in 2 parts horizontally and 2 parts vertically then [X] = 2, [Y] = 2 (-g 2,2), and [NUM_PROC] = X*Y = 4.

The **--inprefix** options define the prefix of your images. When dealing with classical image types (png, pgm, jpg, tif...), your image needs to be divided in tiles before running the program. Their names must be [imgprefix]-[id].[imgtype], with [id] the tile number.

To divide an image into tiles, you can use ImageMagick toolbox, which works as follow (assuming your input image is lena.pgm):
```
convert -crop 2x2@ lena.pgm tile-%d.pgm
```


This will create four files names tile-0.pgm, tile-1.pgm, tile-2.pgm and tile-3.pgm. To run the program, you can run

```
mpirun -np 4 ./area -g 2,2 --inprefix tile --intype pgm --lambda 100
```

That will give you four output files names *out-0.pgm*, *out-1.pgm*, *out-2.pgm*, *out-3.pgm*

If you are dealing with fits or h5 file, you do not need to divide the data manually. You can directly run the program on the full data set. You need to specify the dataset to be read in the h5 file with the option **--dataset**. If the file is lena.h5,and the dataset "data", hence the program call must be:

```
mpirun -np 4 ./area -g 2,2 --inprefix lena --intype h5 --lambda 100 --dataset data
```


The **--intype** option specifies the type of your dataset. If you use h5 dataset, use the option -b to specify the dynamic range of the dataset to be processed

The **--lambda** option specifies the threshold to be applied for the area opening. If you are running the DAP and pattern spectra programs,
then you need to change the lvec.txt file, which details the list of threshold to use to perform the multi-scale analysis. The first value
is the number of thresholds to apply, and the following numbers are the values of the thresholds (need to start with 0).

You can also specify the prefix of the output images, with the **--outprefix** option (by default is "out"). In the pattern spectrum, the results are written in a file called *pattern.txt*
It should be also possible to change the output type of the image by specifying the **--outtype** option (by default is the sane as the input).

There is a verbosity option to display some information while the program is running. 
**--verbosity** (**-v**) [verbosity]: (OPTIONAL) Add verbose output
							"off": No verbose output, except for warnings and errors (default)
							"info": Addtionnal information from proc 0
							"debug": Follow-up of every steps in each process
							"all": Equivalent to debug

## Scaling performance
On a 1.4 by 1.4 Gpixels 8-bit grayscale image of Haiti, decreasing image size: 

<img src="https://github.com/sgazagnes/myimages/blob/main/haitispeedup.png" width="350" />

On a Haiti 16-bit data set constant tile size, increasing image size:

<img src="https://github.com/sgazagnes/myimages/blob/main/Haiti16exp.png" width="350" />

On a Napoli 8 by 8 Gpixels 16-bit data, decreasing image size: 

<img src="https://github.com/sgazagnes/myimages/blob/main/napoli.png" width="350" />

Memory Usage on  CPU 0:

<img src="https://github.com/sgazagnes/myimages/blob/main/memusage.png" width="350" />


## Version

Version 1.0

## Authors

Simon Gazagnes, Michael H.F. Wilkinson


	
