# Fitting 2D histogram with linear function
## Brief Explanation

With this library, you can fit 2D histograom with a linear function. This fit method is based on principal component analysis.


## Requirement
- ROOT 6
- Eigen 3
- C++ 11

## How to Build

Most probably, you already installed ROOT 6. So, here we start from introducing how to install Eigen 3.

At first, do git clone wherever you want to install.
```
$ git clone https://gitlab.com/libeigen/eigen.git
```
And after going to directory of Eigen, type following commands.
```
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=${PATH_TO_INSTALL_DIR} ..
$ make
$ make install
```
You can specify install directory by `${PATH_TO_INSTALL_DIR}`. When this is omitted, install directory is usually `/usr/local/include`.

After installing Eigen 3, go to directory of this library `PrincipalFit`. And then,

```
$ mkdir build
$ cd build
$ cmake -DEIGEN_DIR=${EIGEN_INCLUDE_DIR}/eigen3 ..
```
Please specify the path to include directory of Eigen3 by ${EIGEN_INCLUDE_DIR}. (In case you install Eigen3 at `/usr/local/include/`, cmake would automatically search for install directory of Eigen3.)

Then, type
```
$ make
```

## How to Use

Start ROOT at first and type command,
```
root[0] .x ${PrincipleFit_DIR}/macro/load.C
(int) 0
(int) 0
```

Then, draw some 2D histogram. Like
```
root[1] tree->Draw("var1:var2>>h1(100, -1, 1, 200, -5, 5)", "", "colz")
```
To fit this histogram,
```
root[2] pca_fit("h1", -0.5, 0.5, 2, 2)
```
The 2D histogram should be fitted by 2 linear functions. Since the name of the histogram is "h1" in above example, the first argument should be "h1". The argument meanings are
```c++
pca_fit(const char* histname, minX, maxX, minY, maxY)
```

You should select the fit region by the arguments  `minX`, `maxX`, `minY`, `maxY`.

