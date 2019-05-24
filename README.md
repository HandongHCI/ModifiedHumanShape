Modified MPII's Human Shape
=====

Please refer the  discription and source code of the original Human Shape from their [website](http://humanshape.mpi-inf.mpg.de/) and [GitHub](https://github.com/leonid-pishchulin/humanshape).

Related publication: **Pishchulin, L., Wuhrer, S., Helten, T., Theobalt, C., and Schiele, B. (2015). Building Statistical Shape Spaces for 3D Human Modeling. _ArXiv_** ([see](http://arxiv.org/abs/1503.05860))


Installation guideline for the modified MPII's Human Shape
---

1. The original MPII's Human Shape algorithm uses a well-known optimization library, L-BFGS-B, which was initially(?) developed based on Fortran language. To use Fortran language directly in Matlab, we need 'Intel Visual Fortran Compiler' which is not for free (price: 400 USD per year). Instead of using the L-BFGS-B included in MPII's source code, I used a different L-BFGS-B library which was converted to C language and shared through [Matlab FileExchange](https://nl.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb--l-bfgs-b--mex-wrapper) and [Github](https://github.com/stephenbeckr/L-BFGS-B-C) by Stephen Becker. Do not download the source code from the Matlab File Exchange, which is an old version. Download the L-BFGS-B from his GitHub (https://github.com/stephenbeckr/L-BFGS-B-C). By the way, the L-BFGS-B is included in the modified Human Shape, so you don't actually need to get it.


1. Clone or download the modified MPII's Human Shape from my [GitHub](https://github.com/HandongHCI/humanshape).

1. Install `Microsoft Visual Studio Pro 2017`. It should be the `Professional` version, not the `Community` version. When you install the Visual Studio Pro 2017, `VC++ 2017 toolset` and `Windows 10 SDK` should be checked.

1. Run `compile_mex.m` file which is located in the folder `\external\lbfgsb_C\`, which is Stephen Becker's work.

1. In your Matlab, go to the folder same to where `demo.m` locates. Then, copy and paste the following code in the Matlab command line.

    ```
    mex -output shapepose.mexw64 -Ishapemodel\lib\nr\ -Ishapemodel\lib\include\ "shapemodel\shapepose.cpp" "shapemodel\Show.cpp" "shapemodel\NMath.cpp" "shapemodel\NRBM.cpp" "shapemodel\paramMap.cpp" "shapemodel\CTMesh-30DOF.cpp"
    ```

1. Copy and paste the following code in the Matlab command line, too.
    ```
    mex -largeArrayDims -output rigidAlign.mexw64 -Imatlabroot\extern\include\ "evaluation\statQuality\rigidAlign.cpp" "evaluation\statQuality\GeneralizedProcrustes.cpp" -Lmatlabroot\extern\lib\win64\microsoft\ -llibmwblas.lib -llibmwlapack.lib
    ```

Notes
---

1. Get necessary models and fitted meshes from http://humanshape.mpi-inf.mpg.de/

1. `fitting\NRD.m` file has been revised to use Stephen Becker's `LBFGSB_C`.

1. In `fitting\getOptionsOptimizer.m`, `options.UseParallel = 1;` is added to run the Matlab optimization with Parallel Computing Toolbox. If you don't have this toolbox, please delete or comment this line, or make the value 0 (`options.UseParallel = 0;`).