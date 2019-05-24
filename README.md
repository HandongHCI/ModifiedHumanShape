Modified MPII's Human Shape
=====

Please refer the  discription and source code of the original Human Shape from their [website (http://humanshape.mpi-inf.mpg.de/)](http://humanshape.mpi-inf.mpg.de/) and [GitHub (https://github.com/leonid-pishchulin/humanshape)](https://github.com/leonid-pishchulin/humanshape).

Related publication: **Pishchulin, L., Wuhrer, S., Helten, T., Theobalt, C., and Schiele, B. (2015). Building Statistical Shape Spaces for 3D Human Modeling In _ArXiv_** ([download](http://arxiv.org/abs/1503.05860))


Installation guideline for the modified MPII's Human Shape
---

1. The original MPII's Human Shape algorithm uses a well-known optimization library, L-BFGS-B, which was initially(?) developed based on Fortran language. To use Fortran language directly in Matlab, we need 'Intel Visual Fortran Compiler' which is not for free (price: 400 USD per year). Instead of using the L-BFGS-B included in MPII's source code, I used a different L-BFGS-B library which was converted to C language and shared through [Matlab FileExchange (https://nl.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb--l-bfgs-b--mex-wrapper)](https://nl.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb--l-bfgs-b--mex-wrapper) and [Github (https://github.com/stephenbeckr/L-BFGS-B-C)](https://github.com/stephenbeckr/L-BFGS-B-C) by Stephen Becker. Also, he made a `Mex wrapper` which is to create a mex file to use the C code in Matlab programming. By the way, do not download the source code from the Matlab File Exchange. That is an old version. Download the L-BFGS-B from his GitHub (https://github.com/stephenbeckr/L-BFGS-B-C). But, this `Mex wrapper` is included in my GitHub, so you don't actually need to download it.


1. Clone or download the modified MPII's Human Shape from this GitHub.

1. Install `Microsoft Visual Studio Pro 2017`. It should be the `Professional` version, not the `Community` version. When you install the Visual Studio Pro 2017, `VC++ 2017 toolset` and `Windows 10 SDK` should be checked.

1. Run `compile_mex.m` file which is located in the folder `\external\lbfgsb_C\`.

1. In your Matlab, go to the folder same to where `demo.m` locates.

1. Copy and paste the following code in the Matlab command line.

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
