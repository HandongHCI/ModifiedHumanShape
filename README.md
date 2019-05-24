Modified MPII's Human Shape
=====

Please refer the  discription and source code of the original Human Shape from their [website (http://humanshape.mpi-inf.mpg.de/)](http://humanshape.mpi-inf.mpg.de/) and [GitHub (https://github.com/leonid-pishchulin/humanshape)](https://github.com/leonid-pishchulin/humanshape).

Related publication:

Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobalt and Bernt Schiele  
Building Statistical Shape Spaces for 3D Human Modeling
In _ArXiv_, March 2015


Installation Guideline for the modified MPII’s Human Shape
---

1. The original MPII’s Human Shape algorithm uses a well-known optimization library, L-BFGS-B which was originally(?) developed based on Fortran language. To use Fortran language directly in Matlab, we need ‘Intel Visual Fortran Compiler’ which is not for free (price: 400 USD per year). Instead of using the L-BFGS-B included in MPII's source code, I used a different L-BFGS-B library which was converted to C language and shared through [Matlab FileExchange (https://nl.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb--l-bfgs-b--mex-wrapper)](https://nl.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb--l-bfgs-b--mex-wrapper) and [Github (https://github.com/stephenbeckr/L-BFGS-B-C)](https://github.com/stephenbeckr/L-BFGS-B-C) by Stephen Becker. Also, he made a ‘Mex wrapper’ which is to create a mex file to use the C code in Matlab programming. By the way, do not download the source code from the Matlab File Exchange. That is an old version. Download the L-BFGS-B from his GitHub (https://github.com/stephenbeckr/L-BFGS-B-C). But, this ‘Mex wrapper’ is included in my GitHub, so you don’t actually need to download it.



1. Set `MATLAB_HOME` in following make files:

    ```
    external/lbfgsb-for-matlab/Makefile
    shapemodel/Makefile
    evaluation/statQuality/align.mk
    evaluation/statQuality/evaluation.mk
    ```
2. Switch to the top level directory of the source code, issue `make` from the command line

Getting the models
---

1. Download the models
```
    wget http://datasets.d2.mpi-inf.mpg.de/humanshape/caesar.zip
    wget http://datasets.d2.mpi-inf.mpg.de/humanshape/caesar-norm-wsx.zip
    wget http://datasets.d2.mpi-inf.mpg.de/humanshape/caesar-norm-nh.zip
```

2. Unzip the models
```
    unzip caesar.zip && rm -f caesar.zip
    unzip caesar-norm-wsx.zip && rm -f caesar-norm-wsx.zip
    unzip caesar-norm-nh.zip && rm -f caesar-norm-nh.zip
```

Getting the fitted meshes
---

1. Download the fitted meshes
```
    wget http://datasets.d2.mpi-inf.mpg.de/humanshape/caesar-fitted-meshes.zip
    wget http://datasets.d2.mpi-inf.mpg.de/humanshape/caesar-norm-wsx-fitted-meshes.zip
    wget http://datasets.d2.mpi-inf.mpg.de/humanshape/caesar-norm-nh-fitted-meshes.zip
```

2. Unzip the fitted meshes
```
    unzip caesar-fitted-meshes.zip && rm -f caesar-fitted-meshes.zip
    unzip caesar-norm-wsx-fitted-meshes.zip && rm -f caesar-norm-wsx-fitted-meshes.zip
    unzip caesar-norm-nh-fitted-meshes.zip && rm -f caesar-norm-nh-fitted-meshes.zip
```

Running
---
1. Start matlab
2. Edit file `fitting/expParams.m`
   	1) point `p.rootDir` to the full path to the source code directory
   	2) point `p.modelInDir` to the model directory, e.g. `caesar/`
3. Run `demo`

TODO
---
Add evaluation code for 

1. Per-vertex mean fitting accuracy
1. Total fitting accuracy
1. Compactness, generalization and specificity

If you have any questions, send an email to leonid@mpi-inf.mpg.de with a topic "humanshape".
