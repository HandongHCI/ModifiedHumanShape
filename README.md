Modified MPII's Human Shape (updated: January 2018)
=====

Please refer the  discription and source code of the original Human Shape from MPII's [website](http://humanshape.mpi-inf.mpg.de/) and [GitHub](https://github.com/leonid-pishchulin/humanshape).

Related publication: **Pishchulin, L., Wuhrer, S., Helten, T., Theobalt, C., and Schiele, B. (2015). Building Statistical Shape Spaces for 3D Human Modeling. _ArXiv_** ([read](http://arxiv.org/abs/1503.05860))


Installation guideline for the modified MPII's Human Shape
---

1. The original MPII's Human Shape algorithm uses a well-known optimization library, L-BFGS-B, which was initially(?) developed based on Fortran language. To use Fortran language directly in Matlab, we need 'Intel Visual Fortran Compiler' which is not for free (price: 400 USD per year). Thanks to [Dr. Stephen Becker](http://amath.colorado.edu/faculty/becker/), a professor at University of Colorado at Boulder, I could use the L-BFGS-B algorithm without purchasing the Intel Visual Fortran Compiler. Instead of using the L-BFGS-B included in MPII's source code, I used Dr. Becker's C-based L-BFGS-B library (with the mex compiler for Matlab), which is shared through [Matlab FileExchange](https://nl.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb--l-bfgs-b--mex-wrapper) and [Github](https://github.com/stephenbeckr/L-BFGS-B-C). Dr. Becker's L-BFGS-B is included in my modified Human Shape, so you don't need to get it.


1. Clone or download the **modified MPII's Human Shape** from my [GitHub](https://github.com/HandongHCI/humanshape).

1. Install `Microsoft Visual Studio Pro 2017`. It should be the `Professional` version, not the `Community` version. When you install the Visual Studio Pro 2017, `VC++ 2017 toolset` and `Windows 10 SDK` should be checked.

1. In Matlab, run `compile_mex.m` located in the folder `\external\lbfgsb_C\Matlab`, which is Dr. Becker's work.

1. In Matlab, go to the folder same to where `demo.m` locates. Then, copy and paste the following code in the Matlab command line.

    ```
    mex -output shapepose.mexw64 -Ishapemodel\lib\nr\ -Ishapemodel\lib\include\ "shapemodel\shapepose.cpp" "shapemodel\Show.cpp" "shapemodel\NMath.cpp" "shapemodel\NRBM.cpp" "shapemodel\paramMap.cpp" "shapemodel\CTMesh-30DOF.cpp"
    ```

1. Copy and paste the following code in the Matlab command line, too.
    ```
    mex -largeArrayDims -output rigidAlign.mexw64 -Imatlabroot\extern\include\ "evaluation\statQuality\rigidAlign.cpp" "evaluation\statQuality\GeneralizedProcrustes.cpp" -Lmatlabroot\extern\lib\win64\microsoft\ -llibmwblas.lib -llibmwlapack.lib
    ```

    ```
    mex -largeArrayDims -output ErrorEvaluation.mexw64 -Imatlabroot\extern\include\ "evaluation\statQuality\ErrorEvaluation.cpp" "evaluation\statQuality\GaussVector.cpp" "evaluation\statQuality\patternRecognitionPCA.cpp" "evaluation\statQuality\UnsupervisedLearning.cpp" "evaluation\statQuality\Mle.cpp" -Lmatlabroot\extern\lib\win64\microsoft\ -llibmwblas.lib -llibmwlapack.lib
    ```

1. Get necessary models (approx. 1.8 GB) from http://humanshape.mpi-inf.mpg.de/. These are necessary to run the code. Put the unzipped models in `\experiments\models\`, then edit `p.modelInDir` variable in `fitting\expParams.m` to point one of the models (e.g., `models/caesar`, but it is already set in the modification of Human Shape).

1. Run `demo.m`.

Notes
---

1. `fitting\NRD.m` file has been revised to use Dr. Becker's `L-BFGS-B`.

1. In `fitting\getOptionsOptimizer.m`, `options.UseParallel = 1;` is added to run the Matlab optimization with Parallel Computing Toolbox in order to increase the calculation speed. If you don't have this toolbox, please delete or comment this line, or make the value 0 (`options.UseParallel = 0;`).


To-do
---

1. Change of the template model (fullbody, also hand) with different landmark sets.


Q&A
---
Dr. Wonsup Lee (W (dot) Lee (at) Handong (dot) edu) at Handong Global University
GitHub: https://handonghci.github.io/WonsupLee/
ResearchGate: https://www.researchgate.net/profile/Wonsup_Lee