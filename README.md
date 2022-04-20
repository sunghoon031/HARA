# HARA: A Hierarchical Approach for Robust Rotation Averaging

[Paper](https://arxiv.org/abs/2111.08831), [Supplementary material](https://github.com/seonghun-lee/seonghun-lee.github.io/blob/master/pdf/SupplementaryMaterial_HARA_A_Hierarchical_Approach_for_Robust_Rotation_Averaging.pdf)

In this repository, we provide the implementation of HARA. If you use our code, please cite it as follows:

````
@article{hara,
author = {Seong Hun Lee and Javier Civera},
title = {HARA: A Hierarchical Approach for Robust Rotation Averaging},
journal   = {CoRR},
volume    = {abs/2111.08831},
year      = {2021},
}
````

### Quick start
Run `Test_HARA.m` to try it on a synthetic data WITHOUT using the number of inlier matches.

### Main functions:
1. `CreateSyntheticData.m`: Generate a synthetic dataset, as described in the main paper.
2. `RunHARA.m`: Run HARA without using the number of inlier matches.
3. `RunHARA_usingNumberOfInlierMatches.m`: Run HARA using the number of inlier matches. Note that we only provide the function, without a sample dataset. This function is quite similar to `RunHARA.m`, so it shouldn't be too difficult to use on your own dataset.
