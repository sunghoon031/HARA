# HARA: A Hierarchical Approach for Robust Rotation Averaging

[Paper](https://arxiv.org/abs/2111.08831), [Video](https://www.youtube.com/watch?v=oAR-LMStRS4), [Supplementary material](https://github.com/seonghun-lee/seonghun-lee.github.io/blob/master/pdf/SupplementaryMaterial_HARA_A_Hierarchical_Approach_for_Robust_Rotation_Averaging.pdf)

In this repository, we provide the implementation of HARA. If you use our code, please cite it as follows:

````
@InProceedings{Lee_2022_CVPR,
    author    = {Lee, Seong Hun and Civera, Javier},
    title     = {{HARA}: A Hierarchical Approach for Robust Rotation Averaging},
    booktitle = {Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR)},
    month     = {June},
    year      = {2022},
    pages     = {15777--15786}
}
````
### Update (January 8th 2023)

[This recent change in the code](https://github.com/sunghoon031/HARA/commit/00d88296a6be1e2693d4f1e50397f84bad21003b) leads to a significant speedup of the local optimization step compared to the version used for the CVPR paper. We thank Chitturi Sidhartha, the first author of a 3DV paper titled ['It Is All In The Weights: Robust Rotation Averaging Revisited'](https://ieeexplore.ieee.org/document/9665962), for the discussion that led to this finding.

### Quick start
Run `Test_HARA.m` to try it on a synthetic data WITHOUT using the number of inlier matches.

### Main functions:
1. `CreateSyntheticData.m`: Generate a synthetic dataset, as described in the main paper.
2. `RunHARA.m`: Run HARA without using the number of inlier matches.
3. `RunHARA_usingNumberOfInlierMatches.m`: Run HARA using the number of inlier matches. Note that we only provide the function, without the test script or a sample dataset. This function is quite similar to `RunHARA.m`, so it shouldn't be too difficult to use on your own dataset.
