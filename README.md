# k-εCC (Two-equation turbulence model with Spalart-Shur Curvature Correction Term)

Like other turbulence models, the k-ε doesn't fit for curvature flow problems with heavily streamline curvature, as the anisotropic part 
of Reynolds Stress Tensor is not resolved. The assesment of curvature correction term in the k-ε turbulence model was made based in 
the implementation of the kOmegaSSTCC that can be found in [kOmegaSSTCC](https://github.com/ancolli/kOmegaSSTCC). There are only minors 
changes to the code to fit the k-ε model equations.

