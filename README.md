# SingularityFree

**FourDimensionalProjection.wl**, **FourDimensionalProjection2.wl**, and **LocalRing.wl** are three Mathematica packages for constructing **singularity-free bases** of Master Integrals. These packages implement algorithms to find Feynman integral bases where the reduction of any integral through integration-by-parts (IBP) relations does not produce coefficients that are singular as the dimensional regularization parameter $\epsilon \to 0$ (where $d = d_0 - 2\epsilon$, with $d_0$ being a positive integer representing the spacetime dimension).

## Background

Standard IBP reduction methods typically yield Feynman integral bases where the reduction of some integrals gives rise to coefficients that become singular as $\epsilon \to 0$. These singular coefficients can also appear in scattering amplitudes, obscuring their structure and making their evaluation more complicated. This repository provides tools to construct bases that are free from such singularities.

## Algorithms

These packages implement two general algorithms for constructing singularity-free bases, as described in the accompanying paper "Singularity-Free Feynman Integral Bases":

1. **Sequential $d=4$ IBP reduction**: Constructs a basis iteratively by projecting onto the finite part of the set of IBP relations. This algorithm has two different implementations involving either step 6 or step 6'.

2. **Gaussian elimination in a local ring**: Performs elimination within a local ring that forbids division by $\epsilon$ while permitting division by polynomials in $\epsilon$ that are finite at $\epsilon=0$.

## Package Contents

- **FourDimensionalProjection.wl**: Implementation of the sequential $d=4$ IBP reduction algorithm using step 6
- **FourDimensionalProjection2.wl**: Alternative implementation of the sequential $d=4$ IBP reduction algorithm using step 6'
- **LocalRing.wl**: Implementation of the Gaussian elimination in local ring algorithm

## Usage

These Mathematica packages are designed to work with existing IBP reduction workflows to produce singularity-free master integral bases that facilitate cleaner amplitude calculations and avoid spurious singularities in intermediate steps. The packages integrate easily with [NeatIBP](https://github.com/yzhphy/NeatIBP) for streamlined IBP reduction workflows.

### Main Functions

#### FourDimensionalProjection.wl
```mathematica
(* Load the package *)
Get["FourDimensionalProjection.wl"]

(* Main function call - using step 6 *)
result = finiteCoefficients[ibps, G, eps]
```

**Inputs:**
- `ibps`: List of integration-by-parts relations. It is assumed that all the relations are linearly independent, *i.e.* the corresponding matrix is of full rank.
- `G`: Symbol representing the Feynman integrals (e.g., `G[a1, a2, ...]` where `ai` are propagator powers)
- `eps`: Dimensional regularization parameter (typically $\epsilon$)

**Output:** List of master integrals forming a singularity-free basis

#### FourDimensionalProjection2.wl
```mathematica
(* Load the package *)
Get["FourDimensionalProjection2.wl"]

(* Main function call - using step 6' *)
result = finiteCoefficients[ibps, G, eps]

(* Additional helper functions available: *)
pMatrix = permutationMatrix[pVector]

{reduced, inverse} = gaussReductionMatrix[matrix]
```

**Main Function Inputs:** Same as FourDimensionalProjection.wl
- `ibps`: List of integration-by-parts relations  
- `G`: Symbol representing the Feynman integrals
- `eps`: Dimensional regularization parameter

**Additional Functions:**
- `permutationMatrix[p]`: represents the permutation matrix given by permutation vector `p` as a sparse array(compatibility function for older Mathematica versions, before version 13.1).
    - `p`: is list of numbers from $1$ to $n$, representing the permutation of $n$ elements.
- `gaussReductionMatrix[matrix]`: Performs Gaussian-Jordan elimination and returns both the reduced matrix and the (pseudo)inverse matrix.
    - `matrix` is the matrix to be reduced to row echelon,
    - `reduced` is the row echelon form of `matrix`,
    - `inverse` is the (pseudo)inverse of `matrix`.

**Output:** List of master integrals forming a singularity-free basis (alternative implementation)

#### LocalRing.wl
```mathematica
(* Load the package *)
Get["LocalRing.wl"]

(* Main functions for Gaussian elimination in local ring *)
{reducedMatrix, columnIndex} = GaussianElimination[matrix, "PivotStrategy" -> "OriginalIfPossible"]
```

**Inputs for GaussianElimination:**
- `matrix`: matrix to be reduced
- `"PivotStrategy"`: either `"OriginalIfPossible"` (default) or `"Markowitz"`

**Output:** Matrices reduced using local ring operations that avoid division by `eps`

## Reference

For detailed background and algorithm descriptions, please refer to the accompanying paper "Singularity-Free Feynman Integral Bases" by S. De Angelis, D. Kosower, R. Ma, Z. Wu, and Y. Zhang.

## Copyright

© 2025 S. De Angelis, D. Kosower, R. Ma, Z. Wu, and Y. Zhang

This software is provided for academic and research purposes. If you use this software in your research, please cite the accompanying paper "Singularity-Free Feynman Integral Bases".