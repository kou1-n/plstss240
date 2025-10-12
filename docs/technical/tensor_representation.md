# Tensor Representation

This project stores stress and strain tensors in two different ways depending on the context.

## Internal calculations

Most subroutines work with second-order tensors as `3×3` matrices:

```
dimension str(3,3), sig(3,3)
dimension ctens(3,3,3,3)
```

Here `str` and `sig` hold the strain and stress tensors, while `ctens` is the fourth-order tangent modulus. Calculations such as element stiffness assembly or the constitutive update (`stress.f`, `elastc.f`, etc.) use these array shapes directly.

## Output and storage

When stresses or strains are saved per element, the program converts each `3×3` matrix to a six-component vector (Voigt notation). In `postpr.f` this mapping is explicit:

```
sigma(1,nel) = sts(1,1)
sigma(2,nel) = sts(2,1)
...
sigma(6,nel) = sts(3,3)
```

Arrays such as `sigma(6,nel)` and `epsln(6,nel)` therefore contain the tensor components in the order `(xx, yx, yy, zx, zy, zz)`.

## Summary

The finite element solver keeps tensors as `3×3` matrices for its internal updates and only converts them to six-element vectors when writing results. Understanding this mixed representation helps when reading the Fortran sources or interpreting the output files.
