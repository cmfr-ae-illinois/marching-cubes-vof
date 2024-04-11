# Marching cubes

This project contains an implementation of the Marching-Cubes algorithm of [Lorensen and Cline (1987)](https://doi.org/10.1145/37402.37422) for reconstructing a surface from a discrete indicator function.

The code is largely taken from [Paul Bourke's website](https://paulbourke.net/geometry/polygonise/).

Three types of vertex interpolation are implemented:
- `MIDDLE`: the vertex is always positioned in the middle of the intersected edge.
- `LINEAR`: the vertex position is linearily interpolated from the end-point values of the discrete indicator function.
- `MANSON`: using the method proposed by [Manson, Smith, and Schaefer (2011)](http://dx.doi.org/10.1111/j.1467-8659.2011.01869.x).

## Requires 
- [CMake](https://cmake.org/) (version > 2.8)
- [Vofi](https://github.com/VOFTracking/Vofi)

## Compilation of marching-cubes
To compile, follow these steps:
```bash
mkdir build
cd build
VOFI_DIR=<path_to_vofi_library> MC_DIR=$PWD/.. cmake ..
make
```

## Execution
`marching-cubes <CASE> <N> <INTERPOLATION>`

Choose from:

- `CASE` = `PLANE` or `SPHERE` or `ELLIPSOID` or `ORTHOCIRCLE` or `SINWAVE1` or `JET` or `DODECAHEDRON`
- `N` = number of cells in each direction
- `INTERPOLATION` = `MANSON` or `LINEAR` or `MIDDLE`
