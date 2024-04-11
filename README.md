# Marching cubes

## Requires 
- [CMake](https://cmake.org/) (version > 2.8)
- [Vofi](https://github.com/VOFTracking/Vofi)

## Compilation of marching-cubes
To compile, follow these steps:
```bash
git clone https://git.multiflow.org/software-toolbox/marching-cubes.git
cd marching-cubes
VOFI_DIR=<path_to_vofilibrary> MC_DIR=$PWD ./compile
```

## Execution
`./build/marching-cubes <CASE> <N> <INTERPOLATION>`

Choose from:

- CASE = PLANE -- SPHERE -- ELLIPSOID -- ORTHOCIRCLE -- SINWAsVE1 -- JET -- DODECAHEDRON
- N = number of cells in each direction
- INTERPOLATION = MASON -- LINEAR -- MIDDLE
