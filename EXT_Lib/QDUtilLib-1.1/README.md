# QDUtilLib

List of Fortran modules which contain utilities for codes.
All modules contain a testing unit.

It has been tested with:

- gfortran (9.1.0 on linux, 12.1.0 on macos)
- ifort (19.1.3.304 on linux)

## 1) Installation and testing

### a) with a makefile:

To build the library, **libQD_gfortran_opt1_omp1.a**, with the default options (OPt=1, OMP=1, LAPACK=1)
```bash
make lib
```
It creates a library with the folowing name: **libQD_XXX_optY_ompZ.a**
with XXX, the compiler name (like gfortran), Y or Z the value O or 1
For instance, the default library is: **libQD_gfortran_opt1_omp1.a**
The module file (.mod) are in the OBJ/obj__XXX_optY_ompZ directory.

Two options to clean:
```bash
make clean
```

Remove some files, but keep the libraries, **libQD_XXX_optY_ompZ.a**

```bash
make cleanall
```

Remove all files

To test the module, in TESTS directory, run

```bash
./run_tests.sh
```

The tests are running with gfortran and several option combinations:

- OPT=1 or 0: compilation with optimization or without optimization
- OMP=1 or 0: with or without openmp
- LAPACK=1 or 0: with or without blas and lapack libraries

The file, **ALL_Tests.log**, contains a summary of all the tests.

### b) with fpm:

To build the library, **libQDUtilLib.a**, with the Lapack (LAPACK=1)
```bash
fpm build
```
The library is in **build/gfortran_xxxx/QDUtilLib** directory.
To remove Lapack library, the **fpm.tom** file must be edited.

Two options to clean:
```bash
fpm clean
```

Remove the build directory.

Remove some files, but keep the libraries, **libQDUtilLib.a**

```bash
fpm cleanall
```

To test the module:

```bash
fpm test | grep TESTING | grep Number
```

To run an example:

```bash
fpm run AppQDLib 
```

## 2) List of modules

### NumParameters

This module contains:

- **kind** definitions from the intrinsinc ISO_FORTRAN_ENV module
  - Rkind: default real kind (normaly real64)
  - Rk4, Rk8, Rk16 for real32, real64, real128
  - Ikind: default integer kind (usualy int32)
  - ILkind: default long integer kind (usualy int64)
  - Ik4, Ik8 for int32, int64
- Some real (kind=Rkind) numbers : ZERO to TWELVE, HALF, HUNDRED,ONETENTH ..., PI
- Some complexe (kind=Rkind) numbers : EYE (i), CZERO, CONE, CHALF
- out_unit and in_unit: the standard output and input units from the intrinsinc ISO_FORTRAN_ENV module
- print_level: level of printing:     0 minimal, 1 default, 2 large, -1 nothing
To change print_level:
CALL set_print_level(prtlev), where **prtlev** is an integer.

### String

This module contains functions and subroutines to manipulate character string (character(len=)):

- Functions to change to upper to lowercase or the reverse: **TO_lowercase**, **TO_uppercase**
- Functions to convert numbers in string: **TO_string**
  All kinds defined in NumParameters module and logical type are possible
  For real convertion, an optional format (Rformat) can be given.
- Function to read a line from a file define with its unit: **Read_line**
- Function to check is a string is empty: **string_IS_empty**

### Frac

This module contains functions and subroutines to manipulate fractions. 
All functions **elemental** (except the ones with a string). Therefore, one can use those functions with table.

- A type: **Frac_t**, which contains the numerator and denominator
  When the fraction is negative, the numerator has the minus sign.
  When the denominator is **0**, the fraction is simplified as 1/0.
- The operators (+, -, *, /, **, ==, >, <, <=, >=, /=) and the affectation (=) are orverloaded.
- When Frac_t is used for initialization and after numerical operations, the fraction is always reduced  (with the Greatest common divisor).
- One can test if a fraction is an integer with: frac_IS_integer(Frac)
- A fraction can be exported to:
  - to string: **TO_string(frac)**
  - to real (with the default kind=Rkind): **TO_real(frac)**
- The affectation can be done as follows:
  - Frac2 = Frac_t(3,-6) ! stored as -1/2
  - Frac2 = Frac1
  - Frac2 = Int1
  - Frac2 = '1/4' ! using a string of characters
  - where Frac1, Frac2 are fraction, Int1 is an integer
- Examples:
```Fortran
  TYPE(Frac_t) :: Frac1, Frac2
  TYPE(Frac_t), allocatable :: tab_Frac(:)

  Frac1 = '1/-2' ! use the conversion from string to Frac_t
  write(*,*) 'Frac1: ',TO_String(Frac1) ! it give "Frac1: -1/2"
  Frac2 = -2*Frac1 ! here the result is one and it is simplified
  write(*,*) 'Frac2: ',TO_String(Frac2) ! it give "Frac2: 1"
  Frac2 = Frac1**3
  write(*,*) 'Frac2: ',TO_String(Frac2) ! it give "Frac2: -1/8"
  tab_Frac = Frac_t(1,[2,3,4])
  write(*,*) 'tab_Frac: ',(TO_String(tab_Frac(i)) // ' ',i=1,size(tab_Frac)) ! it give "tab_Frac: 1/2 1/3 1/4 "
```

### RW_MatVec

This module contains public subroutines to read and write real and complex (kind=Rkind) vectors and matrices:

- Subroutines to read a matrix from a file: **Read_Mat** or **Read_Vec**
  All real kinds, all complexe kinds and one integer kind from NumParameters module are possible
  **Read_Mat(Mat,nio,nbcol,err)** or **Read_Vec(Vec,nio,nbcol,err)**
  - Mat or Vec: the matrice or the Vector to be readed. Mat or Vec have to be allocated
  - nio: the file unit
  - nbcol: the number of columns to be readed per block. If the number of columns of Mat or Vec is larger than nbcol, then several blocks are readed
  - err: integer to catch error (err=0 => no error)
- Subroutines to write a matrix to a file or string: **Write_Mat** or **Write_Vec**
  All real and complexe kinds from NumParameters module are possible. For complex number, it is not possible to write into a string.
  **Write_Mat(Mat,nio,nbcol,Rformat,info,iprint)** or **Write_Mat(Mat,string,nbcol,Rformat,info,iprint)** or **Write_Vec(Vec,nio,nbcol,Rformat,info,iprint)** or **Write_Vec(Vec,string,nbcol,Rformat,info,iprint)**
  - Mat or Vec: the matrice or the Vector to be printed. Mat or Vec have to be allocated
  - nio: the file unit
  - string: a character string with len long enough. It is better an allocatable character length [character (len=:), allocatable :: string]
  - nbcol: the number of columns to be printed per block. If the number of columns of Mat or Vec is larger than nbcol, then several blocks are printed
  - Rformat: (optional) format for real and complex. If Rformat is not present, then the default formats (RMatIO_format for real and CMatIO_format for complex) are used.
  - iprint: (optional) the Mat or Vec are printed if iprint is not present or if iprint=0
- The module contains the two default formats which are public:
  - RMatIO_format: (f18.10)
  - CMatIO_format: '(',f15.7,',',f15.7,')'

### Matrix

This module contains public functions and subroutines to perform some lienar algebra operations with the default real kind (kind=Rkind). It can use some LAPACK subroutines:

- Function to compute the inverse of a matrix:     **inv_OF_Mat_TO**
- Function to solve a linear system of equation:   **LinearSys_Solve**
- Function to compute the determinant of a matrix: **Det_OF**
- Function to initialyze a real identity matrix:   **Identity_Mat**

### Diago

This module contains public subroutines to digonlize a matrix with the default real kind (kind=Rkind). It can use some LAPACK subroutines:

- Subroutine to diagnalize a real matrix: **diagonalization**
