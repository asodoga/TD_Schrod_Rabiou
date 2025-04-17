# Quantum Molecular Dynamics Simulator (Hagedorn Wavepacket Propagation)

**A Fortran-based simulation tool for quantum molecular dynamics using generalized Hagedorn wavepackets with multi-basis support**

## Table of Contents
- [Key Features](#key-features)
- [Theoretical Background](#theoretical-background)
- [Installation](#installation)
- [Compilation](#compilation)
- [Running Simulations](#running-simulations)
- [Examples](#examples)
- [License](#license)

## Key Features

### Advanced Wavepacket Propagation
- Generalized Hagedorn basis implementation (multi-dimensional)
- Support for alternative primitive bases (Fourier, etc.)
- Three-step propagation scheme (non-variational approach)

### Technical Specifications
- High-performance Fortran 95+ core

### Simulation Capabilities
- Flexible propagator options (Hagedorn/standard)
- Customizable basis parameters

## Theoretical Background

This code implements:
- Time-dependent Schrödinger equation solver
- Hagedorn's wavepacket formulation extended to:
  - Multi-dimensional systems
  - Generalized basis sets
- Specialized propagation algorithm (3-step process)

## Installation

### Dependencies
| Package        | Minimum Version | Notes                     |
|--------------  |-----------------|---------------------------|
| gfortran       | 9.5             | Fortran compiler          |
| gcc            | 9.0             | C compiler                |
| CMake          | 3.12            | Build system              |
| GNUplot        | 5.4             | Visualization             |
| QuantumModelLib| -               | Potential library         |
| QDUtilLib      | -               | Math/IO utilities         |

### Setup
```bash
# Create required directories
mkdir -p obj EXT_LIB

# Clone and build dependencies
cd EXT_LIB
git clone https://github.com/lauvergn/QuantumModelLib.git
git clone https://github.com/lauvergn/QDUtilLib.git
# Follow build instructions in each repository
```

### Compilation

Edit default_comp.sh:
Set library paths in Makefile
Configure compilation flags
For the default compilation in the bash thes folowing :
``` bash
./default.comp nb nb**2

```
The script performs:

Source compilation via Make
Output directory creation
Generation of TD_SCHROD.x executable

### Configuration File Structure (default.comp)
# ====== Potential Section ======
&potential
    pot_name  = 'model'
    option= number
    adiabatic = f/t
    nsurf = number
    ndim = number
   /

# ====== Basis Configuration ======
&basis_nd name='dp' nb_basis=3 /
  &basis_nd name='basis type' nb =$nb nq=$nq Q0  SCALEQ  Imp_k  alpha
  &basis_nd name= 'el' nb=1 /

# ====== Initial WP Settings ======
&defGWP ndim  Elecindex  Coef/
  &defWP0  sigma=   Beta=  Qeq=    imp_k gamma /

# ====== Propagation Settings ======
&prop t0=0.0  tf=60. delta_t=0.1 eps=1.e-20
      max_iter=500   propa_name='hagedorn'  propa_name2='taylor' Beta=t  P=t renorm=f/

PROPAGATOR name            # Hagedorn or Standard
TIME_STEP=delta_t           # in atomic unit (ua)
TOTAL_TIME=tf               # Total propagation time


MIT License
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!    Copyright 2025 Rabiou ISSA [1,*]
!      with contributions of
!        David Lauvergnat [2]
!        Komi SODOGA [1]
!        
!
![1]:   Laboratoire de Physique des Matériaux et des Composants à semi-conducteurs (LPMCS),Université de Lomé
![2]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![*]: issagoudouya@gmail.com  
