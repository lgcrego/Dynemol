# Dynemol
DynEMol: tools for studying Dynamics of Electrons in Molecules

🧪 DynEMol: Dynamics of Electrons in Molecules
DynEMol is an open-source suite for simulating nonadiabatic excited-state dynamics and molecular systems using a hybrid semiempirical approach. It has been in continuous development since 2013 and is actively maintained. If you have questions or encounter issues, feel free to open an issue on this GitHub repository or contact the maintainer directly at luis.gc.rego [at] gmail.com.

🔧 What DynEMol Does

DynEMol combines two computational strategies:

Extended Hückel theory for the electronic structure

Classical force fields (AMBER/OPLS-compatible) for nuclear motion

This allows it to model a wide range of systems with nonadiabatic dynamics, including:

Molecules in vacuum or solution

Molecule–solid and solid–liquid interfaces

Bulk materials and hybrid perovskites

DynEMol supports:

Excited-state nonadiabatic dynamics (Ehrenfest and fewest-switches)

Classical molecular dynamics

Geometry optimization

Force-field and parameter optimization

Time-dependent electron and hole population analysis

Electronic and vibronic density of states (DOS, VDOS)

Pre- and post-processing tools for simulations and analysis

💻 Code & Compilation
DynEMol is primarily written in Modern Fortran (2003/2008/2018) and is designed to work best with Intel compilers:

Use make-ifx for Intel oneAPI (LLVM-based compilers)

Use make-ifort for legacy Intel compilers (pre-2022)

Other compilation targets (e.g., for debugging or profiling) are documented in the makefile.

🖥️ Parallelization & Branches
The repository is divided into two main branches:

master — optimized for HPC clusters with MPI + OpenMP support and some GPU acceleration

SingleNode — optimized for single-node systems, parallelized with OpenMP and supports the same GPU routines

Both branches include a manipulate/ folder containing standalone utilities for:

Editing input structures (copy, rotate, translate, delete atoms)

Calculating autocorrelation functions, spatial distributions

Generating plots from simulation output

Subdirectories such as general.util/ and short.time-VDOS.util/ include additional utilities for graphing and post-processing.

📁 Input & File Formats
DynEMol supports a range of input formats:

Coordinates: .pdb

Force fields: AMBER/OPLS-compatible .psf, .prm, .itp, .top

Configuration files: card.inpt, velocity.inpt, opt_eht_parms.inpt

Input file specifications can be found in:

card_file_formats

IO_file_formats

IO_file_MM-formats

cost_file_formats

See the tutorial-Dynemol_code.pdf for detailed examples and instructions.

📘 Examples in Literature
DynEMol has been used in peer-reviewed studies across a variety of systems, including:

Photochemistry in solution

Charge-transfer at molecule–semiconductor interfaces

Hybrid organic–inorganic perovskites

Chiral nanowires and electronic transport

Excited-state proton-coupled dynamics

A list of selected publications will be included in this repository.

🚀 Getting Started
To run a basic simulation:

Set your environment:

export DYNEMOLDIR=/path/to/dynemol

export DYNEMOLWORKDIR=$(pwd)

Place your input files (card.inpt, input.pdb, etc.) in $DYNEMOLWORKDIR

Run:

$DYNEMOLDIR/dynemol

Outputs will be saved in folders such as dyn.trunk/, dos.trunk/, MO.trunk/, etc., and refreshed at each execution.
