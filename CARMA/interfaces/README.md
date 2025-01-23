# PARMA: Interfaces for running CARMA without a broader model
Developed by Parker Case (parker.a.case@nasa.gov)
Current version: v1.0 (2024/08/20)
## Introduction
While the testing suite in CARMA is useful for making sure code changes haven't
altered the basic processes, it is cumbersome to change the experiments and
look at the output.

Here we've developed PARMA, a python-wrapped interface for CARMA. We leverage
f2py, which compiles fortran in a way that allows fortran subroutines to be
called by Python as functions.

## List of interfaces and related python drivers
Below is a list of Fotran interfaces and their related python drivers:
- carma_column_sulfate.F90: the simplest column model with only pure sulfate
    - parma_column_sulfate.py: the python interface for running the sulfate
        column model. The main method of this file allows you to run this
        model with `python parma_column_sulfate.py`.
    - parma_column_sulfate_analysis.py: an example for reading the output of
        this model.
- carma_column_dust.F90: a column with a pure sulfate group and a mixed group
    with a dust element and a sulfate element.
    - parma_column_dust.py: the python interface for running the dust
        column model. The main method of this file allows you to run this
        model with `python parma_column_dust.py`.
    - parma_column_dust_analysis.py: an example for reading the output of
        this model.
- carma_box.F90: A box model interface for a single group, single element.
    - parma_box.py: Simple python script to drive carma_box.F90 based on
                     stratospheric aerosol.
    - parma_box.ipynb: The same is parma_box.py, in notebook form.

## Environment
PARMA has been tested using gfortran and python 3.9. If you are on NASA NCCS
resources, you can replicate my environment by loading the following modules:

`module load python/GEOSpyD/Min4.11.0_py3.9`
`module load comp/gcc/9.3.0`
`module load lib/mkl/20.0.0.166`

Compilers besides gfortran are probably useable, as are other versions of
python 3.x, but have not been tested.

## Compiling PARMA
You can compile PARMA using the normal CARMA make structure. For more details,
see `make-carma.csh` and the `Makefile` in the CARMA directory. To use this to
compile carma_box.F90 and prepare the associated python scripts, run the
following command:

`./make-carma.csh parma_box`

This will compile carma_box.F90 using f2py (as well as compiling all of CARMA).
Additionally, this will copy parma_box.py, parma_box.ipynb, and the quite of
python CARMA tools into the build directory.

Now you can navigate to the `build/carma/` directory. From there, you can
run the standalone PARMA script by running:

`python parma_box.py`

You can also use your favorite notebooks tools (Jupyter notebooks, for example)
to open parma_box.ipynb.
