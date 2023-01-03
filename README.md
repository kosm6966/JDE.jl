# JDE.jl
[![DOI](https://zenodo.org/badge/539056584.svg)](https://zenodo.org/badge/latestdoi/539056584)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

The **JDE.jl** software written in [Julia](https://julialang.org) numerically calculates the EPR spectrum from transitions between sublevels of the coupled triplet pair quintet state from parallel and bridged TIPS-BP1'-like dimers that undergo singlet fission. The code was written with an emphasis on clarity to be useful to scientists of various backgrounds and customizable for various dimer systems. The population calculations are a distinctive feature of this software. See (DOI: XXXXXX) for a complete explanation of our underlying theory. Features include:

* Compute the EPR spectrum for the triplet pair from singlet fission.
  * Calculate the full TT hamiltonian using the diabatic or adiabatic basis.
  * Calculate resonances for a given excitation energy.
  * Calculate intensities from the hamiltonian based on the JDE model.
  * Calculate single-orientation spectra.
  * Calculate the spectrum from a glassy sample whose orientation distribution is given in an input file or generate a random distribution for specified quadrants.
  * Vary spectroscopic parameters `g`, `J`, `D`, `E`, `β`, `X`, `Θ`, `Φ`. For parallel dimers, set the bridging angle, `β = 180`.

* Calculate populations for the triplet pair from singlet fission using the JDE model.

* Generate figures for single orientation spectra that identify where transitions occur along with quintet sublevel populations.

* Calculate population and entropy distributions for quintet sublevels with respect to a dimer's orientation in the field.


## Installation
If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). The following usage instructions assume that Julia has been added to your path. The JDE software works with Julia v1.7.3.

## Usage
Download the package and all contents. Then, open a terminal and run the package from the cloned directory:
```bash
git clone git@github.com:joeleaves/JDE.git
cd JDE
julia --project=@. -e 'import Pkg; Pkg.instantiate()' # Install JDE's dependencies
julia -e 'import Pkg; Pkg.add(["Plots", "PyPlot"])' # Install postprocessing tools
```
If you installed JDE this way, you always have to start Julia with the `--project`
flag set to your local JDE clone, e.g.,
```bash
julia --project=@.
```

Type "?" for help. A list of exported function names are given in the JDE software help page. Main functions also have help pages that give detailed descriptions and provide examples for calling the functions,
```julia
help?> JDE
```


 For examples of how to run some of the main functions, see the .jl files in the test folder. These can be run from the command line as,
```bash
julia --project=@. test/SAtest.jl
```

SAtest.jl provides an example of how to fit parameters using simulated annealing. TIPSBP1test-data.jl and TIPSBP1test-figs.jl generate the data and figures published in the paper (DOI: XXXXXX). To run an interactive script use the flag,
```bash
julia -i --project=@. test/TIPSBP1test-figs.jl
```
which loads a julia REPL that must be exited using,
```julia
julia> exit()
```

Most functions take optional inputs and are run using the julia REPL,
```julia
julia> figure3de()
julia> figure3de(1280.,0.,0.)
```

## Referencing
If you use the JDE software in your own research or write a paper using results obtained with the help of the JDE software, please cite the following article:
```bibtex
@article{dill2022,
  title={Entangled Spin-polarized Excitons from Singlet Fission in a Rigid Dimer},
  author={Dill, Ryan D. and Smyser, Kori E. and Rugg, Brandon K. and Damrauer, Niels H. and Eaves, Joel D.},
  journal={Nature Communications},
  volume={X},
  number={X},
  pages={XX},
  year={2022},
  doi={10.XXXX}
}
```

In addition, you can refer to the JDE software directly as
```bibtex
@misc{smyser2022,
  title={JDE.jl: A nonadiabatic transition theory for the triplet pair from singlet fission},
  author={Smyser, Kori E. and Eaves, Joel D.},
  year={2022},
  month={10},
  howpublished={\url{https://github.com/joeleaves/JDE}},
  doi={10.5281/zenodo.7464138}
}
```

## Authors
The JDE software was written by Kori E. Smyser and is based on the JDE model, designed by Kori E. Smyser and [Joel D. Eaves](https://www.colorado.edu/lab/eavesgroup), who maintains this repository.

## Acknowledgments
Funding was provided by the United States Department of Energy, Office of Basic Energy Sciences (ERW7404). This work also utilized resources from the University of Colorado Boulder Research Computing Group, which is supported by the National Science Foundation (awards ACI-1532235 and ACI-1532236), the University of Colorado Boulder, and Colorado State University.

## License and contributing
The JDE software is licensed under the Apache-2.0 license (see [LICENSE.md](LICENSE.md)).

*Disclaimer:* Using the software and the scripts in this repository is at your own risk. The software is freeware and it comes without warranty. Therefore, the Authors are not responsible for the loss of data, time, bad results, or anything else deriving by the use of the software, data or any contents of this repository. We do not make any warranty, express or implied, or assume any legal liability or responsibility for the accuracy, completeness, or usefulness of the software, scipts, data, or results, or represent that its use would not infringe privately owned rights.

*Note on performance:* Julia uses just-in-time compilation to transform its source code to native, optimized machine code at the *time of execution* and caches the compiled methods for further use. That means that the first execution of a Julia method is typically slow, with subsequent runs being much faster. 

