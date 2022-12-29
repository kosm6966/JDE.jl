module TIPSBP1

# Load Package Dependencies
import CSV
import DataFrames
using Plots
using DelimitedFiles
using LinearAlgebra
using Parameters
import Dates
using Optim

include("parameters.jl")
include("makefigures.jl")
include("spinmats.jl")
include("functions.jl")
include("hamiltonians.jl")
include("calculatespectrum.jl")
include("makespectrum.jl")
include("makepopdist.jl")
include("makeweber.jl")
include("simulatedannealing.jl")


export BridgedDimer
export Experiment
export TestData
export Spectrum
export Polarization
export Populations
export quintetspectrum
export initializetest
export polarplot
export weberplot
export makerandp0
export run_optim
export figure2
export figure3a
export figure3de


end