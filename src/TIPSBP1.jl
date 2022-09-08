module TIPSBP1

# Load Package Dependencies
import CSV
import DataFrames
using Plots
import DelimitedFiles
using LinearAlgebra
using Parameters
import Dates

include("parameters.jl")
include("preparetest.jl")
include("spinmats.jl")
include("functions.jl")
include("hamiltonians.jl")
include("calculatespectrum.jl")
include("makespectrum.jl")
include("makepopdist.jl")
include("makeweber.jl")


export BridgedDimer
export Experiment
export Spectrum
export Polarization
export Populations
export quintetspectrum
export initializetest
export polarplot
export weberplot


end