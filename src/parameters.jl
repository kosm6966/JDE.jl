"""

    Polarization()

    Polarization is a structure that is sent to the `polarplot()` function to calculate a 
plot of the polarization of the quintet state (`plottype = :Polarization`), as a function of
the dimer's orientation with respect to the field.

# Examples
```julia-repl
julia> pp = Polarization();
```
"""
@with_kw struct Polarization a::Bool = true end

"""

    Populations(M::Int64)

    Populations is a structure that is sent to the `polarplot()` function to calculate a 
plot of the populations for the specified M-sublevel (`plottype = :Populations`), as a function of
the dimer's orientation with respect to the field.

# Examples
```julia-repl
julia> pp = Populations(M = -2)
Populations
  M: Int64 -2
```
"""
@with_kw struct Populations
    M::Int64 = 0
end

struct Gaussian a::Bool end
struct Lorentzian a::Bool end
struct Diabatic a::Bool end
struct Adiabatic a::Bool end
"""

    Spectrum(nOutput::Int64,fwhm::Float64,line::Symbol,basis::Symbol)

Spectrum is a structure that is sent to the `quintetspectrum()` function to calculate  
the quintet spectrum in the diabatic basis. 

`nOutput::Int64` specifies whether the output of the function `S` is the total quintet spectrum (`1`); 
is separated into two spectra (`2`), `S(|M|) = [1↔2, 0↔1]`; or is seperated into four 
spectra (`4`), `S(M) = [-2↔-1, -1↔0, 0↔+1, +1↔+2]`.

`multiOut::Bool` specifies whether the output of the function `S` is integrated over the specified
orientations (`false`, default) or has a third dimension in which each S[:,:,j] corresponds to the 
spectrum for the given j-th orientation (`true`).

`fwhm::Float64` gives the linewidth fo the line function.

`line::Symbol` gives the line function, either `:Gaussian` or `:Lorentzian`.

`basis::Symbol` gives the chosen basis, either `:Adiabatic` or `:Diabatic`.

# Examples
```julia-repl
julia> spec = Spectrum()
Spectrum
  nOutput: Int64 1
  multiOut: Bool false
  fwhm: Float64 0.6
  line: Gaussian
  basis: Diabatic
```

```julia-repl
julia> spec = Spectrum(fwhm=0.8,line=:Lorentzian)
Spectrum
  nOutput: Int64 1
  multiOut: Bool false
  fwhm: Float64 0.8
  line: Lorentzian
  basis: Diabatic
```
"""
@with_kw struct Spectrum
    nOutput::Int64 = 1
    multiOut::Bool = false
    fwhm::Float64 = 0.6
    line = :Gaussian
    basis = :Diabatic
    function Spectrum(nOutput,multiOut,fwhm,line,basis)
        if nOutput != 1 && nOutput != 2 && nOutput != 4
                    error("nOutput must be 1, 2 or 4")
        end
        if line == :Gaussian && basis == :Diabatic
            new(nOutput,multiOut,fwhm,Gaussian(true),Diabatic(true))
        elseif line == :Lorentzian && basis == :Diabatic
            new(nOutput,multiOut,fwhm,Lorentzian(true),Diabatic(true))
        elseif line == :Gaussian && basis == :Adiabatic
            new(nOutput,multiOut,fwhm,Gaussian(true),Adiabatic(true))
        elseif line == :Lorentzian && basis == :Adiabatic
            new(nOutput,multiOut,fwhm,Lorentzian(true),Adiabatic(true))
        else error("Structure not called properly, check line and basis!")
        end
    end
end

"""

    Experiment(distribution::String,θ::{Float64,Vector{Float64}},w::{Float64,Vector{Float64}},
    field::Vector{Float64},ν::Float64`)

Experiment is a structure that is sent to the `quintetspectrum()` function to calculate  
the quintet powder spectrum. The default call to `quintetspectrum()` uses the values
nOctants and nPoints to calculate nPoints random rotations, projected on a nOctants of 
a sphere. The default likewise uses the values b0_left, b0_right, and nField to 
generate a vector of field points over which the spectrum is calculated.

`distribution::String` gives the path and filename to the orientation distribution of a 
powder spectrum. The file is comma separated and has the header `theta,phi,weight`. The
angles are given in degree.

`θ::{Float64,Vector{Float64}}` gives the values of theta in degrees. Alternate to reading from file.

`ϕ::{Float64,Vector{Float64}}` gives the values of phi in degrees. Alternate to reading from file.

`w::{Float64,Vector{Float64}}` gives the weights. Alternate to reading from file.

`field::Vector{Float64}` gives the field values.

`ν::Float64` gives the microwave excitation frequency in GHz.

# Examples
To calculate a powder spectrum where the orientations are specified in a file and the user
passes the field points as a vector called `b0vec`:
```julia-repl
julia> experiment = Experiment(distribution="/path/to/filename.txt",field=b0vec);
```
To calculate the spectrum for B||z, B||x, and B||y:
```julia-repl
julia> experiment = Experiment(θ=[0.,90.,90.],ϕ=[0.,0.,90.],field=b0vec)
```
"""
@with_kw struct Experiment 
    # Experiment parameters
    distribution = true
    θ = true
    ϕ = true
    w = true
    field = true

    ν::Float64 = 9.7310     # excitation frequency, GHz

    function Experiment(distribution,θ,ϕ,w,field,ν)
        if typeof(distribution) == String
            powder = CSV.read(distribution, DataFrames.DataFrame)
            θ = powder.theta
            ϕ = powder.phi
            w = (powder.weight ./maximum(powder.weight))
        elseif length(w) != length(θ)
            w = ones(Float64,length(θ))
        end
        ν = ν*1000
        new(distribution,θ,ϕ,w,field,ν)
    end
end

"""

    BridgedDimer{Float64}(d,e,β,x,Θ,Φ,g,j)

BridgedDimer is a structure that holds all spectroscopic parameters for the bridged dimer
system that is specified by the single parameter β, like TIPS-BP1'. The default values
are the best fit parameters for TIPS-BP1'.

`d::Float64` gives the axial zfs interaction.

`e::Float64` gives the rhombic zfs interaction.

`β::Float64` gives the angle `β` in degree as defined in Fig. 3c.

`x::Float64` gives the inter-chromophore anisotropic interaction.

`Θ::Float64` gives the angle `Θ` (capital theta) that defines the inter-chromophore 
anisotropic interaction in degree.

`Φ::Float64` gives the angle `Φ` (capital phi) that defines the inter-chromophore 
anisotropic interaction in degree.

# Examples
```julia-repl
julia> dimer = BridgedDimer()
BridgedDimer
  d: Float64 1322.0
  e: Float64 0.0
  β: Float64 111.0
  x: Float64 59.0
  Θ: Float64 90.0
  Φ: Float64 180.0
  g: Float64 2.003067
  j: Float64 -20000.0
```
"""
@with_kw struct BridgedDimer @deftype Float64
    # Optimized parameters
    d = 1322.         # axial ZFS, MHz
    e = 0.            # rhombic ZFS, MHz
    β = 111.
    
    x = 59.           # anisotropic exchange, MHz
    Θ = 90.	          # for X
	Φ = 180.		  # for X

    g = 2.003067      # g-tensor

    j = -20000.       # isotropic exchange, MHz
end