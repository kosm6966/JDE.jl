"""

    TestData()

TestData is a structure that opens data sets to reproduce the TIPS-BP1' EPR spectrum.

# Examples
```julia-repl
julia> test = TestData()
TestData
    prompt: Vector{Float64} 0.0
    hankel: Vector{Float64} 0.0
    residual: Vector{Float64} 0.0
    x: Vector{Float64} 0.0
```
"""
@with_kw struct TestData @deftype Vector{Float64}
    prompt = zeros(Float64,1)
    hankel = zeros(Float64,1)
    residual = zeros(Float64,1)
    x = zeros(Float64,1)
    aspect_ratio::Float64 = 1.
    function TestData(prompt,hankel,residual,x,aspect_ratio)
        datapath = "./test/testdata/Peak Normed - Residual.txt"
        A = DelimitedFiles.readdlm(datapath, '\t')
        x = A[:,1]
        prompt = A[:,2]
        residual = A[:,4]
        hankel = A[:,3]

        h=16.251
        w=48.838
        aspect_ratio = h/w/(2/abs(x[end]-x[1]))

        new(prompt,hankel,residual,x,aspect_ratio)
    end
end 

function initializetest(D=1322.,X=59.,beta=111.)
    dimer = BridgedDimer(d=D,x=X,β=beta)
    test = TestData()
    experiment = Experiment(distribution="./test/testdata/powder_distribution.csv",field=test.x)
    return test, experiment, dimer
end

function saveascsv(A,filename)
    datapath = joinpath(@__DIR__,"..", filename)
    open(datapath,"w") do f
        DelimitedFiles.writedlm(f,A,",")
end
    return nothing
end

"""

    figure2([D,X,β]; A)

figure2 is a function that reproduces the TIPS-BP1' EPR spectra from Figure 2. The 
default values for the spectroscopic parameters are the best fit values for TIPS-BP1'.
The paramters D, X, and β are optional parameters that default to the TIPS-BP1' values
if not specified. The function outputs a csv text file called "Year-Month-Day_Figure2.txt"
where the first column is the field values (mT), the second is the calculated spectrum 
in 2a, the third 2b, and the fourth 2c. A is an optional inputs for amplitude.

# Examples
```julia-repl
julia> figure2()

julia> figure2(1280.,0.,0.)

```
"""
function figure2(D=1322.,X=59.,beta=111.; A=0.8)
    test,experiment,dimer = initializetest(D,X,beta)
    S1 = quintetspectrum(dimer,Spectrum(basis=:Adiabatic),experiment);
    S2 = quintetspectrum(dimer,Spectrum(nOutput=2),experiment);
    data = [test.x S1*A S2*A]
    filename="Figure2"
    mkpath("data/")
    saveascsv(data,string("data/",Dates.today(),"_",filename,".txt"))
end

"""

    figure3a([D,X,β])

figure3a is a function that reproduces the TIPS-BP1' data from Figure 3a. The 
default values for the spectroscopic parameters are the best fit values for 
TIPS-BP1'. The paramters D, X, and β are optional parameters that default to 
the TIPS-BP1' values if not specified. The function outputs five csv files, 
detailed below.

# Output
Year-Month-Day_Figure3a-Hankel.txt: column one is the field in mT, columns 2-4 
    are the M = 0 ↔ ± 1 transitions for B||z, B||x, B||y, respectively.
Year-Month-Day_Figure3a-Residual.txt: column one is the field in mT, columns 2-4 
    are the M = ± 1 ↔ ± 2 transitions for B||z, B||x, B||y, respectively.
Year-Month-Day_Figure3a-Eigenvalues-Bi.txt: column one is the field in mT, columns
    2-6 are the eigenvalues of the quintet M-levels for B||i, where i is x, y, or z.
Year-Month-Day_Figure3a-Populations-Bi.txt: column one is the field in mT, column 
    two is the energy in GHz, column three is the population for the four quintet 
    transitions when the field is along the dimer axes i = x, y, or z. Arrows in 
    the figure are drawn in the figure from the states with low to high population 
    at a given resonance.

# Examples
```julia-repl
julia> figure3a()

julia> figure3a(1280.,0.,0.)

```
"""
function figure3a(D=1322.,X=59.,beta=111.)
    test,experiment,dimer = initializetest(D,X,beta)
    experiment = Experiment(θ=[0.,90.,90.],ϕ=[0.,0.,90.],field=test.x);

    S = quintetspectrum(dimer,Spectrum(nOutput=2,multiOut=true),experiment);
    data = [test.x S[:,:,1]]
    filename="Figure3a-Residual"
    mkpath("data/")
    saveascsv(data,string("data/",Dates.today(),"_",filename,".txt"))
    data = [test.x S[:,:,2]]
    filename="Figure3a-Hankel"
    saveascsv(data,string("data/",Dates.today(),"_",filename,".txt"))

	En, Populations, Signal = calculateweber(dimer,experiment);
    En ./= 1000
    Signal[2] ./= 1000

    data = [test.x En[:,:,1]']
    filename="Figure3a-Eigenvalues-Bz"
    saveascsv(data,string("data/",Dates.today(),"_",filename,".txt"))
    data = [test.x En[:,:,2]']
    filename="Figure3a-Eigenvalues-Bx"
    saveascsv(data,string("data/",Dates.today(),"_",filename,".txt"))
    data = [test.x En[:,:,3]']
    filename="Figure3a-Eigenvalues-By"
    saveascsv(data,string("data/",Dates.today(),"_",filename,".txt"))

    En, Populations, Signal = calculateweber(dimer,experiment);
    En ./= 1000
    data = [reshape(Signal[1][:, :, 1],:,1) reshape(Signal[2][:, :, 1],:,1) reshape(Populations[:, :, 1],:,1)]
    filename="Figure3a-Populations-Bz"
    saveascsv(data,string("data/",Dates.today(),"_",filename,".txt"))
    data = [reshape(Signal[1][:, :, 2],:,1) reshape(Signal[2][:, :, 2],:,1) reshape(Populations[:, :, 2],:,1)]
    filename="Figure3a-Populations-Bx"
    saveascsv(data,string("data/",Dates.today(),"_",filename,".txt"))
    data = [reshape(Signal[1][:, :, 3],:,1) reshape(Signal[2][:, :, 3],:,1) reshape(Populations[:, :, 3],:,1)]
    filename="Figure3a-Populations-By"
    saveascsv(data,string("data/",Dates.today(),"_",filename,".txt"))
end

"""

    figure3de([D,X,β])

figure3de is a function that reproduces the TIPS-BP1' data from Figures 3d and 3e. 
The default values for the spectroscopic parameters are the best fit values for 
TIPS-BP1'. The paramters D, X, and β are optional parameters that default to 
the TIPS-BP1' values if not specified. The function outputs six csv text files, 
detailed below.

# Output
Year-Month-Day_Figure3d.txt: The first row of this matrix gives values of theta
    in degrees. The first column gives values of phi in degrees. The inner matrix
    values are the entropy for the corresponding (phi,theta).
Year-Month-Day_Figure3e_i.txt: The first row of this matrix gives values of theta
    in degrees. The first column gives values of phi in degrees. The inner matrix
    values are the populations for the corresponding (phi,theta) and the M-level,
    i = -2, -1, 0, 1, 2.

# Examples
```julia-repl
julia> figure3de()

julia> figure3de(1280.,0.,0.)

```
"""
function figure3de(D=1322.,X=59.,beta=111.)
    _,_,dimer = initializetest(D,X,beta)
    polarplot(dimer,Polarization())
    ϕ=[0.:0.5:360.;]
    θ=[LinRange( 0., 90., Int(floor(0.5*(length(ϕ)+1))));]
    Populations = JDE.calculatepopdist(dimer,θ,ϕ)
    fn(M) = string("Figure3e_",M-3)
    mkpath("data/")
    for ii = 1:5
        data = [[0.;θ]'; [ϕ Populations[:,:,ii]]]
        filename=fn(ii)
        saveascsv(data,string("data/",Dates.today(),"_",filename,".txt"))
    end

    Entropy = calculateentropy(θ,ϕ,Populations)
    data = [[0.;θ]'; [ϕ Entropy]]
    filename="Figure3d"
    saveascsv(data,string("data/",Dates.today(),"_",filename,".txt"))
end
