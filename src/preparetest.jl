"""

    TestData()

TestData is a structure that opens data sets to reproduce the TIPS-BP1' EPR spectrum.

# Examples
```julia-repl
julia> test = TestData()
TestData
    promp: Vector{Float64} 0.0
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
    function TestData(prompt,hankel,residual,x)
        datapath = "./test/testdata/Peak Normed - Residual.txt"
        A = DelimitedFiles.readdlm(datapath, '\t')
        x = A[:,1]
        prompt = A[:,2]
        prompt ./= maximum(abs.(prompt))
        residual = A[:,4]
        residual ./= maximum(abs.(residual))
        hankel = A[:,3] + A[:,4]
        hankel ./= maximum(abs.(hankel))
        new(prompt,hankel,residual,x)
    end
end 

function initializetest()
    dimer = BridgedDimer()
    test = TestData()
    experiment = Experiment(distribution="./test/testdata/powder_distribution.csv",field=test.x);
    return test, experiment, dimer
end

function saveascsv(A,filename)
    datapath = joinpath(@__DIR__,"..", filename)
    open(datapath,"w") do f
        DelimitedFiles.writedlm(f,A,",")
end
    return nothing
end

function figure2()
    test,experiment,dimer = initializetest()
    S1 = quintetspectrum(dimer,Spectrum(basis=:Adiabatic),experiment);
    S2 = quintetspectrum(dimer,Spectrum(nOutput=2),experiment);
    data = [test.x S1*0.8 S2*0.8]
    filename="Figure2"
    saveascsv(data,string("Data/",Dates.today(),"_",filename,".txt"))
end

function figure3a()
    test,experiment,dimer = initializetest()
    experiment = Experiment(θ=[0.,90.,90.],ϕ=[0.,0.,90.],field=test.x);

    S = quintetspectrum(dimer,Spectrum(nOutput=2,multiOut=true),experiment);
    data = [test.x S[:,:,1]]
    filename="Figure2a-Residual"
    saveascsv(data,string("Data/",Dates.today(),"_",filename,".txt"))
    data = [test.x S[:,:,2]]
    filename="Figure2a-Hankel"
    saveascsv(data,string("Data/",Dates.today(),"_",filename,".txt"))

	En, Populations, Signal = calculateweber(dimer,experiment);
    En ./= 1000
    Signal[2] ./= 1000

    data = [test.x'; En[:,:,1]]
    filename="Figure2a-Eigenvalues-Bz"
    saveascsv(data,string("Data/",Dates.today(),"_",filename,".txt"))
    data = [test.x'; En[:,:,2]]
    filename="Figure2a-Eigenvalues-Bx"
    saveascsv(data,string("Data/",Dates.today(),"_",filename,".txt"))
    data = [test.x'; En[:,:,3]]
    filename="Figure2a-Eigenvalues-By"
    saveascsv(data,string("Data/",Dates.today(),"_",filename,".txt"))

    En, Populations, Signal = calculateweber(dimer,experiment);
    En ./= 1000
    data = [reshape(Signal[1][:, :, 1],:,1) reshape(Signal[2][:, :, 1],:,1) reshape(Populations[:, :, 1],:,1)]
    filename="Figure2a-Populations-Bz"
    saveascsv(data,string("Data/",Dates.today(),"_",filename,".txt"))
    data = [reshape(Signal[1][:, :, 2],:,1) reshape(Signal[2][:, :, 2],:,1) reshape(Populations[:, :, 2],:,1)]
    filename="Figure2a-Populations-Bx"
    saveascsv(data,string("Data/",Dates.today(),"_",filename,".txt"))
    data = [reshape(Signal[1][:, :, 3],:,1) reshape(Signal[2][:, :, 3],:,1) reshape(Populations[:, :, 3],:,1)]
    filename="Figure2a-Populations-By"
    saveascsv(data,string("Data/",Dates.today(),"_",filename,".txt"))
end

function figure3de()
    _,_,dimer = initializetest()
    polarplot(dimer,Polarization())
    ϕ=[0.:0.5:360.;]
    θ=[LinRange( 0., 90., Int(floor(0.5*(length(ϕ)+1))));]
    Populations = TIPSBP1.calculatepopdist(dimer,θ,ϕ)
    fn(M) = string("Figure3e_",M-3)
    for ii = 1:5
        data = [[0.;θ]'; [ϕ Populations[:,:,ii]]]
        filename=fn(ii)
        saveascsv(data,string("Data/",Dates.today(),"_",filename,".txt"))
    end

    Entropy = calculateentropy(θ,ϕ,Populations)
    data = [[0.;θ]'; [ϕ Entropy]]
    filename="Figure3d"
    saveascsv(data,string("Data/",Dates.today(),"_",filename,".txt"))
end
