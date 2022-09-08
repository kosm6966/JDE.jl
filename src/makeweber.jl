function weberplot(dimer,spec,exp)
    S = quintetspectrum(dimer,spec,exp);
    l=["z" "x" "y"]
	En, Populations, Signal = calculateweber(dimer,exp);
    En ./= 1000
    Signal[2] ./= 1000
    
    i = 1
    fig1 = plot(exp.field,En[:,:,i]',linewidth=2, xaxis=("B₀ (mT)"),title = l[i],yaxis=("Energy (GHz)"),label=false,linecolor=:black)
    scatter!(Signal[1][:,1:2,i],Signal[2][:,1:2,i], xaxis=("B₀ (mT)"),yaxis=("Energy (GHz)"),markersize=20*(Populations[:,1:2,i]).^0.5,label=false)
    plot!(exp.field,S[:,i,1]*10 .- 45,label="⁵TT₁↔⁵TT₂",l=2)
    plot!(exp.field,S[:,i,2]*10 .- 45,label="⁵TT₀↔⁵TT₁",l=2)
    i = 2
    fig2 = plot(exp.field,En[:,:,i]',linewidth=2, xaxis=("B₀ (mT)"),title = l[i],yaxis=("Energy (GHz)"),label=false,linecolor=:black)
    scatter!(Signal[1][:,1:2,i],Signal[2][:,1:2,i], xaxis=("B₀ (mT)"),yaxis=("Energy (GHz)"),markersize=20*(Populations[:,1:2,i]).^0.5,label=false)
    plot!(exp.field,S[:,i,1]*10 .- 45,label="⁵TT₁↔⁵TT₂",l=2)
    plot!(exp.field,S[:,i,2]*10 .- 45,label="⁵TT₀↔⁵TT₁",l=2)
    i = 3
    fig3 = plot(exp.field,En[:,:,i]',linewidth=2, xaxis=("B₀ (mT)"),title = l[i],yaxis=("Energy (GHz)"),label=false,linecolor=:black)
    scatter!(Signal[1][:,1:2,i],Signal[2][:,1:2,i], xaxis=("B₀ (mT)"),yaxis=("Energy (GHz)"),markersize=20*(Populations[:,1:2,i]).^0.5,label=false)
    plot!(exp.field,S[:,i,1]*10 .- 45,label="⁵TT₁↔⁵TT₂",l=2)
    plot!(exp.field,S[:,i,2]*10 .- 45,label="⁵TT₀↔⁵TT₁",l=2)
    return fig1, fig2, fig3
end

function calculateweber(dimer,exp)
    basis = Diabatic(true)
    b0 = exp.field
    θ = exp.θ
    ϕ = exp.ϕ
    
    ## Compute HJ
	HJ = zeros(Float64,9,9)
	SA,SB,S = spin9()
    isoexchange!(HJ,S,dimer.j)
    S[3] .*= dimer.g*9.274009994e-24/6.62607015e-34*1e-9
 
	##ZFS##
	R = zeros(Float64,3,3)
    Tmp3 = zeros(Float64,3,3)
	DA = zeros(Float64,3,3)
	DB = zeros(Float64,3,3)
    zfstensor!(R,Tmp3,DA,DB,dimer)

	##X##
	X = zeros(Float64,3,3)
    anisotensor!(X,dimer)

	H0 = zeros(Complex{Float64},9,9)
    H = zeros(Complex{Float64},9,9)
    Tmp9 = zeros(Complex{Float64},9,9)
    
    ## Preallocate arrays
	H0 = zeros(Complex{Float64},9,9)
	H = zeros(Complex{Float64},9,9)
    Tmp9 = zeros(Complex{Float64},9,9)
    index = Array{CartesianIndex{2},2}(undef,4,1)
    En = zeros(Float64,5,length(b0))
	dE = zeros(Float64,4,length(b0))
    
    Diags = zeros(Float64,5,length(b0),length(θ))
    SSS = zeros(Float64,4,length(θ))
    bSignal = zeros(Float64,4,2,length(θ))
    enSignal = zeros(Float64,4,2,length(θ))
    Populations = zeros(Float64,4,2,length(θ))

	for i in eachindex(θ)
        aniso!(H0,R,Tmp3,Tmp9,DA,DB,X,SA,SB,θ[i],ϕ[i])
        H0 .+= HJ

        ## Add the zeeman field to the full hamiltonian 
        enfield(En,H,H0,S,b0,basis)
        Diags[:,:,i] = En

        ## Compute the resonances and intensities
        findRes!(index,dE,En,SSS,b0,exp,i)
        findPop!(Populations,bSignal,enSignal,H,H0,b0,S,index,i)
	end
	return Diags, Populations, (bSignal, enSignal)
end

function findPop!(Populations,Res,En,H,H0,b0,S,index,i)
    for j in 1:4
        H .= H0 .+ b0[index[j][2]].*S[3]
        #[Q(+1)-Q(+2)... 8-9
        # Q( 0)-Q(+1)... 7-8
        # Q(-1)-Q( 0)... 6-7
        # Q(-2)-Q(-1)]   5-6
        rho = real(abs.( H[1,5:9] ).^2)
        rho ./= sum(rho)
        Populations[j,1,i] = rho[5-j]
        Populations[j,2,i] = rho[6-j]

        En[j,1,i] = real( H[9-j,9-j] )
        En[j,2,i] = real(H[10-j,10-j] )

        Res[j,1,i] = b0[index[j][2]]
        Res[j,2,i] = b0[index[j][2]]
    end
    return nothing
end