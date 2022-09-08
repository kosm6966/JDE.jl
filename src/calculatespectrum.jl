"""
    quintethamiltonian(dimer,exp,basis)
	
Compute the hamitlonian for the quintet pair and return the resonances and intensities of four transitions. 

"""
function quintethamiltonian(dimer,exp,basis)
    b0 = exp.field
    θ = exp.θ
    ϕ = exp.ϕ
    
    ## Compute HJ
	HJ = zeros(Float64,9,9)
	SA,SB,S = spin9()
    isoexchange!(HJ,S,dimer.j)
    S[3] .*= dimer.g*9.274009994e-24/6.62607015e-34*1e-9
 
	## Compute ZFS tensors
	R = zeros(Float64,3,3)
    Tmp3 = zeros(Float64,3,3)
	DA = zeros(Float64,3,3)
	DB = zeros(Float64,3,3)
    zfstensor!(R,Tmp3,DA,DB,dimer)

	## Compute X tensor
	X = zeros(Float64,3,3)
    anisotensor!(X,dimer)

    ## Preallocate arrays
	H0 = zeros(Complex{Float64},9,9)
	H = zeros(Complex{Float64},9,9)
    Tmp9 = zeros(Complex{Float64},9,9)
    index = Array{CartesianIndex{2},2}(undef,4,1)
    En = zeros(Float64,5,length(b0))
	dE = zeros(Float64,4,length(b0))
    Signal = zeros(Float64,4,length(θ))
    Intensity = zeros(Float64,4,length(θ))

	for i in eachindex(θ)
        ## Compute the anisotropic hamiltonian for a given orientation
        aniso!(H0,R,Tmp3,Tmp9,DA,DB,X,SA,SB,θ[i],ϕ[i])
        H0 .+= HJ

        ## Add the zeeman field to the full hamiltonian 
        enfield(En,H,H0,S,b0,basis)

        ## Compute the resonances and intensities
        findRes!(index,dE,En,Signal,b0,exp,i)
        findInt!(Tmp9,Intensity,H,H0,S,b0,index,i,basis)
    end
	
	return real(Signal), real(Intensity)
end

function enfield(En,H,H0,S,b0,basis::Diabatic)
    for j in eachindex(b0)
        H .= H0 .+ b0[j].*S[3]
        En[:,j] .= real(diag(H[5:9,5:9]))
    end
end

function enfield(En,H,H0,S,b0,basis::Adiabatic)
    for j in eachindex(b0)
        H .= H0 .+ b0[j].*S[3]
        H .= (H .+ H') .* 0.5
        En[:,j] .= real(eigvals(H[5:9,5:9]))
    end
end

function findRes!(index,dE,En,Signal,b0,exp,i)
    for k in eachindex(b0), j in 1:4
            #[Q(+1)-Q(+2)...
            # Q( 0)-Q(+1)...
            # Q(-1)-Q( 0)...
            # Q(-2)-Q(-1)]
        dE[j,k] = abs.(En[5-j,k] - En[6-j,k])
    end

    index .= argmin( abs.( dE .- exp.ν), dims=2)

    for j in 1:4
        Signal[j,i] = b0[index[j][2]]
    end
end

function findInt!(Tmp9,Intensity,H,H0,S,b0,index,i,basis::Diabatic)
    for j in 1:4
        H .= H0 .+ b0[index[j][2]].*S[3]
        #[Q(+1)-Q(+2)... 8-9
        # Q( 0)-Q(+1)... 7-8
        # Q(-1)-Q( 0)... 6-7
        # Q(-2)-Q(-1)]   5-6
		dP = real(abs.( H[1,9-j] ).^2 .- abs.( H[1,10-j] ).^2)
        Mu2 = real(abs( S[1][9-j,10-j] )^2)
        Intensity[j,i] = dP*Mu2
    end
    return nothing
end

function findInt!(Tmp9,Intensity,H,H0,S,b0,index,i,basis::Adiabatic)
    Tmp9 .*= 0.
    Tmp9[1,1] = 1.
    for j in 1:4
        H .= H0 .+ b0[index[j][2]].*S[3]
        H .= (H .+ H') .* 0.5
        Tmp9[5:9,5:9] = eigvecs(H[5:9,5:9])
        Tmp9[2:4,2:4] = eigvecs(H[2:4,2:4])
        conjtransUAU!(Tmp9,H)
        #[Q(+1)-Q(+2)... 8-9
        # Q( 0)-Q(+1)... 7-8
        # Q(-1)-Q( 0)... 6-7
        # Q(-2)-Q(-1)]   5-6
		dP = real(abs.( H[1,9-j] ).^2 .- abs.( H[1,10-j] ).^2)
        H .= S[1]
        conjtransUAU!(Tmp9,H)
        Mu2 = real(abs( H[9-j,10-j] )^2)
        Intensity[j,i] = dP*Mu2
    end
    return nothing
end