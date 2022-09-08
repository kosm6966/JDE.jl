"""
    isoexchange!(HJ,S,j)
	
Compute the inter-chromophore isotropic interaction hamiltonian `HJ` for the coupled triplet pair.

"""
function isoexchange!(HJ,S,j)
	HJ .= zero(HJ)
	for i = 1:3
			HJ .+= S[i]*S[i] # S² = Sx² + Sy² + Sz²
	end
	HJ .-=  diagm(repeat([4.],9)) # S² - 4
	HJ .*= 0.5*j # J/2(S² - 4)

    return nothing
end

"""
    zfstensor!(R,Tmp3,DA,DB,dimer::BrigedDimer)

Compute	the intra-chromophore anisotropic interaction tensors DA and DB for bridged chromophores. β is the angle between them. See TIPS-BP1' paper.

"""
function zfstensor!(R,Tmp3,DA,DB,dimer::BridgedDimer)
	Tmp3 .= zero(Tmp3)
    onethird = 0.3333333333333333
    Tmp3[1,1] = -dimer.d*onethird + dimer.e
	Tmp3[2,2] = -dimer.d*onethird - dimer.e
	Tmp3[3,3] = 2. * dimer.d*onethird
    β = deg2rad(180. - dimer.β)*0.5
	rmat_active!(R, 0., -β, 0.)
	rotate!(DA,Tmp3,R)
	rmat_active!(R, 0., β, 0.)
	rotate!(DB,Tmp3,R)

    return nothing
end

"""
    anisotensor!(X,dimer) Computes the inter-chromophore anisotropic interaction tensor X, as defined in the dimer frame.

"""
function anisotensor!(X,dimer)
	sθ = sin(deg2rad(dimer.Θ))
	cθ = cos(deg2rad(dimer.Θ))
	sϕ = sin(deg2rad(dimer.Φ))
	cϕ = cos(deg2rad(dimer.Φ))
	X[1,1] = 1. - (3. *sθ^2 *cϕ^2)
	X[2,1] = -3. *sθ^2 *sϕ *cϕ
	X[3,1] = -3. *sθ *cθ *cϕ
	X[1,2] = X[2,1]
	X[2,2] = 1. - (3. *sθ^2 *sϕ^2)
	X[3,2] = -3. *sθ *cθ *sϕ
	X[1,3] = X[3,1]
	X[2,3] = X[3,2]
	X[3,3] = 1. - (3. *cθ^2)
	X .*= dimer.x

    return nothing
end

"""
    aniso!(H0,R,Tmp3,Tmp9,DA,DB,X,SA,SB,fθ,fϕ) 

Computes the anisotropic interaction hamiltonian `H0 = HA + HB + HX` for a given orientation `(θ,ϕ)` with respect to the field.

"""
function aniso!(H0,R,Tmp3,Tmp9,DA,DB,X,SA,SB,θ,ϕ)
    H0 .= zero(H0)
    rmat_passive!(R,deg2rad(ϕ),deg2rad(θ),0.0)

    rotate!(Tmp3,DA,R)
    for j = 1:3, i = 1:3
        Tmp9 .= SA[i]*SA[j]
        H0 .+= Tmp3[i,j] .* Tmp9
    end
    rotate!(Tmp3,DB,R)
    for j = 1:3, i = 1:3
        Tmp9 .= SB[i]*SB[j]
        H0 .+= Tmp3[i,j] .* Tmp9
    end
    rotate!(Tmp3,X,R)
    for j = 1:3, i = 1:3
        Tmp9 .= SA[i]*SB[j]
        H0 .+= Tmp3[i,j] .* Tmp9
    end

    return nothing
end

