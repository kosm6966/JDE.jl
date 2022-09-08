function polarplot(dimer,plottype::Populations;ϕ::Vector{Float64}=[0.:0.5:360.;],θ::Vector{Float64}=[LinRange( 0., 90., Int(floor(0.5*(length(ϕ)+1))));])

    Populations = calculatepopdist(dimer,θ,ϕ)

    θ .= deg2rad.(θ)
    ϕ .= deg2rad.(ϕ)

    fig = drawpopplot(θ,ϕ,Populations,plottype.M)

	return fig
end

function polarplot(dimer,plottype::Polarization;ϕ::Vector{Float64}=[0.:0.5:360.;],θ::Vector{Float64}=[LinRange( 0., 90., Int(floor(0.5*(length(ϕ)+1))));])

    Populations = calculatepopdist(dimer,θ,ϕ)

    θ .= deg2rad.(θ)
    ϕ .= deg2rad.(ϕ)

    fig = drawpolarplot(θ,ϕ,Populations)

	return fig
end

function calculatepopdist(dimer,θ,ϕ)
    SA,SB,_ = spin9()

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
    Tmp9 = zeros(Complex{Float64},9,9)

    Populations = Array{Float64,3}(undef,length(ϕ),length(θ),5)

	for i in eachindex(ϕ)
		for j in eachindex(θ)
			aniso!(H0,R,Tmp3,Tmp9,DA,DB,X,SA,SB,θ[j],ϕ[i])
			for k in 1:5
				# Q(-2),Q(-1),Q( 0),Q(+1),Q(+2)
		        Populations[i,j,k] = real(abs( H0[1,4+k] )^2)
			end
			Populations[i,j,:] = Populations[i,j,:]/sum(Populations[i,j,:])
		end
	end
	return real(Populations)
end

function drawpopplot(theta,phi,pop,m)
    n = m+3
    a = "Population, M= "
    # Make polar plot
    if m > 0
        a = string(a,"+",m)
    else a = string(a,m); end
    pyplot()
    fig = contourf(phi,theta, transpose(pop[:,:,n]),xticks=(pi.*[0,1/2,1,3/2],round.(Int,180 .*[0,1/2,1,3/2])), yticks=(pi.*[0,1/6,1/3,1/2],round.(Int,180 .*[0,1/6,1/3,1/2])),proj=:polar, grid=true, levels=LinRange(0.,1.,11),  alpha=0.75,clims=(0.,1.),bg_inside=nothing,color=cgrad(:Spectral_11, categorical = true,  rev=true),title=a)
    return fig
end

function drawpolarplot(theta,phi,pop)
    Entropy = calculateentropy(theta,phi,pop)
    pyplot()
    fig = contourf(phi,theta, transpose(Entropy),xticks=(pi.*[0,1/2,1,3/2],round.(Int,180 .*[0,1/2,1,3/2])), yticks=(pi.*[0,1/6,1/3,1/2],round.(Int,180 .*[0,1/6,1/3,1/2])),proj=:polar, grid=true, levels=LinRange(0.,1.,11),  alpha=0.75,clims=(0.,1.),bg_inside=nothing,color=cgrad(:Spectral_11, categorical = true,  rev=true),title="Polarization")
    return fig
end

function calculateentropy(theta,phi,pop)
    Entropy = ones(length(phi),length(theta))
    pop .+= 1e-20
    for i in eachindex(phi)
        for j in eachindex(theta)
            for k in 1:5
                Entropy[i,j] += pop[i,j,k]*log2(pop[i,j,k])*0.43067655807339306
            end
        end
    end
    return Entropy
end