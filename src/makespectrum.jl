function quintetspectrum(dimer,spec,exp)

    Signal, Intensity = quintethamiltonian(dimer,exp,spec.basis)

    S = makeline(Signal,Intensity,exp,spec)

	return S
end

function makeline(Signal,Intensity,exp,spec)
    b0 = exp.field
    f(idx) = linefunc(spec.line,spec.fwhm,b0,idx)

    if spec.multiOut
        S = zeros(Float64,length(b0),length(Signal[1,:]),4)
        for j in eachindex(Signal[1,:])
            for k = 1:4
                index = argmin(abs.(b0.-Signal[k,j]))
                S[:,j,k] .+= Intensity[k,j].*f(b0[index]).*exp.w[j]
            end
        end
        n = maximum(abs.(sum(S,dims=3)))
        if spec.nOutput == 2
            S[:,:,1] .+= S[:,:,4]
            S[:,:,2] .+= S[:,:,3]
            S = S[:,:,1:2]
        elseif spec.nOutput == 1
            S = sum(S,dims=3)
            S = S[:,:]
        end
    else
        S = zeros(Float64,length(b0),4)
        for j in eachindex(Signal[1,:])
            for k = 1:4
                index = argmin(abs.(b0.-Signal[k,j]))
                S[:,k] .+= Intensity[k,j].*f(b0[index]).*exp.w[j]
            end
        end
        n = maximum(abs.(sum(S,dims=2)))
        if spec.nOutput == 2
            S[:,1] .+= S[:,4]
            S[:,2] .+= S[:,3]
            S = S[:,1:2]
        elseif spec.nOutput == 1
            S = sum(S,dims=2)
        end
    end

    S ./= n

    return S
end

function linefunc(l::Gaussian,fwhm,b0,bn)
    one_Γ = 1. /(fwhm*0.8493218002880191)
    ls = (0.7978845608028654*one_Γ * exp.(-2. *((b0.-bn).*one_Γ).^2. ))
    return ls
end

function linefunc(l::Lorentzian,fwhm,b0,bn)
    one_Γ = 1. /(fwhm*0.5773502691896258)
    ls = (0.3675525969478614*one_Γ ./(1 .+ 1.3333333333333333*((b0.-bn).*one_Γ).^2 ))
    return ls
end
