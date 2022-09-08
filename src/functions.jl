function conjtransUAU!(U,A)
    A .= U'*A
    A .= A*U
    return nothing
end

function UtransAU!(U,A)
    A .= transpose(U)*A
    A .= A*U
    return nothing
end

function rmat_passive!(R,α,β,γ)
    cα = cos(α)
    sα = sin(α)
    cβ = cos(β)
    sβ = sin(β)
    cγ = cos(γ)
    sγ = sin(γ)

    R[1,1] =  cα*cβ*cγ-sα*sγ
    R[2,1] = -cα*cβ*sγ-sα*cγ
    R[3,1] =  cα*sβ
    R[1,2] =  sα*cβ*cγ+cα*sγ
    R[2,2] = -sα*cβ*sγ+cα*cγ
    R[3,2] =  sα*sβ
    R[1,3] = -sβ*cγ
    R[2,3] =  sβ*sγ
    R[3,3] =  cβ
    return nothing
end

function rmat_active!(R,α,β,γ)
    cα = cos(α)
    sα = sin(α)
    cβ = cos(β)
    sβ = sin(β)
    cγ = cos(γ)
    sγ = sin(γ)

    R[1,1] =  cα*cβ*cγ-sα*sγ
    R[2,1] =  sα*cβ*cγ+cα*sγ
    R[3,1] = -sβ*cγ
    R[1,2] = -cα*cβ*sγ-sα*cγ
    R[2,2] = -sα*cβ*sγ+cα*cγ
    R[3,2] =  sβ*sγ
    R[1,3] =  cα*sβ
    R[2,3] =  sα*sβ
    R[3,3] =  cβ
    return nothing
end

function rotate!(H,H0,R)
    H .= R*H0*transpose(R)
    return nothing
end