"""
`cgmatrix!(C)` populates a preallocated matrix with Clebsch-Gordan
                coefficients for S = 0-2.
                Outputs:
                C - A matrix containing the required Clebsch-Gordan coefficient.
                        m1m2 are ordered down rows, ex: m1m2 = -1-1; -10; -11; 0-1;...
                        s increases across columns, ex: s = 0, 1, 1, 1, 2, 2, 2, 2, 2,...
                        m increases across columns, ex: m = 0, -1, 0, 1, -2, -1, 0, 1,...
"""
function cgmatrix!(C)
    C .= zero(C)
    one_sqrt2 = 1. /sqrt(2.)
    one_sqrt3 = 1. /sqrt(3.)
    one_sqrt6 = 1. /sqrt(6.)
    sqrt2_3   = sqrt(2. /3.)
    C[3,1] =  one_sqrt3
    C[5,1] = -one_sqrt3
    C[7,1] =  one_sqrt3
    C[2,2] = -one_sqrt2
    C[4,2] =  one_sqrt2
    C[3,3] = -one_sqrt2
    C[7,3] =  one_sqrt2
    C[6,4] = -one_sqrt2
    C[8,4] =  one_sqrt2
    C[1,5] =  1.0
    C[2,6] =  one_sqrt2
    C[4,6] =  one_sqrt2
    C[3,7] =  one_sqrt6
    C[5,7] =  sqrt2_3
    C[7,7] =  one_sqrt6
    C[6,8] =  one_sqrt2
    C[8,8] =  one_sqrt2
    C[9,9] =  1.0
    return nothing
end

"""
`spin9()` generates two vectors SA and SB where SA = [sAx, sAy, sAz] and SB = [sBx, sBy, sBz]

"""
function spin9()
    C = zeros(Complex{Float64},9,9)
    cgmatrix!(C)
    sx = zeros(Complex{Float64},3,3)
    sy = zeros(Complex{Float64},3,3)
    sz = zeros(Complex{Float64},3,3)
    I3 = Array{Complex{Float64},2}(I(3))

    sx[2,1] = 1.0
    sx[1,2] = 1.0
    sx[3,2] = 1.0
    sx[2,3] = 1.0
    sx ./=sqrt(2.)

    sy[2,1] = -1.0
    sy[1,2] = 1.0
    sy[3,2] = -1.0
    sy[2,3] = 1.0
    sy .*= im/sqrt(2.)

    sz[1,1] = -1.0
    sz[3,3] = 1.0

    Asx = kron(sx,I3)
    UtransAU!(C,Asx)
    Asy = kron(sy,I3)
    UtransAU!(C,Asy)
    Asz = kron(sz,I3)
    UtransAU!(C,Asz)

    SA = [Asx,Asy,Asz]

    Bsx = kron(I3,sx)
    UtransAU!(C,Bsx)
    Bsy = kron(I3,sy)
    UtransAU!(C,Bsy)
    Bsz = kron(I3,sz)
    UtransAU!(C,Bsz)

    SB = [Bsx,Bsy,Bsz]
    S = [SA[1]+SB[1], SA[2]+SB[2], SA[3]+SB[3]]

    return SA, SB, S
end