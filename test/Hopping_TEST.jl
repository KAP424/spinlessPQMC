push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/spinlessPQMC/Hopping/")
using DelimitedFiles
using BenchmarkTools
using KAPDQMC_tV
using LinearAlgebra
using Random

function main()
    rng=MersenneTwister(time_ns())

    t=1;   Lattice="HoneyComb60"    
    U=1;     Δt=0.05;     Θ=5.15;
    BatchSize=5;

    L=3
    site=[L,L]

    model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"H0")
    # println(model.nodes)

    s=Initial_s(model,rng)
    path="C:/Users/admin/Desktop/JuliaDQMC/code/spinlessPQMC/test/"

    s=phy_update(path,model,s,2,true)

    println(norm(model.UV[:,:,1]-model.UV[:,:,1]'))

    # Half
    indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
    # println(indexA)

    # HalfHalf
    indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))
    # println(indexB)

    ss=[s[:,:,:],s[:,:,:]]
    λ=0.5
    Nλ=2
    Sweeps=2

    # ss=ctrl_SCEEicr(path,model,indexA,indexB,Sweeps,λ,Nλ,ss,true)
end

main()
# println(@btime main())
# NO MKL 241.254 ms (204009 allocations: 15.66 MiB)
# With MKL 782.518 ms (204056 allocations: 15.87 MiB)
# First 304.746 ms (537237 allocations: 457.84 MiB)
# Secord 251.443 ms (204507 allocations: 20.56 MiB)
# 959.025 ms (810702 allocations: 63.29 MiB)
# 168.357 s (3793634 allocations: 866.99 MiB)

# -----------------------------------------------
# t=1;   Lattice="HoneyComb60"    
# U=0;     Δt=0.02;     Θ=0.0;
# BatchSize=5;
# L=9
# site=[L,L]

# rng=MersenneTwister(2)
# model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")
# s=Initial_s(model,rng)

# println(size(model.nnidx))
# println(length(model.nnidx))
# println(model.nnidx[1])


# G=Gτ(model,s,div(model.Nt,2))

# tmpN=Vector{Float64}(undef,model.Ns)
# tmpNN=Matrix{Float64}(undef,model.Ns,model.Ns)
# E,V,R0,R1=phy_measure(tmpN,tmpNN,model,G,div(model.Nt,2),s)  
# println("E: ",E,"  V: ",V)
# println("R0: ",R0,"\nR1: ",R1)
# -----------------------------------------------


# -----------------------------------------------
# BR=model.HalfeK*model.Pt

# BL=model.Pt'*model.HalfeKinv


# println(norm(BL*BR-I(div(model.Ns,2))))

# -----------------------------------------------


# s=phy_update(path,model,s,30,true)



# ss=ctrl_EEicr(path,model,indexA,Sweeps,λ,Nλ,ss,true)

# ----------------------------------------------------------------------------------------
# print(norm(Gt-G0))
# for _ in 1:1000
#     i=rand(1:size(s)[1])
#     j=rand(1:size(s)[2])
#     k=rand(1:size(s)[3])

#     s[i,j,k]=-s[i,j,k]

#     if Poss(model,s)<0
#         println("NNNNN")
#     else
#         print("-")
#     end
# end                 
# ----------------------------------------------------------------------------------------

# uv=[-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]

# uv'*[0 1 ;1 0]*uv

# lt=1
# j=1
# V=zeros(Float64,model.Ns,model.Ns)
# for i in 1:size(s)[2]
#     x,y=model.nnidx[i,j]
#     V[x,y]=V[y,x]=s[lt+1,i,j]
# end
# V
# E=diag(model.UV[j,:,:]*V*model.UV[j,:,:]')

# model.UV[j,:,:]'*diagm(E)*model.UV[j,:,:]-V
# V
# E

# VV=zeros(model.Ns)
# UV=zeros(model.Ns,model.Ns)
# for i in 1:size(s)[2]
#     x,y=model.nnidx[i,j]
#     UV[x,x]=UV[x,y]=UV[y,x]=-2^0.5/2
#     UV[y,y]=2^0.5/2
#     VV[x]=-s[lt+1,i,j]
#     VV[y]=s[lt+1,i,j]
# end

# E=diag(UV*V*UV')

# norm(E+VV)
# ----------------------------------------------------------------------------------------



# using LinearAlgebra

# A=rand(Float64,10,2)
# B=rand(Float64,2,10)

# r1=det(I(10)+A*B)

# r2=1+dot(A[:,1],B[1,:])+dot(A[:,2],B[2,:])
# r2=1+dot(A[:,1]+A[:,2],B[1,:]+B[2,:])

# r2=det(I(2)+B*A)


# function dot22(A,B)
#     return dot(A[:,1],B[1,:])+dot(A[:,2],B[2,:])+dot(A[:,1],B[2,:])*dot(A[:,2],B[1,:])
# end


# rho1=inv(I(10)+A*B)

# rho2=I(10)-A*inv(I(2)+B*A)*B

# norm(rho1-rho2)