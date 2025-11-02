push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/spinlessPQMC/Hopping/")
using DelimitedFiles

using KAPDQMC_spinless_H
using LinearAlgebra
using Random


function main()
    rng=MersenneTwister(time_ns())

    t=1;   Lattice="HoneyComb60"    
    U=5;     Δt=0.05;     Θ=3.0;
    BatchSize=5;

    L=3
    site=[L,L]

    model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")
    println(model.nodes)

    s=Initial_s(model,rng)

    path="C:/Users/admin/Desktop/JuliaDQMC/code/spinlessPQMC/test/"
    s=phy_update(path,model,s,1,true)

end


# main()

# -----------------------------------------------
t=1;   Lattice="HoneyComb60"    
U=0;     Δt=0.02;     Θ=2.0;
BatchSize=5;
L=9
site=[L,L]

rng=MersenneTwister(2)
model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"H0")
s=Initial_s(model,rng)

println(model.α)

G=Gτ(model,s,div(model.Nt,2))


E,V,R0,R1=phy_measure(model,G,div(model.Nt,2),s)  
println("E: ",E,"  V: ",V)
println("R0: ",R0,"\nR1: ",R1)
# -----------------------------------------------


# -----------------------------------------------
BR=model.HalfeK*model.Pt

BL=model.Pt'*model.HalfeKinv


println(norm(BL*BR-I(div(model.Ns,2))))

# -----------------------------------------------


# s=phy_update(path,model,s,30,true)


# # # Half
# indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
# println(indexA)

# # # HalfHalf
# indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))
# println(indexB)

# ss=[s[:,:,:],s[:,:,:]]
# λ=0.5
# Nλ=2
# Sweeps=10

# ss=ctrl_SCEEicr(path,model,indexA,indexB,Sweeps,λ,Nλ,ss,true)
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









