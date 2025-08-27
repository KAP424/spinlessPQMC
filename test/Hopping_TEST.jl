push!(LOAD_PATH,"E:/桌面/JuliaDQMC/code/spinlessPQMC/Hopping/")
using DelimitedFiles

using KAPDQMC_spinless_M
using LinearAlgebra
using Random
rng=MersenneTwister(1)

t=1;   Lattice="HoneyComb"    
U=5;     Δt=0.05;     Θ=2.0;
BatchSize=10;
  

L=3
site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")


s=Initial_s(model,rng)
# G0=Gτ(model,s,div(model.Nt,2))
# Gt,G0,Gt0,G0t=G4(model,s,div(model.Nt,2),2)

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

path="E:/桌面/JuliaDQMC/code/spinlessPQMC/test/"
s=phy_update(path,model,s,3,false)
s=phy_update(path,model,s,100,true)

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

