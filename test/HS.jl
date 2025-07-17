push!(LOAD_PATH,"E:/桌面/JuliaDQMC/code/spinlessPQMC/source/")
using DelimitedFiles

using KAPDQMC_spinless
using LinearAlgebra
using Random


t=1;   Lattice="HoneyComb"    
V=0;     Δt=0.1;     Θ=1.2;
BatchSize=10;
  

L=12
site=[L,L]

model=Hubbard_Para(t,V,Lattice,site,Δt,Θ,BatchSize,"V")

D=zeros(model.Ns)

# # Half
# indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
# # println(indexA)

# # HalfHalf
# indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))
# # println(indexB)
# II=I(length(indexA[:]))

# rng=MersenneTwister(Int(floor(time_ns())))


# s=Initial_s(model,rng)
# G01=Gτ(model,s,div(model.Nt,2))
# s=Initial_s(model,rng)
# G02=Gτ(model,s,div(model.Nt,2))

# norm(G01-G02)

# gup1=GroverMatrix(G01[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
# gup2=GroverMatrix(G02[indexA[:],indexA[:]],G01[indexA[:],indexA[:]])

# Gdn1=G01'[:,:]
# Gdn2=G02'[:,:]
# for i in 1:model.Ns
#     for j in 1:model.Ns
#         Gdn1[i,j]=(-1)^(i+j)*Gdn1[i,j]
#         Gdn2[i,j]=(-1)^(i+j)*Gdn2[i,j]
#     end
# end


# gdn=GroverMatrix( Gdn1[indexA[:],indexA[:]] , Gdn2[indexA[:],indexA[:]] )

# U=zeros(size(II))
# for i in 1:size(U)[1]
#     # print((-1)^(i+1))
#     U[i,i]=(-1)^(i+1)
# end


# norm(U*gup2'*U-gdn)

# norm(det(gup1)-conj(det(gdn)))

# norm(det(gup1)-det(gup2))
# norm(G01[indexA[:],indexA[:]]*G02[indexA[:],indexA[:]]-G02[indexA[:],indexA[:]]*G01[indexA[:],indexA[:]])
# norm(gup1-transpose(gup2))

