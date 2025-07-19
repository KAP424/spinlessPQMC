push!(LOAD_PATH,"E:/桌面/JuliaDQMC/code/spinlessPQMC/source/")
using DelimitedFiles

using KAPDQMC_spinless
using LinearAlgebra
using Random
rng=MersenneTwister(1)

t=1;   Lattice="SQUARE"    
U=1;     Δt=0.1;     Θ=0.3;
BatchSize=10;
  

L=4
site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")

# --------------------------------
# TEST for nn2idx
# for i in 1:18
#     print(i)
#     println(nn2idx(Lattice,site,i))
# end

# K=K_Matrix(Lattice,site)
# print(K)
# --------------------------------


# --------------------------------
# TEST for Green function
s=Initial_s(model,rng)

# τ=model.Nt
# G=Gτ(model,s,div(model.Nt,2))
# Gt,G0,Gt0,G0t=G4(model,s,τ,div(model.Nt,2))
# Gt0_,G0t_=G12FF(model,s,τ,div(model.Nt,2))
# println(norm(Gt0-Gt0_),',',norm(G0t-G0t_))
# Gt_=Gτ(model,s,τ)
# G0_=Gτ(model,s,div(model.Nt,2))
# println(norm(Gt-Gt_),',',norm(G0-G0_))
# --------------------------------


# --------------------------------
# TEST for phy_update
path="E:/桌面/JuliaDQMC/code/spinlessPQMC/test/"
s=phy_update(path,model,1,s)

lt=2
x=1
y=1
G=Gτ(model,s,lt)

xidx=2*x-1
yidx=findall(model.K[xidx,:].!=0)[y]
Δ=diagm(exp.( 2*model.α.*[-s[lt,x,y],s[lt,x,y]] ))-I(2)
subidx=[xidx,yidx]
# r=1+tr(Δ*(I(2)-G[subidx,subidx]) )
r=det( I(2)+Δ*(I(2)-G[subidx,subidx]) )
G=G-(G[:,subidx]*Δ*(I(model.Ns)-G)[subidx,:])./( 1+tr(Δ*(I(2)-G[subidx,subidx]) ) )

ss=s[:,:,:]
ss[lt,x,y]=-ss[lt,x,y]

GG=Gτ(model,ss,lt)
println(findall(model.K[xidx,:].!=0))
println(r-Poss(model,ss)/Poss(model,s))
println(norm(G-GG))

A=diagm([0,1,0,-1,0])
B=reshape(collect(1:25),5,5)
A*B

collect(1:2:4)
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

