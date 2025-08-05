push!(LOAD_PATH,"E:/桌面/JuliaDQMC/code/spinlessPQMC/CDW/")
using DelimitedFiles

using KAPDQMC_spinless_CDW
using LinearAlgebra
using Random
rng=MersenneTwister(1)

t=1;   Lattice="HoneyComb"    
U=8;     Δt=0.1;     Θ=0.3;
BatchSize=10;
  

L=6
site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")
nn2idx(Lattice,site,1)
# for x in 1:size(s)[2]
#     xidx=2*x-1
#     println(findall(model.K[xidx,:].!=0))
# end

# findall(model.K[1,:].!=0)

# ------------------------------------------------------------------------
# TEST for nn2idx
# for i in 1:18
#     print(i)
#     println(nn2idx(Lattice,site,i))
# end

# K=K_Matrix(Lattice,site)
# print(K)
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# TEST for Green function
s=Initial_s(model,rng)

# τ=model.Nt
# # τ=1
# G=Gτ(model,s,div(model.Nt,2))
# Gt,G0,Gt0,G0t=G4(model,s,τ,div(model.Nt,2))
# Gt0_,G0t_=G12FF(model,s,τ,div(model.Nt,2))
# println(norm(Gt0-Gt0_),',',norm(G0t-G0t_))
# Gt_=Gτ(model,s,τ)
# G0_=Gτ(model,s,div(model.Nt,2))
# println(norm(Gt-Gt_),',',norm(G0-G0_))
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# TEST for phy_update
path="E:/桌面/JuliaDQMC/code/spinlessPQMC/test/"
s=phy_update(path,model,s,1,true)


# lt=1
# x=1
# y=1
# G=Gτ(model,s,lt)

# ss=s[:,:,:]
# ss[lt,x,y]=-ss[lt,x,y]

# xidx=2*x-1
# yidx=findall(model.K[xidx,:].!=0)[y]
# Δ=diagm(exp.( 2*model.α.*[-s[lt,x,y],s[lt,x,y]] ))-I(2)
# subidx=[xidx,yidx]
# r=I(2)+Δ*(I(2)-G[subidx,subidx])
# println(det(r)-Poss(model,ss)/Poss(model,s))

# GG=Gτ(model,ss,lt)
# ((G-GG)./(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:])))
# ------------------------------------------------------------------------
# TEST for phy_measure
# G=Gτ(model,s,div(model.Nt,2))
# tmp=phy_measure(model,G,div(model.Nt,2),s).*sign(1.1)
# tmp.*sign(1)


# Ek=model.t*sum(model.K.*G)
# Ev=R0=R1=0
# for x in 1:div(model.Ns,2)
#     xidx=2*x-1
#     nnidx=findall(model.K[xidx,:].!=0)
#     for y in eachindex(nnidx)
#         Ev+=(1-G[xidx,xidx])*(1-G[nnidx[y],nnidx[y]])-G[xidx,nnidx[y]]*G[nnidx[y],xidx]-1/4
#     end
# end
# Ev


# for rx in 1:model.site[1]
#     for ry in 1:model.site[2]
#         tmp=0
#         for ix in 1:model.site[1]
#             for iy in 1:model.site[2]
#                 idx1=2*( mod(iy-1,site[2])*model.site[1]+mod1(ix,site[1]) )-1
#                 idx2=2*( mod(iy+ry-2,site[2])*model.site[1]+mod1(ix+rx,site[1]) )-1
#                 tmp+=(1-G[idx1,idx1])*(1-G[idx2,idx2])-G[idx1,idx2]*G[idx2,idx1]
#                 tmp+=(1-G[idx1+1,idx1+1])*(1-G[idx2+1,idx2+1])-G[idx1+1,idx2+1]*G[idx2+1,idx1+1]
#                 tmp-=(1-G[idx1+1,idx1+1])*(1-G[idx2,idx2])-G[idx1+1,idx2]*G[idx2,idx1+1]
#                 tmp-=(1-G[idx1,idx1])*(1-G[idx2+1,idx2+1])-G[idx1,idx2+1]*G[idx2+1,idx1]
#             end
#         end
#         tmp/=prod(model.site)
#         R0+=tmp
#         R1+=cos(π/model.site[1]*(rx+ry))*tmp
#     end
# end
# R0
# R1
# 1-R1/R0
# ------------------------------------------------------------------------


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

