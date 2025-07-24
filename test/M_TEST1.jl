push!(LOAD_PATH,"E:/桌面/JuliaDQMC/code/spinlessPQMC/Majorana/")
using DelimitedFiles

using KAPDQMC_spinless_M
using LinearAlgebra
using Random
rng=MersenneTwister(1)

t=1;   Lattice="HoneyComb"    
U=8;     Δt=0.1;     Θ=0.5;
BatchSize=10;
  

L=3
site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")


# ------------------------------------------------------------------------
# TEST for Green function
s=Initial_s(model,rng)


# τ=model.Nt
# τ=1
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
# s=phy_update(path,model,s,1,false)

lt=1
k=1

G=Gτ(model,s,lt)

ss=s[:,:]
ss[lt,k]=-ss[lt,k]

x,y=model.nnidx[k].I
Δ=[0 1im*s[lt,k]; -1im*s[lt,k] 0]
E,V=eigen(Δ) 
Δ=V*diagm(exp.(model.α*E))*V'- I(2)
subidx=[x,y]
r=I(2)+Δ*(I(2)-G[subidx,subidx])
detR=abs2(det(r))

GG=Gτ(model,ss,lt)
G=G-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
# ((G-GG)-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:])))


println(norm(G-GG))
println(detR-Poss(model,ss)/Poss(model,s))

subidx
D=model.K[:,:]
for k in 1:size(s)[2]
    x,y=model.nnidx[k].I
    D[x,y]*=s[lt,x]
    D[y,x]*=-s[lt,x]
end

DD=model.K[:,:]
for k in 1:size(s)[2]
    x,y=model.nnidx[k].I
    DD[x,y]*=ss[lt,x]
    DD[y,x]*=-ss[lt,x]
end

findall((s-ss).!=0)
(DD-D)[subidx,subidx]
findall((DD-D).!=0)
DD-D
# ------------------------------------------------------------------------

# dim=3
# A=zeros(ComplexF64,dim,dim)

# A[1,2]=1im
# A[2,1]=-1im
# E,V=eigen(A)
# V*diagm(exp.(E))*V'-I(dim)


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

index=findall(UpperTriangular(model.K).!=0)

s=zeros(model.Nt,length(index))

model.K[index[1]]
x,y=index[1].I
x
reverse(idx)
idx