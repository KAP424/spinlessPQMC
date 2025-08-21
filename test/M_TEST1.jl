push!(LOAD_PATH,"E:/桌面/JuliaDQMC/code/spinlessPQMC/Majorana/")
using DelimitedFiles

using KAPDQMC_spinless_M
using LinearAlgebra
using Random
rng=MersenneTwister(time_ns())

t=1;   Lattice="HoneyComb"    
U=0;     Δt=0.05;     Θ=0.1;
BatchSize=10;
  

L=3
site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"H0")
s=Initial_s(model,rng)

path="E:/桌面/JuliaDQMC/code/spinlessPQMC/test/"
# s=phy_update(path,model,s,1,true)

G0=Gτ(model,s,div(model.Nt,2))

H0=model.t*1im*UpperTriangular(model.K)
H0=(H0+H0')/4

Ek=2*(real(sum(H0.*G0)))

println(Ek)


println((sum(real.(model.K.*(G0))))/2)


println((sum(imag.(model.K.*(G0))))/2)

# -model.U*sum(G.*conj.(G).*model.K)

# norm(G+transpose(G))

# G[1,3]+G[3,1]

# ------------------------------------------------------------------------
# TEST for Green function
# println(size(s))

# V=zeros(ComplexF64,3,model.Ns,model.Ns)
# lt=1
# for i in 1:size(s)[2]
#     for j in 1:size(s)[3]
#         x,y=model.nnidx[i,j]
#         V[j,x,y]=s[lt,i,j]*1im/2
#         V[j,y,x]=-s[lt,i,j]*1im/2
#     end
# end
# ------------------------------------------------------------------------
# TEST for Trotter V error
# Vall=zeros(ComplexF64,model.Ns,model.Ns)
# for i in 1:size(s)[2]
#     for j in 1:size(s)[3]
#         x,y=model.nnidx[i,j]
#         Vall[x,y]=s[lt,i,j]*1im/2
#         Vall[y,x]=-s[lt,i,j]*1im/2
#     end
# end

# norm(Vall-V[1,:,:]-V[2,:,:]-V[3,:,:])

# expV1=I(model.Ns)
# expV2=I(model.Ns)
# expV3=I(model.Ns)
# for i in (1,2,3)
#     E,U=eigen(V[i,:,:])
#     expV1*=U*diagm(exp.(model.α*E))*U'
# end
# for i in (1,3,2)
#     E,U=eigen(V[i,:,:])
#     expV2*=U*diagm(exp.(model.α*E))*U'
# end
# for i in (2,1,3)
#     E,U=eigen(V[i,:,:])
#     expV3*=U*diagm(exp.(model.α*E))*U'
# end
# norm(expV1-expV2)
# norm(expV1-expV3)
# norm(expV2-expV3)

# norm(V[1,:,:]*V[2,:,:]-V[2,:,:]*V[1,:,:])


# E,U=eigen(Vall)
# expVall=U*diagm(exp.(model.α*E))*U'

# norm(expV-expVall)
# ------------------------------------------------------------------------
# TEST for general UV
# UV=zeros(ComplexF64,3,model.Ns,model.Ns)
# for i in 1:size(UV)[1]
#     _,UV[i,:,:]=eigen(V[i,:,:])
#     UV[i,:,:]=UV[i,:,:]'
# end 

# for lt in  1:size(s)[1]
#     for j in 1:size(s)[3]
#         for i in 1:size(s)[2]
#             x,y=model.nnidx[i,j]
#             V[j,x,y]=s[lt,i,j]*1im/2
#             V[j,y,x]=-s[lt,i,j]*1im/2
#         end
#         tmp=UV[j,:,:]*V[j,:,:]*UV[j,:,:]'
#         if norm(tmp-diagm(diag(tmp)))>1e-4
#             println(j)
#         end
#     end
# end
# ------------------------------------------------------------------------
# TEST for update
# s=Initial_s(model,rng)
# lt=1
# i=1
# j=3

# G=Gτ(model,s,lt)
# x,y=model.nnidx[i,j]
# subidx=[x,y]


# Δ=[0 -1im*s[lt,i,j]; 1im*s[lt,i,j] 0]
# E,U=eigen(Δ)
# Δ=U*diagm(exp.(model.α*E))*U'-I(2)
# r=I(2)+Δ*(I(2)-G[subidx,subidx])

# G1=G-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))

# ss=copy(s)
# ss[lt,i,j]=-ss[lt,i,j]
# findall((ss-s).!=0)
# GG=Gτ(model,ss,lt)

# norm(G1-GG)

# V=zeros(ComplexF64,3,model.Ns,model.Ns)
# for i in 1:size(s)[2]
#     for j in 1:size(s)[3]
#         x,y=model.nnidx[i,j]
#         V[j,x,y]=s[lt,i,j]*1im/2
#         V[j,y,x]=-s[lt,i,j]*1im/2
#     end
# end
# findall(V[j,:,:].!=0)
# findall(VV.!=0)



# lt=1
# i=1
# j=2
# x,y=model.nnidx[i,j]
# VV=V[:,:,:]
# VV[j,x,y]=-VV[j,x,y]
# VV[j,y,x]=-VV[j,y,x]

# norm(V[j,:,:]*VV[j,:,:]-VV[j,:,:]*V[j,:,:])
# findall(VV.!=V)

# ΔΔ=diag(model.UV[j,:,:]*(VV[j,:,:]-V[j,:,:])*model.UV[j,:,:]')
# ΔΔ=model.UV[j,:,:]'*diagm(exp.(model.α*ΔΔ))*model.UV[j,:,:]-I(model.Ns)

# idx=findall(abs.(ΔΔ).>1e-3)
# ΔΔ[subidx,subidx]
# rr=I(model.Ns)+ΔΔ*(I(model.Ns)-G)
# G2=G-(G/rr*ΔΔ*((I(model.Ns)-G)))

# norm(G2-G1)
# norm(G2-GG)

# BR=model.eK*model.Pt
# BR_=model.eK*model.Pt
# for j in 1:size(V)[1]
#     E=diag(model.UV[j,:,:]*V[j,:,:]*model.UV[j,:,:]')
#     BR=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:]*BR

#     E=diag(model.UV[j,:,:]*VV[j,:,:]*model.UV[j,:,:]')
#     BR_=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:]*BR_
# end
# norm( (I(model.Ns)+ΔΔ)*BR-BR_ )



# ΔΔ=(V[j,:,:]-VV)
# E,U=eigen(ΔΔ)
# ΔΔ=U*diagm(exp.(model.α*E))*U'-I(model.Ns)
# subidx=[x,y]
# ΔΔ[subidx,subidx]


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
# # TEST for phy_update 正方向

# lt=0

# G=model.eK*Gτ(model,s,lt)*model.eKinv

# V=zeros(ComplexF64,3,model.Ns,model.Ns)

# for j in size(s)[3]:-1:1
# # j=3
#     for i in 1:size(s)[2]
#         x,y=model.nnidx[i,j]
#         V[j,x,y]=s[lt,i,j]*1im/2
#         V[j,y,x]=-s[lt,i,j]*1im/2
#     end
#     E=diag(model.UV[j,:,:]*V[j,:,:]*model.UV[j,:,:]')
#     G=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:] *G* model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:]
    
#     for i in 1:size(s)[2]
#     # i=2
#         x,y=model.nnidx[i,j]
#         subidx=[x,y]

#         Δ=[0 -1im*s[lt,i,j]; 1im*s[lt,i,j] 0]
#         E,U=eigen(Δ)
#         Δ=U*diagm(exp.(model.α*E))*U'-I(2)
#         r=I(2)+Δ*(I(2)-G[subidx,subidx])

#         detR=abs2(det(r))

#         println(detR)
#         if rand(rng)<detR
#             G-=(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
#             s[lt,i,j]=-s[lt,i,j]
        
#             ss=copy(s)
#             GG=model.eK*Gτ(model,ss,lt)*model.eKinv
#             VV=zeros(ComplexF64,3,model.Ns,model.Ns)
#             for jj in size(s)[3]:-1:j
#                 for ii in 1:size(s)[2]
#                     x,y=model.nnidx[ii,jj]
#                     VV[jj,x,y]=ss[lt,ii,jj]*1im/2
#                     VV[jj,y,x]=-ss[lt,ii,jj]*1im/2
#                 end
#                 E=diag(model.UV[jj,:,:]*VV[jj,:,:]*model.UV[jj,:,:]')
#                 GG=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *GG* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
#             end

#             if(norm(G-GG)>1e-4)
#                 println(j," error: ",norm(G-GG),"  ",i)
#             end
#         end

#     end
# end


# # ------------------------------------------------------------------------
# REVERSE
# lt=model.Nt

# G=Gτ(model,s,lt)

# V=zeros(ComplexF64,3,model.Ns,model.Ns)

# for j in 1:size(s)[3]
# # j=1
#     for i in 1:size(s)[2]
#     # i=2
#         x,y=model.nnidx[i,j]
#         subidx=[x,y]

#         Δ=[0 -1im*s[lt,i,j]; 1im*s[lt,i,j] 0]
#         E,U=eigen(Δ)
#         Δ=U*diagm(exp.(model.α*E))*U'-I(2)
#         r=I(2)+Δ*(I(2)-G[subidx,subidx])

#         detR=abs2(det(r))

#         println(detR)
#         if rand(rng)<detR
#             G-=(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
#             s[lt,i,j]=-s[lt,i,j]
        
#             ss=copy(s)
#             GG=Gτ(model,ss,lt)
#             VV=zeros(ComplexF64,3,model.Ns,model.Ns)
#             for jj in 1:j-1
#                 for ii in 1:size(s)[2]
#                     x,y=model.nnidx[ii,jj]
#                     VV[jj,x,y]=ss[lt,ii,jj]*1im/2
#                     VV[jj,y,x]=-ss[lt,ii,jj]*1im/2
#                 end
#                 E=diag(model.UV[jj,:,:]*VV[jj,:,:]*model.UV[jj,:,:]')
#                 GG=model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:] *GG* model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:]
#             end

#             if(norm(G-GG)>1e-4)
#                 println(j," error: ",norm(G-GG),"  ",i)
#             end
#         end

#     end

#     for i in 1:size(s)[2]
#         x,y=model.nnidx[i,j]
#         V[j,x,y]=s[lt,i,j]*1im/2
#         V[j,y,x]=-s[lt,i,j]*1im/2
#     end
#     E=diag(model.UV[j,:,:]*V[j,:,:]*model.UV[j,:,:]')
#     G=model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:] *G* model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:]
# end

# # ------------------------------------------------------------------------
