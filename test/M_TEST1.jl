push!(LOAD_PATH,"E:/桌面/JuliaDQMC/code/spinlessPQMC/Majorana/")
using DelimitedFiles

using KAPDQMC_spinless_M
using LinearAlgebra
using Random
rng=MersenneTwister(1)

t=1;   Lattice="HoneyComb"    
U=8;     Δt=0.05;     Θ=0.1;
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
model.nnidx

for k in 1:size(s)[2]
    x,y=model.nnidx[k].I
    println(x," ",y)
end

# ------------------------------------------------------------------------
# TEST for phy_update
path="E:/桌面/JuliaDQMC/code/spinlessPQMC/test/"
s=phy_update(path,model,s,1,false)

i=1
k=2

G=Gτ(model,s,i)

ss=s[:,:]
ss[i,k]=-ss[i,k]
GG=Gτ(model,ss,i)

subidx=collect(model.nnidx[k].I)
Δ=[0 -1im*s[i,k]; 1im*s[i,k] 0]
E,V=eigen(Δ) 
Δ=V*diagm(exp.(model.α.*E))*V'- I(2)
r=I(2)+Δ*(I(2)-G[subidx,subidx])
detR=abs2(det(r))

norm(Δ/r-I(2)/r*Δ)

G1=G-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
norm((G-GG)-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:])))

println(norm(G1-GG))
println(detR-Poss(model,ss)/Poss(model,s))
# ------------------------------------------------------------------------

D=zeros(ComplexF64,model.Ns,model.Ns)
for k in 1:size(s)[2]
    x,y=model.nnidx[k].I
    D[x,y]+=s[i,k]*1im/2
    D[y,x]-=s[i,k]*1im/2
end
D


DD=zeros(ComplexF64,model.Ns,model.Ns)
for k in 1:size(s)[2]
    x,y=model.nnidx[k].I
    DD[x,y]=ss[i,k]*1im/2
    DD[y,x]=-ss[i,k]*1im/2
end
DD
findall((s-ss).!=0)
(DD-D)[subidx,subidx]
findall((DD-D).!=0)
E,V=eigen(DD-D)
Δ1 =V*diagm(exp.(model.α*E))*V'-I(model.Ns)
norm(Δ1[subidx,subidx]-Δ)
r1=I(model.Ns)+Δ1*(I(model.Ns)-G)
norm(Δ1/r1-I(model.Ns)/r1*Δ1)
G2=G-(G/r1*Δ1*((I(model.Ns)-G)))
norm(GG-G2)
norm(G2-G1)
abs2(det(r1))-detR


# ------------------------------------------------------------------------

# BL,BR=Gτ(model,s,i)
# BL_,BR_=Gτ(model,ss,i)
# norm(BL-BL_)
# ΔΔ=zeros(ComplexF64,model.Ns,model.Ns)
# ΔΔ[subidx,subidx]=Δ
# norm( (I(model.Ns)+ΔΔ)*BR -BR_   )

# norm(ΔΔ-Δ1)

E,V=eigen(D)
EV1=V*diagm(exp.(model.α*E))*V'
E,V=eigen(DD)
EV2=V*diagm(exp.(model.α*E))*V'
norm(EV1*EV2-EV2*EV1)

ΔΔ =EV2/EV1-I(model.Ns)
rr=I(model.Ns)+ΔΔ*(I(model.Ns)-G)

norm((I(model.Ns)+ΔΔ)*EV1-EV2)

G3=G-(G/rr*ΔΔ*((I(model.Ns)-G)))

norm(G3-GG)



ΔΔ
for iii in 1:size(ΔΔ)[1]
    for jjj in 1:size(ΔΔ)[2]
        if abs(ΔΔ[iii,jjj])<1e-5
            ΔΔ[iii,jjj]=0
        end
    end
end
ΔΔ
