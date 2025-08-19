push!(LOAD_PATH,"E:/桌面/JuliaDQMC/code/spinlessPQMC/Hopping/")
using DelimitedFiles

using KAPDQMC_spinless_M
using LinearAlgebra
using Random
rng=MersenneTwister(time_ns())

t=1;   Lattice="HoneyComb"    
U=8;     Δt=0.05;     Θ=0.2;
BatchSize=10;
  

L=3
site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"H0")


s=Initial_s(model,rng)
G0=Gτ(model,s,div(model.Nt,2))


path="E:/桌面/JuliaDQMC/code/spinlessPQMC/test/"
s=phy_update(path,model,s,1,true)


