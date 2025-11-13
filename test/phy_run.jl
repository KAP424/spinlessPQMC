push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/spinlessPQMC/Hopping/")
using DelimitedFiles
using KAPDQMC_tV
using LinearAlgebra
using Random


t=1;   Lattice="HoneyComb60"    
U=1;     Δt=0.05;     Θ=10.0;
BatchSize=5;
  

L=6
site=[L,L]

path="C:/Users/admin/Desktop/JuliaDQMC/code/spinlessPQMC/test/"

println("Threads: ",Threads.nthreads())

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"H0")

print(Threads.nthreads())

Threads.@threads for i in 1:Threads.nthreads()  
    rng=MersenneTwister(time_ns()+Threads.threadid())
    s=Initial_s(model,rng)
    s=phy_update(path,model,s,10,false)
    s=phy_update(path,model,s,50,true)
end

