push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/spinlessPQMC/Hopping/")
using DelimitedFiles
using KAPDQMC_tV
using LinearAlgebra
using Random


t=1;   Lattice="HoneyComb120"    
U=1;     Δt=0.025;     Θ=4.0;
BatchSize=5;

println("Threads: ",Threads.nthreads())
path="C:/Users/admin/Desktop/JuliaDQMC/code/spinlessPQMC/test/Nicr/"
# 3,6,9,12,15
for L in [6,9,12,15]
    println("L = ",L)
    site=[L,L]

    model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")
    println(model.nodes)
    Threads.@threads for i in 1:Threads.nthreads()  
        rng=MersenneTwister(time_ns()+Threads.threadid())

        # Half
        indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
        # println(indexA)

        # HalfHalf
        indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))
        # println(indexB)

        ss=[Initial_s(model,rng),Initial_s(model,rng)]
        λ=0.
        Nλ=1
        Sweeps=20

        ss=ctrl_SCEEicr(path,model,indexA,indexB,5,λ,Nλ,ss,false)
        ss=ctrl_SCEEicr(path,model,indexA,indexB,Sweeps,λ,Nλ,ss,true)
    end
end

