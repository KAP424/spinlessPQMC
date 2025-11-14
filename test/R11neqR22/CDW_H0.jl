push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/spinlessPQMC/Hopping/")
using DelimitedFiles
using BenchmarkTools
using KAPDQMC_tV
using LinearAlgebra
using Random

Lattice="HoneyComb120"    

Θ=collect(0:0.1:6.0)
dΘ=Θ[2]-Θ[1]

L= [6]

path="./"
Initial="V"

R=zeros(length(Θ),5)
R[:,1]=Θ

for i in eachindex(L)

    site=[L[i],L[i]]
    model=Hubbard_Para(1,0,Lattice,site,0.1,0,10,"V")
    rng=MersenneTwister(time_ns()+Threads.threadid())
    s=Initial_s(model,rng)

    K=K_Matrix(Lattice,site)
    dt=0.05
    E,V=LAPACK.syevd!('V', 'L',K[:,:])
    eK=V*Diagonal(exp.(-dt.*E))*V'

    Ns=size(K)[1]
    ns=div(Ns, 2)

    G = Array{Float64}(undef, Ns, Ns)
    BL = Array{Float64}(undef, ns, Ns)
    BR = Array{Float64}(undef, Ns, ns)
    tmpN = Vector{Float64}(undef, Ns)
    tmpNn = Matrix{Float64}(undef, Ns, ns)
    tmpNN = Matrix{Float64}(undef, Ns, Ns)
    tmpnN = Matrix{Float64}(undef, ns, Ns)
    tmpnn= Matrix{Float64}(undef, ns, ns)
    tau = Vector{Float64}(undef, ns)
    ipiv = Vector{LAPACK.BlasInt}(undef, ns)


    Pt = zeros(Float64, Ns, div(Ns, 2))  # 预分配 Pt
    if Initial=="H0"
        KK=K[:,:]
        μ=1e-5
        if occursin("HoneyComb", Lattice)
            KK+=μ*Diagonal(repeat([-1, 1], div(Ns, 2)))
        elseif Lattice=="SQUARE"
            for i in 1:Ns
                x,y=i_xy(Lattice,site,i)
                KK[i,i]+=μ*(-1)^(x+y)
            end
        end
        E,V=LAPACK.syevd!('V', 'L',K[:,:])
        Pt=V[:,1:div(Ns,2)]
    elseif Initial=="V" 
        if occursin("HoneyComb", Lattice)
            for i in 1:Int(Ns/2)
                Pt[i*2,i]=1
            end
        else
            count=1
            for i in 1:Ns
                x,y=i_xy(Lattice,site,i)
                if (x+y)%2==1
                    Pt[i,count]=1
                    count+=1
                end
            end
        end
    end

    BL.=Pt'
    BR.=Pt

    mul!(tmpnn,BL,BR)
    LAPACK.getrf!(tmpnn,ipiv)
    LAPACK.getri!(tmpnn, ipiv)
    mul!(tmpNn,BR,tmpnn)
    mul!(G, tmpNn,BL)
    lmul!(-1.0,G)
    for i in diagind(G)
        G[i]+=1
    end

    tmp=phy_measure(tmpN,tmpNN,model,G,div(model.Nt,2),s)[3]
    R[1,2:5]=tmp

    for j in 2:length(Θ)
        for iii in 1:div(dΘ,dt)
            BL.=BL*eK
            BR.=eK*BR
            
            LAPACK.gerqf!(BL, tau)
            LAPACK.orgrq!(BL, tau, ns)

            LAPACK.geqrf!(BR, tau)
            LAPACK.orgqr!(BR, tau, ns)
        end
        mul!(tmpnn,BL,BR)
        LAPACK.getrf!(tmpnn,ipiv)
        LAPACK.getri!(tmpnn, ipiv)
        mul!(tmpNn,BR,tmpnn)
        mul!(G, tmpNn,BL)
        lmul!(-1.0,G)
        for i in diagind(G)
            G[i]+=1
        end
        tmp=phy_measure(tmpN,tmpNN,model,G,div(model.Nt,2),s)[3]
        R[j,2:5]=tmp    
    end

    open("$(path)R_CDW_H0.csv", "a") do io
        lock(io)
        writedlm(io, R, ',')
        unlock(io)
    end
    
end