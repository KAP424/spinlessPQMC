# For particle hole symmetric of t-V model
# attractive-U and repulsive-U get the same S_2

function ctrl_SCEEicr(path::String,model::_Hubbard_Para,indexA::Vector{Int64},indexB::Vector{Int64},Sweeps::Int64,λ::Float64,Nλ::Int64,ss::Vector{Array{UInt8,3}},record)
    global LOCK=ReentrantLock()
    Ns=model.Ns
    ns=div(Ns, 2)
    NN=length(model.nodes)
    tau = Vector{Float64}(undef, ns)

    Θidx=div(NN,2)+1

    uv=[-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]
    name = if model.Lattice=="SQUARE" "□" 
    elseif model.Lattice=="HoneyComb60" "HC60" 
    elseif model.Lattice=="HoneyComb120" "HC120" 
    else error("Lattice: $(model.Lattice) is not allowed !") end  

    file="$(path)/SCEEicr$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)N$(Nλ)BS$(model.BatchSize).csv"
    
    atexit() do
        if record
            lock(LOCK) do
                open(file, "a") do io
                    writedlm(io, O', ',')
                end
            end
        end
    end
    
    rng=MersenneTwister(Threads.threadid()+time_ns())

    G1=G4Workspace(
        Matrix{Float64}(undef ,Ns, Ns),
        Matrix{Float64}(undef ,Ns, Ns),
        Matrix{Float64}(undef ,Ns, Ns),
        Matrix{Float64}(undef ,Ns, Ns)      )
    G2=G4Workspace(
        Matrix{Float64}(undef ,Ns, Ns),
        Matrix{Float64}(undef ,Ns, Ns),
        Matrix{Float64}(undef ,Ns, Ns),
        Matrix{Float64}(undef ,Ns, Ns)      )

    GMA=GMWorkspace(
        indexA,
        0.0,
        Matrix{Float64}(undef ,length(indexA),length(indexA)),
        Matrix{Float64}(undef ,length(indexA),2),
        Matrix{Float64}(undef ,2,length(indexA)),
        Matrix{Float64}(undef ,2,2),
        Matrix{Float64}(undef, length(indexA), length(indexA)),
        Matrix{Float64}(undef, length(indexA), 2),    
        Matrix{Float64}(undef, 2, length(indexA)),
        Vector{LAPACK.BlasInt}(undef, length(indexA))     )
    GMB=GMWorkspace(
        indexB,
        0.0,
        Matrix{Float64}(undef ,length(indexB),length(indexB)),
        Matrix{Float64}(undef ,length(indexB),2),
        Matrix{Float64}(undef ,2,length(indexB)),
        Matrix{Float64}(undef ,2,2),
        Matrix{Float64}(undef, length(indexB), length(indexB)),
        Matrix{Float64}(undef, length(indexB), 2),    
        Matrix{Float64}(undef, 2, length(indexB)),
        Vector{LAPACK.BlasInt}(undef, length(indexB))     )
    tmpG=tmpSCEEWorkspace(
        Matrix{Float64}(undef, 2,2),
        Matrix{Float64}(undef, 2,2),
        Vector{Float64}(undef, Ns),
        Vector{Float64}(undef, Ns),
        Matrix{Float64}(undef, Ns, Ns),
        Matrix{Float64}(undef, Ns, Ns),
        Matrix{Float64}(undef, Ns, ns),
        Matrix{Float64}(undef, ns, Ns),
        Matrix{Float64}(undef, ns, ns),
        Matrix{Float64}(undef, Ns, 2),
        Matrix{Float64}(undef, 2, Ns),
        Matrix{Float64}(undef, 2, 2),
        Vector{Float64}(undef, 2),
        Vector{LAPACK.BlasInt}(undef, ns),
        [-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]
    )

    tmpO=0.0
    counter=0
    O=zeros(Float64,Sweeps+1)
    O[1]=λ

    BMs1=Array{Float64}(undef,Ns,Ns,NN-1)  # Number_of_BM*Ns*Ns
    BMs2=Array{Float64}(undef,Ns,Ns,NN-1)  # Number_of_BM*Ns*Ns
    BMsinv1=Array{Float64}(undef,Ns,Ns,NN-1)  # Number_of_BM*Ns*Ns
    BMsinv2=Array{Float64}(undef,Ns,Ns,NN-1)  # Number_of_BM*Ns*Ns

    for idx in 1:NN-1
        BM_F!(tmpG,view(BMs1,:, : , idx),model,ss[1],idx)
        BM_F!(tmpG,view(BMs2,:,:,idx),model,ss[2],idx)
        BMinv_F!(tmpG,view(BMsinv1,:,:,idx),model,ss[1],idx)
        BMinv_F!(tmpG,view(BMsinv2,:,:,idx),model,ss[2],idx)
        @assert norm(view(BMs1,:,:,idx)*view(BMsinv1,:,:,idx)-I(Ns))<1e-8 "BM1 inv error at idx=$idx"
    end

    BLMs1=Array{Float64}(undef,ns,Ns,NN)
    BRMs1=Array{Float64}(undef,Ns,ns,NN)
    transpose!(view(BLMs1,:,:,NN) , model.Pt)
    copyto!(view(BRMs1,:,:,1) , model.Pt)
    
    BLMs2=Array{Float64}(undef,ns,Ns,NN)
    BRMs2=Array{Float64}(undef,Ns,ns,NN)
    transpose!(view(BLMs2,:,:,NN) , model.Pt)
    copyto!(view(BRMs2,:,:,1) , model.Pt)


    # 没办法优化BL和BR的初始化，只能先全部算出来
    for i in 1:NN-1
        mul!(tmpG.nN,view(BLMs1,:,:,NN-i+1),view(BMs1,:,:,NN-i))
        LAPACK.gerqf!(tmpG.nN, tau)
        LAPACK.orgrq!(tmpG.nN, tau, ns)
        copyto!(view(BLMs1,:,:,NN-i) , tmpG.nN)
        
        mul!(tmpG.Nn, view(BMs1,:,:,i), view(BRMs1,:,:,i))
        LAPACK.geqrf!(tmpG.Nn, tau)
        LAPACK.orgqr!(tmpG.Nn, tau, ns)
        copyto!(view(BRMs1,:,:,i+1) , tmpG.Nn)
        # ---------------------------------------------------------------
        mul!(tmpG.nN,view(BLMs2,:,:,NN-i+1),view(BMs2,:,:,NN-i))
        LAPACK.gerqf!(tmpG.nN, tau)
        LAPACK.orgrq!(tmpG.nN, tau, ns)
        copyto!(view(BLMs2,:,:,NN-i) , tmpG.nN)

        mul!(tmpG.Nn, view(BMs2,:,:,i), view(BRMs2,:,:,i))
        LAPACK.geqrf!(tmpG.Nn, tau)
        LAPACK.orgqr!(tmpG.Nn, tau, ns)
        copyto!(view(BRMs2,:,:,i+1) , tmpG.Nn)

    end

    G4!(tmpG,G1,model.nodes,1,BLMs1,BRMs1,BMs1,BMsinv1)
    G4!(tmpG,G2,model.nodes,1,BLMs2,BRMs2,BMs2,BMsinv2)
    GroverMatrix!(GMA.gmInv,view(G1.O,indexA,indexA),view(G2.O,indexA,indexA))
    detg_A=abs(det(GMA.gmInv))
    LAPACK.getrf!(GMA.gmInv,GMA.ipiv)
    LAPACK.getri!(GMA.gmInv, GMA.ipiv)
    GroverMatrix!(GMB.gmInv,view(G1.O,indexB,indexB),view(G2.O,indexB,indexB))
    detg_B=abs(det(GMB.gmInv))
    LAPACK.getrf!(GMB.gmInv,GMB.ipiv)
    LAPACK.getri!(GMB.gmInv, GMB.ipiv)
    idx=1
    for loop in 1:Sweeps
        # println("\n ====== Sweep $loop / $Sweeps ======")
        for lt in 1:model.Nt
            #####################################################################
            # println("\n WrapTime check at lt=$lt")
            if lt-1 != div(model.Nt,2)
                Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2))
                Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2))
                    
                if norm(G1.t-Gt1_)+norm(G2.t-Gt2_)+norm(G1.tO-Gt01_)+norm(G2.tO-Gt02_)+norm(G1.Ot-G0t1_)+norm(G2.Ot-G0t2_)>1e-3
                    println( norm(G1.t-Gt1_),' ',norm(G2.t-Gt2_),'\n',norm(G1.O-G1.O_),' ',norm(G2.O-G2.O_),'\n',norm(G1.tO-G1.tO_),'\n',norm(G2.tO-G2.tO_),'\n',norm(G1.Ot-G1.Ot_),'\n',norm(G2.Ot-G2.Ot_) )
                    error("$lt : WrapTime")
                end
            end
            #####################################################################

            WrapK!(tmpG.NN,G1,model.eK,model.eKinv)
            WrapK!(tmpG.NN,G2,model.eK,model.eKinv)
            
            for j in 3:-1:1
                for i in axes(ss[1],2)
                    x,y=model.nnidx[i,j]
                    tmpG.N[x]=model.η[ss[1][lt,i,j]]
                    tmpG.N[y]=-model.η[ss[1][lt,i,j]]
                    tmpG.N_[x]=model.η[ss[2][lt,i,j]]
                    tmpG.N_[y]=-model.η[ss[2][lt,i,j]]
                end
                tmpG.N.= exp.(model.α.*tmpG.N)
                tmpG.N_.= exp.(model.α.*tmpG.N_)
                
                WrapV!(tmpG,G1.tO,view(model.UV,:,:,j),"L")
                WrapV!(tmpG,G1.t,view(model.UV,:,:,j),"B")
                WrapV!(tmpG,G1.Ot,view(model.UV,:,:,j),"R")

                copyto!(tmpG.N, tmpG.N_)
                WrapV!(tmpG,G2.tO,view(model.UV,:,:,j),"L")
                WrapV!(tmpG,G2.t,view(model.UV,:,:,j),"B")
                WrapV!(tmpG,G2.Ot,view(model.UV,:,:,j),"R")

                # update
                SCEE_LayerUpdate(tmpG,j,rng,model,G1,G2,GMA,GMB,view(ss[1],lt,:,j),view(ss[2],lt,:,j),λ)
            end

            ##------------------------------------------------------------------------
            tmpO+=(detg_A/detg_B)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
            
            if  any(model.nodes .== lt) 
                idx+=1
                BM_F!(tmpG,view(BMs1,:,:,idx-1),model,ss[1],idx-1)
                BMinv_F!(tmpG,view(BMsinv1,:,:,idx-1),model,ss[1],idx-1)
                BM_F!(tmpG,view(BMs2,:,:,idx-1),model,ss[2],idx-1)
                BMinv_F!(tmpG,view(BMsinv2,:,:,idx-1),model,ss[2],idx-1)
                for i in idx:max(Θidx,idx)
                    # println("update BR i=",i)
                    mul!(tmpG.Nn, view(BMs1,:,:,i-1), view(BRMs1,:,:,i-1))
                    LAPACK.geqrf!(tmpG.Nn,tau)
                    LAPACK.orgqr!(tmpG.Nn, tau, ns)
                    copyto!(view(BRMs1,:,:,i) , tmpG.Nn)
                    # ---------------------------------------------------------------
                    mul!(tmpG.Nn, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpG.Nn,tau)
                    LAPACK.orgqr!(tmpG.Nn, tau, ns)
                    copyto!(view(BRMs2,:,:,i) , tmpG.Nn)
                end

                for i in idx-1:-1:min(Θidx,idx)
                    # println("update BL i=",i)
                    mul!(tmpG.nN,view(BLMs1,:,:,i+1),view(BMs1,:,:,i))
                    LAPACK.gerqf!(tmpG.nN,tau)
                    LAPACK.orgrq!(tmpG.nN, tau, ns)
                    copyto!(view(BLMs1,:,:,i) , tmpG.nN)
                    # ---------------------------------------------------------------
                    mul!(tmpG.nN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    LAPACK.gerqf!(tmpG.nN,tau)
                    LAPACK.orgrq!(tmpG.nN, tau, ns)
                    copyto!(view(BLMs2,:,:,i) , tmpG.nN)
                end
                G4!(tmpG,G1,model.nodes,idx,BLMs1,BRMs1,BMs1,BMsinv1,"Forward")
                G4!(tmpG,G2,model.nodes,idx,BLMs2,BRMs2,BMs2,BMsinv2,"Forward")
                GroverMatrix!(GMA.gmInv,view(G1.O,indexA,indexA),view(G2.O,indexA,indexA))
                detg_A=abs(det(GMA.gmInv))
                LAPACK.getrf!(GMA.gmInv,GMA.ipiv)
                LAPACK.getri!(GMA.gmInv, GMA.ipiv)
                GroverMatrix!(GMB.gmInv,view(G1.O,indexB,indexB),view(G2.O,indexB,indexB))
                detg_B=abs(det(GMB.gmInv))
                LAPACK.getrf!(GMB.gmInv,GMB.ipiv)
                LAPACK.getri!(GMB.gmInv, GMB.ipiv)
            end
        end
        
        # println("\n #####################################################################")

        for lt in model.Nt:-1:1
            
            #####################################################################
            if lt-1 != div(model.Nt,2)
                Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt,div(model.Nt,2))
                Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt,div(model.Nt,2))
                    
                if norm(G1.t-Gt1_)+norm(G2.t-Gt2_)+norm(G1.tO-Gt01_)+norm(G2.tO-Gt02_)+norm(G1.Ot-G0t1_)+norm(G2.Ot-G0t2_)>1e-3
                    println( norm(G1.t-Gt1_),' ',norm(G2.t-Gt2_),'\n',norm(G1.O-G1.O_),' ',norm(G2.O-G2.O_),'\n',norm(G1.tO-G1.tO_),'\n',norm(G2.tO-G2.tO_),'\n',norm(G1.Ot-G1.Ot_),'\n',norm(G2.Ot-G2.Ot_) )
                    error("$lt : WrapTime")
                end
            end
            #####################################################################

            for j in 1:3
                # update
                SCEE_LayerUpdate(tmpG,j,rng,model,G1,G2,GMA,GMB,view(ss[1],lt,:,j),view(ss[2],lt,:,j),λ)

                for i in axes(ss[1],2)
                    x,y=model.nnidx[i,j]
                    tmpG.N[x]=model.η[ss[1][lt,i,j]]
                    tmpG.N[y]=-model.η[ss[1][lt,i,j]]
                    tmpG.N_[x]=model.η[ss[2][lt,i,j]]
                    tmpG.N_[y]=-model.η[ss[2][lt,i,j]]
                end
                tmpG.N.= exp.(.-model.α.*tmpG.N)
                tmpG.N_.= exp.(.-model.α.*tmpG.N_)

                WrapV!(tmpG,G1.tO,view(model.UV,:,:,j),1)
                WrapV!(tmpG,G1.t,view(model.UV,:,:,j),3)
                WrapV!(tmpG,G1.Ot,view(model.UV,:,:,j),2)

                copyto!(tmpG.N, tmpG.N_)

                WrapV!(tmpG,G2.tO,view(model.UV,:,:,j),1)
                WrapV!(tmpG,G2.t,view(model.UV,:,:,j),3)
                WrapV!(tmpG,G2.Ot,view(model.UV,:,:,j),2)
            end

            ##------------------------------------------------------------------------
            tmpO+=(detg_A/detg_B)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
            if  any(model.nodes.== (lt-1)) 
                idx-=1
                BM_F!(tmpG,view(BMs1,:,:,idx),model,ss[1],idx)
                BM_F!(tmpG,view(BMs2,:,:,idx),model,ss[2],idx)
                BMinv_F!(tmpG,view(BMsinv1,:,:,idx),model,ss[1],idx)
                BMinv_F!(tmpG,view(BMsinv2,:,:,idx),model,ss[2],idx)
                for i in idx:-1:min(Θidx,idx)
                    # println("update BL i=",i)
                    mul!(tmpG.nN,view(BLMs1,:,:,i+1),view(BMs1,:,:,i))
                    LAPACK.gerqf!(tmpG.nN,tau)
                    LAPACK.orgrq!(tmpG.nN, tau, ns)
                    copyto!(view(BLMs1,:,:,i) , tmpG.nN)

                    mul!(tmpG.nN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    LAPACK.gerqf!(tmpG.nN,tau)
                    LAPACK.orgrq!(tmpG.nN, tau, ns)
                    copyto!(view(BLMs2,:,:,i) , tmpG.nN)
                end
                for i in idx+1:max(Θidx,idx)
                    # println("update BR i=",i)
                    mul!(tmpG.nN, view(BMs1,:,:,i-1), view(BRMs1,:,:,i-1))
                    LAPACK.geqrf!(tmpG.nN,tau)
                    LAPACK.orgqr!(tmpG.nN, tau, ns)
                    copyto!(view(BRMs1,:,:,i) , tmpG.nN)

                    mul!(tmpG.nN, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpG.nN,tau)
                    LAPACK.orgqr!(tmpG.nN, tau, ns)
                    copyto!(view(BRMs2,:,:,i) , tmpG.nN)
                end
                G4!(tmpG,G1,model.nodes,idx,BLMs1,BRMs1,BMs1,BMsinv1,"Backward")
                G4!(tmpG,G2,model.nodes,idx,BLMs2,BRMs2,BMs2,BMsinv2,"Backward")
                GroverMatrix!(GMA.gmInv,view(G1.O,indexA,indexA),view(G2.O,indexA,indexA))
                detg_A=abs(det(GMA.gmInv))
                LAPACK.getrf!(GMA.gmInv,GMA.ipiv)
                LAPACK.getri!(GMA.gmInv, GMA.ipiv)
                GroverMatrix!(GMB.gmInv,view(G1.O,indexB,indexB),view(G2.O,indexB,indexB))
                detg_B=abs(det(GMB.gmInv))
                LAPACK.getrf!(GMB.gmInv,GMB.ipiv)
                LAPACK.getri!(GMB.gmInv, GMB.ipiv)
            else
                WrapK!(tmpG.NN,G1,model.eKinv,model.eK)
                WrapK!(tmpG.NN,G2,model.eKinv,model.eK)
            end
        end


        O[loop+1]=tmpO/counter
        tmpO=counter=0
    end
    return ss
end

"""
    Update G4 by overwriting Gt,G0,Gt0,G0t with r and subidx.
        G0 = G0 + G0t(:,subidx) ⋅ r ⋅ Gt0[subidx,:]
        Gt0 = Gt0 + r ⋅ Gt[subidx,:] ⋅ G0t(:,subidx)

        G0t = G0t - G0t(:,subidx) ⋅ r ⋅ (I-Gt[subidx,:])
        Gt = Gt - r ⋅ Gt[subidx,:] ⋅ r ⋅ (I-Gt[subidx,:])
    with r ≡ inv(r) ⋅ Δ
    ------------------------------------------------------------------------------
"""
function G4update!(tmpG,G,subidx::Vector{Int64})

    mul!(tmpG.zN,r,view(G.tO,subidx,:))   # useful for GΘ,GτΘ
    mul!(tmpG.NN, view(G.Ot,:,subidx),tmpG.zN)
    axpy!(1, tmpG.NN, G.O)
    mul!(tmpG.NN, view(G.t,:,subidx),tmpG.zN)
    axpy!(1, tmpG.NN, G.t0)

    mul!(tmpG.zN,r,view(G.t,subidx,:))
    lmul!(-1.0,tmpG.zN)
    axpy!(1.0,r,view(tmpG.zN,:,subidx))   # useful for GΘτ,Gτ
    mul!(tmpG.NN, view(G.Ot,:,subidx),tmpG.zN)
    axpy!(-1.0, tmpG.NN, G.Ot)
    mul!(tmpG.NN, view(G.t,:,subidx),tmpG.zN)
    axpy!(-1.0, tmpG.NN, G.t)
end

"""
    Update gmInv with a,b,Tau.
        gmInv = gmInv - gmInv ⋅ a ⋅ inv(Tau) ⋅ b
    Universal for s1 and s2.
    ------------------------------------------------------------------------------
"""
function GMupdate!(GM)
    mul!(GM.A2, GM.gmInv,GM.a )
    mul!(GM.zA,inv(GM.Tau),GM.b)
    mul!(GM.AA, GM.A2, GM.zA)
    axpy!(-1.0, GM.AA, GM.gmInv)
end

"""
    Return det(Tau)
    Update s1 and overwrite a , b , Tau.
        a = G0t1(:,subidx) ⋅ r
        b = Gt01[subidx,:] ⋅ (2G02-I) ⋅ gmInv
    with r ≡ inv(r) ⋅ Δ
    Warning : G02 here !!!  Gt01,G0t1
    ------------------------------------------------------------------------------
"""
function get_abTau1!(GM,G,G0,subidx::Vector{Int64})
    copyto!(GM.AA, view(G0,GM.idx,GM.idx))
    lmul!(2.0, GM.AA)
    for i in diagind(GM.AA)
        GM.AA[i] -= 1
    end
    mul!(GM.zA,view(G.tO,subidx,GM.idx),GM.AA)
    mul!(GM.b, GM.zA, GM.gmInv)
    mul!(GM.a,view(G.Ot,GM.idx,subidx),G.r)
    mul!(GM.Tau, GM.b, GM.a)
    GM.Tau[1,1]+=1; GM.Tau[2,2]+=1;
    return det(GM.Tau)
end

"""
    Return det(Tau)
    Update s2 and overwrite a , b , Tau.
        a = (2G1.O-I) ⋅ G2.Ot(:,subidx)
        b = r ⋅ G2.tO[subidx,:] ⋅ gmInv
    with r ≡ inv(r) ⋅ Δ
    Warning : G1.O here !!!  G2.tO,G2.Ot
    ------------------------------------------------------------------------------
"""
function get_abTau2!(GM,G,G0,subidx::Vector{Int64})
    copyto!(GM.AA , view(G0,GM.idx,GM.idx))
    lmul!(2.0, GM.AA)
    for i in diagind(GM.AA)
        GM.AA[i] -= 1
    end
    mul!(GM.a,view(G.Ot,GM.idx,subidx),G.r)
    mul!(GM.b,view(G.tO,subidx,GM.idx),GM.AA)
    mul!(GM.Tau, GM.b, GM.gmInv)
    GM.Tau[1,1]+=1; GM.Tau[2,2]+=1;
    return det(GM.Tau)
end

"""
    Overwrite G according to eK and eKinv , with option mid
        Forward Wrap :      WrapK!(tmpNN,Gt,Gt0,G0t,eK,eKinv)
            Gt = eK ⋅ Gt ⋅ eKinv
            Gt0 = eK ⋅ Gt0
            G0t = G0t ⋅ eKinv
        Inverse Wrap:       WrapK!(tmpNN,Gt,Gt0,G0t,eKinv,eK)    
            Gt = eKinv ⋅ Gt ⋅ eK
            Gt0 = eKinv ⋅ Gt0
            G0t = G0t ⋅ eK
    Only wrap Kinetic part forward direction  
    ------------------------------------------------------------------------------
"""
function WrapK!(tmpNN::Matrix{Float64},G::G4Workspace,eK::Matrix{Float64},eKinv::Matrix{Float64})
    mul!(tmpNN,G.t,eKinv)
    mul!(G.t,eK,tmpNN)
    mul!(tmpNN, eK, G.tO)
    copyto!(G.tO, tmpNN)
    mul!(tmpNN, G.Ot,eKinv)
    copyto!(G.Ot, tmpNN)
end

"""
    Overwrite G according to UV , D and option LR
        LR=1 : Only Left:   G = (UV * D * UV') * G
        LR=2 : Only Right   G = G * (UV * D * UV')'
        LR=3 : Both Side and D will be changed to 1/D !!!   G = (UV * D * UV') * G * (UV * inv(D) * UV')'
    Only wrap interaction part 
    ------------------------------------------------------------------------------
"""


function SCEE_LayerUpdate(tmpG,j,rng,model,G1,G2,GMA,GMB,s1,s2,λ)
    
    for i in eachindex(s1)
        x,y=model.nnidx[i,j]
        subidx=[x,y]
        # update ss[1]
        begin
            sx = rand(rng,  model.samplers_dict[s1[i]])
            p=get_r!(tmpG,model.α,model.η[sx]-model.η[s1[i]],subidx,G1.t)
            p*=model.γ[sx]/model.γ[s1[i]]

            detTau_A=get_abTau1!(GMA,G1,G2.O,subidx)
            detTau_B=get_abTau1!(GMB,G2,G1.O,subidx)

            @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
            if rand(rng)<p
                detg_A*=detTau_A
                detg_B*=detTau_B

                GMupdate!(GMA)
                GMupdate!(GMB)
                G4update!(tmpG,G1,subidx)

                s1[i]=sx
            end
        end

        # update ss[2]
        begin
            sx = rand(rng,  model.samplers_dict[s2[i]])
            p=get_r!(tmpG,model.α,model.η[sx]- model.η[s2[i]],subidx,G2.t)
            p*=model.γ[sx]/model.γ[s2[i]]
            
            detTau_A=get_abTau2!(GMA,G2,G1.O,subidx)
            detTau_B=get_abTau2!(GMB,G2,G1.O,subidx)
            @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
            if rand(rng)<p
                detg_A*=detTau_A
                detg_B*=detTau_B

                GMupdate!(GMA)
                GMupdate!(GMB)
                G4update!(tmpG,G2,subidx)
                s2[i]=sx
            end
        end
    end

    return s1,s2,detg_A,detg_B
end