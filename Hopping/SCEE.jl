# For particle hole symmetric of t-V model
# attractive-U and repulsive-U get the same S_2

mutable struct G4Workspace
    t::Matrix{Float64}
    O::Matrix{Float64}
    tO::Matrix{Float64}
    Ot::Matrix{Float64}
end

mutable struct GMWorkspace
    detg::Float64
    gmInv::Matrix{Float64}
    a::Matrix{Float64}
    b::Matrix{Float64}
    Tau::Matrix{Float64}
end

mutable struct tmpEEWorkspace
    NN::Matrix{Float64}
    2N::Matrix{Float64}
    N2::Matrix{Float64}
end

mutable struct tmpGMWorkspace
    # tmpA2,tmp2A,tmp22,tmpAA

    # tmpNN,tmp2N

    tmpN::Vector{Float64}
    tmpNN::Matrix{Float64}
    tmp2N::Matrix{Float64}
    tmpN2::Matrix{Float64}
    tmp22::Matrix{Float64}
    tmp2::Vector{Float64}
    r::Matrix{Float64}
    Δ::Matrix{Float64}
    tmpAA::Matrix{Float64}
    tmpBB::Matrix{Float64}
    tmpA2::Matrix{Float64}
    tmpB2::Matrix{Float64}
    tmp2A::Matrix{Float64}
    tmp2B::Matrix{Float64}
end

function ctrl_SCEEicr(path::String,model::_Hubbard_Para,indexA::Vector{Int64},indexB::Vector{Int64},Sweeps::Int64,λ::Float64,Nλ::Int64,ss::Vector{Array{UInt8,3}},record)
    global LOCK=ReentrantLock()
    Ns=model.Ns
    ns=div(Ns, 2)
    NN=length(model.nodes)
    tau = Vector{Float64}(undef, ns)
    ipiv = Vector{LAPACK.BlasInt}(undef, ns)
    ipivA = Vector{LAPACK.BlasInt}(undef, length(indexA))
    ipivB = Vector{LAPACK.BlasInt}(undef, length(indexB))
    

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
    elements = (1, 2, 3, 4)
    samplers_dict = Dict{UInt8, Random.Sampler}()
    for excluded in elements
        allowed = [i for i in elements if i != excluded]
        samplers_dict[excluded] = Random.Sampler(rng, allowed)
    end

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
        0.0,
        Matrix{Float64}(undef ,length(indexA),length(indexA)),
        Matrix{Float64}(undef ,length(indexA),2),
        Matrix{Float64}(undef ,2,length(indexA)),
        Matrix{Float64}(undef ,2,2)     )
    GMB=GMWorkspace(
        0.0,
        Matrix{Float64}(undef ,length(indexB),length(indexB)),
        Matrix{Float64}(undef ,length(indexB),2),
        Matrix{Float64}(undef ,2,length(indexB)),
        Matrix{Float64}(undef ,2,2)     )
    tmpG=tmpGWorkspace(
        Matrix{Float64}(undef, 2,2),
        Matrix{Float64}(undef, 2,2),
        Vector{Float64}(undef, Ns),
        Matrix{Float64}(undef, Ns, Ns),
        Matrix{Float64}(undef, Ns, ns),
        Matrix{Float64}(undef, ns, Ns),
        Matrix{Float64}(undef, ns, ns),
        Matrix{Float64}(undef, Ns, 2),
        Matrix{Float64}(undef, 2, Ns),
        Matrix{Float64}(undef, 2,2),
        Vector{Float64}(undef,2),
        Vector{LAPACK.BlasInt}(undef, ns),
        [-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]    )
    tmpA=tmpEEWorkspace(
        Matrix{Float64}(undef, length(indexA), length(indexA)),
        Matrix{Float64}(undef, 2, length(indexA)),
        Matrix{Float64}(undef, length(indexA), 2)    )
    tmpA=tmpEEWorkspace(
        Matrix{Float64}(undef, length(indexB), length(indexB)),
        Matrix{Float64}(undef, 2, length(indexB)),
        Matrix{Float64}(undef, length(indexB), 2)    
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
        # @assert norm(view(BMs1,:,:,idx)*view(BMsinv1,:,:,idx)-I(Ns))<1e-8 "BM1 inv error at idx=$idx"
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
    LAPACK.getrf!(GMA.gmInv,ipivA)
    LAPACK.getri!(GMA.gmInv, ipivA)
    GroverMatrix!(GMB.gmInv,view(G1.O,indexB,indexB),view(G2.O,indexB,indexB))
    detg_B=abs(det(GMB.gmInv))
    LAPACK.getrf!(GMB.gmInv,ipivB)
    LAPACK.getri!(GMB.gmInv, ipivB)
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
                for i in 1:Ns
                    x,y=model.nnidx[i,j]
                    tmpG.N[x]=model.η[ss[1][lt,i,j]]
                    tmpG.N[y]=-model.η[ss[1][lt,i,j]]
                    tmpG.N_[x]=model.η[ss[2][lt,i,j]]
                    tmpG.N_[y]=-model.η[ss[2][lt,i,j]]
                end
                tmpG.N.= exp.(model.α.*tmpG.N)
                tmpG.N_.= exp.(model.α.*tmpG.N_)
                
                WrapV!(tmpG,G1.tO,view(model.UV,:,:,j),1)
                WrapV!(tmpG,G2.tO,.N_,view(model.UV,:,:,j),1)

                WrapV!(tmpNN,G1.t,tmpN,view(model.UV,:,:,j),3)
                WrapV!(tmpNN,G2.t,tmpN_,view(model.UV,:,:,j),3)

                WrapV!(tmpNN,G1.Ot,tmpN,view(model.UV,:,:,j),2)
                WrapV!(tmpNN,G2.Ot,tmpN_,view(model.UV,:,:,j),2)

                # update
                for i in 1:Ns
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]

                    # # update ss[1]
                    # begin
                    #     sx = rand(rng,  samplers_dict[ss[1][lt,i,j]])
                    #     p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,model.η[sx]-model.η[ss[1][lt,i,j]],subidx,G1.t)
                    #     p*=model.γ[sx]/model.γ[ss[1][lt,i,j]]

                    #     detTau_A=get_abTau1!(tmpAA,tmp2A,a_A,b_A,Tau_A,indexA,subidx,r,G2.O,G1.tO,G1.Ot,gmInv_A)
                    #     detTau_B=get_abTau1!(tmpBB,tmp2B,a_B,b_B,Tau_B,indexB,subidx,r,G2.O,G1.tO,G1.Ot,gmInv_B)

                    #     @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
                    #     if rand(rng)<p
                    #         detg_A*=detTau_A
                    #         detg_B*=detTau_B
        
                    #         GMupdate!(tmpA2,tmp2A,tmp22,tmpAA,a_A,b_A,Tau_A,gmInv_A)
                    #         GMupdate!(tmpB2,tmp2B,tmp22,tmpBB,a_B,b_B,Tau_B,gmInv_B)
                    #         G4update!(tmpNN,tmp2N,subidx,r,G1.t,G1.O,G1.tO,G1.Ot)

                    #         ss[1][lt,i,j]=sx
                    #         #####################################################################
                    #         print('-')
                    #         if lt==div(model.Nt,2)+1
                    #             G1.t_,G1.O_,G1.Ot_,G1.tO_=G4(model,ss[1],lt-1,div(model.Nt,2))
                    #         else
                    #             G1.t_,G1.O_,G1.tO_,G1.Ot_=G4(model,ss[1],lt-1,div(model.Nt,2))
                    #         end
                    #         G1.t_=model.eK*G1.t_*model.eKinv
                    #         G1.tO_=model.eK*G1.tO_
                    #         G1.Ot_=G1.Ot_*model.eKinv
                    #         GM_A_=GroverMatrix(G1.O_[indexA[:],indexA[:]],G2.O[indexA[:],indexA[:]])
                    #         gmInv_A_=inv(GM_A_)
                    #         GM_B_=GroverMatrix(G1.O_[indexB[:],indexB[:]],G2.O[indexB[:],indexB[:]])
                    #         gmInv_B_=inv(GM_B_)
                    #         detg_A_=det(GM_A_)
                    #         detg_B_=det(GM_B_)

                    #         for jj in size(ss[1])[3]:-1:j
                    #             E=zeros(Ns)
                    #             for ii in 1:size(ss[1])[2]
                    #                 x,y=model.nnidx[ii,jj]
                    #                 E[x]=model.η[ss[1][lt,ii,jj]]
                    #                 E[y]=-model.η[ss[1][lt,ii,jj]]
                    #             end
                    #             G1.t_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]' *G1.t_* model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                    #             G1.tO_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]'*G1.tO_
                    #             G1.Ot_=G1.Ot_*model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                    #         end
        
                    #         if norm(G1.t-G1.t_)+norm(G1.O-G1.O_)+norm(G1.tO-G1.tO_)+norm(G1.Ot-G1.Ot_)+
                    #            norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                    #             println('\n',norm(G1.t-G1.t_),'\n',norm(G1.O-G1.O_),'\n',norm(G1.tO-G1.tO_),'\n',norm(G1.Ot-G1.Ot_))
                    #             println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
                    #             error("s1:  $lt  $j:,,,asdasdasd")
                    #         end
                    #         ####################################################################
                    #     end
                    # end

                    # # update ss[2]
                    # begin
                    #     sx = rand(rng,  samplers_dict[ss[2][lt,i,j]])
                    #     p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,model.η[sx]- model.η[ss[2][lt,i,j]],subidx,G2.t)
                    #     p*=model.γ[sx]/model.γ[ss[2][lt,i,j]]
                        
                    #     detTau_A=get_abTau2!(tmpAA,tmp2A,a_A,b_A,Tau_A,indexA,subidx,r,G1.O,G2.tO,G2.Ot,gmInv_A)
                    #     detTau_B=get_abTau2!(tmpBB,tmp2B,a_B,b_B,Tau_B,indexB,subidx,r,G1.O,G2.tO,G2.Ot,gmInv_B)

                    #     @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
                    #     if rand(rng)<p
                    #         detg_A*=detTau_A
                    #         detg_B*=detTau_B
        
                    #         GMupdate!(tmpA2,tmp2A,tmp22,tmpAA,a_A,b_A,Tau_A,gmInv_A)
                    #         GMupdate!(tmpB2,tmp2B,tmp22,tmpBB,a_B,b_B,Tau_B,gmInv_B)
                    #         G4update!(tmpNN,tmp2N,subidx,r,G2.t,G2.O,G2.tO,G2.Ot)
        
                    #         ss[2][lt,i,j]=sx
                    #         #####################################################################
                    #         # print('*')
                    #         # if lt==div(model.Nt,2)+1
                    #         #     G2.t_,G2.O_,G2.Ot_,G2.tO_=G4(model,ss[2],lt-1,div(model.Nt,2))
                    #         # else
                    #         #     G2.t_,G2.O_,G2.tO_,G2.Ot_=G4(model,ss[2],lt-1,div(model.Nt,2))
                    #         # end
                    #         # G2.t_=model.eK*G2.t_*model.eKinv
                    #         # G2.tO_=model.eK*G2.tO_
                    #         # G2.Ot_=G2.Ot_*model.eKinv
                    #         # GM_A_=GroverMatrix(G1.O[indexA[:],indexA[:]],G2.O_[indexA[:],indexA[:]])
                    #         # gmInv_A_=inv(GM_A_)
                    #         # GM_B_=GroverMatrix(G1.O[indexB[:],indexB[:]],G2.O_[indexB[:],indexB[:]])
                    #         # gmInv_B_=inv(GM_B_)
                    #         # detg_A_=det(GM_A_)
                    #         # detg_B_=det(GM_B_)

                    #         # for jj in size(ss[1])[3]:-1:j
                    #         #     E=zeros(Ns)
                    #         #     for ii in 1:size(ss[1])[2]
                    #         #         x,y=model.nnidx[ii,jj]
                    #         #         E[x]=model.η[ss[2][lt,ii,jj]]
                    #         #         E[y]=-model.η[ss[2][lt,ii,jj]]
                    #         #     end
                    #         #     G2.t_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]' *G2.t_* model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                    #         #     G2.tO_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]'*G2.tO_
                    #         #     G2.Ot_=G2.Ot_*model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                    #         # end
        
                    #         # if norm(G2.t-G2.t_)+norm(G2.O-G2.O_)+norm(G2.tO-G2.tO_)+norm(G2.Ot-G2.Ot_)+
                    #         #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                    #         #     println('\n',norm(G2.t-G2.t_),'\n',norm(G2.O-G2.O_),'\n',norm(G2.tO-G2.tO_),'\n',norm(G2.Ot-G2.Ot_))
                    #         #     println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
                    #         #     error("s2:  $lt  $x:,,,asdasdasd")
                    #         # end
                    #         #####################################################################

                    #     end
                    # end
                    
                end

            end

            ##------------------------------------------------------------------------
            tmpO+=(detg_A/detg_B)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
            
            if  any(model.nodes .== lt) 
                idx+=1
                BM_F!(tmpN,tmpNN,view(BMs1,:,:,idx-1),model,ss[1],idx-1)
                BMinv_F!(tmpN,tmpNN,view(BMsinv1,:,:,idx-1),model,ss[1],idx-1)
                BM_F!(tmpN,tmpNN,view(BMs2,:,:,idx-1),model,ss[2],idx-1)
                BMinv_F!(tmpN,tmpNN,view(BMsinv2,:,:,idx-1),model,ss[2],idx-1)
                for i in idx:max(Θidx,idx)
                    # println("update BR i=",i)
                    mul!(tmpNn, view(BMs1,:,:,i-1), view(BRMs1,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    copyto!(view(BRMs1,:,:,i) , tmpNn)
                    # ---------------------------------------------------------------
                    mul!(tmpNn, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    copyto!(view(BRMs2,:,:,i) , tmpNn)
                end

                for i in idx-1:-1:min(Θidx,idx)
                    # println("update BL i=",i)
                    mul!(tmpnN,view(BLMs1,:,:,i+1),view(BMs1,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau, ns)
                    copyto!(view(BLMs1,:,:,i) , tmpnN)
                    # ---------------------------------------------------------------
                    mul!(tmpnN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau, ns)
                    copyto!(view(BLMs2,:,:,i) , tmpnN)
                end
                G4!(II,tmpnn,tmpNn,tmpNN,tmpNN2,ipiv,G1.t,G1.O,G1.tO,G1.Ot,model.nodes,idx,BLMs1,BRMs1,BMs1,BMsinv1,"Forward")
                G4!(II,tmpnn,tmpNn,tmpNN,tmpNN2,ipiv,G2.t,G2.O,G2.tO,G2.Ot,model.nodes,idx,BLMs2,BRMs2,BMs2,BMsinv2,"Forward")
                GroverMatrix!(gmInv_A,view(G1.O,indexA,indexA),view(G2.O,indexA,indexA))
                detg_A=abs(det(gmInv_A))
                LAPACK.getrf!(gmInv_A,ipivA)
                LAPACK.getri!(gmInv_A, ipivA)
                GroverMatrix!(gmInv_B,view(G1.O,indexB,indexB),view(G2.O,indexB,indexB))
                detg_B=abs(det(gmInv_B))
                LAPACK.getrf!(gmInv_B,ipivB)
                LAPACK.getri!(gmInv_B, ipivB)
            end
        end
        
        # println("\n #####################################################################")

        for lt in model.Nt:-1:1
            
            #####################################################################
            # if lt-1 != div(model.Nt,2)
            #     G1.t_,G1.O_,G1.tO_,G1.Ot_=G4(model,ss[1],lt,div(model.Nt,2))
            #     G2.t_,G2.O_,G2.tO_,G2.Ot_=G4(model,ss[2],lt,div(model.Nt,2))
                    
            #     if norm(G1.t-G1.t_)+norm(G2.t-G2.t_)+norm(G1.tO-G1.tO_)+norm(G2.tO-G2.tO_)+norm(G1.Ot-G1.Ot_)+norm(G2.Ot-G2.Ot_)>1e-3
            #         println( norm(G1.t-G1.t_),'\n',norm(G2.t-G2.t_),'\n',norm(G1.tO-G1.tO_),'\n',norm(G2.tO-G2.tO_),'\n',norm(G1.Ot-G1.Ot_),'\n',norm(G2.Ot-G2.Ot_) )
            #         error("$lt : WrapTime")
            #     end
            # end
            #####################################################################

            for j in 1:3
                # update
                for i in 1:Ns
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]

                    # update ss[1]
                    begin
                        sx = rand(rng,  samplers_dict[ss[1][lt,i,j]])
                        p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,model.η[sx]-model.η[ss[1][lt,i,j]],subidx,G1.t)
                        p*=model.γ[sx]/model.γ[ss[1][lt,i,j]]

                        detTau_A=get_abTau1!(tmpAA,tmp2A,a_A,b_A,Tau_A,indexA,subidx,r,G2.O,G1.tO,G1.Ot,gmInv_A)
                        detTau_B=get_abTau1!(tmpBB,tmp2B,a_B,b_B,Tau_B,indexB,subidx,r,G2.O,G1.tO,G1.Ot,gmInv_B)

                        @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
                        if rand(rng)<p
                            detg_A*=detTau_A
                            detg_B*=detTau_B
        
                            GMupdate!(tmpA2,tmp2A,tmp22,tmpAA,a_A,b_A,Tau_A,gmInv_A)
                            GMupdate!(tmpB2,tmp2B,tmp22,tmpBB,a_B,b_B,Tau_B,gmInv_B)
                            G4update!(tmpNN,tmp2N,subidx,r,G1.t,G1.O,G1.tO,G1.Ot)

                            ss[1][lt,i,j]=sx
                            #####################################################################
                            print('-')
                            if lt==div(model.Nt,2)+1
                                G1.t_,G1.O_,G1.Ot_,G1.tO_=G4(model,ss[1],lt-1,div(model.Nt,2))
                            else
                                G1.t_,G1.O_,G1.tO_,G1.Ot_=G4(model,ss[1],lt-1,div(model.Nt,2))
                            end
                            G1.t_=model.eK*G1.t_*model.eKinv
                            G1.tO_=model.eK*G1.tO_
                            G1.Ot_=G1.Ot_*model.eKinv
                            GM_A_=GroverMatrix(G1.O_[indexA[:],indexA[:]],G2.O[indexA[:],indexA[:]])
                            gmInv_A_=inv(GM_A_)
                            GM_B_=GroverMatrix(G1.O_[indexB[:],indexB[:]],G2.O[indexB[:],indexB[:]])
                            gmInv_B_=inv(GM_B_)
                            detg_A_=det(GM_A_)
                            # detg_B_=det(GM_B_)

                            # for jj in size(ss[1])[3]:-1:j
                            #     E=zeros(Ns)
                            #     for ii in 1:size(ss[1])[2]
                            #         x,y=model.nnidx[ii,jj]
                            #         E[x]=model.η[ss[1][lt,ii,jj]]
                            #         E[y]=-model.η[ss[1][lt,ii,jj]]
                            #     end
                            #     G1.t_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]' *G1.t_* model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                            #     G1.tO_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]'*G1.tO_
                            #     G1.Ot_=G1.Ot_*model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                            # end
        
                            # if norm(G1.t-G1.t_)+norm(G1.O-G1.O_)+norm(G1.tO-G1.tO_)+norm(G1.Ot-G1.Ot_)+
                            #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                            #     println('\n',norm(G1.t-G1.t_),'\n',norm(G1.O-G1.O_),'\n',norm(G1.tO-G1.tO_),'\n',norm(G1.Ot-G1.Ot_))
                            #     println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
                            #     error("s1:  $lt  $j:,,,asdasdasd")
                            # end
                            ####################################################################
                        end
                    end

                    # update ss[2]
                    begin
                        sx = rand(rng,  samplers_dict[ss[2][lt,i,j]])
                        p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,model.η[sx]- model.η[ss[2][lt,i,j]],subidx,G2.t)
                        p*=model.γ[sx]/model.γ[ss[2][lt,i,j]]
                        
                        detTau_A=get_abTau2!(tmpAA,tmp2A,a_A,b_A,Tau_A,indexA,subidx,r,G1.O,G2.tO,G2.Ot,gmInv_A)
                        detTau_B=get_abTau2!(tmpBB,tmp2B,a_B,b_B,Tau_B,indexB,subidx,r,G1.O,G2.tO,G2.Ot,gmInv_B)

                        @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
                        if rand(rng)<p
                            detg_A*=detTau_A
                            detg_B*=detTau_B
        
                            GMupdate!(tmpA2,tmp2A,tmp22,tmpAA,a_A,b_A,Tau_A,gmInv_A)
                            GMupdate!(tmpB2,tmp2B,tmp22,tmpBB,a_B,b_B,Tau_B,gmInv_B)
                            G4update!(tmpNN,tmp2N,subidx,r,G2.t,G2.O,G2.tO,G2.Ot)
        
                            ss[2][lt,i,j]=sx
                            #####################################################################
                            # print('*')
                            # if lt==div(model.Nt,2)+1
                            #     G2.t_,G2.O_,G2.Ot_,G2.tO_=G4(model,ss[2],lt-1,div(model.Nt,2))
                            # else
                            #     G2.t_,G2.O_,G2.tO_,G2.Ot_=G4(model,ss[2],lt-1,div(model.Nt,2))
                            # end
                            # G2.t_=model.eK*G2.t_*model.eKinv
                            # G2.tO_=model.eK*G2.tO_
                            # G2.Ot_=G2.Ot_*model.eKinv
                            # GM_A_=GroverMatrix(G1.O[indexA[:],indexA[:]],G2.O_[indexA[:],indexA[:]])
                            # gmInv_A_=inv(GM_A_)
                            # GM_B_=GroverMatrix(G1.O[indexB[:],indexB[:]],G2.O_[indexB[:],indexB[:]])
                            # gmInv_B_=inv(GM_B_)
                            # detg_A_=det(GM_A_)
                            # detg_B_=det(GM_B_)

                            # for jj in size(ss[1])[3]:-1:j
                            #     E=zeros(Ns)
                            #     for ii in 1:size(ss[1])[2]
                            #         x,y=model.nnidx[ii,jj]
                            #         E[x]=model.η[ss[2][lt,ii,jj]]
                            #         E[y]=-model.η[ss[2][lt,ii,jj]]
                            #     end
                            #     G2.t_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]' *G2.t_* model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                            #     G2.tO_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]'*G2.tO_
                            #     G2.Ot_=G2.Ot_*model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                            # end
        
                            # if norm(G2.t-G2.t_)+norm(G2.O-G2.O_)+norm(G2.tO-G2.tO_)+norm(G2.Ot-G2.Ot_)+
                            #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                            #     println('\n',norm(G2.t-G2.t_),'\n',norm(G2.O-G2.O_),'\n',norm(G2.tO-G2.tO_),'\n',norm(G2.Ot-G2.Ot_))
                            #     println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
                            #     error("s2:  $lt  $x:,,,asdasdasd")
                            # end
                            #####################################################################

                        end
                    end

                end

                for i in 1:Ns
                    x,y=model.nnidx[i,j]
                    tmpN[x]=model.η[ss[1][lt,i,j]]
                    tmpN[y]=-model.η[ss[1][lt,i,j]]
                    tmpN_[x]=model.η[ss[2][lt,i,j]]
                    tmpN_[y]=-model.η[ss[2][lt,i,j]]
                end
                tmpN.= exp.(.-model.α.*tmpN)
                tmpN_.= exp.(.-model.α.*tmpN_)

                WrapV!(tmpNN,G1.tO,tmpN,view(model.UV,:,:,j),1)
                WrapV!(tmpNN,G2.tO,tmpN_,view(model.UV,:,:,j),1)

                WrapV!(tmpNN,G1.t,tmpN,view(model.UV,:,:,j),3)
                WrapV!(tmpNN,G2.t,tmpN_,view(model.UV,:,:,j),3)

                WrapV!(tmpNN,G1.Ot,tmpN,view(model.UV,:,:,j),2)
                WrapV!(tmpNN,G2.Ot,tmpN_,view(model.UV,:,:,j),2)

            end

            ##------------------------------------------------------------------------
            tmpO+=(detg_A/detg_B)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
            if  any(model.nodes.== (lt-1)) 
                idx-=1
                BM_F!(tmpN,tmpNN,view(BMs1,:,:,idx),model,ss[1],idx)
                BM_F!(tmpN,tmpNN,view(BMs2,:,:,idx),model,ss[2],idx)
                BMinv_F!(tmpN,tmpNN,view(BMsinv1,:,:,idx),model,ss[1],idx)
                BMinv_F!(tmpN,tmpNN,view(BMsinv2,:,:,idx),model,ss[2],idx)
                for i in idx:-1:min(Θidx,idx)
                    # println("update BL i=",i)
                    mul!(tmpnN,view(BLMs1,:,:,i+1),view(BMs1,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau, ns)
                    copyto!(view(BLMs1,:,:,i) , tmpnN)

                    mul!(tmpnN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau, ns)
                    copyto!(view(BLMs2,:,:,i) , tmpnN)
                end
                for i in idx+1:max(Θidx,idx)
                    # println("update BR i=",i)
                    mul!(tmpNn, view(BMs1,:,:,i-1), view(BRMs1,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    copyto!(view(BRMs1,:,:,i) , tmpNn)

                    mul!(tmpNn, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    copyto!(view(BRMs2,:,:,i) , tmpNn)
                end
                G4!(II,tmpnn,tmpNn,tmpNN,tmpNN2,ipiv,G1.t,G1.O,G1.tO,G1.Ot,model.nodes,idx,BLMs1,BRMs1,BMs1,BMsinv1,"Backward")
                G4!(II,tmpnn,tmpNn,tmpNN,tmpNN2,ipiv,G2.t,G2.O,G2.tO,G2.Ot,model.nodes,idx,BLMs2,BRMs2,BMs2,BMsinv2,"Backward")
                GroverMatrix!(gmInv_A,view(G1.O,indexA,indexA),view(G2.O,indexA,indexA))
                detg_A=abs(det(gmInv_A))
                LAPACK.getrf!(gmInv_A,ipivA)
                LAPACK.getri!(gmInv_A, ipivA)
                GroverMatrix!(gmInv_B,view(G1.O,indexB,indexB),view(G2.O,indexB,indexB))
                detg_B=abs(det(gmInv_B))
                LAPACK.getrf!(gmInv_B,ipivB)
                LAPACK.getri!(gmInv_B, ipivB)
            else
                WrapK!(tmpNN,G1.t,G1.tO,G1.Ot,model.eKinv,model.eK)
                WrapK!(tmpNN,G2.t,G2.tO,G2.Ot,model.eKinv,model.eK)
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
function G4update!(tmpNN::Matrix{Float64},tmp2N::Matrix{Float64},subidx::Vector{Int64},r::Matrix{Float64},Gt::Matrix{Float64},G0::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64})

    mul!(tmp2N,r,view(Gt0,subidx,:))   # useful for GΘ,GτΘ
    mul!(tmpNN, view(G0t,:,subidx),tmp2N)
    axpy!(1, tmpNN, G0)
    mul!(tmpNN, view(Gt,:,subidx),tmp2N)
    axpy!(1, tmpNN, Gt0)

    mul!(tmp2N,r,view(Gt,subidx,:))
    lmul!(-1.0,tmp2N)
    axpy!(1.0,r,view(tmp2N,:,subidx))   # useful for GΘτ,Gτ
    mul!(tmpNN, view(G0t,:,subidx),tmp2N)
    axpy!(-1.0, tmpNN, G0t)
    mul!(tmpNN, view(Gt,:,subidx),tmp2N)
    axpy!(-1.0, tmpNN, Gt)
end

"""
    Update gmInv with a,b,Tau.
        gmInv = gmInv - gmInv ⋅ a ⋅ inv(Tau) ⋅ b
    Universal for s1 and s2.
    ------------------------------------------------------------------------------
"""
function GMupdate!(tmpA2::Matrix{Float64},tmp2A::Matrix{Float64},tmp22::Matrix{Float64},tmpAA::Matrix{Float64},a::Matrix{Float64},b::Matrix{Float64},Tau::Matrix{Float64},gmInv::Matrix{Float64})
    mul!(tmpA2, gmInv,a )
    inv22!(tmp22,Tau)
    mul!(tmp2A,tmp22,b)
    mul!(tmpAA, tmpA2, tmp2A)
    axpy!(-1.0, tmpAA, gmInv)
end

"""
    Return det(Tau)
    Update s1 and overwrite a , b , Tau.
        a = (2G2.O-I) ⋅ G1.Ot(:,subidx)
        b = r ⋅ G1.tO[subidx,:] ⋅ gmInv
    with r ≡ inv(r) ⋅ Δ
    Warning : G2.O here !!!  G1.tO,G1.Ot
    ------------------------------------------------------------------------------
"""
function get_abTau1!(tmpA,GMA,G4A,index::Vector{Int64},subidx::Vector{Int64},r::Matrix{Float64})
    copyto!(tmpAA, view(G0,index,index))
    lmul!(2.0, tmpAA)
    for i in diagind(tmpAA)
        tmpAA[i] -= 1
    end
    mul!(tmp2A,view(Gt0,subidx,index),tmpAA)
    mul!(b, tmp2A, gmInv)
    mul!(a,view(G0t,index,subidx),r)
    mul!(Tau, b, a)
    Tau[1,1]+=1; Tau[2,2]+=1;
    return det(Tau)
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
function get_abTau2!(tmpAA::Matrix{Float64},tmp2A::Matrix{Float64},a::Matrix{Float64},b::Matrix{Float64},Tau::Matrix{Float64},index::Vector{Int64},subidx::Vector{Int64},r::Matrix{Float64},G0::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64},gmInv::Matrix{Float64})
    copyto!(tmpAA , view(G0,index,index))
    lmul!(2.0, tmpAA)
    for i in diagind(tmpAA)
        tmpAA[i] -= 1
    end
    mul!(a,tmpAA,view(G0t,index,subidx))
    mul!(tmp2A,r,view(Gt0,subidx,index))
    mul!(b,tmp2A,gmInv)
    mul!(Tau, b, a)
    Tau[1,1]+=1; Tau[2,2]+=1;
    return det(Tau)
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
function WrapK!(tmpNN::Matrix{Float64},G4::G4Workspace,eK::Matrix{Float64},eKinv::Matrix{Float64})
    mul!(tmpNN,G4.t,eKinv)
    mul!(G4.t,eK,tmpNN)
    
    mul!(tmpNN, eK, G4.t0)
    copyto!(G4.t0, tmpNN)
    mul!(tmpNN, G4.Ot,eKinv)
    copyto!(G4.Ot, tmpNN)
end

"""
    Overwrite G according to UV , D and option LR
        LR=1 : Only Left:   G = (UV * D * UV') * G
        LR=2 : Only Right   G = G * (UV * D * UV')'
        LR=3 : Both Side and D will be changed to 1/D !!!   G = (UV * D * UV') * G * (UV * inv(D) * UV')'
    Only wrap interaction part 
    ------------------------------------------------------------------------------
"""


function SCEE_LayerUpdate(j,rng,samplers_dict,model,s1,s2,G1.t,G1.O,G1.tO,G1.Ot,G2.t,G2.O,G2.tO,G2.Ot,gmInv_A,detg_A,gmInv_B,detg_B,λ)
    
    for i in 1:Ns
        x,y=model.nnidx[i,j]
        subidx=[x,y]
        # update ss[1]
        begin
            sx = rand(rng,  samplers_dict[s1[i]])
            p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,model.η[sx]-model.η[s1[i]],subidx,G1.t)
            p*=model.γ[sx]/model.γ[s1[i]]

            detTau_A=get_abTau1!(tmpAA,tmp2A,a_A,b_A,Tau_A,indexA,subidx,r,G2.O,G1.tO,G1.Ot,gmInv_A)
            detTau_B=get_abTau1!(tmpBB,tmp2B,a_B,b_B,Tau_B,indexB,subidx,r,G2.O,G1.tO,G1.Ot,gmInv_B)

            @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
            if rand(rng)<p
                detg_A*=detTau_A
                detg_B*=detTau_B

                GMupdate!(tmpA2,tmp2A,tmp22,tmpAA,a_A,b_A,Tau_A,gmInv_A)
                GMupdate!(tmpB2,tmp2B,tmp22,tmpBB,a_B,b_B,Tau_B,gmInv_B)
                G4update!(tmpNN,tmp2N,subidx,r,G1.t,G1.O,G1.tO,G1.Ot)

                s1[i]=sx
            end
        end

        # update ss[2]
        begin
            sx = rand(rng,  samplers_dict[s2[i]])
            p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,model.η[sx]- model.η[s2[i]],subidx,G2.t)
            p*=model.γ[sx]/model.γ[s2[i]]
            
            detTau_A=get_abTau2!(tmpAA,tmp2A,a_A,b_A,Tau_A,indexA,subidx,r,G1.O,G2.tO,G2.Ot,gmInv_A)
            detTau_B=get_abTau2!(tmpBB,tmp2B,a_B,b_B,Tau_B,indexB,subidx,r,G1.O,G2.tO,G2.Ot,gmInv_B)

            @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
            if rand(rng)<p
                detg_A*=detTau_A
                detg_B*=detTau_B

                GMupdate!(tmpA2,tmp2A,tmp22,tmpAA,a_A,b_A,Tau_A,gmInv_A)
                GMupdate!(tmpB2,tmp2B,tmp22,tmpBB,a_B,b_B,Tau_B,gmInv_B)
                G4update!(tmpNN,tmp2N,subidx,r,G2.t,G2.O,G2.tO,G2.Ot)

                s2[i]=sx
            end
        end
    end

    return s1,s2,detg_A,detg_B
end