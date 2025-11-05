# For particle hole symmetric of t-V model
# attractive-U and repulsive-U get the same S_2

function ctrl_SCEEicr(path::String,model::_Hubbard_Para,indexA::Vector{Int64},indexB::Vector{Int64},Sweeps::Int64,λ::Float64,Nλ::Int64,ss::Vector{Array{Int8,3}},record)
    Ns=model.Ns
    ns=div(model.Ns, 2)
    NN=length(model.nodes)
    tau = Vector{Float64}(undef, ns)
    ipiv = Vector{LAPACK.BlasInt}(undef, ns)
    ipivA = Vector{LAPACK.BlasInt}(undef, length(indexA))
    ipivB = Vector{LAPACK.BlasInt}(undef, length(indexB))
    II=Diagonal(ones(Float64,Ns))

    Θidx=div(NN,2)+1

    uv=[-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]
    name = if model.Lattice=="SQUARE" "□" 
    elseif model.Lattice=="HoneyComb60" "HC" 
    elseif model.Lattice=="HoneyComb120" "HC120" 
    else error("Lattice: $(model.Lattice) is not allowed !") end  

    file="$(path)SCEEicr$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)N$(Nλ)BS$(model.BatchSize).csv"
    
    atexit() do
        if record
            open(file, "a") do io
                lock(io)
                writedlm(io, O', ',')
                unlock(io)
            end
        end
        writedlm("$(path)ss/SS$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)λ$(Int(round(Nλ*λ))).csv", [ss[1] ss[2]],",")
    end
    
    
    rng=MersenneTwister(Threads.threadid()+time_ns())

    Gt1= Matrix{Float64}(undef ,model.Ns, model.Ns)
    Gt2= Matrix{Float64}(undef ,model.Ns, model.Ns)
    G01= Matrix{Float64}(undef ,model.Ns, model.Ns)
    G02= Matrix{Float64}(undef ,model.Ns, model.Ns)
    Gt01= Matrix{Float64}(undef ,model.Ns, model.Ns)
    Gt02= Matrix{Float64}(undef ,model.Ns, model.Ns)
    G0t1= Matrix{Float64}(undef ,model.Ns, model.Ns)
    G0t2= Matrix{Float64}(undef ,model.Ns, model.Ns)
    gmInv_A=Matrix{Float64}(undef ,length(indexA),length(indexA))
    gmInv_B=Matrix{Float64}(undef ,length(indexB),length(indexB))
    detg_A=detg_B=0 

    b_A= Matrix{Float64}(undef ,2,length(indexA))
    a_A= Matrix{Float64}(undef ,length(indexA),2)
    Tau_A= Matrix{Float64}(undef ,2,2)
    b_B= Matrix{Float64}(undef ,2,length(indexB))
    a_B= Matrix{Float64}(undef ,length(indexB),2)
    Tau_B= Matrix{Float64}(undef ,2,2)

    # 预分配临时数组
    tmpN = Vector{Float64}(undef, Ns)
    tmpN_ = Vector{Float64}(undef, Ns)
    tmpNN = Matrix{Float64}(undef, Ns, Ns)
    tmpNN2 = Matrix{Float64}(undef, Ns, Ns)
    tmpNn = Matrix{Float64}(undef, Ns, ns)
    tmpnN = Matrix{Float64}(undef, ns, Ns)
    tmpnn = Matrix{Float64}(undef, ns, ns)
    tmpAA = Matrix{Float64}(undef ,length(indexA),length(indexA))
    tmpBB = Matrix{Float64}(undef ,length(indexB),length(indexB))

    tmp2N = Matrix{Float64}(undef, 2, Ns)
    tmp2A= Matrix{Float64}(undef ,2,length(indexA))
    tmp2B= Matrix{Float64}(undef ,2,length(indexB))
    tmpA2= Matrix{Float64}(undef ,length(indexA),2)
    tmpB2= Matrix{Float64}(undef ,length(indexB),2)
    tmp22 = Matrix{Float64}(undef, 2,2)
    tmp2 = Vector{Float64}(undef,2)

    r = Matrix{Float64}(undef, 2,2)
    Δ = Matrix{Float64}(undef, 2,2)

    tmpO=0.0
    counter=0
    O=zeros(Float64,Sweeps+1)
    O[1]=λ

    BMs1=Array{Float64}(undef,model.Ns,model.Ns,NN-1)  # Number_of_BM*Ns*Ns
    BMs2=Array{Float64}(undef,model.Ns,model.Ns,NN-1)  # Number_of_BM*Ns*Ns
    BMsinv1=Array{Float64}(undef,model.Ns,model.Ns,NN-1)  # Number_of_BM*Ns*Ns
    BMsinv2=Array{Float64}(undef,model.Ns,model.Ns,NN-1)  # Number_of_BM*Ns*Ns

    for idx in 1:NN-1
        BM_F!(view(BMs1,:, : , idx),model,ss[1],idx)
        BM_F!(view(BMs2,:,:,idx),model,ss[2],idx)
        BMinv_F!(view(BMsinv1,:,:,idx),model,ss[1],idx)
        BMinv_F!(view(BMsinv2,:,:,idx),model,ss[2],idx)
        @assert norm(view(BMs1,:,:,idx)*view(BMsinv1,:,:,idx)-I(model.Ns))<1e-8 "BM1 inv error at idx=$idx"
    end

    BLMs1=Array{Float64}(undef,ns,model.Ns,NN)
    BRMs1=Array{Float64}(undef,model.Ns,ns,NN)
    view(BLMs1,:,:,NN) .= model.Pt'
    view(BRMs1,:,:,1) .= model.Pt
    
    BLMs2=Array{Float64}(undef,ns,model.Ns,NN)
    BRMs2=Array{Float64}(undef,model.Ns,ns,NN)
    view(BLMs2,:,:,NN) .= model.Pt'
    view(BRMs2,:,:,1) .= model.Pt

    # 没办法优化BL和BR的初始化，只能先全部算出来
    for i in 1:NN-1
        mul!(tmpnN,view(BLMs1,:,:,NN-i+1),view(BMs1,:,:,NN-i))
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau, ns)
        view(BLMs1,:,:,NN-i) .= tmpnN
        
        mul!(tmpNn, view(BMs1,:,:,i), view(BRMs1,:,:,i))
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau, ns)
        view(BRMs1,:,:,i+1) .= tmpNn
        # ---------------------------------------------------------------
        mul!(tmpnN,view(BLMs2,:,:,NN-i+1),view(BMs2,:,:,NN-i))
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau, ns)
        view(BLMs2,:,:,NN-i) .= tmpnN

        mul!(tmpNn, view(BMs2,:,:,i), view(BRMs2,:,:,i))
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau, ns)
        view(BRMs2,:,:,i+1) .= tmpNn

    end

    for loop in 1:Sweeps
        println("\n ====== Sweep $loop / $Sweeps ======")
        G4!(II,tmpnn,tmpNn,tmpNN,tmpNN2,ipiv,Gt1,G01,Gt01,G0t1,model.nodes,1,BLMs1,BRMs1,BMs1,BMsinv1)
        G4!(II,tmpnn,tmpNn,tmpNN,tmpNN2,ipiv,Gt2,G02,Gt02,G0t2,model.nodes,1,BLMs2,BRMs2,BMs2,BMsinv2)
        GroverMatrix!(gmInv_A,view(G01,indexA,indexA),view(G02,indexA,indexA))
        detg_A=abs(det(gmInv_A))
        LAPACK.getrf!(gmInv_A,ipivA)
        LAPACK.getri!(gmInv_A, ipivA)
        GroverMatrix!(gmInv_B,view(G01,indexB,indexB),view(G02,indexB,indexB))
        detg_B=abs(det(gmInv_B))
        LAPACK.getrf!(gmInv_B,ipivB)
        LAPACK.getri!(gmInv_B, ipivB)
        idx=1

        for lt in 1:model.Nt
            #####################################################################
            # # println("\n WrapTime check at lt=$lt")
            # if lt-1 != div(model.Nt,2)
            #     Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2))
            #     Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2))
                    
            #     if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
            #         println( norm(Gt1-Gt1_),' ',norm(Gt2-Gt2_),'\n',norm(G01-G01_),' ',norm(G02-G02_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
            #         error("$lt : WrapTime")
            #     end
            # end
            #####################################################################

            WrapK!(tmpNN,Gt1,Gt01,G0t1,model.eK,model.eKinv)
            WrapK!(tmpNN,Gt2,Gt02,G0t2,model.eK,model.eKinv)
            
            for j in 3:-1:1
                for i in 1:ns
                    x,y=model.nnidx[i,j]
                    tmpN[x]=ss[1][lt,i,j]
                    tmpN[y]=-ss[1][lt,i,j]
                    tmpN_[x]=ss[2][lt,i,j]
                    tmpN_[y]=-ss[2][lt,i,j]
                end
                tmpN.= exp.(model.α.*tmpN)
                tmpN_.= exp.(model.α.*tmpN_)
                
                WrapV!(tmpNN,Gt01,tmpN,view(model.UV,:,:,j),1)
                WrapV!(tmpNN,Gt02,tmpN_,view(model.UV,:,:,j),1)

                WrapV!(tmpNN,Gt1,tmpN,view(model.UV,:,:,j),3)
                WrapV!(tmpNN,Gt2,tmpN_,view(model.UV,:,:,j),3)

                WrapV!(tmpNN,G0t1,tmpN,view(model.UV,:,:,j),2)
                WrapV!(tmpNN,G0t2,tmpN_,view(model.UV,:,:,j),2)

                # update
                for i in 1:ns
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]

                    # update ss[1]
                    begin
                        p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,ss[1][lt,i,j],subidx,Gt1)

                        detTau_A=get_abTau1!(tmpAA,tmp2A,a_A,b_A,Tau_A,indexA,subidx,r,G02,Gt01,G0t1,gmInv_A)
                        detTau_B=get_abTau1!(tmpBB,tmp2B,a_B,b_B,Tau_B,indexB,subidx,r,G02,Gt01,G0t1,gmInv_B)

                        # println("detTau_A: ",detTau_A," detTau_B: ",detTau_B," p:",p)
                        @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
                        if rand(rng)<p
                            detg_A*=detTau_A
                            detg_B*=detTau_B
        
                            GMupdate!(tmpA2,tmp2A,tmp22,tmpAA,a_A,b_A,Tau_A,gmInv_A)
                            GMupdate!(tmpB2,tmp2B,tmp22,tmpBB,a_B,b_B,Tau_B,gmInv_B)
                            G4update!(tmpNN,tmp2N,subidx,r,Gt1,G01,Gt01,G0t1)

                            ss[1][lt,i,j]=-ss[1][lt,i,j]
                            #####################################################################
                            # print('-')
                            # if lt==div(model.Nt,2)+1
                            #     Gt1_,G01_,G0t1_,Gt01_=G4(model,ss[1],lt-1,div(model.Nt,2))
                            # else
                            #     Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2))
                            # end
                            # Gt1_=model.eK*Gt1_*model.eKinv
                            # Gt01_=model.eK*Gt01_
                            # G0t1_=G0t1_*model.eKinv
                            # GM_A_=GroverMatrix(G01_[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
                            # gmInv_A_=inv(GM_A_)
                            # GM_B_=GroverMatrix(G01_[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
                            # gmInv_B_=inv(GM_B_)
                            # detg_A_=det(GM_A_)
                            # detg_B_=det(GM_B_)

                            # for jj in size(ss[1])[3]:-1:j
                            #     E=zeros(model.Ns)
                            #     for ii in 1:size(ss[1])[2]
                            #         x,y=model.nnidx[ii,jj]
                            #         E[x]=ss[1][lt,ii,jj]
                            #         E[y]=-ss[1][lt,ii,jj]
                            #     end
                            #     Gt1_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]' *Gt1_* model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                            #     Gt01_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]'*Gt01_
                            #     G0t1_=G0t1_*model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                            # end
        
                            # if norm(Gt1-Gt1_)+norm(G01-G01_)+norm(Gt01-Gt01_)+norm(G0t1-G0t1_)+
                            #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                            #     println('\n',norm(Gt1-Gt1_),'\n',norm(G01-G01_),'\n',norm(Gt01-Gt01_),'\n',norm(G0t1-G0t1_))
                            #     println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
                            #     error("s1:  $lt  $j:,,,asdasdasd")
                            # end
                            ####################################################################
                        end
                    end

                    # update ss[2]
                    begin
                        p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,ss[2][lt,i,j],subidx,Gt2)

                        detTau_A=get_abTau2!(tmpAA,tmp2A,a_A,b_A,Tau_A,indexA,subidx,r,G01,Gt02,G0t2,gmInv_A)
                        detTau_B=get_abTau2!(tmpBB,tmp2B,a_B,b_B,Tau_B,indexB,subidx,r,G01,Gt02,G0t2,gmInv_B)

                        @fastmath p*=detTau_A^λ * detTau_B^(1-λ)
                        if rand(rng)<p
                            detg_A*=detTau_A
                            detg_B*=detTau_B
        
                            GMupdate!(tmpA2,tmp2A,tmp22,tmpAA,a_A,b_A,Tau_A,gmInv_A)
                            GMupdate!(tmpB2,tmp2B,tmp22,tmpBB,a_B,b_B,Tau_B,gmInv_B)
                            G4update!(tmpNN,tmp2N,subidx,r,Gt2,G02,Gt02,G0t2)
        
                            ss[2][lt,i,j]=-ss[2][lt,i,j]
                            #####################################################################
                            # print('*')
                            # if lt==div(model.Nt,2)+1
                            #     Gt2_,G02_,G0t2_,Gt02_=G4(model,ss[2],lt-1,div(model.Nt,2))
                            # else
                            #     Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2))
                            # end
                            # Gt2_=model.eK*Gt2_*model.eKinv
                            # Gt02_=model.eK*Gt02_
                            # G0t2_=G0t2_*model.eKinv
                            # GM_A_=GroverMatrix(G01[indexA[:],indexA[:]],G02_[indexA[:],indexA[:]])
                            # gmInv_A_=inv(GM_A_)
                            # GM_B_=GroverMatrix(G01[indexB[:],indexB[:]],G02_[indexB[:],indexB[:]])
                            # gmInv_B_=inv(GM_B_)
                            # detg_A_=det(GM_A_)
                            # detg_B_=det(GM_B_)

                            # for jj in size(ss[1])[3]:-1:j
                            #     E=zeros(model.Ns)
                            #     for ii in 1:size(ss[1])[2]
                            #         x,y=model.nnidx[ii,jj]
                            #         E[x]=ss[2][lt,ii,jj]
                            #         E[y]=-ss[2][lt,ii,jj]
                            #     end
                            #     Gt2_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]' *Gt2_* model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                            #     Gt02_=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]'*Gt02_
                            #     G0t2_=G0t2_*model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                            # end
        
                            # if norm(Gt2-Gt2_)+norm(G02-G02_)+norm(Gt02-Gt02_)+norm(G0t2-G0t2_)+
                            #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                            #     println('\n',norm(Gt2-Gt2_),'\n',norm(G02-G02_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t2-G0t2_))
                            #     println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
                            #     error("s2:  $lt  $x:,,,asdasdasd")
                            # end
                            #####################################################################

                        end
                    end
                    
                end

            end
            ##------------------------------------------------------------------------
            tmpO+=(detg_A/detg_B)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
            if  any(model.nodes .== lt) 
                idx+=1
                BM_F!(view(BMs1,:,:,idx-1),model,ss[1],idx-1)
                BMinv_F!(view(BMsinv1,:,:,idx-1),model,ss[1],idx-1)
                BM_F!(view(BMs2,:,:,idx-1),model,ss[2],idx-1)
                BMinv_F!(view(BMsinv2,:,:,idx-1),model,ss[2],idx-1)
                for i in idx:max(Θidx,idx)
                    # println("update BR i=",i)
                    mul!(tmpNn, view(BMs1,:,:,i-1), view(BRMs1,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BRMs1,:,:,i) .= tmpNn
                    # ---------------------------------------------------------------
                    mul!(tmpNn, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BRMs2,:,:,i) .= tmpNn
                end

                for i in idx-1:-1:min(Θidx,idx)
                    # println("update BL i=",i)
                    mul!(tmpnN,view(BLMs1,:,:,i+1),view(BMs1,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau, ns)
                    view(BLMs1,:,:,i) .= tmpnN
                    # ---------------------------------------------------------------
                    mul!(tmpnN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau, ns)
                    view(BLMs2,:,:,i) .= tmpnN
                end
                G4!(II,tmpnn,tmpNn,tmpNN,tmpNN2,ipiv,Gt1,G01,Gt01,G0t1,model.nodes,idx,BLMs1,BRMs1,BMs1,BMsinv1,true)
                G4!(II,tmpnn,tmpNn,tmpNN,tmpNN2,ipiv,Gt2,G02,Gt02,G0t2,model.nodes,idx,BLMs2,BRMs2,BMs2,BMsinv2,true)
                GroverMatrix!(gmInv_A,view(G01,indexA,indexA),view(G02,indexA,indexA))
                detg_A=abs(det(gmInv_A))
                LAPACK.getrf!(gmInv_A,ipivA)
                LAPACK.getri!(gmInv_A, ipivA)
                GroverMatrix!(gmInv_B,view(G01,indexB,indexB),view(G02,indexB,indexB))
                detg_B=abs(det(gmInv_B))
                LAPACK.getrf!(gmInv_B,ipivB)
                LAPACK.getri!(gmInv_B, ipivB)
            end
        end
        
        # println("\n #####################################################################")

        for lt in model.Nt:-1:1
            
            #####################################################################
            # if lt-1 != div(model.Nt,2)
            #     Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt,div(model.Nt,2))
            #     Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt,div(model.Nt,2))
                    
            #     if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
            #         println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
            #         error("$lt : WrapTime")
            #     end
            # end
            #####################################################################

            for j in 1:3
                # update
                for i in 1:ns
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]

                    # update ss[1]
                    begin
                        p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,ss[1][lt,i,j],subidx,Gt1)

                        detTau_A=get_abTau1!(tmpAA,tmp2A,a_A,b_A,Tau_A,indexA,subidx,r,G02,Gt01,G0t1,gmInv_A)
                        detTau_B=get_abTau1!(tmpBB,tmp2B,a_B,b_B,Tau_B,indexB,subidx,r,G02,Gt01,G0t1,gmInv_B)

                        # println("detTau_A: ",detTau_A," detTau_B: ",detTau_B," p:",p)
                        @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
                        if rand(rng)<p
                            detg_A*=detTau_A
                            detg_B*=detTau_B
        
                            GMupdate!(tmpA2,tmp2A,tmp22,tmpAA,a_A,b_A,Tau_A,gmInv_A)
                            GMupdate!(tmpB2,tmp2B,tmp22,tmpBB,a_B,b_B,Tau_B,gmInv_B)
                            G4update!(tmpNN,tmp2N,subidx,r,Gt1,G01,Gt01,G0t1)

                            ss[1][lt,i,j]=-ss[1][lt,i,j]
                        end
                    end

                    # update ss[2]
                    begin
                        p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,ss[2][lt,i,j],subidx,Gt2)

                        detTau_A=get_abTau2!(tmpAA,tmp2A,a_A,b_A,Tau_A,indexA,subidx,r,G01,Gt02,G0t2,gmInv_A)
                        detTau_B=get_abTau2!(tmpBB,tmp2B,a_B,b_B,Tau_B,indexB,subidx,r,G01,Gt02,G0t2,gmInv_B)

                        @fastmath p*=detTau_A^λ * detTau_B^(1-λ)
                        if rand(rng)<p
                            detg_A*=detTau_A
                            detg_B*=detTau_B
        
                            GMupdate!(tmpA2,tmp2A,tmp22,tmpAA,a_A,b_A,Tau_A,gmInv_A)
                            GMupdate!(tmpB2,tmp2B,tmp22,tmpBB,a_B,b_B,Tau_B,gmInv_B)
                            G4update!(tmpNN,tmp2N,subidx,r,Gt2,G02,Gt02,G0t2)
        
                            ss[2][lt,i,j]=-ss[2][lt,i,j]
                        end
                    end
                end

                for i in 1:ns
                    x,y=model.nnidx[i,j]
                    tmpN[x]=ss[1][lt,i,j]
                    tmpN[y]=-ss[1][lt,i,j]
                    tmpN_[x]=ss[2][lt,i,j]
                    tmpN_[y]=-ss[2][lt,i,j]
                end
                tmpN.= exp.(.-model.α.*tmpN)
                tmpN_.= exp.(.-model.α.*tmpN_)

                WrapV!(tmpNN,Gt01,tmpN,view(model.UV,:,:,j),1)
                WrapV!(tmpNN,Gt02,tmpN_,view(model.UV,:,:,j),1)

                WrapV!(tmpNN,Gt1,tmpN,view(model.UV,:,:,j),3)
                WrapV!(tmpNN,Gt2,tmpN_,view(model.UV,:,:,j),3)

                WrapV!(tmpNN,G0t1,tmpN,view(model.UV,:,:,j),2)
                WrapV!(tmpNN,G0t2,tmpN_,view(model.UV,:,:,j),2)

            end

            WrapK!(tmpNN,Gt1,Gt01,G0t1,model.eKinv,model.eK)
            WrapK!(tmpNN,Gt2,Gt02,G0t2,model.eKinv,model.eK)

            ##------------------------------------------------------------------------
            tmpO+=(detg_A/detg_B)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
            if  any(model.nodes.== (lt-1)) 
                idx-=1
                BM_F!(view(BMs1,:,:,idx),model,ss[1],idx)
                BM_F!(view(BMs2,:,:,idx),model,ss[2],idx)
                BMinv_F!(view(BMsinv1,:,:,idx),model,ss[1],idx)
                BMinv_F!(view(BMsinv2,:,:,idx),model,ss[2],idx)
                for i in idx:-1:min(Θidx,idx)
                    # println("update BL i=",i)
                    mul!(tmpnN,view(BLMs1,:,:,i+1),view(BMs1,:,:,i))
                    tmpNn.=tmpnN'
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BLMs1,:,:,i) .= tmpNn'

                    mul!(tmpnN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    tmpNn.=tmpnN'
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BLMs2,:,:,i) .= tmpNn'
                end
                for i in idx+1:max(Θidx,idx)
                    # println("update BR i=",i)
                    mul!(tmpNn, view(BMs1,:,:,i-1), view(BRMs1,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BRMs1,:,:,i) .= tmpNn

                    mul!(tmpNn, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BRMs2,:,:,i) .= tmpNn
                end
                G4!(II,tmpnn,tmpNn,tmpNN,tmpNN2,ipiv,Gt1,G01,Gt01,G0t1,model.nodes,idx,BLMs1,BRMs1,BMs1,BMsinv1,false)
                G4!(II,tmpnn,tmpNn,tmpNN,tmpNN2,ipiv,Gt2,G02,Gt02,G0t2,model.nodes,idx,BLMs2,BRMs2,BMs2,BMsinv2,false)
                GroverMatrix!(gmInv_A,view(G01,indexA,indexA),view(G02,indexA,indexA))
                detg_A=abs(det(gmInv_A))
                LAPACK.getrf!(gmInv_A,ipivA)
                LAPACK.getri!(gmInv_A, ipivA)
                GroverMatrix!(gmInv_B,view(G01,indexB,indexB),view(G02,indexB,indexB))
                detg_B=abs(det(gmInv_B))
                LAPACK.getrf!(gmInv_B,ipivB)
                LAPACK.getri!(gmInv_B, ipivB)
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
        a = (2G02-I) ⋅ G0t1(:,subidx)
        b = r ⋅ Gt01[subidx,:] ⋅ gmInv
    with r ≡ inv(r) ⋅ Δ
    Warning : G02 here !!!  Gt01,G0t1
    ------------------------------------------------------------------------------
"""
function get_abTau1!(tmpAA::Matrix{Float64},tmp2A::Matrix{Float64},a::Matrix{Float64},b::Matrix{Float64},Tau::Matrix{Float64},index::Vector{Int64},subidx::Vector{Int64},r::Matrix{Float64},G0::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64},gmInv::Matrix{Float64})
    tmpAA .= view(G0,index,index)
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
        a = (2G01-I) ⋅ G0t2(:,subidx)
        b = r ⋅ Gt02[subidx,:] ⋅ gmInv
    with r ≡ inv(r) ⋅ Δ
    Warning : G01 here !!!  Gt02,G0t2
    ------------------------------------------------------------------------------
"""
function get_abTau2!(tmpAA::Matrix{Float64},tmp2A::Matrix{Float64},a::Matrix{Float64},b::Matrix{Float64},Tau::Matrix{Float64},index::Vector{Int64},subidx::Vector{Int64},r::Matrix{Float64},G0::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64},gmInv::Matrix{Float64})
    tmpAA .= view(G0,index,index)
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
function WrapK!(tmpNN::Matrix{Float64},Gt::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64},eK::Matrix{Float64},eKinv::Matrix{Float64})
    mul!(tmpNN,Gt,eKinv)
    mul!(Gt,eK,tmpNN)
    
    mul!(tmpNN, eK, Gt0)
    copyto!(Gt0, tmpNN)
    mul!(tmpNN, G0t,eKinv)
    copyto!(G0t, tmpNN)
end

"""
    Overwrite G according to UV , D and option LR
        LR=1 : Only Left:   G = (UV * D * UV') * G
        LR=2 : Only Right   G = G * (UV * D * UV')'
        LR=3 : Both Side and D will be changed to 1/D !!!   G = (UV * D * UV') * G * (UV * inv(D) * UV')'
    Only wrap interaction part 
    ------------------------------------------------------------------------------
"""
function WrapV!(tmpNN::Matrix{Float64},G::Matrix{Float64},D::Vector{Float64},UV::SubArray{Float64, 2, Array{Float64, 3}},LR::Int64)
    if LR==1
        mul!(tmpNN,UV',G)
        mul!(G,Diagonal(D),tmpNN)
        mul!(tmpNN,UV,G)
        copyto!(G, tmpNN)
    elseif LR==2
        mul!(tmpNN, G , UV)
        mul!(G, tmpNN , Diagonal(D))
        mul!(tmpNN, G , UV')
        copyto!(G, tmpNN)
    else
        mul!(tmpNN,UV',G)
        mul!(G,tmpNN,UV)
        mul!(tmpNN,Diagonal(D),G)
        D.= 1 ./D
        mul!(G,tmpNN,Diagonal(D))
        mul!(tmpNN,UV,G)
        mul!(G,tmpNN,UV')
    end
end


