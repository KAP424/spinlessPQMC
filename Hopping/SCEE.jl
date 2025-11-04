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
    tmpNn = Matrix{Float64}(undef, Ns, ns)
    tmpnN = Matrix{Float64}(undef, ns, Ns)
    tmpAA = Matrix{Float64}(undef ,length(indexA),length(indexA))
    tmpBB = Matrix{Float64}(undef ,length(indexB),length(indexB))

    tmp2N = Matrix{Float64}(undef, 2, Ns)
    tmpN2 = Matrix{Float64}(undef, Ns,2)
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
        # println("\n ====== Sweep $loop / $Sweeps ======")
        for lt in 1:model.Nt
            if  any(model.nodes.==(lt-1)) 
                idx= (lt==1) ? 2 : findfirst(model.nodes .== (lt-1))
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

                for i in idx:-1:min(Θidx,idx)-1
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
                idx=findfirst(model.nodes .== (lt-1))
                G4!(Gt1,G01,Gt01,G0t1,model.nodes,idx,BLMs1,BRMs1,BMs1,BMsinv1)
                G4!(Gt2,G02,Gt02,G0t2,model.nodes,idx,BLMs2,BRMs2,BMs2,BMsinv2)
                GroverMatrix!(gmInv_A,view(G01,indexA,indexA),view(G02,indexA,indexA))
                detg_A=abs(det(gmInv_A))
                LAPACK.getrf!(gmInv_A,ipivA)
                LAPACK.getri!(gmInv_A, ipivA)
                GroverMatrix!(gmInv_B,view(G01,indexB,indexB),view(G02,indexB,indexB))
                detg_B=abs(det(gmInv_B))
                LAPACK.getrf!(gmInv_B,ipivB)
                LAPACK.getri!(gmInv_B, ipivB)
            end

            #####################################################################
            # Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2))
            # Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2))
                
            # if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
            #     println( norm(Gt1-Gt1_),' ',norm(Gt2-Gt2_),'\n',norm(G01-G01_),' ',norm(G02-G02_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
            #     error("$lt : WrapTime")
            # end
            #####################################################################

            Forward_WrapK!(Gt1,G01,Gt01,G0t1,model.eK,model.eKinv,lt==div(model.Nt,2)+1)
            Forward_WrapK!(Gt2,G02,Gt02,G0t2,model.eK,model.eKinv,lt==div(model.Nt,2)+1)
            
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
                
                WrapV!(Gt01,tmpN,view(model.UV,:,:,j),1)
                WrapV!(Gt02,tmpN_,view(model.UV,:,:,j),1)

                WrapV!(Gt1,tmpN,view(model.UV,:,:,j),3)
                WrapV!(Gt2,tmpN_,view(model.UV,:,:,j),3)

                WrapV!(G0t1,tmpN,view(model.UV,:,:,j),2)
                WrapV!(G0t2,tmpN_,view(model.UV,:,:,j),2)

                # update
                for i in 1:ns
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]

                    # update ss[1]
                    begin
                        p=get_r!(r,model.α,ss[1][lt,i,j],subidx,Gt1)

                        get_abTau1(a_A,b_A,Tau_A,indexA,subidx,r,G02,Gt01,G0t1,gmInv_A)
                        get_abTau1(a_B,b_B,Tau_B,indexB,subidx,r,G02,Gt01,G0t1,gmInv_B)
                        detTau_A=det(Tau_A)
                        detTau_B=det(Tau_B)

                        # println("detTau_A: ",detTau_A," detTau_B: ",detTau_B," p:",p)
                        @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
                        if rand(rng)<p
                            detg_A*=detTau_A
                            detg_B*=detTau_B
        
                            GMupdate!(a_A,b_A,Tau_A,gmInv_A)
                            GMupdate!(a_B,b_B,Tau_B,gmInv_B)
                            G4update(subidx,r,Gt1,G01,Gt01,G0t1)

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
                        p=get_r!(r,model.α,ss[2][lt,i,j],subidx,Gt2)

                        get_abTau2(a_A,b_A,Tau_A,indexA,subidx,r,G01,Gt02,G0t2,gmInv_A)
                        get_abTau2(a_B,b_B,Tau_B,indexB,subidx,r,G01,Gt02,G0t2,gmInv_B)
                        detTau_A=det(Tau_A)
                        detTau_B=det(Tau_B)

                        @fastmath p*=detTau_A^λ * detTau_B^(1-λ)
                        if rand(rng)<p
                            detg_A*=detTau_A
                            detg_B*=detTau_B
        
                            GMupdate!(a_A,b_A,Tau_A,gmInv_A)
                            GMupdate!(a_B,b_B,Tau_B,gmInv_B)
                            G4update(subidx,r,Gt2,G02,Gt02,G0t2)
        
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
        end
        # println("\n #####################################################################")

        # for lt in model.Nt-1:-1:1
        #     if mod1(model.Nt-lt,model.WrapTime)==1 || lt==div(model.Nt,2)
        #         Gt1,G01,Gt01,G0t1=G4(model,ss[1],lt,div(model.Nt,2))
        #         Gt2,G02,Gt02,G0t2=G4(model,ss[2],lt,div(model.Nt,2))

        #         GM_A=GroverMatrix(G01[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
        #         gmInv_A=inv(GM_A)
        #         GM_B=GroverMatrix(G01[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
        #         gmInv_B=inv(GM_B)
        #         detg_A=det(GM_A)
        #         detg_B=det(GM_B)
        #     else
        #         #####################################################################
        #         # Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt,div(model.Nt,2))
        #         # Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt,div(model.Nt,2))
                    
        #         # if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
        #         #     println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
        #         #     error("$lt : WrapTime")
        #         # end
        #         #####################################################################
        #     end

        #     for j in 1:size(ss[1])[3]
        #         for i in 1:size(ss[1])[2]
        #             # update
        #             x,y=model.nnidx[i,j]
        #             subidx=[x,y]

        #             # 更新s1
        #             E=[-2*ss[1][lt,i,j] , 2*ss[1][lt,i,j]]
        #             Δ1=uv*diagm(exp.(model.α.*E))*uv'-I(2)
        #             r1=I(2)+Δ1*(I(2)-Gt1[subidx,subidx])

        #             b_A=(Gt01[subidx,indexA[:]]) *(2*G02[indexA[:],indexA[:]]-IA)*gmInv_A
        #             a_A=G0t1[indexA[:],subidx]
        #             Tau_A=b_A*a_A
                    
        #             b_B=(Gt01[subidx,indexB[:]]) *(2*G02[indexB[:],indexB[:]]-IB)*gmInv_B
        #             a_B=G0t1[indexB[:],subidx]
        #             Tau_B=b_B*a_B

        #             p=det(r1+Δ1*Tau_A)^λ*det(r1+Δ1*Tau_B)^(1-λ)

        #             if p<0
        #                 println("Negative Sign: $(p)")
        #             end

        #             if rand(rng)<p
        #                 rho_A=inv(r1+Δ1*Tau_A)*Δ1
        #                 gmInv_A-=gmInv_A* a_A * rho_A * b_A
        #                 detg_A*=det(I(2)+inv(r1)*Δ1*Tau_A)
    
        #                 rho_B=inv(r1+Δ1*Tau_B)*Δ1
        #                 gmInv_B-=gmInv_B* ( a_B * rho_B * b_B)
        #                 detg_B*=det(I(2)+inv(r1)*Δ1*Tau_B)

        #                 G01+= (G0t1[:,subidx] /r1*Δ1 *(Gt01[subidx,:]))
        #                 Gt01+=(Gt1[:,subidx] /r1*Δ1 *(Gt01[subidx,:]))
        #                 G0t1-=(G0t1[:,subidx] /r1*Δ1  * (II-Gt1)[subidx,:]  )
        #                 Gt1-= (Gt1[:,subidx] /r1*Δ1 *((II-Gt1)[subidx,:]) )         
        #                 ss[1][lt,i,j]=-ss[1][lt,i,j]
        #                 #####################################################################
        #                 # print('-')
        #                 # if lt==div(model.Nt,2)+1
        #                 #     Gt1_,G01_,G0t1_,Gt01_=G4(model,ss[1],lt-1,div(model.Nt,2))
        #                 # else
        #                 #     Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2))
        #                 # end
        #                 # Gt1_=model.eK*Gt1_*model.eKinv
        #                 # Gt01_=model.eK*Gt01_
        #                 # G0t1_=G0t1_*model.eKinv
        #                 # GM_A_=GroverMatrix(G01_[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
        #                 # gmInv_A_=inv(GM_A_)
        #                 # GM_B_=GroverMatrix(G01_[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
        #                 # gmInv_B_=inv(GM_B_)
        #                 # detg_A_=det(GM_A_)
        #                 # detg_B_=det(GM_B_)

        #                 # for jj in size(ss[1])[3]:-1:j
        #                 #     E=zeros(model.Ns)
        #                 #     for ii in 1:size(ss[1])[2]
        #                 #         x,y=model.nnidx[ii,jj]
        #                 #         E[x]=ss[1][lt,ii,jj]
        #                 #         E[y]=-ss[1][lt,ii,jj]
        #                 #     end
        #                 #     Gt1_=model.UV[:,:,jj]'*diagm(exp.(model.α*E))*model.UV[:,:,jj] *Gt1_* model.UV[:,:,jj]'*diagm(exp.(-model.α*E))*model.UV[:,:,jj]
        #                 #     Gt01_=model.UV[:,:,jj]'*diagm(exp.(model.α*E))*model.UV[:,:,jj]*Gt01_
        #                 #     G0t1_=G0t1_*model.UV[:,:,jj]'*diagm(exp.(-model.α*E))*model.UV[:,:,jj]
        #                 # end
    
        #                 # if norm(Gt1-Gt1_)+norm(G01-G01_)+norm(Gt01-Gt01_)+norm(G0t1-G0t1_)+
        #                 #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
        #                 #     println('\n',norm(Gt1-Gt1_),'\n',norm(G01-G01_),'\n',norm(Gt01-Gt01_),'\n',norm(G0t1-G0t1_))
        #                 #     println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
        #                 #     error("s1:  $lt  $j:,,,asdasdasd")
        #                 # end
        #                 #####################################################################
        #             end

        #             # 更新s2
        #             E=[-2*ss[2][lt,i,j] , 2*ss[2][lt,i,j]]
        #             Δ2=uv*diagm(exp.(model.α.*E))*uv'-I(2)
        #             r2=I(2)+Δ2*(I(2)-Gt2[subidx,subidx])

        #             b_A=Gt02[subidx,indexA[:]]*gmInv_A
        #             a_A=(2*G01[indexA[:],indexA[:]]-IA)*G0t2[indexA[:],subidx]
        #             Tau_A=b_A*a_A

        #             b_B=Gt02[subidx,indexB[:]]*gmInv_B
        #             a_B=(2*G01[indexB[:],indexB[:]]-IB)*G0t2[indexB[:],subidx]
        #             Tau_B=b_B*a_B
                    
        #             p=det(r2+Δ2*Tau_A)^λ*det(r2+Δ2*Tau_B)^(1-λ)

        #             if p<0
        #                 println("Negative Sign: $(p)")
        #             end

        #             if rand(rng)<p
        #                 rho_A=inv(r2+Δ2*Tau_A)*Δ2
        #                 gmInv_A-=gmInv_A* a_A * rho_A * b_A
        #                 detg_A*=det(I(2)+inv(r2)*Δ2*Tau_A) 
    
        #                 rho_B=inv(r2+Δ2*Tau_B)*Δ2
        #                 gmInv_B-=gmInv_B* a_B * rho_B * b_B
        #                 detg_B*=det(I(2)+inv(r2)*Δ2*Tau_B)

        #                 G02+= (G0t2[:,subidx] /r2*Δ2 *(Gt02[subidx,:]))
        #                 Gt02+=(Gt2[:,subidx] /r2*Δ2 *(Gt02[subidx,:]))
        #                 G0t2-=(G0t2[:,subidx] /r2*Δ2  * (II-Gt2)[subidx,:]  )
        #                 Gt2-= (Gt2[:,subidx] /r2*Δ2 *((II-Gt2)[subidx,:]) )         
        #                 ss[2][lt,i,j]=-ss[2][lt,i,j]
        #                 #####################################################################
        #                 # print('*')
        #                 # if lt==div(model.Nt,2)+1
        #                 #     Gt2_,G02_,G0t2_,Gt02_=G4(model,ss[2],lt-1,div(model.Nt,2))
        #                 # else
        #                 #     Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2))
        #                 # end
        #                 # Gt2_=model.eK*Gt2_*model.eKinv
        #                 # Gt02_=model.eK*Gt02_
        #                 # G0t2_=G0t2_*model.eKinv
        #                 # GM_A_=GroverMatrix(G01[indexA[:],indexA[:]],G02_[indexA[:],indexA[:]])
        #                 # gmInv_A_=inv(GM_A_)
        #                 # GM_B_=GroverMatrix(G01[indexB[:],indexB[:]],G02_[indexB[:],indexB[:]])
        #                 # gmInv_B_=inv(GM_B_)
        #                 # detg_A_=det(GM_A_)
        #                 # detg_B_=det(GM_B_)

        #                 # for jj in size(ss[1])[3]:-1:j
        #                 #     E=zeros(model.Ns)
        #                 #     for ii in 1:size(ss[1])[2]
        #                 #         x,y=model.nnidx[ii,jj]
        #                 #         E[x]=ss[2][lt,ii,jj]
        #                 #         E[y]=-ss[2][lt,ii,jj]
        #                 #     end
        #                 #     Gt2_=model.UV[:,:,jj]'*diagm(exp.(model.α*E))*model.UV[:,:,jj] *Gt2_* model.UV[:,:,jj]'*diagm(exp.(-model.α*E))*model.UV[:,:,jj]
        #                 #     Gt02_=model.UV[:,:,jj]'*diagm(exp.(model.α*E))*model.UV[:,:,jj]*Gt02_
        #                 #     G0t2_=G0t2_*model.UV[:,:,jj]'*diagm(exp.(-model.α*E))*model.UV[:,:,jj]
        #                 # end
    
        #                 # if norm(Gt2-Gt2_)+norm(G02-G02_)+norm(Gt02-Gt02_)+norm(G0t2-G0t2_)+
        #                 #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
        #                 #     println('\n',norm(Gt2-Gt2_),'\n',norm(G02-G02_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t2-G0t2_))
        #                 #     println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
        #                 #     error("s2:  $lt  $x:,,,asdasdasd")
        #                 # end
        #                 #####################################################################
        #             end

        #         end
        #         # 开始卸载eV eK
        #         E1=zeros(model.Ns)
        #         E2=zeros(model.Ns)
        #         for i in 1:size(ss[1])[2]
        #             x,y=model.nnidx[i,j]
        #             E1[x]=ss[1][lt,i,j]
        #             E1[y]=-ss[1][lt,i,j]
        #             E2[x]=ss[2][lt,i,j]
        #             E2[y]=-ss[2][lt,i,j]
        #         end

        #         Gt1=model.UV[j,:,:]'*diagm(exp.(-model.α*E1))*model.UV[j,:,:] *Gt1* model.UV[j,:,:]'*diagm(exp.(model.α*E1))*model.UV[j,:,:]
        #         Gt2=model.UV[j,:,:]'*diagm(exp.(-model.α*E2))*model.UV[j,:,:] *Gt2* model.UV[j,:,:]'*diagm(exp.(model.α*E2))*model.UV[j,:,:]
        #         Gt01=model.UV[j,:,:]'*diagm(exp.(-model.α*E1))*model.UV[j,:,:] *Gt01
        #         Gt02=model.UV[j,:,:]'*diagm(exp.(-model.α*E2))*model.UV[j,:,:] *Gt02
        #         G0t1=G0t1* model.UV[j,:,:]'*diagm(exp.(model.α*E1))*model.UV[j,:,:]
        #         G0t2=G0t2* model.UV[j,:,:]'*diagm(exp.(model.α*E2))*model.UV[j,:,:]
        #     end
        #     Gt1=model.eKinv*Gt1*model.eK
        #     Gt2=model.eKinv*Gt2*model.eK
        #     Gt01=model.eKinv*Gt01
        #     Gt02=model.eKinv*Gt02
        #     G0t1=G0t1*model.eK
        #     G0t2=G0t2*model.eK

        #     ##------------------------------------------------------------------------
        #     tmpO+=(detg_A/detg_B)^(1/Nλ)
        #     counter+=1
        #     ##------------------------------------------------------------------------
        # end

        O[loop+1]=tmpO/counter
        tmpO=counter=0
    end
    return ss
end
# function inv22!(A,B)
#     detB=det(B)
#     A[1,1]=B[2,2]/detB
#     A[1,2]=-B[1,2]/detB
#     A[2,1]=-B[2,1]/detB
#     A[2,2]=B[1,1]/detB
# end

function G4update(subidx::Vector{Int64},r::Matrix{Float64},Gt::Matrix{Float64},G0::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64})
    Ns=size(Gt,1)
    tmpNN = Matrix{Float64}(undef, Ns, Ns)
    tmp2N = Matrix{Float64}(undef, 2, Ns)

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

function GMupdate!(a::Matrix{Float64},b::Matrix{Float64},Tau::Matrix{Float64},gmInv::Matrix{Float64})
    NA=size(gmInv,1)
    tmpA2=Matrix{Float64}(undef, NA, 2)
    tmp2A=Matrix{Float64}(undef, 2, NA)
    tmpAA=Matrix{Float64}(undef, NA, NA)
    tmp22=Matrix{Float64}(undef, 2, 2)

    mul!(tmpA2, gmInv,a )
    inv22!(tmp22,Tau)
    mul!(tmp2A,tmp22,b)
    mul!(tmpAA, tmpA2, tmp2A)
    axpy!(-1.0, tmpAA, gmInv)

end

function get_abTau1(a::Matrix{Float64},b::Matrix{Float64},Tau::Matrix{Float64},index::Vector{Int64},subidx::Vector{Int64},r::Matrix{Float64},G0::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64},gmInv::Matrix{Float64})
    NA=length(index)
    tmpAA=Matrix{Float64}(undef, NA, NA)
    tmp2A=Matrix{Float64}(undef, 2, NA)

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
end

function get_abTau2(a::Matrix{Float64},b::Matrix{Float64},Tau::Matrix{Float64},index::Vector{Int64},subidx::Vector{Int64},r::Matrix{Float64},G0::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64},gmInv::Matrix{Float64})
    NA=length(index)
    tmpAA=Matrix{Float64}(undef, NA, NA)
    tmp2A=Matrix{Float64}(undef, 2, NA)

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
end

function Forward_WrapK!(Gt::Matrix{Float64},G0::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64},eK::Matrix{Float64},eKinv::Matrix{Float64},mid::Bool)
    tmpNN=similar(Gt)

    mul!(tmpNN,Gt,eKinv)
    mul!(Gt,eK,tmpNN)
    
    if mid
        # 顺序更新Θ+1层时候，得到的是Θ层的四Green函数，再将eK,eV往上安装同时更新
        # 而Θ层的四Green函数是和 1 ~ Θ 层同用一个Wrap规则
        # 所以为了让Θ层的四Green函数与Θ-2Θ同用一个Wrap规则，需要将Θ层的四Green函数的最后两个进行倒换
        tmpNN .= G0
        for i in diagind(tmpNN)
            tmpNN[i] -= 1
        end
        mul!(G0t, tmpNN ,eKinv)
        mul!(Gt0, eK, G0)
    else
        mul!(tmpNN, eK, Gt0)
        copyto!(Gt0, tmpNN)
        mul!(tmpNN, G0t,eKinv)
        copyto!(G0t, tmpNN)
        # G0t1=G0t1*model.eKinv*diagm(exp.(-1im*model.α.*D1))
        # Gt01=diagm(exp.(1im*model.α.*D1))*model.eK*Gt01
    end
end

function WrapV!(G::Matrix{Float64},D::Vector{Float64},UV::SubArray{Float64, 2, Array{Float64, 3}},LR::Int64)
    """
    LR=1 : Only Left
    LR=2 : Only Right
    LR=3 : Both Side and D will be changed to 1/D !!!
    """
    tmpNN=similar(G)

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

function get_r!(r::Matrix{Float64},α::Float64,s::Int8,subidx::Vector{Int64},Gt::Matrix{Float64})
    uv=[1.0 1.0; 1.0 -1.0]/sqrt(2)
    tmp2=Vector{Float64}(undef,2)
    Δ=Matrix{Float64}(undef, 2, 2)
    tmp22=Matrix{Float64}(undef, 2, 2)

    tmp2[1]=-2*s;   tmp2[2]=2*s;
    tmp2 .= exp.(α.*tmp2).-1
    mul!(r,uv,Diagonal(tmp2))
    mul!(Δ,r,uv')
    mul!(r,Δ,view(Gt,subidx,subidx))
    axpby!(1.0,Δ, -1.0, r)   # r = I + Δ ⋅ (I - Gt1[subidx,subidx])
    r[1,1]+=1; r[2,2]+=1;
    p=det(r)
    # redefine r=inv(r) ⋅ ̇Δ 
    inv22!(tmp22,r)
    mul!(r,tmp22,Δ)
    return p
end

function Inverse_WrapK!(Gt::Matrix{Float64},G0::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64},eK::Matrix{Float64},eKinv::Matrix{Float64},mid::Bool)
    tmpNN=similar(Gt)

    mul!(tmpNN,Gt,eK)
    mul!(Gt,eKinv,tmpNN)
    
    if mid
        # 顺序更新Θ+1层时候，得到的是Θ层的四Green函数，再将eK,eV往上安装同时更新
        # 而Θ层的四Green函数是和 1 ~ Θ 层同用一个Wrap规则
        # 所以为了让Θ层的四Green函数与Θ-2Θ同用一个Wrap规则，需要将Θ层的四Green函数的最后两个进行倒换
        tmpNN .= G0
        for i in diagind(tmpNN)
            tmpNN[i] -= 1
        end
        mul!(G0t, tmpNN ,eK)
        mul!(Gt0, eKinv, G0)
    else
        mul!(tmpNN, eKinv, Gt0)
        copyto!(Gt0, tmpNN)
        mul!(tmpNN, G0t,eK)
        copyto!(G0t, tmpNN)
        # G0t1=G0t1*model.eKinv*diagm(exp.(-1im*model.α.*D1))
        # Gt01=diagm(exp.(1im*model.α.*D1))*model.eK*Gt01
    end
end
