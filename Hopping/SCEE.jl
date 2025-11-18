# For particle hole symmetric of t-V model
# attractive-U and repulsive-U get the same S_2

function ctrl_SCEEicr(path::String,model::Hubbard_Para_,indexA::Vector{Int64},indexB::Vector{Int64},Sweeps::Int64,λ::Float64,Nλ::Int64,ss::Vector{Array{UInt8,3}},record)
    global LOCK=ReentrantLock()
    ERROR=1e-6
    WrapErr = Matrix{Float64}(undef, model.Ns, model.Ns)

    NN=length(model.nodes)
    Θidx=div(NN,2)+1

    UPD = UpdateBuffer()
    SCEE=SCEEBuffer(model.Ns)
    A=AreaBuffer(indexA)
    B=AreaBuffer(indexB)
    G1=G4Buffer(model.Ns,NN)
    G2=G4Buffer(model.Ns,NN)

    name = if model.Lattice=="SQUARE" "□" 
    elseif model.Lattice=="HoneyComb60" "HC60" 
    elseif model.Lattice=="HoneyComb120" "HC120" 
    else error("Lattice: $(model.Lattice) is not allowed !") end  
    file="$(path)/SCEEicr$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)N$(Nλ)BS$(model.BatchSize).csv"
    rng=MersenneTwister(Threads.threadid()+time_ns())
    
    atexit() do
        if record
            lock(LOCK) do
                open(file, "a") do io
                    writedlm(io, O', ',')
                end
            end
        end
    end
    
    

    Gt1,G01,Gt01,G0t1,BLMs1,BRMs1,BMs1,BMsinv1 =
        G1.Gt, G1.G0, G1.Gt0, G1.G0t, G1.BLMs, G1.BRMs, G1.BMs, G1.BMinvs
    Gt2,G02,Gt02,G0t2,BLMs2,BRMs2,BMs2,BMsinv2 =
        G2.Gt, G2.G0, G2.Gt0, G2.G0t, G2.BLMs, G2.BRMs, G2.BMs, G2.BMinvs

    # 预分配临时数组
    tmpN = SCEE.N
    tmpN_ = SCEE.N_
    tmpNN = SCEE.NN
    tmpNn = SCEE.Nn
    tmpnN = SCEE.nN
    tau = SCEE.tau

    counter=0
    O=zeros(Float64,Sweeps+1)
    O[1]=λ
    tmpO=0.0

    for idx in 1:NN-1
        BM_F!(tmpN,tmpNN,view(BMs1,:, : , idx),model,ss[1],idx)
        BM_F!(tmpN,tmpNN,view(BMs2,:,:,idx),model,ss[2],idx)
        BMinv_F!(tmpN,tmpNN,view(BMsinv1,:,:,idx),model,ss[1],idx)
        BMinv_F!(tmpN,tmpNN,view(BMsinv2,:,:,idx),model,ss[2],idx)
        # @assert norm(view(BMs1,:,:,idx)*view(BMsinv1,:,:,idx)-I(Ns))<1e-8 "BM1 inv error at idx=$idx"
    end

    transpose!(view(BLMs1,:,:,NN) , model.Pt)
    copyto!(view(BRMs1,:,:,1) , model.Pt)
    
    transpose!(view(BLMs2,:,:,NN) , model.Pt)
    copyto!(view(BRMs2,:,:,1) , model.Pt)


    # 没办法优化BL和BR的初始化，只能先全部算出来
    for i in 1:NN-1
        mul!(tmpnN,view(BLMs1,:,:,NN-i+1),view(BMs1,:,:,NN-i))
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau)
        copyto!(view(BLMs1,:,:,NN-i) , tmpnN)
        
        mul!(tmpNn, view(BMs1,:,:,i), view(BRMs1,:,:,i))
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau)
        copyto!(view(BRMs1,:,:,i+1) , tmpNn)
        # ---------------------------------------------------------------
        mul!(tmpnN,view(BLMs2,:,:,NN-i+1),view(BMs2,:,:,NN-i))
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau)
        copyto!(view(BLMs2,:,:,NN-i) , tmpnN)

        mul!(tmpNn, view(BMs2,:,:,i), view(BRMs2,:,:,i))
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau)
        copyto!(view(BRMs2,:,:,i+1) , tmpNn)

    end

    
    idx=1
    get_ABGM!(G1,G2,A,B,SCEE,model.nodes,idx,"Forward")
    for loop in 1:Sweeps
        # println("\n ====== Sweep $loop / $Sweeps ======")
        for lt in 1:model.Nt
            #####################################################################
            println("\n WrapTime check at lt=$lt")
            Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2),"Forward")
            Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2),"Forward")
                
            if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>ERROR
                println( norm(Gt1-Gt1_),' ',norm(Gt2-Gt2_),'\n',norm(G01-G01_),' ',norm(G02-G02_),'\n',norm(Gt01-Gt01_),' ',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),' ',norm(G0t2-G0t2_) )
                error("$lt : WrapTime")
            end
            #####################################################################

            WrapK!(tmpNN,G1,model.eK,model.eKinv)
            WrapK!(tmpNN,G2,model.eK,model.eKinv)
            
            for j in reverse(axes(ss[1],2))
                for i in axes(ss[1],1)
                    x,y=model.nnidx[i,j]
                    tmpN[x]=model.η[ss[1][i,j,lt]]
                    tmpN[y]=-model.η[ss[1][i,j,lt]]
                    tmpN_[x]=model.η[ss[2][i,j,lt]]
                    tmpN_[y]=-model.η[ss[2][i,j,lt]]
                end
                tmpN.= exp.(tmpN)
                tmpN_.= exp.(tmpN_)
                
                WrapV!(tmpNN,Gt01,tmpN,view(model.UV,:,:,j),1)
                WrapV!(tmpNN,Gt02,tmpN_,view(model.UV,:,:,j),1)

                WrapV!(tmpNN,Gt1,tmpN,view(model.UV,:,:,j),3)
                WrapV!(tmpNN,Gt2,tmpN_,view(model.UV,:,:,j),3)

                WrapV!(tmpNN,G0t1,tmpN,view(model.UV,:,:,j),2)
                WrapV!(tmpNN,G0t2,tmpN_,view(model.UV,:,:,j),2)

                # update
                UpdateSCEELayer!(rng,j,view(ss[1],:,j,lt),view(ss[2],:,j,lt),G1,G2,A,B,model,UPD,SCEE,λ)
                #####################################################################
                    print('-')
                    Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2),"Forward")
                    Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2),"Forward")
                    Gt1_=model.eK*Gt1_*model.eKinv
                    Gt01_=model.eK*Gt01_
                    G0t1_=G0t1_*model.eKinv
                    Gt2_=model.eK*Gt2_*model.eKinv
                    Gt02_=model.eK*Gt02_
                    G0t2_=G0t2_*model.eKinv

                    GM_A_=GroverMatrix(G01_[indexA[:],indexA[:]],G02_[indexA[:],indexA[:]])
                    gmInv_A_=inv(GM_A_)
                    GM_B_=GroverMatrix(G01_[indexB[:],indexB[:]],G02_[indexB[:],indexB[:]])
                    gmInv_B_=inv(GM_B_)
                    detg_A_=det(GM_A_)
                    detg_B_=det(GM_B_)

                    for jj in 3:-1:j
                        E=zeros(model.Ns)
                        E_=zeros(model.Ns)
                        for ii in 1:size(ss[1])[1]
                            x,y=model.nnidx[ii,jj]
                            E[x]=model.η[ss[1][ii,jj,lt]]
                            E[y]=-model.η[ss[1][ii,jj,lt]]
                            E_[x]=model.η[ss[2][ii,jj,lt]]
                            E_[y]=-model.η[ss[2][ii,jj,lt]]
                        end
                        Gt1_=model.UV[:,:,jj]*Diagonal(exp.(E))*model.UV[:,:,jj] *Gt1_* model.UV[:,:,jj]*Diagonal(exp.(-E))*model.UV[:,:,jj]
                        Gt01_=model.UV[:,:,jj]*Diagonal(exp.(E))*model.UV[:,:,jj]*Gt01_
                        G0t1_=G0t1_*model.UV[:,:,jj]*Diagonal(exp.(-E))*model.UV[:,:,jj]
                        Gt2_=model.UV[:,:,jj]*Diagonal(exp.(E_))*model.UV[:,:,jj] *Gt2_* model.UV[:,:,jj]*Diagonal(exp.(-E_))*model.UV[:,:,jj]
                        Gt02_=model.UV[:,:,jj]*Diagonal(exp.(E_))*model.UV[:,:,jj]*Gt02_
                        G0t2_=G0t2_*model.UV[:,:,jj]*Diagonal(exp.(-E_))*model.UV[:,:,jj]
                    end

                    if norm(Gt1-Gt1_)+norm(G01-G01_)+norm(Gt01-Gt01_)+norm(G0t1-G0t1_)+
                        norm(Gt2-Gt2_)+norm(G02-G02_)+norm(Gt02-Gt02_)+norm(G0t2-G0t2_)+
                    norm(gmInv_A_-A.gmInv)+norm(B.gmInv-gmInv_B_)+abs(A.detg-detg_A_)+abs(B.detg-detg_B_)>ERROR

                        println('\n',norm(Gt1-Gt1_),'\n',norm(G01-G01_),'\n',norm(Gt01-Gt01_),'\n',norm(G0t1-G0t1_))
                        println('\n',norm(Gt2-Gt2_),'\n',norm(G02-G02_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t2-G0t2_))
                        println(norm(gmInv_A_-A.gmInv)," ",norm(B.gmInv-gmInv_B_)," ",abs(A.detg-detg_A_)," ",abs(B.detg-detg_B_))
                        error("s1:  $lt  $j:,,,asdasdasd")
                    end
                ######################################################################
            end

            ##------------------------------------------------------------------------
            tmpO+=(A.detg/B.detg)^(1/Nλ)
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
                    LAPACK.orgqr!(tmpNn, tau)
                    copyto!(view(BRMs1,:,:,i) , tmpNn)
                    # ---------------------------------------------------------------
                    mul!(tmpNn, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau)
                    copyto!(view(BRMs2,:,:,i) , tmpNn)
                end

                for i in idx-1:-1:min(Θidx,idx)
                    # println("update BL i=",i)
                    mul!(tmpnN,view(BLMs1,:,:,i+1),view(BMs1,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau)
                    copyto!(view(BLMs1,:,:,i) , tmpnN)
                    # ---------------------------------------------------------------
                    mul!(tmpnN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau)
                    copyto!(view(BLMs2,:,:,i) , tmpnN)
                end
                
                #####################################################################
                if lt != div(model.Nt,2)
                    copyto!(WrapErr, Gt1)
                    axpy!(1.0, Gt2, WrapErr)
                    axpy!(1.0, G01, WrapErr)
                    axpy!(1.0, G02, WrapErr)
                    axpy!(1.0, Gt01, WrapErr)
                    axpy!(1.0, Gt02, WrapErr)
                    axpy!(1.0, G0t1, WrapErr)
                    axpy!(1.0, G0t2, WrapErr)
                end
                #####################################################################
                get_ABGM!(G1,G2,A,B,SCEE,model.nodes,idx,"Forward")
                #####################################################################
                if lt != div(model.Nt,2)
                    axpy!(-1.0, Gt1, WrapErr)
                    axpy!(-1.0, Gt2, WrapErr)
                    axpy!(-1.0, G01, WrapErr)
                    axpy!(-1.0, G02, WrapErr)
                    axpy!(-1.0, Gt01, WrapErr)
                    axpy!(-1.0, Gt02, WrapErr)
                    axpy!(-1.0, G0t1, WrapErr)
                    axpy!(-1.0, G0t2, WrapErr)
                    tmp=norm(WrapErr)
                    if tmp>ERROR
                        println("Forward WrapTime error for at lt=$lt : $tmp")
                    end
                end
                #####################################################################
            end
        end
        
        # println("\n #####################################################################")

        for lt in model.Nt:-1:1
            
            #####################################################################
            Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt,div(model.Nt,2),"Backward")
            Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt,div(model.Nt,2),"Backward")
                
            if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>ERROR
                println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
                error("$lt : WrapTime")
            end
            #####################################################################

            for j in axes(ss[1],2)
                # update
                UpdateSCEELayer!(rng,j,view(ss[1],:,j,lt),view(ss[2],:,j,lt),G1,G2,A,B,model,UPD,SCEE,λ)
                for i in axes(ss[1],1)
                    x,y=model.nnidx[i,j]
                    tmpN[x]=model.η[ss[1][i,j,lt]]
                    tmpN[y]=-model.η[ss[1][i,j,lt]]
                    tmpN_[x]=model.η[ss[2][i,j,lt]]
                    tmpN_[y]=-model.η[ss[2][i,j,lt]]
                end
                tmpN.= exp.(.-tmpN)
                tmpN_.= exp.(.-tmpN_)

                WrapV!(tmpNN,Gt01,tmpN,view(model.UV,:,:,j),1)
                WrapV!(tmpNN,Gt02,tmpN_,view(model.UV,:,:,j),1)

                WrapV!(tmpNN,Gt1,tmpN,view(model.UV,:,:,j),3)
                WrapV!(tmpNN,Gt2,tmpN_,view(model.UV,:,:,j),3)

                WrapV!(tmpNN,G0t1,tmpN,view(model.UV,:,:,j),2)
                WrapV!(tmpNN,G0t2,tmpN_,view(model.UV,:,:,j),2)

            end

            WrapK!(tmpNN,G1,model.eKinv,model.eK)
            WrapK!(tmpNN,G2,model.eKinv,model.eK)

            ##------------------------------------------------------------------------
            tmpO+=(A.detg/B.detg)^(1/Nλ)
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
                    LAPACK.orgrq!(tmpnN, tau)
                    copyto!(view(BLMs1,:,:,i) , tmpnN)

                    mul!(tmpnN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau)
                    copyto!(view(BLMs2,:,:,i) , tmpnN)
                end
                for i in idx+1:max(Θidx,idx)
                    # println("update BR i=",i)
                    mul!(tmpNn, view(BMs1,:,:,i-1), view(BRMs1,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau)
                    copyto!(view(BRMs1,:,:,i) , tmpNn)

                    mul!(tmpNn, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau)
                    copyto!(view(BRMs2,:,:,i) , tmpNn)
                end

                #####################################################################
                    if lt-1 != div(model.Nt,2)
                        copyto!(WrapErr, Gt1)
                        axpy!(1.0, Gt2, WrapErr)
                        axpy!(1.0, G01, WrapErr)
                        axpy!(1.0, G02, WrapErr)
                        axpy!(1.0, Gt01, WrapErr)
                        axpy!(1.0, Gt02, WrapErr)
                        axpy!(1.0, G0t1, WrapErr)
                        axpy!(1.0, G0t2, WrapErr)
                    end
                #####################################################################

                get_ABGM!(G1,G2,A,B,SCEE,model.nodes,idx,"Backward")
                
                #####################################################################
                    if lt-1 != div(model.Nt,2)
                        axpy!(-1.0, Gt1, WrapErr)
                        axpy!(-1.0, Gt2, WrapErr)
                        axpy!(-1.0, G01, WrapErr)
                        axpy!(-1.0, G02, WrapErr)
                        axpy!(-1.0, Gt01, WrapErr)
                        axpy!(-1.0, Gt02, WrapErr)
                        axpy!(-1.0, G0t1, WrapErr)
                        axpy!(-1.0, G0t2, WrapErr)
                        tmp=norm(WrapErr)
                        if tmp>ERROR
                            println("Backward WrapTime error for at lt=$(lt-1) : $tmp")
                        end
                    end
                #####################################################################
            end
        end


        O[loop+1]=tmpO/counter
        tmpO=counter=0
    end
    return ss
end

function get_ABGM!(G1::G4Buffer_,G2::G4Buffer_,A::AreaBuffer_,B::AreaBuffer_,SCEE::SCEEBuffer_,nodes,idx,direction::String="Backward")
    G4!(SCEE,G1,nodes,idx,direction)
    G4!(SCEE,G2,nodes,idx,direction)
    GroverMatrix!(A.gmInv,view(G1.G0,A.index,A.index),view(G2.G0,A.index,A.index))
    A.detg=det(A.gmInv)
    LAPACK.getrf!(A.gmInv, A.ipiv)
    LAPACK.getri!(A.gmInv, A.ipiv)

    GroverMatrix!(B.gmInv,view(G1.G0,B.index,B.index),view(G2.G0,B.index,B.index))
    B.detg=det(B.gmInv)
    LAPACK.getrf!(B.gmInv, B.ipiv)
    LAPACK.getri!(B.gmInv, B.ipiv)
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
function G4update!(SCEE::SCEEBuffer_,UPD::UpdateBuffer_,G::G4Buffer_)
    mul!(SCEE.zN,UPD.r,view(G.Gt0,UPD.subidx,:))   # useful for GΘ,GτΘ
    mul!(SCEE.NN, view(G.G0t,:,UPD.subidx),SCEE.zN)
    axpy!(1.0, SCEE.NN, G.G0)
    mul!(SCEE.NN, view(G.Gt,:,UPD.subidx),SCEE.zN)
    axpy!(1.0, SCEE.NN, G.Gt0)

    mul!(SCEE.zN,UPD.r,view(G.Gt,UPD.subidx,:))
    lmul!(-1.0,SCEE.zN)
    axpy!(1.0,UPD.r,view(SCEE.zN,:,UPD.subidx))   # useful for GΘτ,Gτ
    mul!(SCEE.NN, view(G.G0t,:,UPD.subidx),SCEE.zN)
    axpy!(-1.0, SCEE.NN, G.G0t)
    mul!(SCEE.NN, view(G.Gt,:,UPD.subidx),SCEE.zN)
    axpy!(-1.0, SCEE.NN, G.Gt)
end

"""
    Update gmInv with a,b,Tau.
        gmInv = gmInv - gmInv ⋅ a ⋅ inv(Tau) ⋅ b
    Universal for s1 and s2.
    ------------------------------------------------------------------------------
"""
function GMupdate!(A::AreaBuffer_)
    mul!(A.N2, A.gmInv,A.a )
    # inv22!(A.Tau)
    A.Tau .= inv(A.Tau)
    mul!(A.zN,A.Tau,A.b)
    mul!(A.NN, A.N2, A.zN)
    axpy!(-1.0, A.NN, A.gmInv)
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
function get_abTau1!(A::AreaBuffer_,UPD::UpdateBuffer_,G0::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64})
    copyto!(A.NN, view(G0,A.index,A.index))
    lmul!(2.0, A.NN)
    for i in diagind(A.NN)
        A.NN[i] -= 1
    end
    mul!(A.zN,view(Gt0,UPD.subidx,A.index),A.NN)
    mul!(A.b, A.zN, A.gmInv)
    mul!(A.a,view(G0t,A.index,UPD.subidx),UPD.r)
    mul!(A.Tau, A.b, A.a)
    A.Tau[1,1]+=1; A.Tau[2,2]+=1;
    ans=det(A.Tau)
    @assert ans>0 "detTau1 sign error !"    
    return ans
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
function get_abTau2!(A::AreaBuffer_,UPD::UpdateBuffer_,G0::Matrix{Float64},Gt0::Matrix{Float64},G0t::Matrix{Float64})
    copyto!(A.NN , view(G0,A.index,A.index))
    lmul!(2.0, A.NN)
    for i in diagind(A.NN)
        A.NN[i] -= 1
    end
    mul!(A.a,A.NN,view(G0t,A.index,UPD.subidx))
    mul!(A.zN,UPD.r,view(Gt0,UPD.subidx,A.index))
    mul!(A.b,A.zN,A.gmInv)
    mul!(A.Tau, A.b, A.a)
    A.Tau[1,1]+=1; A.Tau[2,2]+=1;
    ans=det(A.Tau)
    @assert ans>0 "detTau2 sign error !"    
    return ans
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
function WrapK!(tmpNN::Matrix{Float64},G::G4Buffer_,eK::Matrix{Float64},eKinv::Matrix{Float64})
    mul!(tmpNN,G.Gt,eKinv)
    mul!(G.Gt,eK,tmpNN)
    
    mul!(tmpNN, eK, G.Gt0)
    copyto!(G.Gt0, tmpNN)
    mul!(tmpNN, G.G0t,eKinv)
    copyto!(G.G0t, tmpNN)
end

function GroverMatrix!(GM::Matrix{Float64},G1::SubArray{Float64, 2, Matrix{Float64}, Tuple{Vector{Int64}, Vector{Int64}}, false},G2::SubArray{Float64, 2, Matrix{Float64}, Tuple{Vector{Int64}, Vector{Int64}}, false})
    mul!(GM,G1,G2)
    lmul!(2.0, GM)
    axpy!(-1.0, G1, GM)
    axpy!(-1.0, G2, GM)
    for i in diagind(GM)
        GM[i] += 1.0
    end
end


function G4!(SCEE::SCEEBuffer_,G::G4Buffer_,nodes::Vector{Int64},idx::Int64,direction="Forward")
    Θidx=div(length(nodes),2)+1
    BLMs, BRMs, BMs, BMinvs, Gt, G0, Gt0, G0t =
        G.BLMs, G.BRMs, G.BMs, G.BMinvs, G.Gt, G.G0, G.Gt0, G.G0t
    II, tmpnn, tmpNn, tmpNN, tmpNN_, ipiv =
        SCEE.II, SCEE.nn, SCEE.Nn, SCEE.NN, SCEE.NN_, SCEE.ipiv

    get_G!(tmpnn,tmpNn,ipiv,view(BLMs,:,:,idx),view(BRMs,:,:,idx),Gt)
    if idx==Θidx
        G0 .= Gt
        if direction=="Forward"
            Gt0.= Gt
            G0t.= Gt .- II 
        elseif direction=="Backward"
            Gt0.= Gt .- II
            G0t.= Gt
        end
    else
        get_G!(tmpnn,tmpNn,ipiv,view(BLMs,:,:,Θidx),view(BRMs,:,:,Θidx),G0)
    
        Gt0 .= II
        G0t .= II
        if idx<Θidx
            for j in idx:Θidx-1
                if j==idx
                    tmpNN_ .= Gt
                else
                    get_G!(tmpnn,tmpNn,ipiv,view(BLMs,:,:,j),view(BRMs,:,:,j),tmpNN_)
                end
                mul!(tmpNN,tmpNN_, G0t)
                mul!(G0t, view(BMs,:,:,j), tmpNN)
                tmpNN .= II .- tmpNN_
                mul!(tmpNN_,Gt0, tmpNN)
                mul!(Gt0, tmpNN_, view(BMinvs,:,:,j))
                
            end
            lmul!(-1.0, Gt0)
        else
            for j in Θidx:idx-1
                if j==Θidx
                    tmpNN_ .= G0
                else
                    get_G!(tmpnn,tmpNn,ipiv,view(BLMs,:,:,j),view(BRMs,:,:,j),tmpNN_)
                end
                mul!(tmpNN, tmpNN_, Gt0)
                mul!(Gt0, view(BMs,:,:,j), tmpNN)
                tmpNN .= II .- tmpNN_
                mul!(tmpNN_, G0t, tmpNN)
                mul!(G0t, tmpNN_,view(BMinvs,:,:,j))
            end
            lmul!(-1.0, G0t)
        end        
    end
end

function UpdateSCEELayer!(rng,j,s1,s2,G1::G4Buffer_,G2::G4Buffer_,A::AreaBuffer_,B::AreaBuffer_,model::Hubbard_Para_,UPD::UpdateBuffer_,SCEE::SCEEBuffer_,λ)
    for i in axes(s1,1)
        x,y=model.nnidx[i,j]
        UPD.subidx=[x,y]

        # update s1
        begin
            sx = rand(rng,  model.samplers_dict[s1[i]])
            p=get_r!(UPD,model.η[sx]- model.η[s1[i]],G1.Gt)
            p*=model.γ[sx]/model.γ[s1[i]]

            detTau_A=get_abTau1!(A,UPD,G2.G0,G1.Gt0,G1.G0t)
            detTau_B=get_abTau1!(B,UPD,G2.G0,G1.Gt0,G1.G0t)

            @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
            if rand(rng)<p
                A.detg*=detTau_A
                B.detg*=detTau_B

                GMupdate!(A)
                GMupdate!(B)
                G4update!(SCEE,UPD,G1)
                s1[i]=sx
            end
        end

        # update ss[2]
        begin
            sx = rand(rng,  model.samplers_dict[s2[i]])
            p=get_r!(UPD,model.η[sx]- model.η[s2[i]],G2.Gt)
            p*=model.γ[sx]/model.γ[s2[i]]

            detTau_A=get_abTau2!(A,UPD,G1.G0,G2.Gt0,G2.G0t)
            detTau_B=get_abTau2!(B,UPD,G1.G0,G2.Gt0,G2.G0t)

            @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
            if rand(rng)<p
                A.detg*=detTau_A
                B.detg*=detTau_B

                GMupdate!(A)
                GMupdate!(B)
                G4update!(SCEE,UPD,G2)
                s2[i]=sx
            end
        end
    end
end


