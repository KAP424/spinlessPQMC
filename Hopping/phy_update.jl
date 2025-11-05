# Trotter e^V1 e^V2 e^V3 e^K

function phy_update(path::String,model::_Hubbard_Para,s::Array{Int8,3},Sweeps::Int64,record::Bool)
    Ns=model.Ns
    ns=div(model.Ns, 2)
    NN=length(model.nodes)
    tau = Vector{Float64}(undef, ns)
    ipiv = Vector{LAPACK.BlasInt}(undef, ns)

    uv=[-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]

    name = if model.Lattice=="SQUARE" "□" 
    elseif model.Lattice=="HoneyComb60" "HC" 
    elseif model.Lattice=="HoneyComb120" "HC120" 
    else error("Lattice: $(model.Lattice) is not allowed !") end  

    file="$(path)H_phy$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv"

    rng=MersenneTwister(Threads.threadid()+time_ns())
    
    Ek=Ev=0.0
    R0=zeros(Float64,4)
    R1=zeros(Float64,4)
    counter=0

    G = Matrix{Float64}(undef ,model.Ns, model.Ns)

    # 预分配 BL 和 BR
    BLs = Array{Float64}(undef, ns, model.Ns,NN)
    BRs = Array{Float64}(undef, model.Ns, ns,NN)

    # 预分配临时数组
    tmpN = Vector{Float64}(undef, Ns)
    tmpNN = Matrix{Float64}(undef, Ns, Ns)
    BM = Matrix{Float64}(undef, Ns, Ns)
    tmpNn = Matrix{Float64}(undef, Ns, ns)
    tmpnn = Matrix{Float64}(undef, ns, ns)
    tmpnN = Matrix{Float64}(undef, ns, Ns)

    tmp2N = Matrix{Float64}(undef, 2, Ns)
    tmp22 = Matrix{Float64}(undef, 2,2)
    tmp2 = Vector{Float64}(undef,2)

    r = Matrix{Float64}(undef, 2,2)
    Δ = Matrix{Float64}(undef, 2,2)


    Θidx=div(length(model.nodes),2)+1

    view(BRs,:,:,1) .= model.Pt
    view(BLs,:,:,NN) .= model.Pt'

    for idx in NN-1:-1:1
        BM_F!(BM,model, s, idx)
        mul!(tmpnN,view(BLs,:,:,idx+1), BM)
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau, ns)
        copyto!(view(BLs,:,:,idx), tmpnN)
        # view(BLs,:,:,idx) .= Matrix(qr!(tmpNn).Q)'
    end

    for loop in 1:Sweeps
        # println("\n Sweep: $loop ")
        get_G!(tmpnn,tmpNn,ipiv,view(BLs,:,:,1), view(BRs,:,:,1),G)
        idx=1
        for lt in 1:model.Nt
            #####################################################################
            # println(lt)
            # if norm(G-Gτ(model,s,lt-1))>1e-5
            #     println("\n Sweep: $loop ")
            #     error("Wrap-$(lt)   :   $(norm(G-Gτ(model,s,lt-1)))")
            # end
            #####################################################################

            mul!(tmpNN,G,model.eKinv)
            mul!(G,model.eK,tmpNN)
            # G=model.eK*G*model.eKinv
            
            for j in 3:-1:1
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    tmpN[x]=s[lt,i,j]
                    tmpN[y]=-s[lt,i,j]
                end
                tmpN.= exp.(model.α.*tmpN)

                WrapV!(tmpNN,G,tmpN,view(model.UV,:,:,j),3)

                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]

                    p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,s[lt,i,j],subidx,G)
                    if p<-1e-3
                        println("Negative Sign: $(p)")
                    end
                    #####################################################################
                    # ss=copy(s)
                    # ss[lt,i,j]=-ss[lt,i,j]
                    # dassda=p-Poss(model,ss)/Poss(model,s)
                    # if abs(dassda)>1e-5
                    #     error("Poss error lt-$(lt) No.$(j): $(dassda)")
                    # end
                    #####################################################################

                    if rand(rng)<p
                        Gupdate!(tmpNN,tmp2N,subidx,r,G)
                        s[lt,i,j]=-s[lt,i,j]
                        #####################################################################
                        # print("*")
                        # GG=model.eK*Gτ(model,s,lt-1)*model.eKinv
                        # for jj in 3:-1:j
                        #     E=zeros(model.Ns)
                        #     for ii in 1:size(s)[2]
                        #         x,y=model.nnidx[ii,jj]
                        #         E[x]=s[lt,ii,jj]
                        #         E[y]=-s[lt,ii,jj]
                        #     end
                        #     GG=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]' *GG* model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                        # end
                        # if(norm(G-GG)>1e-5)
                        #     println("loop=$(loop) lt=$(lt) j=$(j)")
                        #     error(j," update error: ",norm(G-GG),"  lt=",lt)
                        # end
                        #####################################################################
                    end
                end
            end

            if record && abs(idx-Θidx)<=1
                tmp=phy_measure(tmpN,tmpNN,model,G,lt,s)
                Ek+=tmp[1]
                Ev+=tmp[2]
                axpy!(1,tmp[3],R0)
                axpy!(1,tmp[4],R1)
                counter+=1
            end

            if any(model.nodes.== lt)
                idx+=1
                BM_F!(BM,model, s, idx - 1)
                mul!(tmpNn, BM, view(BRs,:,:,idx-1))
                LAPACK.geqrf!(tmpNn, tau)
                LAPACK.orgqr!(tmpNn, tau, ns)
                copyto!(view(BRs,:,:,idx), tmpNn)
                
                tmpNN .= G

                get_G!(tmpnn,tmpNn,ipiv,view(BLs,:,:,idx), view(BRs,:,:,idx),G)

                #####################################################################
                axpy!(-1.0, G, tmpNN)  
                if norm(tmpNN)>1e-8
                    println("Warning for Batchsize Wrap Error : $(norm(tmpNN))")
                end
                #####################################################################

            end

        end


        for lt in model.Nt:-1:1
            
            #####################################################################
            # if norm(G-Gτ(model,s,lt))>1e-4
            #     error("Wrap-$(lt)   :   $(norm(G-Gτ(model,s,lt+1)))")
            # end
            #####################################################################

            for j in 1:size(s)[3]
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]
                    
                    p=get_r!(uv,tmp2,Δ,tmp22,r,model.α,s[lt,i,j],subidx,G)
                    if p<-1e-3
                        println("Negative Sign: $(p)")
                    end
                    if rand(rng)<p
                        Gupdate!(tmpNN,tmp2N,subidx,r,G)
                        s[lt,i,j]=-s[lt,i,j]
                    end
                end
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    tmpN[x]=s[lt,i,j]
                    tmpN[y]=-s[lt,i,j]
                end
                tmpN.=exp.(.-model.α.*tmpN)
                WrapV!(tmpNN,G,tmpN,view(model.UV,:,:,j),3)
            end
            mul!(tmpNN,model.eKinv,G)
            mul!(G,tmpNN,model.eK)
            
            if record && abs(idx-Θidx)<=1
                tmp=phy_measure(tmpN,tmpNN,model,G,lt-1,s)
                Ek+=tmp[1]
                Ev+=tmp[2]
                axpy!(1,tmp[3],R0)
                axpy!(1,tmp[4],R1)
                counter+=1
            end

            if any(model.nodes.== (lt-1))
                idx-=1
                BM_F!(BM,model, s, idx)
                mul!(tmpnN,view(BLs,:,:,idx+1),BM)
                LAPACK.gerqf!(tmpnN, tau)
                LAPACK.orgrq!(tmpnN, tau, ns)
                view(BLs,:,:,idx).=tmpnN
                # BL .= Matrix(qr(( BL * BM )').Q)'

                get_G!(tmpnn,tmpNn,ipiv,view(BLs,:,:,idx), view(BRs,:,:,idx),G)
            end
        end

        if record
            open(file, "a") do io
                lock(io)
                writedlm(io,vcat([Ek,Ev],R0,R1 )'./counter,',')
                unlock(io)
            end
            Ek=Ev=0.0
            fill!(R0,0.0)
            fill!(R1,0.0)
            counter=0
        end
    end
    return s
end 

"""
    No Return. Overwrite G 
        G = I - BR ⋅ inv(BL ⋅ BR) ⋅ BL 
    ------------------------------------------------------------------------------
"""
function get_G!(tmpnn,tmpNn,ipiv,BL,BR,G)
    mul!(tmpnn, BL,BR)
    LAPACK.getrf!(tmpnn,ipiv)
    LAPACK.getri!(tmpnn, ipiv)
    mul!(tmpNn, BR, tmpnn)
    mul!(G, tmpNn, BL)
    lmul!(-1.0,G)
    for i in diagind(G)
        G[i]+=1
    end
end

function phy_measure(tmpN,tmpNN,model,G,lt,s)
    """
    (Ek,Ev,R0,R1)    
    """
    G0=G[:,:]
    tmp=zeros(Float64,4)
    
    if lt>model.Nt/2
        for t in lt:-1:div(model.Nt,2)+1
            for j in 1:size(s)[3]
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    tmpN[x]=s[t,i,j]
                    tmpN[y]=-s[t,i,j]
                end
                tmpN.=exp.(.-model.α.*tmpN)

                WrapV!(tmpNN,G0,tmpN,view(model.UV,:,:,j),3)
                # G0=model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:] *G0* model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:]
            end
            mul!(tmpNN,model.eKinv,G0)
            mul!(G0,tmpNN,model.eK)
            # G0= model.eKinv*G0*model.eK
        end
    else
        for t in lt+1:div(model.Nt,2)
            mul!(tmpNN,G0,model.eKinv)
            mul!(G0,model.eK,tmpNN)
            # G0=model.eK*G0*model.eKinv
            for j in 3:-1:1
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    tmpN[x]=s[t,i,j]
                    tmpN[y]=-s[t,i,j]
                end
                tmpN.= exp.(model.α.*tmpN)

                WrapV!(tmpNN,G0,tmpN,view(model.UV,:,:,j),3)
                # G0=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:] *G0* model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:]
            end
        end
    end
    #####################################################################
    # if norm(G0-Gτ(model,s,div(model.Nt,2)))>1e-3
    #     error("record error lt=$(lt) : $(norm(G0-Gτ(model,s,div(model.Nt,2))))")
    # end
    #####################################################################
    mul!(tmpNN,model.HalfeK,G0)
    mul!(G0,tmpNN,model.HalfeKinv)
    # G0=model.HalfeK* G0 *model.HalfeKinv

    Ek=model.t*sum(model.K.*G0)
    Ev=0.0
    for k in 1:length(model.nnidx)
        x,y=model.nnidx[k]
        Ev+=(1-G0[x,x])*(1-G0[y,y])-G0[x,y]*G0[y,x]-1/4
    end
    Ev*=model.U

    R0=zeros(Float64,4)
    R1=zeros(Float64,4)
    if occursin("HoneyComb", model.Lattice)
        for rx in 1:model.site[1]
            for ry in 1:model.site[2]
                fill!(tmp,0.0)
                for ix in 1:model.site[1]
                    for iy in 1:model.site[2]
                        idx1=xy_i(model.Lattice,model.site,ix,iy)-1
                        idx2=xy_i(model.Lattice,model.site,mod1(ix+rx,model.site[1]),mod1(iy+ry,model.site[2]))-1
                        if idx1==idx2
                            tmp[1]+=(1-G0[idx1,idx1])
                            tmp[2]+=(1-G0[idx1+1,idx1+1])
                            # tmp[1]+=(1-G0[idx1,idx1]) * (1-G0[idx2,idx2]) - G0[idx1,idx2]*G0[idx2,idx1]
                            # tmp[2]+=(1-G0[idx1+1,idx1+1]) * (1-G0[idx2+1,idx2+1]) - G0[idx1+1,idx2+1]*G[idx2+1,idx1+1]
                        else
                            tmp[1]+=(1-G0[idx1,idx1]) * (1-G0[idx2,idx2]) - G0[idx1,idx2]*G0[idx2,idx1]
                            tmp[2]+=(1-G0[idx1+1,idx1+1]) * (1-G0[idx2+1,idx2+1]) - G0[idx1+1,idx2+1]*G[idx2+1,idx1+1]
                        end
                        tmp[3]+=(1-G0[idx1+1,idx1+1])*(1-G0[idx2,idx2])-G0[idx1+1,idx2]*G0[idx2,idx1+1]
                        tmp[4]+=(1-G0[idx1,idx1])*(1-G0[idx2+1,idx2+1])-G0[idx1,idx2+1]*G0[idx2+1,idx1]
                    end
                end
                axpy!(1,tmp,R0)
                axpy!(cos(2*π/model.site[1]*rx+2*π/model.site[2]*ry)/2,tmp,R1)
            end
        end
        lmul!(4/model.Ns^2,R0)
        lmul!(4/model.Ns^2,R1)
    elseif model.Lattice=="SQUARE"
        for rx in 1:model.site[1]
            for ry in 1:model.site[2]
                tmp=0
                for ix in 1:model.site[1]
                    for iy in 1:model.site[2]
                        idx1=ix+(iy-1)*model.site[1]
                        idx2=mod1(rx+ix,model.site[1])+mod((ry+iy-1),model.site[2])*model.site[1]
                        
                        tmp+=(1-G0[idx1,idx1])*(1-G0[idx2,idx2])-G0[idx1,idx2]*G0[idx2,idx1]
                    end
                end
                tmp/=prod(model.site)
                R0+=tmp*cos(π*(rx+ry))
                R1+=cos(π*(rx+ry)+2*π/model.site[1]*rx+2*π/model.site[2]*ry )*tmp
            end
        end
    end
    # 1-R1/R0
    return Ek,Ev,R0,R1
end

"""
    No Return. Overwrite A with inv(B)
    ------------------------------------------------------------------------------
"""
function inv22!(A,B)
    detB=det(B)
    A[1,1]=B[2,2]/detB
    A[1,2]=-B[1,2]/detB
    A[2,1]=-B[2,1]/detB
    A[2,2]=B[1,1]/detB
end

"""
    Return p=det(r). Overwrite 
        Δ = uv ⋅ diag( exp(α⋅[-2s,2s]) - I ) ⋅ uvᵀ 
        r = I + Δ ⋅ (I - Gt[subidx,subidx])
        r ≡ inv(r) ⋅ ̇Δ .
    ------------------------------------------------------------------------------
"""
function get_r!(uv::Matrix{Float64},tmp2::Vector{Float64},Δ::Matrix{Float64},tmp22::Matrix{Float64},r::Matrix{Float64},α::Float64,s::Int8,subidx::Vector{Int64},Gt::Matrix{Float64})
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

"""
    No Return. Overwrite G = G - G · inv(r) ⋅ Δ · (I-G)
    ------------------------------------------------------------------------------
"""
function Gupdate!(tmpNN::Matrix{Float64},tmp2N::Matrix{Float64},subidx::Vector{Int64},r::Matrix{Float64},G::Matrix{Float64})
    mul!(tmp2N,r,view(G,subidx,:))
    lmul!(-1.0,tmp2N)
    axpy!(1.0,r,view(tmp2N,:,subidx))   # useful for GΘτ,Gτ
    mul!(tmpNN, view(G,:,subidx),tmp2N)
    axpy!(-1.0, tmpNN, G)
end

"""
    Only for short imiginary time debug
    ------------------------------------------------------------------------------
""" 
function Poss(model,s)
    BR=model.Pt[:,:]
    for lt in 1:model.Nt
        BR=model.eK*BR
        for j in size(s)[3]:-1:1
            E=zeros(model.Ns)
            for i in 1:size(s)[2]
                x,y=model.nnidx[i,j]
                E[x]=s[lt,i,j]
                E[y]=-s[lt,i,j]
            end
            BR=model.UV[:,:,j]*Diagonal(exp.(model.α*E))*model.UV[:,:,j]*BR
        end
    end
    BR=model.Pt'*BR
    return det(BR)
end

