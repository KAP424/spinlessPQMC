# Trotter e^V1 e^V2 e^V3 e^K

function phy_update(path::String,model::_Hubbard_Para,s::Array{UInt8,3},Sweeps::Int64,record::Bool)
    global LOCK=ReentrantLock()
    Ns=model.Ns
    ns=div(model.Ns, 2)
    NN=length(model.nodes)
    Θidx=div(NN,2)+1
    name = if model.Lattice=="SQUARE" "□" 
    elseif model.Lattice=="HoneyComb60" "HC" 
    elseif model.Lattice=="HoneyComb120" "HC120" 
    else error("Lattice: $(model.Lattice) is not allowed !") end  
    file="$(path)/H_phy$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv"
    rng=MersenneTwister(Threads.threadid()+time_ns())

    tau = Vector{Float64}(undef, ns)

    Ek=Ev=0.0
    R0=zeros(Float64,4)
    R1=zeros(Float64,4)
    counter=0

    G = Matrix{Float64}(undef ,model.Ns, model.Ns)

    # 预分配 BL 和 BR
    BLs = Array{Float64}(undef, ns, model.Ns,NN)
    BRs = Array{Float64}(undef, model.Ns, ns,NN)
    BM = Matrix{Float64}(undef, Ns, Ns)

    # 预分配临时数组
    tmpG=tmpPhyWorkspace(
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
        [-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]
    )

    copyto!(view(BRs,:,:,1) , model.Pt)
    transpose!(view(BLs,:,:,NN) , model.Pt)

    for idx in NN-1:-1:1
        BM_F!(tmpG,BM,model, s, idx)
        mul!(tmpG.nN,view(BLs,:,:,idx+1), BM)
        LAPACK.gerqf!(tmpG.nN, tau)
        LAPACK.orgrq!(tmpG.nN, tau, ns)
        copyto!(view(BLs,:,:,idx), tmpG.nN)
        # view(BLs,:,:,idx) .= Matrix(qr!(tmpNn).Q)'
    end

    idx=1
    get_G!(tmpG,view(BLs,:,:,1), view(BRs,:,:,1),G)
    for _ in 1:Sweeps
        # println("\n Sweep: $loop ")
        for lt in 1:model.Nt
            #####################################################################
            # println(lt)
            if norm(G-Gτ(model,s,lt-1))>1e-5
                error("Wrap-$(lt)   :   $(norm(G-Gτ(model,s,lt-1)))")
            end
            #####################################################################

            mul!(tmpG.NN,G,model.eKinv)
            mul!(G,model.eK,tmpG.NN)
            # G=model.eK*G*model.eKinv
            
            for j in 3:-1:1
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    tmpG.N[x]=model.η[s[lt,i,j]]
                    tmpG.N[y]=-model.η[s[lt,i,j]]
                end
                tmpG.N.= exp.(model.α.*tmpG.N)

                WrapV!(tmpG,G,view(model.UV,:,:,j),"B")

                phy_LayerUpdate!(tmpG,j,rng,model,G,view(s,lt,:,j))
            end

            if record && abs(idx-Θidx)<=1
                tmp=phy_measure(tmpG,model,G,lt,s)
                Ek+=tmp[1]
                Ev+=tmp[2]
                axpy!(1,tmp[3],R0)
                axpy!(1,tmp[4],R1)
                counter+=1
            end

            if any(model.nodes.== lt)
                idx+=1
                BM_F!(tmpG,BM,model, s, idx - 1)
                mul!(tmpG.Nn, BM, view(BRs,:,:,idx-1))
                LAPACK.geqrf!(tmpG.Nn, tau)
                LAPACK.orgqr!(tmpG.Nn, tau, ns)
                copyto!(view(BRs,:,:,idx), tmpG.Nn)
                copyto!(tmpG.NN , G)
                get_G!(tmpG,view(BLs,:,:,idx), view(BRs,:,:,idx),G)

                #####################################################################
                axpy!(-1.0, G, tmpG.NN)  
                diff=norm(tmpG.NN)
                if diff>1e-7
                    println("Warning for Batchsize Wrap Error : $(diff)")
                end
                #####################################################################

            end

        end

        for lt in model.Nt:-1:1
            
            #####################################################################
            if norm(G-Gτ(model,s,lt))>1e-4
                error("Wrap-$(lt)   :   $(norm(G-Gτ(model,s,lt+1)))")
            end
            ######################################################################

            for j in 1:size(s)[3]
                phy_LayerUpdate!(tmpG,j,rng,model,G,view(s,lt,:,j))
                for i in axes(s,2)
                    x,y=model.nnidx[i,j]
                    tmpG.N[x]=model.η[s[lt,i,j]]
                    tmpG.N[y]=-model.η[s[lt,i,j]]
                end
                tmpG.N.=exp.(.-model.α.*tmpG.N)
                WrapV!(tmpG,G,view(model.UV,:,:,j),"B")
            end
            mul!(tmpG.NN,model.eKinv,G)
            mul!(G,tmpG.NN,model.eK)
            
            if record && abs(idx-Θidx)<=1
                tmp=phy_measure(tmpG,model,G,lt,s)
                Ek+=tmp[1]
                Ev+=tmp[2]
                axpy!(1,tmp[3],R0)
                axpy!(1,tmp[4],R1)
                counter+=1
            end

            if any(model.nodes.== (lt-1))
                idx-=1
                BM_F!(tmpG,BM,model, s, idx)
                mul!(tmpG.nN,view(BLs,:,:,idx+1),BM)
                LAPACK.gerqf!(tmpG.nN, tau)
                LAPACK.orgrq!(tmpG.nN, tau, ns)
                copyto!(view(BLs,:,:,idx) , tmpG.nN)
                # BL .= Matrix(qr(( BL * BM )').Q)'
                get_G!(tmpG,view(BLs,:,:,idx), view(BRs,:,:,idx),G)
            end
        end

        if record
            lock(LOCK) do
                open(file, "a") do io
                    writedlm(io,vcat([Ek, Ev], R0, R1)' ./ counter, ',')
                end
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
function get_G!(tmpG,BL,BR,G)
    mul!(tmpG.nn, BL,BR)
    LAPACK.getrf!(tmpG.nn,tmpG.ipiv)
    LAPACK.getri!(tmpG.nn, tmpG.ipiv)
    mul!(tmpG.Nn, BR, tmpG.nn)
    mul!(G, tmpG.Nn, BL)
    lmul!(-1.0,G)
    for i in diagind(G)
        G[i]+=1
    end
end

function phy_measure(tmpG,model,G,lt,s)
    """
    (Ek,Ev,R0,R1)    
    """
    G0=G[:,:]
    tmp=zeros(Float64,4)
    R0=zeros(Float64,4)
    R1=zeros(Float64,4)
    
    if lt>model.Nt/2
        for t in lt:-1:div(model.Nt,2)+1
            for j in 1:size(s)[3]
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    tmpG.N[x]=model.η[s[t,i,j]]
                    tmpG.N[y]=-model.η[s[t,i,j]]
                end
                tmpG.N.=exp.(.-model.α.*tmpG.N)

                WrapV!(tmpG,G0,view(model.UV,:,:,j),"B")
                # G0=model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:] *G0* model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:]
            end
            mul!(tmpG.NN,model.eKinv,G0)
            mul!(G0,tmpG.NN,model.eK)
            # G0= model.eKinv*G0*model.eK
        end
    else
        for t in lt+1:div(model.Nt,2)
            mul!(tmpG.NN,G0,model.eKinv)
            mul!(G0,model.eK,tmpG.NN)
            # G0=model.eK*G0*model.eKinv
            for j in 3:-1:1
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    tmpG.N[x]=model.η[s[t,i,j]]
                    tmpG.N[y]=-model.η[s[t,i,j]]
                end
                tmpG.N.= exp.(model.α.*tmpG.N)
                WrapV!(tmpG,G0,view(model.UV,:,:,j),"B")
                # G0=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:] *G0* model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:]
            end
        end
    end
    #####################################################################
    # if norm(G0-Gτ(model,s,div(model.Nt,2)))>1e-3
    #     error("record error lt=$(lt) : $(norm(G0-Gτ(model,s,div(model.Nt,2))))")
    # end
    #####################################################################
    mul!(tmpG.NN,model.HalfeK,G0)
    mul!(G0,tmpG.NN,model.HalfeKinv)
    G0=model.HalfeK* G0 *model.HalfeKinv

    Ek=model.t*sum(model.K.*G0)
    Ev=0.0
    for k in 1:length(model.nnidx)
        x,y=model.nnidx[k]
        Ev+=(1-G0[x,x])*(1-G0[y,y])-G0[x,y]*G0[y,x]
    end
    # Ev*=model.U

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
                            # tmp[2]+=(1-G0[idx1+1,idx1+1]) * (1-G0[idx2+1,idx2+1]) - G0[idx1+1,idx2+1]*G0[idx2+1,idx1+1]
                        else
                            tmp[1]+=(1-G0[idx1,idx1]) * (1-G0[idx2,idx2]) - G0[idx1,idx2]*G0[idx2,idx1]
                            tmp[2]+=(1-G0[idx1+1,idx1+1]) * (1-G0[idx2+1,idx2+1]) - G0[idx1+1,idx2+1]*G0[idx2+1,idx1+1]
                        end
                        tmp[3]+=(1-G0[idx1+1,idx1+1])*(1-G0[idx2,idx2])-G0[idx1+1,idx2]*G0[idx2,idx1+1]
                        tmp[4]+=(1-G0[idx1,idx1])*(1-G0[idx2+1,idx2+1])-G0[idx1,idx2+1]*G0[idx2+1,idx1]
                    end
                end
                axpy!(1,tmp,R0)
                axpy!(cos(2*π/model.site[1]*rx+2*π/model.site[2]*ry),tmp,R1)
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
        Δ = uv ⋅ diag( exp(α⋅(s′-s)[1,-1]) - I ) ⋅ uvᵀ 
        r = I + Δ ⋅ (I - Gt[subidx,subidx])
        r ≡ inv(r) ⋅ ̇Δ .
    ------------------------------------------------------------------------------
"""
function get_r!(tmpG,α::Float64,Δs::Float64,subidx::Vector{Int64},Gt::Matrix{Float64})
    tmpG.z .= Δs.*[1.0, -1.0]
    tmpG.z .= exp.(α.*tmpG.z).-1
    mul!(tmpG.r,tmpG.uv,Diagonal(tmpG.z))
    mul!(tmpG.Δ,tmpG.r,tmpG.uv')
    mul!(tmpG.r,tmpG.Δ,view(Gt,subidx,subidx))
    axpby!(1.0,tmpG.Δ, -1.0, tmpG.r)   # r = I + Δ ⋅ (I - Gt1[subidx,subidx])
    tmpG.r[1,1]+=1; tmpG.r[2,2]+=1;
    p=det(tmpG.r)
    # redefine r=inv(r) ⋅ ̇Δ 
    inv22!(tmpG.zz,tmpG.r)
    mul!(tmpG.r,tmpG.zz,tmpG.Δ)
    return p
end

"""
    No Return. Overwrite G = G - G · inv(r) ⋅ Δ · (I-G)
    ------------------------------------------------------------------------------
"""
function Gupdate!(tmpG::tmpPhyWorkspace,subidx::Vector{Int64},G::Matrix{Float64})
    mul!(tmpG.zN,tmpG.r,view(G,subidx,:))
    lmul!(-1.0,tmpG.zN)
    axpy!(1.0,tmpG.r,view(tmpG.zN,:,subidx))   # useful for GΘτ,Gτ
    mul!(tmpG.NN, view(G,:,subidx),tmpG.zN)
    axpy!(-1.0, tmpG.NN, G)
end

"""
    Only for short imiginary time debug
    ------------------------------------------------------------------------------
""" 
# function Poss(model,s)
#     E=zeros(model.Ns)
#     BR=model.Pt[:,:]
#     p=1
#     for lt in 1:model.Nt
#         BR=model.eK*BR
#         for j in size(s)[3]:-1:1
#             for i in 1:size(s)[2]
#                 p*=model.γ[s[lt,i,j]]
#                 x,y=model.nnidx[i,j]
#                 E[x]=model.η[s[lt,i,j]]
#                 E[y]=-model.η[s[lt,i,j]]
#             end
#             BR=model.UV[:,:,j]*Diagonal(exp.(model.α*E))*model.UV[:,:,j]*BR
#         end
#     end
#     BR=model.Pt'*BR
#     return det(BR)*p
# end

function phy_LayerUpdate!(tmpG,j,rng,model,G,s)
    for i in axes(model.nnidx,1)
        x,y=model.nnidx[i,j]
        subidx=[x,y]

        sx = rand(rng,  model.samplers_dict[s[i]])
        p=get_r!(tmpG,model.α,model.η[sx]- model.η[s[i]],subidx,G)
        p*=model.γ[sx]/model.γ[s[i]]
        if p<-1e-3
            println("Negative Sign: $(p)")
        end

        if rand(rng)<p
            Gupdate!(tmpG,subidx,G)
            s[i]=sx
        end
    end
end

function WrapV!(tmpG,G::Matrix{Float64},UV::SubArray{Float64, 2, Array{Float64, 3}},LR::String)
    if LR=="L"
        mul!(tmpG.NN,UV',G)
        mul!(G,Diagonal(tmpG.N),tmpG.NN)
        mul!(tmpG.NN,UV,G)
        copyto!(G, tmpG.NN)
    elseif LR=="R"
        mul!(tmpG.NN, G , UV)
        mul!(G, tmpG.NN , Diagonal(tmpG.N))
        mul!(tmpG.NN, G , UV')
        copyto!(G, tmpG.NN)
    else
        mul!(tmpG.NN,UV',G)
        mul!(G,tmpG.NN,UV)
        mul!(tmpG.NN,Diagonal(tmpG.N),G)
        tmpG.N .= 1 ./tmpG.N
        mul!(G,tmpG.NN,Diagonal(tmpG.N))
        mul!(tmpG.NN,UV,G)
        mul!(G,tmpG.NN,UV')
    end
end