# 2d Trotter Decomposition
# using Hopping channel ±1 HS transformation
# Trotter e^V1 e^V2 e^V3 e^K

function BM_F!(BM,model::_Hubbard_Para, s::Array{Int8, 3}, idx::Int64)
    """
    不包头包尾
    """
    Ns=model.Ns
    nodes=model.nodes
    α=model.α
    @assert 0< idx <=length(model.nodes)

    tmpN = Vector{ComplexF64}(undef, Ns)
    tmpNN = Matrix{ComplexF64}(undef, Ns, Ns)

    fill!(tmpNN,0)
    @inbounds for i in diagind(tmpNN)
        tmpNN[i] = 1
    end

    for lt in nodes[idx] + 1:nodes[idx + 1]
        mul!(BM,model.eK,tmpNN)
        for j in 3:-1:1
            for i in 1:div(Ns,2)
                x,y=model.nnidx[i,j]
                tmpN[x]=s[lt,i,j]
                tmpN[y]=-s[lt,i,j]
            end
            tmpN.= exp.(α*tmpN)

            mul!(tmpNN,view(model.UV,:,:,j),BM)
            mul!(BM,Diagonal(tmpN),tmpNN)
            mul!(tmpNN,view(model.UV,:,:,j)',BM)
            copyto!(BM,tmpNN)
        end
    end
end

function BMinv_F!(BM,model::_Hubbard_Para, s::Array{Int8, 3}, idx::Int64)
    """
    不包头包尾
    """
    Ns=model.Ns
    nodes=model.nodes
    α=model.α
    @assert 0< idx <=length(model.nodes)

    tmpN = Vector{ComplexF64}(undef, Ns)
    tmpNN = Matrix{ComplexF64}(undef, Ns, Ns)

    fill!(tmpNN,0)
    @inbounds for i in diagind(tmpNN)
        tmpNN[i] = 1
    end

    for lt in nodes[idx] + 1:nodes[idx + 1]
        mul!(BM,tmpNN,model.eKinv)
        for j in 3:-1:1
            for i in 1:div(Ns,2)
                x,y=model.nnidx[i,j]
                tmpN[x]=s[lt,i,j]
                tmpN[y]=-s[lt,i,j]
            end
            tmpN.= exp.(-α*tmpN)

            mul!(tmpNN,BM,view(model.UV,:,:,j))
            mul!(BM,tmpNN,Diagonal(tmpN))
            mul!(tmpNN,BM,view(model.UV,:,:,j)')
            copyto!(BM,tmpNN)
        end
    end
end

function Initial_s(model::_Hubbard_Para,rng::MersenneTwister)::Array{Int8,3}
    sp=Random.Sampler(rng,[1,-1])
    a,b=size(model.nnidx)
    s=zeros(Int8,model.Nt,a,b)

    for i = 1:size(s)[1]
        for j = 1:size(s)[2]
            for k = 1:size(s)[3]
                s[i,j,k] =rand(rng,sp)
            end
        end  
    end  
    return s
end

"equal time Green function"
function Gτ(model::_Hubbard_Para,s::Array{Int8,3},τ::Int64)::Array{Float64,2}
    BL::Array{Float64,2}=model.Pt'[:,:]
    BR::Array{Float64,2}=model.Pt[:,:]

    E=zeros(model.Ns)
    counter=0
    for lt in model.Nt:-1:τ+1
        for j in 1:size(s)[3]
            fill!(E,0.0)
            for i in 1:size(s)[2]
                x,y=model.nnidx[i,j]
                E[x]=s[lt,i,j]
                E[y]=-s[lt,i,j]
            end
            BL=BL*model.UV[:,:,j]*Diagonal(exp.(model.α*E))*model.UV[:,:,j]'

            #####################################################################
            # V=zeros(Float64,model.Ns,model.Ns)
            # for i in 1:size(s)[2]
            #     x,y=model.nnidx[i,j]
            #     V[x,y]=V[y,x]=s[lt,i,j]
            # end
            # tmp=model.UV[j,:,:]'*diagm(E)*model.UV[j,:,:]
            # if norm(tmp-V)>1e-6
            #     println("diagnose error")
            # end
            #####################################################################
        end
        BL=BL*model.eK
        counter+=1
        if counter==model.BatchSize
            counter=0
            BL=Matrix(qr(BL').Q)'
        end
    end
    counter=0
    for lt in 1:1:τ
        BR=model.eK*BR
        for j in size(s)[3]:-1:1
            fill!(E,0.0)
            for i in 1:size(s)[2]
                x,y=model.nnidx[i,j]
                E[x]=s[lt,i,j]
                E[y]=-s[lt,i,j]
            end
            BR=model.UV[:,:,j]*Diagonal(exp.(model.α*E))*model.UV[:,:,j]' *BR
            #####################################################################
            # V=zeros(Float64,model.Ns,model.Ns)
            # for i in 1:size(s)[2]
            #     x,y=model.nnidx[i,j]
            #     V[x,y]=V[y,x]=s[lt,i,j]
            # end
            # tmp=model.UV[j,:,:]'*diagm(E)*model.UV[j,:,:]
            # if norm(tmp-V)>1e-6
            #     println("diagnose error")
            # end
            #####################################################################
        end
        counter+=1
        if counter==model.BatchSize
            counter=0
            BR=Matrix(qr(BR).Q)
        end
    end

    BL=Matrix(qr(BL').Q)'
    BR=Matrix(qr(BR).Q)

    return I(model.Ns)-BR*inv(BL*BR)*BL
    # return BL,BR
end


"displaced Green function G(τ₁,τ₂)"
function G4(model::_Hubbard_Para,s::Array{Int8,3},τ1::Int64,τ2::Int64)
    if τ1>τ2
        BBs=zeros(Float64,cld(τ1-τ2,model.BatchSize),model.Ns,model.Ns)
        BBsInv=zeros(Float64,size(BBs))
        
        UL=zeros(Float64,1+size(BBs)[1],div(model.Ns,2),model.Ns)
        UR=zeros(Float64,size(UL)[1],model.Ns,div(model.Ns,2))
        G=zeros(Float64,size(UL)[1],model.Ns,model.Ns)

        UL[end,:,:]=model.Pt'[:,:]
        UR[1,:,:]=model.Pt[:,:]
    
        counter=0
        for lt in 1:τ2
            UR[1,:,:]=model.eK*UR[1,:,:]
            for j in size(s)[3]:-1:1
                E=zeros(model.Ns)
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    E[x]=s[lt,i,j]
                    E[y]=-s[lt,i,j]
                end
                UR[1,:,:]=model.UV[j,:,:]'*diagm(exp.(model.α.*E))*model.UV[j,:,:]*UR[1,:,:]
            end

            counter+=1
            if counter==model.BatchSize
                counter=0
                UR[1,:,:]=Matrix(qr(UR[1,:,:]).Q)
            end
        end
        UR[1,:,:]=Matrix(qr(UR[1,:,:]).Q) 
    
        counter=0
        for lt in model.Nt:-1:τ1+1
            for j in 1:size(s)[3]
                E=zeros(model.Ns)
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    E[x]=s[lt,i,j]
                    E[y]=-s[lt,i,j]
                end
                UL[end,:,:]=UL[end,:,:]*model.UV[j,:,:]'*diagm(exp.(model.α.*E))*model.UV[j,:,:]
            end
            UL[end,:,:]=UL[end,:,:]*model.eK

            counter+=1
            if counter==model.BatchSize
                counter=0
                UL[end,:,:]=Matrix(qr(UL[end,:,:]').Q)'
            end
        end
        UL[end,:,:]=Matrix(qr(UL[end,:,:]').Q)'
    
        for lt in 1:size(BBs)[1]-1
            BBs[lt,:,:]=I(model.Ns)
            BBsInv[lt,:,:]=I(model.Ns)
            for lt2 in 1:model.BatchSize
                BBs[lt,:,:]=model.eK*BBs[lt,:,:]
                BBsInv[lt,:,:]=BBsInv[lt,:,:]*model.eKinv
                for j in size(s)[3]:-1:1
                    E=zeros(model.Ns)
                    for i in 1:size(s)[2]
                        x,y=model.nnidx[i,j]
                        E[x]=s[τ2+(lt-1)*model.BatchSize+lt2,i,j]
                        E[y]=-s[τ2+(lt-1)*model.BatchSize+lt2,i,j]
                    end
                    BBs[lt,:,:]=model.UV[j,:,:]'*diagm(exp.(model.α.*E))*model.UV[j,:,:]*BBs[lt,:,:]
                    BBsInv[lt,:,:]=BBsInv[lt,:,:]*model.UV[j,:,:]'*diagm(exp.(-model.α.*E))*model.UV[j,:,:]

                end
            end
        end
    
        BBs[end,:,:]=I(model.Ns)
        BBsInv[end,:,:]=I(model.Ns)
        for lt in τ2+(size(BBs)[1]-1)*model.BatchSize+1:τ1
            BBs[end,:,:]=model.eK*BBs[end,:,:]
            BBsInv[end,:,:]=BBsInv[end,:,:]*model.eKinv
            for j in size(s)[3]:-1:1
                E=zeros(model.Ns)
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    E[x]=s[lt,i,j]
                    E[y]=-s[lt,i,j]
                end
                BBs[end,:,:]=model.UV[j,:,:]'*diagm(exp.(model.α.*E))*model.UV[j,:,:]*BBs[end,:,:]
                BBsInv[end,:,:]=BBsInv[end,:,:]*model.UV[j,:,:]'*diagm(exp.(-model.α.*E))*model.UV[j,:,:]
            end
        end
    
        for i in 1:size(BBs)[1]
            UL[end-i,:,:]=Matrix(qr( (UL[end-i+1,:,:]*BBs[end-i+1,:,:])' ).Q)'
            UR[i+1,:,:]=Matrix(qr(BBs[i,:,:]*UR[i,:,:]).Q)
        end

        for i in 1:size(G)[1]
            G[i,:,:]=I(model.Ns)-UR[i,:,:]*inv(UL[i,:,:]*UR[i,:,:])*UL[i,:,:]

            #####################################################################
            # if i <size(G)[1]
            #     if norm(Gτ(model,s,τ2+(i-1)*model.BatchSize)-G[i,:,:])>1e-3
            #         error("$i Gt:  $(norm(Gτ(model,s,τ2+(i-1)*model.BatchSize)-G[i,:,:]))")
            #     end
            # else
            #     if norm(Gτ(model,s,τ1)-G[i,:,:])>1e-3
            #         error("$i Gt:  $(norm(Gτ(model,s,τ1)-G[i,:,:]))")
            #     end
            # end
            #####################################################################
        end

        G12=I(model.Ns)
        G21=-I(model.Ns)
        for i in 1:size(BBs)[1]
            G12=G12*BBs[end-i+1,:,:]*G[end-i,:,:]
            G21=G21*( I(model.Ns)-G[i,:,:] )*BBsInv[i,:,:]
        end
        
        return G[end,:,:],G[1,:,:],G12,G21
    
    elseif τ1<τ2
        G2,G1,G21,G12=G4(model,s,τ2,τ1)
        return G1,G2,G12,G21
    else
        G=Gτ(model,s,τ1)
        return G,G,-(I(model.Ns)-G),G
    
    end

end


function GroverMatrix(G1::Array{Float64,2},G2::Array{Float64,2})::Array{Float64,2}
    II=I(size(G1)[1])
    return G1*G2+(II-G1)*(II-G2)
end


