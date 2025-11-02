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

    II = Diagonal(ones(Float64, model.Ns)) 
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
    tmpN2 = Matrix{Float64}(undef, Ns,2)
    tmp22 = Matrix{Float64}(undef, 2,2)
    tmp2 = Vector{Float64}(undef,2)

    r = Matrix{Float64}(undef, 2,2)
    Δ = Matrix{Float64}(undef, 2,2)


    Θidx=div(length(model.nodes),2)+1

    view(BRs,:,:,1) .= model.Pt
    view(BLs,:,:,NN) .= model.Pt'

    for idx in NN-1:-1:2
        BM_F!(BM,model, s, idx)
        mul!(tmpnN,view(BLs,:,:,idx+1), BM)
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau, ns)
        copyto!(view(BLs,:,:,idx), tmpnN)
        # view(BLs,:,:,idx) .= Matrix(qr!(tmpNn).Q)'
    end

    for loop in 1:Sweeps
        # println("\n Sweep: $loop ")
        
        BM_F!(BM,model, s, 1)
        mul!(tmpnN,view(BLs,:,:,2), BM)
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau, ns)
        copyto!(view(BLs,:,:,1), tmpnN)

        mul!(tmpnn, view(BLs,:,:,1), view(BRs,:,:,1))
        LAPACK.getrf!(tmpnn,ipiv)
        LAPACK.getri!(tmpnn, ipiv)
        mul!(tmpNn, view(BRs,:,:,1), tmpnn)
        mul!(tmpNN, tmpNn, view(BLs,:,:,1))
        @fastmath G .= II .- tmpNN

        idx=1
        for lt in 1:model.Nt
            if any(view(model.nodes,2:NN) .== (lt - 1))
                idx+=1
                BM_F!(BM,model, s, idx - 1)
                mul!(tmpNn, BM, view(BRs,:,:,idx-1))
                LAPACK.geqrf!(tmpNn, tau)
                LAPACK.orgqr!(tmpNn, tau, ns)
                copyto!(view(BRs,:,:,idx), tmpNn)
                
                mul!(tmpnn, view(BLs,:,:,idx), view(BRs,:,:,idx))
                LAPACK.getrf!(tmpnn,ipiv)
                LAPACK.getri!(tmpnn, ipiv)
                mul!(tmpNn, view(BRs,:,:,idx), tmpnn)
                mul!(tmpNN, tmpNn, view(BLs,:,:,idx))

                # println(norm(G-II+tmpNN))
                @fastmath G .= II .- tmpNN
                axpy!(1.0, G, tmpNN)  
                axpy!(-1.0, II, tmpNN)  
                if norm(tmpNN)>1e-8
                    println("Warning for Batchsize Wrap Error : $(norm(tmpNN))")
                end

            end

            #####################################################################
            # println(lt)
            if norm(G-Gτ(model,s,lt-1))>1e-5
                println("\n Sweep: $loop ")
                error("Wrap-$(lt)   :   $(norm(G-Gτ(model,s,lt-1)))")
            end
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

                mul!(tmpNN,view(model.UV,:,:,j)',G)
                mul!(G,tmpNN,view(model.UV,:,:,j))
                
                mul!(tmpNN,Diagonal(tmpN),G)
                tmpN.= 1 ./tmpN
                mul!(G,tmpNN,Diagonal(tmpN))

                mul!(tmpNN,view(model.UV,:,:,j),G)
                mul!(G,tmpNN,view(model.UV,:,:,j)')
                # G=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:] *G* model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:]

                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]
                    tmp2[1]=-2*s[lt,i,j];   tmp2[2]=2*s[lt,i,j];
                    tmp2 .= exp.(model.α.*tmp2)
                    mul!(tmp22,uv,Diagonal(tmp2))
                    mul!(Δ,tmp22,uv')
                    Δ[1,1]-=1 ; Δ[2,2]-=1;

                    tmp22.= .-view(G,subidx,subidx)
                    tmp22[1,1]+=1;  tmp22[2,2]+=1;
                    mul!(r,Δ,tmp22)
                    r[1,1]+=1; r[2,2]+=1;
                    # r=I(2)+Δ*(I(2)-G[subidx,subidx])
                    
                    detR=det(r)

                    if detR<-1e-3
                        println("Negative Sign: $(detR)")
                    end
                    #####################################################################
                    # ss=copy(s)
                    # ss[lt,i,j]=-ss[lt,i,j]
                    # dassda=detR-Poss(model,ss)/Poss(model,s)
                    # if abs(dassda)>1e-5
                    #     error("Poss error lt-$(lt) No.$(j): $(dassda)")
                    # end
                    #####################################################################

                    if rand(rng)<detR
                        inv22!(tmp22,r)
                        mul!(r,tmp22,Δ)
                        mul!(tmpN2,view(G,:,subidx),r)

                        tmp2N .= .-view(G, subidx , :)
                        tmp2N[1,x] += 1;   tmp2N[2,y] += 1;

                        mul!(tmpNN,tmpN2,tmp2N)
                        axpy!(-1.0,tmpNN,G)
                        # G-=(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                        s[lt,i,j]=-s[lt,i,j]

                        #####################################################################
                        print("*")
                        GG=model.eK*Gτ(model,s,lt-1)*model.eKinv
                        for jj in 3:-1:j
                            E=zeros(model.Ns)
                            for ii in 1:size(s)[2]
                                x,y=model.nnidx[ii,jj]
                                E[x]=s[lt,ii,jj]
                                E[y]=-s[lt,ii,jj]
                            end
                            GG=model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]' *GG* model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]'
                        end
                        if(norm(G-GG)>1e-5)
                            println("loop=$(loop) lt=$(lt) j=$(j)")
                            error(j," update error: ",norm(G-GG),"  lt=",lt)
                        end
                        #####################################################################
                    end
                end
            end

            if record && abs(idx-Θidx)<=1
                tmp=phy_measure(model,G,lt,s)
                Ek+=tmp[1]
                Ev+=tmp[2]
                axpy!(1,tmp[3],R0)
                axpy!(1,tmp[4],R1)
                counter+=1
            end
        end

        BM_F!(BM,model, s, NN-1)
        mul!(tmpNn , BM, view(BRs,:,:,NN-1))
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau, ns)
        view(BRs,:,:,NN) .= tmpNn

        mul!(tmpnn, view(BLs,:,:,NN), view(BRs,:,:,NN))
        LAPACK.getrf!(tmpnn,ipiv)
        LAPACK.getri!(tmpnn, ipiv)
        mul!(tmpNn, view(BRs,:,:,NN), tmpnn)
        mul!(tmpNN, tmpNn, view(BLs,:,:,NN))
        @fastmath G .= II .- tmpNN

        for lt in model.Nt:-1:1
            if any(view(model.nodes,1:NN-1).== lt)
                BM_F!(BM,model, s, idx)
                mul!(tmpnN,view(BLs,:,:,idx+1),BM)
                LAPACK.gerqf!(tmpnN, tau)
                LAPACK.orgrq!(tmpnN, tau, ns)
                view(BLs,:,:,idx).=tmpnN
                # BL .= Matrix(qr(( BL * BM )').Q)'

                mul!(tmpnn, view(BLs,:,:,idx), view(BRs,:,:,idx))
                LAPACK.getrf!(tmpnn,ipiv)
                LAPACK.getri!(tmpnn, ipiv)
                mul!(tmpNn, view(BRs,:,:,idx), tmpnn)
                mul!(tmpNN, tmpNn, view(BLs,:,:,idx))
                @fastmath G .= II .- tmpNN
                idx-=1
            end
            #####################################################################
            # if norm(G-Gτ(model,s,lt))>1e-4
            #     error("Wrap-$(lt)   :   $(norm(G-Gτ(model,s,lt+1)))")
            # end
            #####################################################################

            for j in 1:size(s)[3]
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]
                    tmp2[1]=-2*s[lt,i,j];   tmp2[2]=2*s[lt,i,j];
                    tmp2 .= exp.(model.α.*tmp2)
                    mul!(tmp22,uv,Diagonal(tmp2))
                    mul!(Δ,tmp22,uv')
                    Δ[1,1]-=1 ; Δ[2,2]-=1;

                    tmp22.= .-view(G,subidx,subidx)
                    tmp22[1,1]+=1;  tmp22[2,2]+=1;
                    mul!(r,Δ,tmp22)
                    r[1,1]+=1; r[2,2]+=1;
                    # r=I(2)+Δ*(I(2)-G[subidx,subidx])
                    
                    detR=det(r)

                    if detR<-1e-3
                        println("Negative Sign: $(detR)")
                    end

                    if rand(rng)<detR
                        inv22!(tmp22,r)
                        mul!(r,tmp22,Δ)
                        mul!(tmpN2,view(G,:,subidx),r)

                        tmp2N .= .-view(G, subidx , :)
                        tmp2N[1,x] += 1;   tmp2N[2,y] += 1;

                        mul!(tmpNN,tmpN2,tmp2N)
                        axpy!(-1,tmpNN,G)
                        # G-=(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                        s[lt,i,j]=-s[lt,i,j]
                    
                        #####################################################################
                        # GG=Gτ(model,s,lt)
                        # for jj in 1:j-1
                        #     E=zeros(model.Ns)
                        #     for ii in 1:size(s)[2]
                        #         x,y=model.nnidx[ii,jj]
                        #         E[x]=s[lt,ii,jj]
                        #         E[y]=-s[lt,ii,jj]
                        #     end
                        #     GG=model.UV[:,:,jj]*Diagonal(exp.(-model.α*E))*model.UV[:,:,jj]' *GG* model.UV[:,:,jj]*Diagonal(exp.(model.α*E))*model.UV[:,:,jj]'
                        # end
                        # if(norm(G-GG)>1e-5)
                        #     error(j," inverse update error: ",norm(G-GG),"  ",i)
                        # end
                        #####################################################################
                    end
                end
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    tmpN[x]=s[lt,i,j]
                    tmpN[y]=-s[lt,i,j]
                end
                tmpN.=exp.(model.α.*tmpN)
                
                mul!(tmpNN,view(model.UV,:,:,j)',G)
                mul!(G,tmpNN,view(model.UV,:,:,j))
                
                mul!(tmpNN,G,Diagonal(tmpN))
                tmpN.= 1 ./tmpN
                mul!(G,Diagonal(tmpN),tmpNN)

                mul!(tmpNN,view(model.UV,:,:,j),G)
                mul!(G,tmpNN,view(model.UV,:,:,j)')
                # G=model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:] *G* model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:]
            end
            mul!(tmpNN,model.eKinv,G)
            mul!(G,tmpNN,model.eK)
            # G=model.eKinv*G*model.eK
            

            if record && abs(idx-Θidx)<=1
                tmp=phy_measure(model,G,lt-1,s)
                Ek+=tmp[1]
                Ev+=tmp[2]
                axpy!(1,tmp[3],R0)
                axpy!(1,tmp[4],R1)
                counter+=1
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
            BR=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:]*BR
        end
    end
    BR=model.Pt'*BR
    return det(BR)
end


function phy_measure(model,G,lt,s)
    """
    (Ek,Ev,R0,R1)    
    """
    #####################################################################
    # if norm(G-Gτ(model,s,lt))>1e-3
    #     println("record error lt=$(lt) : $(norm(G-Gτ(model,s,lt)))")
    # end
    #####################################################################

    G0=G[:,:]
    tmpN = Vector{Float64}(undef, model.Ns)
    tmpNN = Matrix{Float64}(undef, model.Ns, model.Ns)
    tmp=zeros(Float64,4)
    
    if lt>model.Nt/2
        for t in lt:-1:div(model.Nt,2)+1
            for j in 1:size(s)[3]
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    tmpN[x]=s[t,i,j]
                    tmpN[y]=-s[t,i,j]
                end
                tmpN.=exp.(model.α.*tmpN)
                mul!(tmpNN,view(model.UV,:,:,j)',G0)
                mul!(G0,tmpNN,view(model.UV,:,:,j))
                
                mul!(tmpNN,G0,Diagonal(tmpN))
                tmpN.= 1 ./tmpN
                mul!(G0,Diagonal(tmpN),tmpNN)

                mul!(tmpNN,view(model.UV,:,:,j),G0)
                mul!(G0,tmpNN,view(model.UV,:,:,j)')
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

                mul!(tmpNN,view(model.UV,:,:,j)',G0)
                mul!(G0,tmpNN,view(model.UV,:,:,j))
                
                mul!(tmpNN,Diagonal(tmpN),G0)
                tmpN .= 1 ./tmpN
                mul!(G0 ,tmpNN,Diagonal(tmpN))

                mul!(tmpNN,view(model.UV,:,:,j),G0 )
                mul!(G0 ,tmpNN,view(model.UV,:,:,j)')
                # G0=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:] *G0* model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:]
            end
        end
    end
    #####################################################################
    if norm(G0-Gτ(model,s,div(model.Nt,2)))>1e-5
        error("record error lt=$(lt) : $(norm(G0-Gτ(model,s,div(model.Nt,2))))")
    end
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


function inv22!(A,B)
    detB=det(B)
    A[1,1]=B[2,2]/detB
    A[1,2]=-B[1,2]/detB
    A[2,1]=-B[2,1]/detB
    A[2,2]=B[1,1]/detB
end