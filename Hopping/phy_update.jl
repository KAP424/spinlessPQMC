# Trotter e^V1 e^V2 e^V3 e^K

function phy_update(path::String,model::_Hubbard_Para,s::Array{Int8,3},Sweeps::Int64,record::Bool)
    uv=[-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]
    if model.Lattice=="SQUARE"
        name="□"
    elseif model.Lattice=="HoneyComb"
        name="HC"
    end
    file="$(path)H_phy$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv"

    rng=MersenneTwister(Threads.threadid())
    R0=R1=Ek=Ev=0
    counter=0

    
    G=zeros(Float64,model.Ns,model.Ns)
    for loop in 1:Sweeps
        for lt in 1:model.Nt
            if mod(lt,model.WrapTime)==1 || lt==div(model.Nt,2)+1
                G=model.eK*Gτ(model,s,lt-1)*model.eKinv
            else
                #####################################################################
                # if norm(G-Gτ(model,s,lt-1))>1e-5
                #     error("Wrap-$(lt)   :   $(norm(G-Gτ(model,s,lt-1)))")
                # end
                #####################################################################
                G=model.eK*G*model.eKinv
            end
            

            for j in size(s)[3]:-1:1
                E=zeros(model.Ns)
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    E[x]=s[lt,i,j]
                    E[y]=-s[lt,i,j]
                end
                G=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:] *G* model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:]

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
                
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]
                    E=[-2*s[lt,i,j] , 2*s[lt,i,j]]
                    Δ=uv*diagm(exp.(model.α.*E))*uv'-I(2)
                    r=I(2)+Δ*(I(2)-G[subidx,subidx])
                    detR=det(r)

                    if detR<0
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
                        G-=(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                        s[lt,i,j]=-s[lt,i,j]
                        #####################################################################
                        # GG=model.eK*Gτ(model,s,lt-1)*model.eKinv
                        # for jj in size(s)[3]:-1:j
                        #     E=zeros(model.Ns)
                        #     for ii in 1:size(s)[2]
                        #         x,y=model.nnidx[ii,jj]
                        #         E[x]=s[lt,ii,jj]
                        #         E[y]=-s[lt,ii,jj]
                        #     end
                        #     GG=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *GG* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        # end
                        # if(norm(G-GG)>1e-4)
                        #     error(j," error: ",norm(G-GG),"  ",i)
                        # end
                        #####################################################################
                    end
                end
            end

            if record && abs(lt-model.Nt/2)<=model.WrapTime
                tmp=phy_measure(model,G,lt,s)
                Ek+=tmp[1]
                Ev+=tmp[2]
                R0+=tmp[3]
                R1+=tmp[4]
                counter+=1
            end
        end

        for lt in model.Nt:-1:1
            if mod(lt,model.WrapTime)==1
                G=Gτ(model,s,lt)
            else
                #####################################################################
                # if norm(G-Gτ(model,s,lt))>1e-6
                #     error("ltltltl")
                # end
                #####################################################################
            end

            for j in 1:size(s)[3]
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]
                    E=[-2*s[lt,i,j] , 2*s[lt,i,j]]
                    Δ=uv'*diagm(exp.(model.α.*E))*uv-I(2)
                    r=I(2)+Δ*(I(2)-G[subidx,subidx])
                    detR=(det(r))
                    if detR<0
                        println("Negative Sign")
                    end
                    #####################################################################
                    # ss=copy(s)
                    # ss[lt,i,j]=-ss[lt,i,j]
                    # dassda=detR-Poss(model,ss)/Poss(model,s)
                    # if abs(dassda)>1e-5
                    #     error("Poss error: $(dassda)")
                    # end
                    #####################################################################
                    if rand(rng)<detR
                        G-=(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                        s[lt,i,j]=-s[lt,i,j]
                    
                        #####################################################################
                        # GG=model.eK*Gτ(model,s,lt-1)*model.eKinv
                        # for jj in size(s)[3]:-1:j
                        #     E=zeros(model.Ns)
                        #     for ii in 1:size(s)[2]
                        #         x,y=model.nnidx[ii,jj]
                        #         E[x]=s[lt,ii,jj]
                        #         E[y]=-s[lt,ii,jj]
                        #     end
                        #     GG=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *GG* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        # end
                        # if(norm(G-GG)>1e-4)
                        #     error(j," error: ",norm(G-GG),"  ",i)
                        # end
                        #####################################################################
                    end
                end
                E=zeros(model.Ns)
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    E[x]=s[lt,i,j]
                    E[y]=-s[lt,i,j]
                end
                G=model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:] *G* model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:]
            end
            G=model.eKinv*G*model.eK

            if record && abs(lt-1-model.Nt/2)<=model.WrapTime
                tmp=phy_measure(model,G,lt-1,s)
                Ek+=tmp[1]
                Ev+=tmp[2]
                R0+=tmp[3]
                R1+=tmp[4]
                counter+=1
            end
        end
        if record
            open(file, "a") do io
                lock(io)
                writedlm(io,[Ek Ev R0 R1]/counter,',')
                unlock(io)
            end
            R0=R1=Ek=Ev=0
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
    G0=G[:,:]
    
    if lt>model.Nt/2
        for t in lt:-1:div(model.Nt,2)+1
            for j in 1:size(s)[3]
                E=zeros(model.Ns)
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    E[x]=s[t,i,j]
                    E[y]=-s[t,i,j]
                end
                G0=model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:] *G0* model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:]
            end
            G0= model.eKinv*G0*model.eK
        end
    else
        for t in lt+1:div(model.Nt,2)
            G0=model.eK*G0*model.eKinv
            for j in size(s)[3]:-1:1
                E=zeros(model.Ns)
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    E[x]=s[t,i,j]
                    E[y]=-s[t,i,j]
                end
                G0=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:] *G0* model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:]
            end
        end
    end
    #####################################################################
    # if norm(G0-Gτ(model,s,div(model.Nt,2)))>1e-6 
    #     error("record error lt=$(lt) : $(norm(G0-Gτ(model,s,div(model.Nt,2))))")
    # end
    #####################################################################

    G0=model.HalfeK* G0 *model.HalfeKinv

    Ek=model.t*sum(model.K.*G0)
    Ev=R0=R1=0
    for k in 1:length(model.nnidx)
        x,y=model.nnidx[k]
        Ev+=(1-G0[x,x])*(1-G0[y,y])-G0[x,y]*G0[y,x]-1/4
    end
    Ev*=model.U
    if model.Lattice=="HoneyComb"
        for rx in 1:model.site[1]
            for ry in 1:model.site[2]
                tmp=0
                for ix in 1:model.site[1]
                    for iy in 1:model.site[2]
                        idx1=2*( (iy-1)*model.site[1]+ix) -1
                        idx2=2*( mod(iy+ry-1,model.site[2])*model.site[1]+mod1(ix+rx,model.site[1]) )-1
                        tmp+=(1-G0[idx1,idx1])*(1-G0[idx2,idx2])-G0[idx1,idx2]*G0[idx2,idx1]
                        tmp+=(1-G0[idx1+1,idx1+1])*(1-G0[idx2+1,idx2+1])-G0[idx1+1,idx2+1]*G[idx2+1,idx1+1]
                        tmp-=(1-G0[idx1+1,idx1+1])*(1-G0[idx2,idx2])-G0[idx1+1,idx2]*G0[idx2,idx1+1]
                        tmp-=(1-G0[idx1,idx1])*(1-G0[idx2+1,idx2+1])-G0[idx1,idx2+1]*G0[idx2+1,idx1]
                    end
                end
                tmp/=prod(model.site)
                R0+=tmp
                R1+=cos(2*π/model.site[1]*rx+2*π/model.site[2]*ry)*tmp
            end
        end
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