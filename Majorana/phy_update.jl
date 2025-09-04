# Trotter e^V1 e^V2 e^V3 e^K

function phy_update(path::String,model::_Hubbard_Para,s::Array{Int8,3},Sweeps::Int64,record::Bool)
    if model.Lattice=="SQUARE"
        name="□"
    elseif model.Lattice=="HoneyComb"
        name="HC"
    end
    file="$(path)M_PHY$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv"

    rng=MersenneTwister(Threads.threadid())
    R0=R1=Ek=Ev=0
    counter=0

    
    G=zeros(Float64,model.Ns,model.Ns)
    for loop in 1:Sweeps
        for lt in 0:model.Nt-1
            println("t:  ",lt)
            if mod(lt,model.WrapTime)==0
                G=model.eK*Gτ(model,s,lt)*model.eKinv
            else
                #####################################################################
                if norm(G-Gτ(model,s,lt))>1e-6 
                    error("$(lt):asdasdsa")
                end
                #####################################################################
                G=model.eK*G*model.eKinv
            end
            
            

            for j in size(s)[3]:-1:1
                V=zeros(ComplexF64,model.Ns,model.Ns)
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    V[x,y]=s[lt+1,i,j]*1im/4
                    V[y,x]=-s[lt+1,i,j]*1im/4
                end
                E=diag(model.UV[j,:,:]*V*model.UV[j,:,:]')
                G=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:] *G* model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:]
                
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]
                    Δ=[0 -1im/2*s[lt+1,i,j]; 1im/2*s[lt+1,i,j] 0]
                    E,U=eigen(Δ)
                    Δ=U*diagm(exp.(model.α*E))*U'-I(2)
                    r=I(2)+Δ*(I(2)-G[subidx,subidx])
                    detR=abs(det(r))

                    # println(detR)
                    if rand(rng)<detR
                        G-=(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                        s[lt+1,i,j]=-s[lt+1,i,j]
                        #####################################################################
                        ss=copy(s)
                        GG=model.eK*Gτ(model,ss,lt)*model.eKinv
                        VV=zeros(ComplexF64,model.Ns,model.Ns)
                        for jj in size(s)[3]:-1:j
                            for ii in 1:size(s)[2]
                                x,y=model.nnidx[ii,jj]
                                VV[x,y]=ss[lt+1,ii,jj]*1im/4
                                VV[y,x]=-ss[lt+1,ii,jj]*1im/4
                            end
                            E=diag(model.UV[jj,:,:]*VV*model.UV[jj,:,:]')
                            GG=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *GG* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        end

                        if(norm(G-GG)>1e-4)
                            println(j," error: ",norm(G-GG),"  ",i)
                        end
                        #####################################################################

                    end

                end
            end

            if record && abs(lt+1-model.Nt/2)<=model.WrapTime
                tmp=phy_measure(model,G,lt+1,s)
                Ek+=tmp[1]
                Ev+=tmp[2]
                R0+=tmp[3]
                R1+=tmp[4]
                counter+=1
            end
        end

        for lt in model.Nt:-1:1
            println("t:  ",lt)
            if mod(lt,model.WrapTime)==1
                G=Gτ(model,s,lt)
            else
                #####################################################################
                if norm(G-Gτ(model,s,lt))>1e-6
                    error("ltltltl")
                end
                #####################################################################
            end

            for j in 1:size(s)[3]
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]
                    Δ=[0 -1im/2*s[lt,i,j]; 1im/2*s[lt,i,j] 0]
                    E,U=eigen(Δ)
                    Δ=U*diagm(exp.(model.α*E))*U'-I(2)
                    r=I(2)+Δ*(I(2)-G[subidx,subidx])
                    detR=abs(det(r))

                    # println(detR)
                    if rand(rng)<detR
                        G-=(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                        s[lt,i,j]=-s[lt,i,j]
                    
                        #####################################################################
                        ss=copy(s)
                        GG=Gτ(model,ss,lt)
                        VV=zeros(ComplexF64,model.Ns,model.Ns)
                        for jj in 1:j-1
                            for ii in 1:size(s)[2]
                                x,y=model.nnidx[ii,jj]
                                VV[x,y]=ss[lt,ii,jj]*1im/4
                                VV[y,x]=-ss[lt,ii,jj]*1im/4
                            end
                            E=diag(model.UV[jj,:,:]*VV*model.UV[jj,:,:]')
                            GG=model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:] *GG* model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:]
                        end
            
                        if(norm(G-GG)>1e-4)
                            println(j," error: ",norm(G-GG),"  ",i)
                        end
                        #####################################################################
                    end
                end
                V=zeros(ComplexF64,model.Ns,model.Ns)
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    V[x,y]=s[lt,i,j]*1im/4
                    V[y,x]=-s[lt,i,j]*1im/4
                end
                E=diag(model.UV[j,:,:]*V*model.UV[j,:,:]')
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

function phy_measure(model,G,lt,s)
    G0=G[:,:]
    
    if lt>model.Nt/2
        for t in lt:-1:div(model.Nt,2)+1
            for j in 1:size(s)[3]
                V=zeros(ComplexF64,model.Ns,model.Ns)
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    V[x,y]=s[t,i,j]*1im/4
                    V[y,x]=-s[t,i,j]*1im/4
                end
                E=diag(model.UV[j,:,:]*V*model.UV[j,:,:]')
                G0=model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:] *G0* model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:]
            end
            G0= model.eKinv*G0*model.eK
        end
    else
        for t in lt+1:div(model.Nt,2)
            G0=model.eK*G0*model.eKinv
            for j in size(s)[3]:-1:1
                V=zeros(ComplexF64,model.Ns,model.Ns)
                for i in 1:size(s)[2]
                    x,y=model.nnidx[i,j]
                    V[x,y]=s[t,i,j]*1im/4
                    V[y,x]=-s[t,i,j]*1im/4
                end
                E=diag(model.UV[j,:,:]*V*model.UV[j,:,:]')
                G0=model.UV[j,:,:]'*diagm(exp.(model.α*E))*model.UV[j,:,:] *G0* model.UV[j,:,:]'*diagm(exp.(-model.α*E))*model.UV[j,:,:]
            end
        end
    end
    #####################################################################
    if norm(G0-Gτ(model,s,div(model.Nt,2)))>1e-6 
        error("record error")
    end
    #####################################################################

    G0=model.HalfeK* G0 *model.HalfeKinv

    Ek=model.t*sqrt(sum(imag.(G0).* model.K))
    Ev=model.U*sqrt(sum(imag.(G0).^2 .*model.K./2 .+1/4))
    R0=R1=0
    if model.Lattice=="HoneyComb"
        for rx in 1:model.site[1]
            for ry in 1:model.site[2]
                tmp=0
                for ix in 1:model.site[1]
                    for iy in 1:model.site[2]
                        idx1=2*( (iy-1)*model.site[1]+ix) -1
                        idx2=2*( mod(iy+ry-1,model.site[2])*model.site[1]+mod1(ix+rx,model.site[1]) )-1
                        tmp+=sqrt(real(G0[idx1,idx1])*real(G0[idx2,idx2])-real(G0[idx1,idx2])*real(G0[idx2,idx1]))
                        tmp+=sqrt(real(G0[idx1+1,idx1+1])*real(G0[idx2+1,idx2+1])-real(G0[idx1+1,idx2+1])*real(G[idx2+1,idx1+1]))
                        tmp-=sqrt(real(G0[idx1+1,idx1+1])*real(G0[idx2,idx2])-real(G0[idx1+1,idx2])*real(G0[idx2,idx1+1]))
                        tmp-=sqrt(real(G0[idx1,idx1])*real(G0[idx2+1,idx2+1])-real(G0[idx1,idx2+1])*real(G0[idx2+1,idx1]))
                    end
                end
                tmp=tmp/prod(model.site)/4
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
                        tmp+=sqrt(real(G0[idx1,idx1])*real(G0[idx2,idx2])-real(G0[idx1,idx2])*real(G0[idx2,idx1]))
                    end
                end
                tmp=tmp/prod(model.site)/4
                R0+=tmp*cos(π*(rx+ry))
                R1+=cos(π*(rx+ry)+2*π/model.site[1]*rx+2*π/model.site[2]*ry )*tmp
            end
        end
    end
    # 1-R1/R0
    return Ek,Ev,R0,R1
end