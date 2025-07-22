function phy_update(path::String,model::_Hubbard_Para,s::Array{Int8,3},Sweeps::Int64,record)
    if model.Lattice=="SQUARE"
        name="□"
    elseif model.Lattice=="HoneyComb"
        name="HC"
    end
    file="$(path)tVPHY$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv"
    
    
    rng=MersenneTwister(Threads.threadid())
    sg=R0=R1=Ek=Ev=0
    counter=0

    G=zeros(Float64,model.Ns,model.Ns)
    for loop in 1:Sweeps
        for lt in 1:model.Nt
            if mod(lt,model.WrapTime)==1
                G=Gτ(model,s,lt)
            else
                D=zeros(model.Ns)
                for x in 1:size(s)[2]
                    xidx=2*x-1
                    nnidx=findall(model.K[xidx,:].!=0)
                    for k in 1:size(s)[3]
                        D[xidx]+=s[lt,x,k]
                        D[nnidx[k]]-=s[lt,x,k]
                    end
                end
                G=diagm(exp.(model.α.*D))*model.eK *G* model.eKinv*diagm(exp.(-model.α.*D))
                
                #####################################################################
                # if norm(G-Gτ(model,s,lt))>1e-4 
                #     error(norm(G-Gτ(model,s,lt))," $(lt) :asd")
                # end
                #####################################################################
            end

            for x in 1:size(s)[2]
                for y in 1:size(s)[3]
                    xidx=2*x-1
                    yidx=findall(model.K[xidx,:].!=0)[y]
                    Δ=diagm(exp.( 2*model.α.*[-s[lt,x,y],s[lt,x,y]] ))-I(2)
                    subidx=[xidx,yidx]
                    r=I(2)+Δ*(I(2)-G[subidx,subidx])
                    detR=det(r)
                    if detR<0
                        println("Warning for negative possibility!")
                    end
                    ####################################################################
                    # ss=s[:,:,:]
                    # ss[lt,x,y]=-ss[lt,x,y]
                    # if abs(Poss(model,ss)/Poss(model,s)-detR)>1e-3
                    #     error("possibility error")
                    # end
                    ####################################################################
                    if rand(rng)<abs(detR)
                        G=G-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                        s[lt,x,y]=-s[lt,x,y]

                        ####################################################################
                        # if norm(G-Gτ(model,s,lt))>1e-4
                        #     error("G update error:",G-Gτ(model,s,lt))
                        # end
                        #####################################################################

                        if record && abs(lt-model.Nt/2)<=model.WrapTime
                            tmp=phy_measure(model,G,lt,s).*sign(detR)
                            Ek+=tmp[1]
                            Ev+=tmp[2]
                            R0+=tmp[3]
                            R1+=tmp[4]
                            counter+=1
                            sg+=sign(detR)
                        end
                    end
                end
            end
        end

        for lt in model.Nt-1:-1:1
            if mod(lt,model.WrapTime)==1
                G=Gτ(model,s,lt)
            else
                D=zeros(model.Ns)
                for x in 1:size(s)[2]
                    xidx=2*x-1
                    nnidx=findall(model.K[xidx,:].!=0)
                    for k in 1:size(s)[3]
                        D[xidx]+=s[lt+1,x,k]
                        D[nnidx[k]]-=s[lt+1,x,k]
                    end
                end
                G=model.eKinv* diagm(exp.(-model.α.*D)) *G* diagm(exp.(model.α.*D))*model.eK
                
                #####################################################################
                # if norm(G-Gτ(model,s,lt))>1e-4
                #     error(norm(G-Gτ(model,s,lt)),"ltltltl")
                # end
                #####################################################################
            end

            for x in 1:size(s)[2]
                for y in 1:size(s)[3]
                    xidx=2*x-1
                    yidx=findall(model.K[xidx,:].!=0)[y]
                    Δ=diagm(exp.( 2*model.α.*[-s[lt,x,y],s[lt,x,y]] ))-I(2)
                    subidx=[xidx,yidx]
                    r=I(2)+Δ*(I(2)-G[subidx,subidx])
                    detR=det(r)
                    if detR<0
                        println("Warning for negative possibility!")
                    end
                    ####################################################################
                    # ss=s[:,:,:]
                    # ss[lt,x,y]=-ss[lt,x,y]
                    # if abs(Poss(model,ss)/Poss(model,s)-detR)>1e-3
                    #     error("possibility error")
                    # end
                    ####################################################################
                    if rand(rng)<abs(detR)
                        G=G-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                        s[lt,x,y]=-s[lt,x,y]

                        ####################################################################
                        # if norm(G-Gτ(model,s,lt))>1e-4
                        #     error("G update error:",G-Gτ(model,s,lt))
                        # end
                        #####################################################################

                        if record && abs(lt-model.Nt/2)<=model.WrapTime
                            tmp=phy_measure(model,G,lt,s).*sign(detR)
                            Ek+=tmp[1]
                            Ev+=tmp[2]
                            R0+=tmp[3]
                            R1+=tmp[4]
                            counter+=1
                            sg+=sign(detR)
                        end
                    end
                end
            end
        end

        if record
            open(file, "a") do io
                lock(io)
                writedlm(io,[sg Ek Ev R0 R1]/counter,',')
                unlock(io)
            end
            sg=R0=R1=Ek=Ev=0
            counter=0
        end
    end
    return s
end 

function Poss(model,s)
    A=model.Pt[:,:]

    for i in 1:model.Nt
        D=zeros(model.Ns)
        for x in 1:size(s)[2]
            xidx=2*x-1
            nnidx=findall(model.K[xidx,:].!=0)
            for k in 1:size(s)[3]
                D[xidx]+=s[i,x,k]
                D[nnidx[k]]-=s[i,x,k]
            end
        end
        A=diagm(exp.(model.α.*D))*model.eK*A
    end
    A=model.Pt'*A
    return det(A)
end

function phy_measure(model,G,lt,s)
    G0=G[:,:]
    if lt>model.Nt/2
        for i in lt:-1:div(model.Nt,2)+1
            D=zeros(model.Ns)
            for x in 1:size(s)[2]
                xidx=2*x-1
                nnidx=findall(model.K[xidx,:].!=0)
                for k in 1:size(s)[3]
                    D[xidx]+=s[i,x,k]
                    D[nnidx[k]]-=s[i,x,k]
                end
            end
            G0= model.eKinv*diagm(exp.(-model.α.*D)) *G0*  diagm(exp.(model.α.*D))*model.eK
        end
    else
        for i in lt+1:div(model.Nt,2)
            D=zeros(model.Ns)
            for x in 1:size(s)[2]
                xidx=2*x-1
                nnidx=findall(model.K[xidx,:].!=0)
                for k in 1:size(s)[3]
                    D[xidx]+=s[i,x,k]
                    D[nnidx[k]]-=s[i,x,k]
                end
            end
            G0=diagm(exp.(model.α.*D))*model.eK *G0* model.eKinv*diagm(exp.(-model.α.*D))  
        end
    end
    #####################################################################
    # if norm(G0-Gτ(model,s,div(model.Nt,2)))>1e-4 
    #     error("record error")
    # end
    #####################################################################

    G0=model.HalfeK* G0 *model.HalfeKinv

    Ek=model.t*sum(model.K.*G0)
    Ev=R0=R1=0
    for x in 1:div(model.Ns,2)
        xidx=2*x-1
        nnidx=findall(model.K[xidx,:].!=0)
        for y in eachindex(nnidx)
            Ev+=(1-G0[xidx,xidx])*(1-G0[nnidx[y],nnidx[y]])-G0[xidx,nnidx[y]]*G0[nnidx[y],xidx]-1/4
        end
    end
    Ev*=model.U

    for rx in 1:model.site[1]
        for ry in 1:model.site[2]
            tmp=0
            for ix in 1:model.site[1]
                for iy in 1:model.site[2]
                    idx1=2*( mod(iy-1,model.site[2])*model.site[1]+mod1(ix,model.site[1]) )-1
                    idx2=2*( mod(iy+ry-2,model.site[2])*model.site[1]+mod1(ix+rx,model.site[1]) )-1
                    tmp+=(1-G0[idx1,idx1])*(1-G0[idx2,idx2])-G0[idx1,idx2]*G0[idx2,idx1]
                    tmp+=(1-G0[idx1+1,idx1+1])*(1-G0[idx2+1,idx2+1])-G0[idx1+1,idx2+1]*G[idx2+1,idx1+1]
                    tmp-=(1-G0[idx1+1,idx1+1])*(1-G0[idx2,idx2])-G0[idx1+1,idx2]*G0[idx2,idx1+1]
                    tmp-=(1-G0[idx1,idx1])*(1-G0[idx2+1,idx2+1])-G0[idx1,idx2+1]*G0[idx2+1,idx1]
                end
            end
            tmp/=prod(model.site)
            R0+=tmp
            R1+=cos(π/model.site[1]*(rx+ry))*tmp
        end
    end

    # 1-R1/R0
    return Ek,Ev,R0,R1
end