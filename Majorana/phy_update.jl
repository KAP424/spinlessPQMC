function phy_update(path::String,model::_Hubbard_Para,s::Array{Int8,3},Sweeps::Int64,record)
    if model.Lattice=="SQUARE"
        name="□"
    elseif model.Lattice=="HoneyComb"
        name="HC"
    end
            # if loop>WarmSweeps
        #     fid = open("$(path)Phy$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv", "a+")
        #     writedlm(fid,[Ek model.U*nn nn mA mB R0 R1 C0 Cmax]/counter,',')
        #     close(fid)
        #     mA=mB=nn=R0=R1=Ek=C0=Cmax=0
        #     counter=0
        # end
    file="$(path)tVPHY$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv"
    atexit() do
        if record
            open(file, "a") do io
                lock(io)
                writedlm(io,[Ek model.U*nn nn mA mB R0 R1 C0 Cmax]/counter,',')
                unlock(io)
            end
        end
        # writedlm("$(path)s/S$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ).csv", s,",")
    end

    rng=MersenneTwister(Threads.threadid())
    mA=mB=nn=R0=R1=Ek=C0=Cmax=0
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
                # if norm(G-Gτ(model,s,lt))>1e-6 
                #     error("$(lt):asd")
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

                    if rand(rng)<detR
                        G=G-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                        s[lt,x,y]=-s[lt,x,y]
                        ####################################################################
                        # if norm(G-Gτ(model,s,lt))>1e-6
                        #     println("G update error")
                        # end
                        #####################################################################
                    end
                end
            end

            if abs(lt-model.Nt/2)<=model.WrapTime
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
                # if norm(G-Gτ(model,s,lt))>1e-6 
                #     error("record error")
                # end
                #####################################################################
                G0=model.HalfeK* G0 *model.HalfeKinv
                ##------------------------------------------------------------------------
                # tmp=Magnetism(model,G0)
                # mA+=tmp[1]
                # mB+=tmp[2]
                # nn+=NN(model,G0)
                # Ek+=EK(model,G0)
                # tmp=CzzofSpin(model,G0)
                # R0+=tmp[1]
                # R1+=tmp[2]
                # C0+=tmp[3]
                # Cmax+=tmp[4]
                # counter+=1
                ##------------------------------------------------------------------------
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
                # if norm(G-Gτ(model,s,lt))>1e-6
                #     error("ltltltl")
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

                    if rand(rng)<detR
                        G=G-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                        s[lt,x,y]=-s[lt,x,y]
                        ####################################################################
                        # if norm(G-Gτ(model,s,lt))>1e-6
                        #     println("G update error")
                        # end
                        #####################################################################
                    end
                end
            end


            if abs(lt-model.Nt/2)<=model.WrapTime
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
                # if norm(G-Gτ(model,s,lt))>1e-6 
                #     error("record error")
                # end
                #####################################################################
                G0=model.HalfeK* G0 *model.HalfeKinv
                ##------------------------------------------------------------------------
                # tmp=Magnetism(model,G0)
                # mA+=tmp[1]
                # mB+=tmp[2]
                # nn+=NN(model,G0)
                # Ek+=EK(model,G0)
                # tmp=CzzofSpin(model,G0)
                # R0+=tmp[1]
                # R1+=tmp[2]
                # C0+=tmp[3]
                # Cmax+=tmp[4]
                # counter+=1
                ##------------------------------------------------------------------------
            end

        end
        # if loop>WarmSweeps
        #     fid = open("$(path)Phy$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv", "a+")
        #     writedlm(fid,[Ek model.U*nn nn mA mB R0 R1 C0 Cmax]/counter,',')
        #     close(fid)
        #     mA=mB=nn=R0=R1=Ek=C0=Cmax=0
        #     counter=0
        # end
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

