function phy_update(path::String,model::_Hubbard_Para,s::Array{Int8,2},Sweeps::Int64,record)
    if model.Lattice=="SQUARE"
        name="□"
    elseif model.Lattice=="HoneyComb"
        name="HC"
    end
    file="$(path)CDW_phy$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv"
    
    
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
                for k in 1:size(s)[2]
                    x,y=model.nnidx[k].I
                    D[x]+=s[lt,k]
                    D[y]-=s[lt,k]
                end
                G=diagm(exp.(model.α.*D))*model.eK *G* model.eKinv*diagm(exp.(-model.α.*D))
                
                #####################################################################
                if norm(G-Gτ(model,s,lt))>1e-4 
                    error(norm(G-Gτ(model,s,lt))," $(lt) :asd")
                end
                #####################################################################
            end

            for k in 1:size(s)[2]
                x,y=model.nnidx[k].I
                Δ=diagm(exp.( 2*model.α.*[-s[lt,k],s[lt,k]] ))-I(2)
                subidx=[x,y]
                r=I(2)+Δ*(I(2)-G[subidx,subidx])
                detR=det(r)
                # if detR<0
                #     print("*")
                # else
                #     print("-")
                # end
                ####################################################################
                # ss=s[:,:]
                # ss[lt,k]=-ss[lt,k]
                # if abs(Poss(model,ss)/Poss(model,s)-detR)>1e-3
                #     println("possibility error:",abs(Poss(model,ss)/Poss(model,s)-detR))
                # end
                ####################################################################
                if rand(rng)<abs(detR)
                    G=G-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                    s[lt,k]=-s[lt,k]
                    # print(detR > 0 ? "+" : "-")

                    ####################################################################
                    if norm(G-Gτ(model,s,lt))>1e-4
                        println("G update error:",norm(G-Gτ(model,s,lt)))
                    end
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

        for lt in model.Nt-1:-1:1
            if mod(lt,model.WrapTime)==1
                G=Gτ(model,s,lt)
            else
                D=zeros(model.Ns)
                for k in 1:size(s)[2]
                    x,y=model.nnidx[k].I
                    D[x]+=s[lt+1,k]
                    D[y]-=s[lt+1,k]
                end
                G=model.eKinv* diagm(exp.(-model.α.*D)) *G* diagm(exp.(model.α.*D))*model.eK
                
                #####################################################################
                # if norm(G-Gτ(model,s,lt))>1e-4
                #     println("reversal wrap error: ",norm(G-Gτ(model,s,lt)))
                # end
                #####################################################################
            end

            for k in 1:size(s)[2]
                x,y=model.nnidx[k].I
                Δ=diagm(exp.( 2*model.α.*[-s[lt,k],s[lt,k]] ))-I(2)
                subidx=[x,y]
                r=I(2)+Δ*(I(2)-G[subidx,subidx])
                detR=det(r)
                # if detR<0
                #     print("*")
                # else
                #     print("-")
                # end
                ####################################################################
                # ss=s[:,:]
                # ss[lt,k]=-ss[lt,k]
                # if abs(Poss(model,ss)/Poss(model,s)-detR)>1e-3
                #     println("possibility error:",abs(Poss(model,ss)/Poss(model,s)-detR))
                # end
                ####################################################################
                if rand(rng)<abs(detR)
                    G=G-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
                    s[lt,k]=-s[lt,k]
                    # print(detR > 0 ? "+" : "-")

                    ####################################################################
                    # if norm(G-Gτ(model,s,lt))>1e-4
                    #     println("G update error:",norm(G-Gτ(model,s,lt)))
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
                        # if detR<0
                        #     println("\n")
                        # end
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
        for k in 1:size(s)[2]
            x,y=model.nnidx[k].I
            D[x]+=s[i,k]
            D[y]-=s[i,k]
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
            for k in 1:size(s)[2]
                x,y=model.nnidx[k].I
                D[x]+=s[i,k]
                D[y]-=s[i,k]
            end
            G0= model.eKinv*diagm(exp.(-model.α.*D)) *G0*  diagm(exp.(model.α.*D))*model.eK
        end
    else
        for i in lt+1:div(model.Nt,2)
            D=zeros(model.Ns)
            for k in 1:size(s)[2]
                x,y=model.nnidx[k].I
                D[x]+=s[i,k]
                D[y]-=s[i,k]
            end
            G0=diagm(exp.(model.α.*D))*model.eK *G0* model.eKinv*diagm(exp.(-model.α.*D))  
        end
    end
    #####################################################################
    # if norm(G0-Gτ(model,s,div(model.Nt,2)))>1e-4 
    #     println("record error:",norm(G0-Gτ(model,s,div(model.Nt,2))))
    # end
    #####################################################################

    G0=model.HalfeK* G0 *model.HalfeKinv

    Ek=model.t*sum(model.K.*G0)
    Ev=R0=R1=0
    for k in 1:length(model.nnidx)
        x,y=model.nnidx[k].I
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