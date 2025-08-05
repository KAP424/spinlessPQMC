function phy_update(path::String,model::_Hubbard_Para,s::Array{Int8,2},Sweeps::Int64,record)
    if model.Lattice=="SQUARE"
        name="□"
    elseif model.Lattice=="HoneyComb"
        name="HC"
    end
    file="$(path)M_PHY$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv"

    rng=MersenneTwister(Threads.threadid())
    sg=R0=R1=Ek=Ev=0
    counter=0

    G=zeros(Float64,model.Ns,model.Ns)
    for loop in 1:Sweeps
        for lt in 1:model.Nt
            if mod(lt,model.WrapTime)==1
                G=Gτ(model,s,lt)
            else
                D=model.K[:,:]
                for k in 1:size(s)[2]
                    x,y=model.nnidx[k].I
                    D[x,y]*=s[lt,k]*1im/2
                    D[y,x]*=-s[lt,k]*1im/2
                end
                E,V=eigen(D)
                G=V*diagm(exp.(model.α.*E))*V'*model.eK *G* model.eKinv*V*diagm(exp.(-model.α.*E))*V'
                
                #####################################################################
                if norm(G-Gτ(model,s,lt))>1e-6 
                    error("$(lt):asdasdsa")
                end
                #####################################################################
            end

            # for x in 1:size(s)[2]
            #     for y in 1:size(s)[3]
            #         xidx=2*x-1
            #         yidx=findall(model.K[xidx,:].!=0)[y]
            #         Δ=[0 -1im*s[lt,x,y]; 1im*s[lt,x,y] 0]
            #         E,V=eigen(Δ) 
            #         Δ=V*diagm(exp.(model.α*E))*V'- I(2)
            #         subidx=[xidx,yidx]
            #         r=I(2)+Δ*(I(2)-G[subidx,subidx])
            #         detR=abs2(det(r))
            #         ####################################################################
            #         ss=s[:,:,:]
            #         ss[lt,x,y]=-ss[lt,x,y]
            #         if abs(Poss(model,ss)/Poss(model,s)-detR)>1e-3
            #             println("possibility error")
            #         end
            #         ####################################################################

            #         if rand(rng)<detR
            #             G=G-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
            #             s[lt,x,y]=-s[lt,x,y]
            #             ####################################################################
            #             if norm(G-Gτ(model,s,lt))>1e-6
            #                 println("G update error")
            #             end
            #             #####################################################################
            #         end
            #     end
            # end

            # if abs(lt-model.Nt/2)<=model.WrapTime
            #     G0=G[:,:]
            #     if lt>model.Nt/2
            #         for i in lt:-1:div(model.Nt,2)+1
            #             D=zeros(model.Ns)
            #             for x in 1:size(s)[2]
            #                 xidx=2*x-1
            #                 nnidx=findall(model.K[xidx,:].!=0)
            #                 for k in 1:size(s)[3]
            #                     D[xidx]+=s[i,x,k]
            #                     D[nnidx[k]]-=s[i,x,k]
            #                 end
            #             end
            #             G0= model.eKinv*diagm(exp.(-model.α.*D)) *G0*  diagm(exp.(model.α.*D))*model.eK
            #         end
            #     else
            #         for i in lt+1:div(model.Nt,2)
            #             D=zeros(model.Ns)
            #             for x in 1:size(s)[2]
            #                 xidx=2*x-1
            #                 nnidx=findall(model.K[xidx,:].!=0)
            #                 for k in 1:size(s)[3]
            #                     D[xidx]+=s[i,x,k]
            #                     D[nnidx[k]]-=s[i,x,k]
            #                 end
            #             end
            #             G0=diagm(exp.(model.α.*D))*model.eK *G0* model.eKinv*diagm(exp.(-model.α.*D))  
            #         end
            #     end
            #     #####################################################################
            #     # if norm(G-Gτ(model,s,lt))>1e-6 
            #     #     error("record error")
            #     # end
            #     #####################################################################
            #     G0=model.HalfeK* G0 *model.HalfeKinv
            #     ##------------------------------------------------------------------------
            #     # tmp=Magnetism(model,G0)
            #     # mA+=tmp[1]
            #     # mB+=tmp[2]
            #     # nn+=NN(model,G0)
            #     # Ek+=EK(model,G0)
            #     # tmp=CzzofSpin(model,G0)
            #     # R0+=tmp[1]
            #     # R1+=tmp[2]
            #     # C0+=tmp[3]
            #     # Cmax+=tmp[4]
            #     # counter+=1
            #     ##------------------------------------------------------------------------
            # end
        end

        for lt in model.Nt-1:-1:1
            if mod(lt,model.WrapTime)==1
                G=Gτ(model,s,lt)
            else
                D=model.K[:,:]
                for k in 1:size(s)[2]
                    x,y=model.nnidx[k].I
                    D[x,y]*=s[lt+1,k]*1im/2
                    D[y,x]*=-s[lt+1,k]*1im/2
                end
                E,V=eigen(D)
                G=model.eKinv* V*diagm(exp.(-model.α.*E))*V' *G* V*diagm(exp.(model.α.*E))*V'*model.eK
                
                #####################################################################
                if norm(G-Gτ(model,s,lt))>1e-6
                    error("ltltltl")
                end
                #####################################################################
            end

        #     for x in 1:size(s)[2]
        #         for y in 1:size(s)[3]
        #             xidx=2*x-1
        #             yidx=findall(model.K[xidx,:].!=0)[y]
        #             Δ=diagm(exp.( 2*model.α.*[-s[lt,x,y],s[lt,x,y]] ))-I(2)
        #             subidx=[xidx,yidx]
        #             r=I(2)+Δ*(I(2)-G[subidx,subidx])
        #             detR=det(r)
        #             if detR<0
        #                 println("Warning for negative possibility!")
        #             end
        #             ####################################################################
        #                 # ss=s[:,:,:]
        #                 # ss[lt,x,y]=-ss[lt,x,y]
        #                 # if abs(Poss(model,ss)/Poss(model,s)-detR)>1e-3
        #                 #     error("possibility error")
        #                 # end
        #             ####################################################################

        #             if rand(rng)<detR
        #                 G=G-(G[:,subidx]/r*Δ*((I(model.Ns)-G)[subidx,:]))
        #                 s[lt,x,y]=-s[lt,x,y]
        #                 ####################################################################
        #                 # if norm(G-Gτ(model,s,lt))>1e-6
        #                 #     println("G update error")
        #                 # end
        #                 #####################################################################
        #             end
        #         end
        #     end


        #     if abs(lt-model.Nt/2)<=model.WrapTime
        #         G0=G[:,:]
        #         if lt>model.Nt/2
        #             for i in lt:-1:div(model.Nt,2)+1
        #                 D=zeros(model.Ns)
        #                 for x in 1:size(s)[2]
        #                     xidx=2*x-1
        #                     nnidx=findall(model.K[xidx,:].!=0)
        #                     for k in 1:size(s)[3]
        #                         D[xidx]+=s[i,x,k]
        #                         D[nnidx[k]]-=s[i,x,k]
        #                     end
        #                 end
        #                 G0= model.eKinv*diagm(exp.(-model.α.*D)) *G0*  diagm(exp.(model.α.*D))*model.eK
        #             end
        #         else
        #             for i in lt+1:div(model.Nt,2)
        #                 D=zeros(model.Ns)
        #                 for x in 1:size(s)[2]
        #                     xidx=2*x-1
        #                     nnidx=findall(model.K[xidx,:].!=0)
        #                     for k in 1:size(s)[3]
        #                         D[xidx]+=s[i,x,k]
        #                         D[nnidx[k]]-=s[i,x,k]
        #                     end
        #                 end
        #                 G0=diagm(exp.(model.α.*D))*model.eK *G0* model.eKinv*diagm(exp.(-model.α.*D))  
        #             end
        #         end
        #         #####################################################################
        #         # if norm(G-Gτ(model,s,lt))>1e-6 
        #         #     error("record error")
        #         # end
        #         #####################################################################
        #         G0=model.HalfeK* G0 *model.HalfeKinv
        #         ##------------------------------------------------------------------------
        #         # tmp=Magnetism(model,G0)
        #         # mA+=tmp[1]
        #         # mB+=tmp[2]
        #         # nn+=NN(model,G0)
        #         # Ek+=EK(model,G0)
        #         # tmp=CzzofSpin(model,G0)
        #         # R0+=tmp[1]
        #         # R1+=tmp[2]
        #         # C0+=tmp[3]
        #         # Cmax+=tmp[4]
        #         # counter+=1
        #         ##------------------------------------------------------------------------
        #     end

        # end
        # if loop>WarmSweeps
        #     fid = open("$(path)Phy$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv", "a+")
        #     writedlm(fid,[Ek model.U*nn nn mA mB R0 R1 C0 Cmax]/counter,',')
        #     close(fid)
        #     mA=mB=nn=R0=R1=Ek=C0=Cmax=0
        #     counter=0
        end
    end
    return s
end 

function Poss(model,s)
    A=model.Pt[:,:]

    for i in 1:model.Nt
        D=zeros(ComplexF64,model.Ns,model.Ns)
        for k in 1:size(s)[2]
            x,y=model.nnidx[k].I
            D[x,y]=s[i,k]*1im/2
            D[y,x]=-s[i,k]*1im/2
        end
        E,V=eigen(D)
        A=V*diagm(exp.(model.α.*E))*V'*model.eK*A
    end
    A=model.Pt'*A
    return abs2(det(A))
end

