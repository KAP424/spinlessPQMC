function phy_update(path::String,model::_Hubbard_Para,WarmSweeps::Int64,Sweeps::Int64,s::Array{UInt8,2})
    if model.Lattice=="SQUARE"
        name="□"
    elseif model.Lattice=="HoneyComb"
        name="HC"
    end
    rng=MersenneTwister(Threads.threadid())
    elements=(1, 2, 3, 4)

    mA=mB=nn=R0=R1=Ek=C0=Cmax=0

    counter=0

    G=zeros(ComplexF64,model.Ns,model.Ns)
    for loop in 1:Sweeps+WarmSweeps
        for lt in 1:model.Nt
            if mod(lt,model.WrapTime)==1
                G=Gτ(model,s,lt)
            else
                D=[model.η[x] for x in s[:,lt]]
                G=diagm(exp.(1im*model.α.*D))*model.eK *G* model.eKinv*diagm(exp.(-1im*model.α.*D))
                
            #####################################################################
                # if norm(G-Gτ(model,s,lt))>1e-6 
                #         error("asd")
                #     end
            #####################################################################
            end

            for x in 1:model.Ns
                sp=Random.Sampler(rng,[i for i in elements if i != s[x,lt]])
                sx=rand(rng,sp)
                Δ=exp(1im*model.α*(model.η[sx]-model.η[s[x,lt]]))-1
                r=1+Δ*(1-G[x,x])

            ####################################################################
                # ss=s[:,:]
                # ss[x,lt]=sx
                # if abs(Poss(model,ss)/Poss(model,s)-abs2(r))>1e-3
                #     error("sada")
                # end
            ####################################################################

                if rand(rng)<model.γ[sx]/model.γ[s[x,lt]]*abs2(r)
                    G=G-Δ/r.*(  G[:,x]    .*  transpose((I(model.Ns)-G)[x,:])   )
                    s[x,lt]=sx
            ####################################################################
                    # if norm(G-Gτ(model,s,lt))>1e-6
                    #     error("asd")
                    # end
            #####################################################################

                end
            end

            if loop>WarmSweeps && abs(lt-model.Nt/2)<=model.WrapTime
                G0=G[:,:]
                if lt>model.Nt/2
                    for i in lt:-1:div(model.Nt,2)+1
                        D=[model.η[x] for x in s[:,i]]
                        G0= model.eKinv*diagm(exp.(-1im*model.α.*D)) *G0*  diagm(exp.(1im*model.α.*D))*model.eK
                    end
                else
                    for i in lt+1:div(model.Nt,2)
                        D=[model.η[x] for x in s[:,i]]
                        G0=diagm(exp.(1im*model.α.*D))*model.eK *G0* model.eKinv*diagm(exp.(-1im*model.α.*D))  
                    end
                end
                #####################################################################
                # if norm(G-Gτ(model,s,lt))>1e-6 
                #     error("asd")
                # end
                #####################################################################
                G0=model.HalfeK* G0 *model.HalfeKinv
            
                ##------------------------------------------------------------------------
                tmp=Magnetism(model,G0)
                mA+=tmp[1]
                mB+=tmp[2]
                nn+=NN(model,G0)
                Ek+=EK(model,G0)
                tmp=CzzofSpin(model,G0)
                R0+=tmp[1]
                R1+=tmp[2]
                C0+=tmp[3]
                Cmax+=tmp[4]
                counter+=1
                ##------------------------------------------------------------------------
            end
        end

        for lt in model.Nt-1:-1:1
            if mod(lt,model.WrapTime)==1
                G=Gτ(model,s,lt)
            else
                D=[model.η[x] for x in s[:,lt+1]]
                G=model.eKinv*diagm(exp.(-1im*model.α.*D)) *G* diagm(exp.(1im*model.α.*D))*model.eK 
                
                #####################################################################
                # if norm(G-Gτ(model,s,lt))>1e-6
                #     error("asd")
                # end
                #####################################################################
            end

            for x in 1:model.Ns
                sp=Random.Sampler(rng,[i for i in elements if i != s[x,lt]])
                sx=rand(rng,sp)
                Δ=exp(1im*model.α*(model.η[sx]-model.η[s[x,lt]]))-1
                r=1+Δ*(1-G[x,x])

                if rand(rng)<model.γ[sx]/model.γ[s[x,lt]]*abs2(r)
                    G=G-Δ/r.*(  G[:,x]    .*  transpose((I(model.Ns)-G)[x,:])   )
                    s[x,lt]=sx
            ####################################################################
                    # if norm(G-Gτ(model,s,lt))>1e-6
                    #     error("asd")
                    # end
            #####################################################################

                end
            end

            if loop>WarmSweeps && abs(lt-model.Nt/2)<=model.WrapTime
                G0=G[:,:]
                if lt>model.Nt/2
                    for i in lt:-1:div(model.Nt,2)+1
                        D=[model.η[x] for x in s[:,i]]
                        G0= model.eKinv*diagm(exp.(-1im*model.α.*D)) *G0*  diagm(exp.(1im*model.α.*D))*model.eK
                    end
                else
                    for i in lt+1:div(model.Nt,2)
                        D=[model.η[x] for x in s[:,i]]
                        G0=diagm(exp.(1im*model.α.*D))*model.eK *G0* model.eKinv*diagm(exp.(-1im*model.α.*D))  
                    end
                end
                #####################################################################
                # if norm(G-Gτ(model,s,lt))>1e-6
                #     error("asd")
                # end
                #####################################################################
                G0=model.HalfeK* G0 *model.HalfeKinv
            
                ##------------------------------------------------------------------------
                tmp=Magnetism(model,G0)
                mA+=tmp[1]
                mB+=tmp[2]
                nn+=NN(model,G0)
                Ek+=EK(model,G0)
                tmp=CzzofSpin(model,G0)
                R0+=tmp[1]
                R1+=tmp[2]
                C0+=tmp[3]
                Cmax+=tmp[4]
                counter+=1
                ##------------------------------------------------------------------------
            end

        end
        if loop>WarmSweeps
            fid = open("$(path)Phy$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv", "a+")
            writedlm(fid,[Ek model.U*nn nn mA mB R0 R1 C0 Cmax]/counter,',')
            close(fid)
            mA=mB=nn=R0=R1=Ek=C0=Cmax=0
            counter=0
        end
    end
    return s
end 

# function Poss(model,s)
#     A=model.Pt[:,:]

#     for i in 1:model.Nt
#         D=[model.η[x] for x in s[:,i]]
#         A=diagm(exp.(1im*model.α.*D))*model.eK*A
#     end
#     A=model.Pt'*A

#     ans=det(A)

#     return abs2(ans)
    
# end
