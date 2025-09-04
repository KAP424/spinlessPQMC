# attractive-U and repulsive-U get the same S_2

function ctrl_SCEEicr(path::String,model::_Hubbard_Para,indexA::Vector{Int64},indexB::Vector{Int64},Sweeps::Int64,λ::Float64,Nλ::Int64,ss::Vector{Array{Int8,3}},record)
    uv=[-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]
    if model.Lattice=="SQUARE"
        name="□"
    elseif model.Lattice=="HoneyComb"
        name="HC"
    end
    file="$(path)SCEEicr$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)N$(Nλ)BS$(model.BatchSize).csv"
    
    atexit() do
        if record
            open(file, "a") do io
                lock(io)
                writedlm(io, O', ',')
                unlock(io)
            end
        end
        writedlm("$(path)ss/SS$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)λ$(Int(round(Nλ*λ))).csv", [ss[1] ss[2]],",")
    end
    
    
    rng=MersenneTwister(Threads.threadid()+round(Int,time()*1000))
    elements=(1, 2, 3, 4)

    Gt1=zeros(ComplexF64,model.Ns,model.Ns)
    Gt2=zeros(ComplexF64,model.Ns,model.Ns)
    G01=zeros(ComplexF64,model.Ns,model.Ns)
    G02=zeros(ComplexF64,model.Ns,model.Ns)
    Gt01=zeros(ComplexF64,model.Ns,model.Ns)
    Gt02=zeros(ComplexF64,model.Ns,model.Ns)
    G0t1=zeros(ComplexF64,model.Ns,model.Ns)
    G0t2=zeros(ComplexF64,model.Ns,model.Ns)
    gmInv_A=zeros(ComplexF64,length(indexA),length(indexA))
    gmInv_B=zeros(ComplexF64,length(indexB),length(indexB))
    detg_A=detg_B=0

    tmpO=0
    counter=0
    O=zeros(Sweeps+1)
    O[1]=λ

    II=I(model.Ns)
    IA=I(length(indexA))
    IB=I(length(indexB))

    for loop in 1:Sweeps
        for lt in 1:model.Nt
            if lt==div(model.Nt,2)+1
                # 更新Θ+1层时候，得到的是Θ层的四Green函数，再将eK,eV往上安装同时更新
                # 而Θ层的四Green函数是和1-Θ层同用一个Wrap规则
                # 所以为了让Θ层的四Green函数与Θ-2Θ同用一个Wrap规则，需要将Θ层的四Green函数的最后两个进行倒换
                Gt1,G01,G0t1,Gt01=G4(model,ss[1],lt-1,div(model.Nt,2))
                Gt2,G02,G0t2,Gt02=G4(model,ss[2],lt-1,div(model.Nt,2))
                Gt1=model.eK*Gt1*model.eKinv
                Gt2=model.eK*Gt2*model.eKinv
                Gt01=model.eK*Gt01
                Gt02=model.eK*Gt02
                G0t1=G0t1*model.eKinv
                G0t2=G0t2*model.eKinv

                GM_A=GroverMatrix(G01[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
                gmInv_A=inv(GM_A)
                GM_B=GroverMatrix(G01[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
                gmInv_B=inv(GM_B)
                detg_A=det(GM_A)
                detg_B=det(GM_B)
            elseif  mod1(lt,model.WrapTime)==1 
                Gt1,G01,Gt01,G0t1=G4(model,ss[1],lt-1,div(model.Nt,2))
                Gt2,G02,Gt02,G0t2=G4(model,ss[2],lt-1,div(model.Nt,2))
                Gt1=model.eK*Gt1*model.eKinv
                Gt2=model.eK*Gt2*model.eKinv
                Gt01=model.eK*Gt01
                Gt02=model.eK*Gt02
                G0t1=G0t1*model.eKinv
                G0t2=G0t2*model.eKinv

                GM_A=GroverMatrix(G01[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
                gmInv_A=inv(GM_A)
                GM_B=GroverMatrix(G01[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
                gmInv_B=inv(GM_B)
                detg_A=det(GM_A)
                detg_B=det(GM_B)
            else
                #####################################################################
                Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2))
                Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2))
                    
                if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
                    println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
                    error("$lt : WrapTime")
                end
                #####################################################################
                Gt1=model.eK*Gt1*model.eKinv
                Gt2=model.eK*Gt2*model.eKinv
                Gt01=model.eK*Gt01
                Gt02=model.eK*Gt02
                G0t1=G0t1*model.eKinv
                G0t2=G0t2*model.eKinv
            end

            for j in size(ss[1])[3]:-1:1
                # 装载eV eK
                E1=zeros(model.Ns)
                E2=zeros(model.Ns)
                for i in 1:size(ss[1])[2]
                    x,y=model.nnidx[i,j]
                    E1[x]=ss[1][lt,i,j]
                    E1[y]=-ss[1][lt,i,j]
                    E2[x]=ss[2][lt,i,j]
                    E2[y]=-ss[2][lt,i,j]
                end
                Gt1=model.UV[j,:,:]'*diagm(exp.(model.α*E1))*model.UV[j,:,:] *Gt1* model.UV[j,:,:]'*diagm(exp.(-model.α*E1))*model.UV[j,:,:]
                Gt2=model.UV[j,:,:]'*diagm(exp.(model.α*E2))*model.UV[j,:,:] *Gt2* model.UV[j,:,:]'*diagm(exp.(-model.α*E2))*model.UV[j,:,:]
                Gt01=model.UV[j,:,:]'*diagm(exp.(model.α*E1))*model.UV[j,:,:]*Gt01
                Gt02=model.UV[j,:,:]'*diagm(exp.(model.α*E2))*model.UV[j,:,:]*Gt02
                G0t1=G0t1*model.UV[j,:,:]'*diagm(exp.(-model.α*E1))*model.UV[j,:,:]
                G0t2=G0t2*model.UV[j,:,:]'*diagm(exp.(-model.α*E2))*model.UV[j,:,:]

                for i in 1:size(ss[1])[2]
                    # update  
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]

                    # 更新s1
                    E=[-2*ss[1][lt,i,j] , 2*ss[1][lt,i,j]]
                    Δ1=uv*diagm(exp.(model.α.*E))*uv'-I(2)
                    r1=I(2)+Δ1*(I(2)-Gt1[subidx,subidx])

                    b_A=(Gt01[subidx,indexA[:]]) *(2*G02[indexA[:],indexA[:]]-IA)*gmInv_A
                    a_A=G0t1[indexA[:],subidx]
                    Tau_A=b_A*a_A
                    
                    b_B=(Gt01[subidx,indexB[:]]) *(2*G02[indexB[:],indexB[:]]-IB)*gmInv_B
                    a_B=G0t1[indexB[:],subidx]
                    Tau_B=b_B*a_B

                    p=det(r1+Δ1*Tau_A)^λ*det(r1+Δ1*Tau_B)^(1-λ)

                    if p<0
                        println("Negative Sign: $(p)")
                    end

                    if rand(rng)<p
                        rho_A=inv(r1+Δ1*Tau_A)*Δ1
                        gmInv_A-=gmInv_A* a_A * rho_A * b_A
                        detg_A*=det(I(2)+inv(r1)*Δ1*Tau_A)
    
                        rho_B=inv(r1+Δ1*Tau_B)*Δ1
                        gmInv_B-=gmInv_B* ( a_B * rho_B * b_B)
                        detg_B*=det(I(2)+inv(r1)*Δ1*Tau_B)

                        G01+= (G0t1[:,subidx] /r1*Δ1 *(Gt01[subidx,:]))
                        Gt01+=(Gt1[:,subidx] /r1*Δ1 *(Gt01[subidx,:]))
                        G0t1-=(G0t1[:,subidx] /r1*Δ1  * (II-Gt1)[subidx,:]  )
                        Gt1-= (Gt1[:,subidx] /r1*Δ1 *((II-Gt1)[subidx,:]) )         
                        ss[1][lt,i,j]=-ss[1][lt,i,j]
                        #####################################################################
                        print('-')
                        if lt==div(model.Nt,2)+1
                            Gt1_,G01_,G0t1_,Gt01_=G4(model,ss[1],lt-1,div(model.Nt,2))
                        else
                            Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2))
                        end
                        Gt1_=model.eK*Gt1_*model.eKinv
                        Gt01_=model.eK*Gt01_
                        G0t1_=G0t1_*model.eKinv
                        GM_A_=GroverMatrix(G01_[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
                        gmInv_A_=inv(GM_A_)
                        GM_B_=GroverMatrix(G01_[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
                        gmInv_B_=inv(GM_B_)
                        detg_A_=det(GM_A_)
                        detg_B_=det(GM_B_)

                        for jj in size(ss[1])[3]:-1:j
                            E=zeros(model.Ns)
                            for ii in 1:size(ss[1])[2]
                                x,y=model.nnidx[ii,jj]
                                E[x]=ss[1][lt,ii,jj]
                                E[y]=-ss[1][lt,ii,jj]
                            end
                            Gt1_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *Gt1_* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                            Gt01_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:]*Gt01_
                            G0t1_=G0t1_*model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        end
    
                        if norm(Gt1-Gt1_)+norm(G01-G01_)+norm(Gt01-Gt01_)+norm(G0t1-G0t1_)+
                           norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                            println('\n',norm(Gt1-Gt1_),'\n',norm(G01-G01_),'\n',norm(Gt01-Gt01_),'\n',norm(G0t1-G0t1_))
                            println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
                            error("s1:  $lt  $j:,,,asdasdasd")
                        end
                        #####################################################################
                    end

                    # 更新s2
                    E=[-2*ss[2][lt,i,j] , 2*ss[2][lt,i,j]]
                    Δ2=uv*diagm(exp.(model.α.*E))*uv'-I(2)
                    r2=I(2)+Δ2*(I(2)-Gt2[subidx,subidx])

                    b_A=Gt02[subidx,indexA[:]]*gmInv_A
                    a_A=(2*G01[indexA[:],indexA[:]]-IA)*G0t2[indexA[:],subidx]
                    Tau_A=b_A*a_A

                    b_B=Gt02[subidx,indexB[:]]*gmInv_B
                    a_B=(2*G01[indexB[:],indexB[:]]-IB)*G0t2[indexB[:],subidx]
                    Tau_B=b_B*a_B
                    
                    p=det(r2+Δ2*Tau_A)^λ*det(r2+Δ2*Tau_B)^(1-λ)

                    if p<0
                        println("Negative Sign: $(p)")
                    end

                    if rand(rng)<p
                        rho_A=inv(r2+Δ2*Tau_A)*Δ2
                        gmInv_A-=gmInv_A* a_A * rho_A * b_A
                        detg_A*=det(I(2)+inv(r2)*Δ2*Tau_A) 
    
                        rho_B=inv(r2+Δ2*Tau_B)*Δ2
                        gmInv_B-=gmInv_B* a_B * rho_B * b_B
                        detg_B*=det(I(2)+inv(r2)*Δ2*Tau_B)

                        G02+= (G0t2[:,subidx] /r2*Δ2 *(Gt02[subidx,:]))
                        Gt02+=(Gt2[:,subidx] /r2*Δ2 *(Gt02[subidx,:]))
                        G0t2-=(G0t2[:,subidx] /r2*Δ2  * (II-Gt2)[subidx,:]  )
                        Gt2-= (Gt2[:,subidx] /r2*Δ2 *((II-Gt2)[subidx,:]) )         
                        ss[2][lt,i,j]=-ss[2][lt,i,j]
                        #####################################################################
                        print('*')
                        if lt==div(model.Nt,2)+1
                            Gt2_,G02_,G0t2_,Gt02_=G4(model,ss[2],lt-1,div(model.Nt,2))
                        else
                            Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2))
                        end
                        Gt2_=model.eK*Gt2_*model.eKinv
                        Gt02_=model.eK*Gt02_
                        G0t2_=G0t2_*model.eKinv
                        GM_A_=GroverMatrix(G01[indexA[:],indexA[:]],G02_[indexA[:],indexA[:]])
                        gmInv_A_=inv(GM_A_)
                        GM_B_=GroverMatrix(G01[indexB[:],indexB[:]],G02_[indexB[:],indexB[:]])
                        gmInv_B_=inv(GM_B_)
                        detg_A_=det(GM_A_)
                        detg_B_=det(GM_B_)

                        for jj in size(ss[1])[3]:-1:j
                            E=zeros(model.Ns)
                            for ii in 1:size(ss[1])[2]
                                x,y=model.nnidx[ii,jj]
                                E[x]=ss[2][lt,ii,jj]
                                E[y]=-ss[2][lt,ii,jj]
                            end
                            Gt2_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *Gt2_* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                            Gt02_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:]*Gt02_
                            G0t2_=G0t2_*model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        end
    
                        if norm(Gt2-Gt2_)+norm(G02-G02_)+norm(Gt02-Gt02_)+norm(G0t2-G0t2_)+
                           norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                            println('\n',norm(Gt2-Gt2_),'\n',norm(G02-G02_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t2-G0t2_))
                            println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
                            error("s2:  $lt  $x:,,,asdasdasd")
                        end
                        #####################################################################
                    end
                end
            end
            ##------------------------------------------------------------------------

            tmpO+=(detg_A/detg_B)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
        end
        println("\n #####################################################################")
        for lt in model.Nt-1:-1:1
            if mod1(model.Nt-lt,model.WrapTime)==1 || lt==div(model.Nt,2)
                Gt1,G01,Gt01,G0t1=G4(model,ss[1],lt,div(model.Nt,2))
                Gt2,G02,Gt02,G0t2=G4(model,ss[2],lt,div(model.Nt,2))

                GM_A=GroverMatrix(G01[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
                gmInv_A=inv(GM_A)
                GM_B=GroverMatrix(G01[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
                gmInv_B=inv(GM_B)
                detg_A=det(GM_A)
                detg_B=det(GM_B)
            else
                #####################################################################
                Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt,div(model.Nt,2))
                Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt,div(model.Nt,2))
                    
                if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
                    println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
                    error("$lt : WrapTime")
                end
                #####################################################################
            end

            for j in 1:size(ss[1])[3]
                for i in 1:size(ss[1])[2]
                    # update
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]

                    # 更新s1
                    E=[-2*ss[1][lt,i,j] , 2*ss[1][lt,i,j]]
                    Δ1=uv*diagm(exp.(model.α.*E))*uv'-I(2)
                    r1=I(2)+Δ1*(I(2)-Gt1[subidx,subidx])

                    b_A=(Gt01[subidx,indexA[:]]) *(2*G02[indexA[:],indexA[:]]-IA)*gmInv_A
                    a_A=G0t1[indexA[:],subidx]
                    Tau_A=b_A*a_A
                    
                    b_B=(Gt01[subidx,indexB[:]]) *(2*G02[indexB[:],indexB[:]]-IB)*gmInv_B
                    a_B=G0t1[indexB[:],subidx]
                    Tau_B=b_B*a_B

                    p=det(r1+Δ1*Tau_A)^λ*det(r1+Δ1*Tau_B)^(1-λ)

                    if p<0
                        println("Negative Sign: $(p)")
                    end

                    if rand(rng)<p
                        rho_A=inv(r1+Δ1*Tau_A)*Δ1
                        gmInv_A-=gmInv_A* a_A * rho_A * b_A
                        detg_A*=det(I(2)+inv(r1)*Δ1*Tau_A)
    
                        rho_B=inv(r1+Δ1*Tau_B)*Δ1
                        gmInv_B-=gmInv_B* ( a_B * rho_B * b_B)
                        detg_B*=det(I(2)+inv(r1)*Δ1*Tau_B)

                        G01+= (G0t1[:,subidx] /r1*Δ1 *(Gt01[subidx,:]))
                        Gt01+=(Gt1[:,subidx] /r1*Δ1 *(Gt01[subidx,:]))
                        G0t1-=(G0t1[:,subidx] /r1*Δ1  * (II-Gt1)[subidx,:]  )
                        Gt1-= (Gt1[:,subidx] /r1*Δ1 *((II-Gt1)[subidx,:]) )         
                        ss[1][lt,i,j]=-ss[1][lt,i,j]
                        #####################################################################
                        print('-')
                        if lt==div(model.Nt,2)+1
                            Gt1_,G01_,G0t1_,Gt01_=G4(model,ss[1],lt-1,div(model.Nt,2))
                        else
                            Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2))
                        end
                        Gt1_=model.eK*Gt1_*model.eKinv
                        Gt01_=model.eK*Gt01_
                        G0t1_=G0t1_*model.eKinv
                        GM_A_=GroverMatrix(G01_[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
                        gmInv_A_=inv(GM_A_)
                        GM_B_=GroverMatrix(G01_[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
                        gmInv_B_=inv(GM_B_)
                        detg_A_=det(GM_A_)
                        detg_B_=det(GM_B_)

                        for jj in size(ss[1])[3]:-1:j
                            E=zeros(model.Ns)
                            for ii in 1:size(ss[1])[2]
                                x,y=model.nnidx[ii,jj]
                                E[x]=ss[1][lt,ii,jj]
                                E[y]=-ss[1][lt,ii,jj]
                            end
                            Gt1_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *Gt1_* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                            Gt01_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:]*Gt01_
                            G0t1_=G0t1_*model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        end
    
                        if norm(Gt1-Gt1_)+norm(G01-G01_)+norm(Gt01-Gt01_)+norm(G0t1-G0t1_)+
                           norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                            println('\n',norm(Gt1-Gt1_),'\n',norm(G01-G01_),'\n',norm(Gt01-Gt01_),'\n',norm(G0t1-G0t1_))
                            println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
                            error("s1:  $lt  $j:,,,asdasdasd")
                        end
                        #####################################################################
                    end

                    # 更新s2
                    E=[-2*ss[2][lt,i,j] , 2*ss[2][lt,i,j]]
                    Δ2=uv*diagm(exp.(model.α.*E))*uv'-I(2)
                    r2=I(2)+Δ2*(I(2)-Gt2[subidx,subidx])

                    b_A=Gt02[subidx,indexA[:]]*gmInv_A
                    a_A=(2*G01[indexA[:],indexA[:]]-IA)*G0t2[indexA[:],subidx]
                    Tau_A=b_A*a_A

                    b_B=Gt02[subidx,indexB[:]]*gmInv_B
                    a_B=(2*G01[indexB[:],indexB[:]]-IB)*G0t2[indexB[:],subidx]
                    Tau_B=b_B*a_B
                    
                    p=det(r2+Δ2*Tau_A)^λ*det(r2+Δ2*Tau_B)^(1-λ)

                    if p<0
                        println("Negative Sign: $(p)")
                    end

                    if rand(rng)<p
                        rho_A=inv(r2+Δ2*Tau_A)*Δ2
                        gmInv_A-=gmInv_A* a_A * rho_A * b_A
                        detg_A*=det(I(2)+inv(r2)*Δ2*Tau_A) 
    
                        rho_B=inv(r2+Δ2*Tau_B)*Δ2
                        gmInv_B-=gmInv_B* a_B * rho_B * b_B
                        detg_B*=det(I(2)+inv(r2)*Δ2*Tau_B)

                        G02+= (G0t2[:,subidx] /r2*Δ2 *(Gt02[subidx,:]))
                        Gt02+=(Gt2[:,subidx] /r2*Δ2 *(Gt02[subidx,:]))
                        G0t2-=(G0t2[:,subidx] /r2*Δ2  * (II-Gt2)[subidx,:]  )
                        Gt2-= (Gt2[:,subidx] /r2*Δ2 *((II-Gt2)[subidx,:]) )         
                        ss[2][lt,i,j]=-ss[2][lt,i,j]
                        #####################################################################
                        print('*')
                        if lt==div(model.Nt,2)+1
                            Gt2_,G02_,G0t2_,Gt02_=G4(model,ss[2],lt-1,div(model.Nt,2))
                        else
                            Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2))
                        end
                        Gt2_=model.eK*Gt2_*model.eKinv
                        Gt02_=model.eK*Gt02_
                        G0t2_=G0t2_*model.eKinv
                        GM_A_=GroverMatrix(G01[indexA[:],indexA[:]],G02_[indexA[:],indexA[:]])
                        gmInv_A_=inv(GM_A_)
                        GM_B_=GroverMatrix(G01[indexB[:],indexB[:]],G02_[indexB[:],indexB[:]])
                        gmInv_B_=inv(GM_B_)
                        detg_A_=det(GM_A_)
                        detg_B_=det(GM_B_)

                        for jj in size(ss[1])[3]:-1:j
                            E=zeros(model.Ns)
                            for ii in 1:size(ss[1])[2]
                                x,y=model.nnidx[ii,jj]
                                E[x]=ss[2][lt,ii,jj]
                                E[y]=-ss[2][lt,ii,jj]
                            end
                            Gt2_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *Gt2_* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                            Gt02_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:]*Gt02_
                            G0t2_=G0t2_*model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        end
    
                        if norm(Gt2-Gt2_)+norm(G02-G02_)+norm(Gt02-Gt02_)+norm(G0t2-G0t2_)+
                           norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                            println('\n',norm(Gt2-Gt2_),'\n',norm(G02-G02_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t2-G0t2_))
                            println(norm(gmInv_A_-gmInv_A)," ",norm(gmInv_B-gmInv_B_)," ",abs(detg_A-detg_A_)," ",abs(detg_B-detg_B_))
                            error("s2:  $lt  $x:,,,asdasdasd")
                        end
                        #####################################################################
                    end

                end
                # 开始卸载eV eK
                E1=zeros(model.Ns)
                E2=zeros(model.Ns)
                for i in 1:size(ss[1])[2]
                    x,y=model.nnidx[i,j]
                    E1[x]=ss[1][lt,i,j]
                    E1[y]=-ss[1][lt,i,j]
                    E2[x]=ss[2][lt,i,j]
                    E2[y]=-ss[2][lt,i,j]
                end

                Gt1=model.UV[j,:,:]'*diagm(exp.(-model.α*E1))*model.UV[j,:,:] *Gt1* model.UV[j,:,:]'*diagm(exp.(model.α*E1))*model.UV[j,:,:]
                Gt2=model.UV[j,:,:]'*diagm(exp.(-model.α*E2))*model.UV[j,:,:] *Gt2* model.UV[j,:,:]'*diagm(exp.(model.α*E2))*model.UV[j,:,:]
                Gt01=model.UV[j,:,:]'*diagm(exp.(-model.α*E1))*model.UV[j,:,:] *Gt01
                Gt02=model.UV[j,:,:]'*diagm(exp.(-model.α*E2))*model.UV[j,:,:] *Gt02
                G0t1=G0t1* model.UV[j,:,:]'*diagm(exp.(model.α*E1))*model.UV[j,:,:]
                G0t2=G0t2* model.UV[j,:,:]'*diagm(exp.(model.α*E2))*model.UV[j,:,:]
            end
            Gt1=model.eKinv*Gt1*model.eK
            Gt2=model.eKinv*Gt2*model.eK
            Gt01=model.eKinv*Gt01
            Gt02=model.eKinv*Gt02
            G0t1=G0t1*model.eK
            G0t2=G0t2*model.eK

            ##------------------------------------------------------------------------
            tmpO+=(detg_A/detg_B)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
        end

        O[loop+1]=tmpO/counter
        tmpO=counter=0
    end
    return ss
end
