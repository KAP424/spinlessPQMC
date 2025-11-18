


function ctrl_EEicr(path::String,model::Hubbard_Para_,index::Vector{Int64},Sweeps::Int64,λ::Float64,Nλ::Int64,ss::Vector{Array{UInt8,3}},record::Bool)::Vector{Array{UInt8,3}}
    uv=[-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]
    if model.Lattice=="SQUARE"
        name="□"
    elseif model.Lattice=="HoneyComb60"
        name="HC60"
    elseif model.Lattice=="HoneyComb120"
        name="HC120"
    end
    file="$(path)$(length(index))_EEicr$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)N$(Nλ)BS$(model.BatchSize).csv"
    

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

    Gt1=zeros(Float64,model.Ns,model.Ns)
    Gt2=zeros(Float64,model.Ns,model.Ns)
    G01=zeros(Float64,model.Ns,model.Ns)
    G02=zeros(Float64,model.Ns,model.Ns)
    Gt01=zeros(Float64,model.Ns,model.Ns)
    Gt02=zeros(Float64,model.Ns,model.Ns)
    G0t1=zeros(Float64,model.Ns,model.Ns)
    G0t2=zeros(Float64,model.Ns,model.Ns)
    gmInv=zeros(Float64,length(index),length(index))
    detg=0

    tmpO=0
    counter=0
    O=zeros(Sweeps+1)
    O[1]=λ

    I1=I(model.Ns)
    I2=I(length(index))

    
    for loop in 1:Sweeps
        for lt in 1:model.Nt
            if lt==div(model.Nt,2)+1
                # 更新Θ+1层时候，得到的是Θ层的四Green函数，再将eK,eV往上安装同时更新
                # 而Θ层的四Green函数是和 1 ~ Θ 层同用一个Wrap规则
                # 所以为了让Θ层的四Green函数与Θ-2Θ同用一个Wrap规则，需要将Θ层的四Green函数的最后两个进行倒换
                Gt1,G01,G0t1,Gt01=G4(model,ss[1],lt-1,div(model.Nt,2))
                Gt2,G02,G0t2,Gt02=G4(model,ss[2],lt-1,div(model.Nt,2))
                Gt1=model.eK*Gt1*model.eKinv
                Gt2=model.eK*Gt2*model.eKinv
                Gt01=model.eK*Gt01
                Gt02=model.eK*Gt02
                G0t1=G0t1*model.eKinv
                G0t2=G0t2*model.eKinv

                GM=GroverMatrix(G01[index[:],index[:]],G02[index[:],index[:]])
                gmInv=inv(GM)
                detg=det(GM)
            elseif  mod1(lt,model.WrapTime)==1 
                Gt1,G01,Gt01,G0t1=G4(model,ss[1],lt-1,div(model.Nt,2))
                Gt2,G02,Gt02,G0t2=G4(model,ss[2],lt-1,div(model.Nt,2))
                Gt1=model.eK*Gt1*model.eKinv
                Gt2=model.eK*Gt2*model.eKinv
                Gt01=model.eK*Gt01
                Gt02=model.eK*Gt02
                G0t1=G0t1*model.eKinv
                G0t2=G0t2*model.eKinv

                GM=GroverMatrix(G01[index[:],index[:]],G02[index[:],index[:]])
                gmInv=inv(GM)
                detg=det(GM)
            else
                #####################################################################
                # Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2))
                # Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2))
                    
                # if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
                #     println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
                #     error("$lt : WrapTime")
                # end
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
                    Δ=uv*diagm(exp.(model.α.*E))*uv'-I(2)
                    r=I(2)+Δ*(I(2)-Gt1[subidx,subidx])

                    b=(Gt01[subidx,index[:]]) *(2*G02[index[:],index[:]]-I2)*gmInv
                    a=G0t1[index[:],subidx]
                    Tau=b*a

                    p=det(r)^(1-λ)*det(r+Δ*Tau)^λ

                    if p<0
                        println("Negative Sign: $(p)")
                    end

                    if rand(rng)<p
                        rho=inv(r+Δ*Tau)*Δ
                        gmInv-=gmInv* a * rho * b
                        detg*=det(I(2)+inv(r)*Δ*Tau)
    
                        G01+= (G0t1[:,subidx] /r*Δ *(Gt01[subidx,:]))
                        Gt01+=(Gt1[:,subidx] /r*Δ *(Gt01[subidx,:]))
                        G0t1-=(G0t1[:,subidx] /r*Δ  * (I1-Gt1)[subidx,:]  )
                        Gt1-= (Gt1[:,subidx] /r*Δ *((I1-Gt1)[subidx,:]) )         
                        ss[1][lt,i,j]=-ss[1][lt,i,j]
                        #####################################################################
                        # print('-')
                        # if lt==div(model.Nt,2)+1
                        #     Gt1_,G01_,G0t1_,Gt01_=G4(model,ss[1],lt-1,div(model.Nt,2))
                        # else
                        #     Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2))
                        # end
                        # Gt1_=model.eK*Gt1_*model.eKinv
                        # Gt01_=model.eK*Gt01_
                        # G0t1_=G0t1_*model.eKinv
                        # GM_=GroverMatrix(G01_[index[:],index[:]],G02[index[:],index[:]])
                        # gmInv_=inv(GM_)
                        # detg_=det(GM_)

                        # for jj in size(ss[1])[3]:-1:j
                        #     E=zeros(model.Ns)
                        #     for ii in 1:size(ss[1])[2]
                        #         x,y=model.nnidx[ii,jj]
                        #         E[x]=ss[1][lt,ii,jj]
                        #         E[y]=-ss[1][lt,ii,jj]
                        #     end
                        #     Gt1_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *Gt1_* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        #     Gt01_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:]*Gt01_
                        #     G0t1_=G0t1_*model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        # end
    
                        # if norm(Gt1-Gt1_)+norm(G01-G01_)+norm(Gt01-Gt01_)+norm(G0t1-G0t1_)+
                        #    norm(gmInv-gmInv_)+abs(detg-detg_)>1e-3
                        #     println('\n',norm(Gt1-Gt1_),'\n',norm(G01-G01_),'\n',norm(Gt01-Gt01_),'\n',norm(G0t1-G0t1_))
                        #     println(norm(gmInv-gmInv_)," ",abs(detg-detg_))
                        #     error("s1:  $lt  $j:,,,asdasdasd")
                        # end
                        ####################################################################
                    end

                    # 更新s2
                    E=[-2*ss[2][lt,i,j] , 2*ss[2][lt,i,j]]
                    Δ=uv*diagm(exp.(model.α.*E))*uv'-I(2)
                    r=I(2)+Δ*(I(2)-Gt2[subidx,subidx])

                    b=Gt02[subidx,index[:]]*gmInv
                    a=(2*G01[index[:],index[:]]-I2)*G0t2[index[:],subidx]
                    Tau=b*a

                    p=det(r)^(1-λ)*det(r+Δ*Tau)^λ

                    if p<0
                        println("Negative Sign: $(p)")
                    end

                    if rand(rng)<p
                        rho=inv(r+Δ*Tau)*Δ
                        gmInv-=gmInv* a * rho * b
                        detg*=det(I(2)+inv(r)*Δ*Tau) 
    
                        G02+= (G0t2[:,subidx] /r*Δ *(Gt02[subidx,:]))
                        Gt02+=(Gt2[:,subidx] /r*Δ *(Gt02[subidx,:]))
                        G0t2-=(G0t2[:,subidx] /r*Δ  * (I1-Gt2)[subidx,:]  )
                        Gt2-= (Gt2[:,subidx] /r*Δ *((I1-Gt2)[subidx,:]) )         
                        ss[2][lt,i,j]=-ss[2][lt,i,j]
                        #####################################################################
                        # print('*')
                        # if lt==div(model.Nt,2)+1
                        #     Gt2_,G02_,G0t2_,Gt02_=G4(model,ss[2],lt-1,div(model.Nt,2))
                        # else
                        #     Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2))
                        # end
                        # Gt2_=model.eK*Gt2_*model.eKinv
                        # Gt02_=model.eK*Gt02_
                        # G0t2_=G0t2_*model.eKinv
                        # GM_=GroverMatrix(G01[index[:],index[:]],G02_[index[:],index[:]])
                        # gmInv_=inv(GM_)
                        # detg_=det(GM_)

                        # for jj in size(ss[1])[3]:-1:j
                        #     E=zeros(model.Ns)
                        #     for ii in 1:size(ss[1])[2]
                        #         x,y=model.nnidx[ii,jj]
                        #         E[x]=ss[2][lt,ii,jj]
                        #         E[y]=-ss[2][lt,ii,jj]
                        #     end
                        #     Gt2_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *Gt2_* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        #     Gt02_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:]*Gt02_
                        #     G0t2_=G0t2_*model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        # end
    
                        # if norm(Gt2-Gt2_)+norm(G02-G02_)+norm(Gt02-Gt02_)+norm(G0t2-G0t2_)+
                        #    norm(gmInv-gmInv_)+abs(detg-detg_)>1e-3
                        #     println('\n',norm(Gt2-Gt2_),'\n',norm(G02-G02_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t2-G0t2_))
                        #     println(norm(gmInv-gmInv_)," ",abs(detg-detg_))
                        #     error("s2:  $lt  $x:,,,asdasdasd")
                        # end
                        #####################################################################
                    end
                end
            end
            ##------------------------------------------------------------------------

            tmpO+=(detg)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
        end

        # println("\n #####################################################################")
        
        for lt in model.Nt-1:-1:1
            if mod1(model.Nt-lt,model.WrapTime)==1 || lt==div(model.Nt,2)
                Gt1,G01,Gt01,G0t1=G4(model,ss[1],lt,div(model.Nt,2))
                Gt2,G02,Gt02,G0t2=G4(model,ss[2],lt,div(model.Nt,2))

                GM=GroverMatrix(G01[index[:],index[:]],G02[index[:],index[:]])
                gmInv=inv(GM)
                detg=det(GM)
            else
                #####################################################################
                # Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt,div(model.Nt,2))
                # Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt,div(model.Nt,2))
                    
                # if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
                #     println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
                #     error("$lt : WrapTime")
                # end
                #####################################################################
            end

            for j in 1:size(ss[1])[3]
                for i in 1:size(ss[1])[2]
                    # update
                    x,y=model.nnidx[i,j]
                    subidx=[x,y]

                    # 更新s1
                    E=[-2*ss[1][lt,i,j] , 2*ss[1][lt,i,j]]
                    Δ=uv*diagm(exp.(model.α.*E))*uv'-I(2)
                    r=I(2)+Δ*(I(2)-Gt1[subidx,subidx])

                    b=(Gt01[subidx,index[:]]) *(2*G02[index[:],index[:]]-I2)*gmInv
                    a=G0t1[index[:],subidx]
                    Tau=b*a
                    
                    p=det(r)^(1-λ)*det(r+Δ*Tau)^λ

                    if p<0
                        println("Negative Sign: $(p)")
                    end

                    if rand(rng)<p
                        rho=inv(r+Δ*Tau)*Δ
                        gmInv-=gmInv* a * rho * b
                        detg*=det(I(2)+inv(r)*Δ*Tau)
    
                        G01+= (G0t1[:,subidx] /r*Δ *(Gt01[subidx,:]))
                        Gt01+=(Gt1[:,subidx] /r*Δ *(Gt01[subidx,:]))
                        G0t1-=(G0t1[:,subidx] /r*Δ  * (I1-Gt1)[subidx,:]  )
                        Gt1-= (Gt1[:,subidx] /r*Δ *((I1-Gt1)[subidx,:]) )         
                        ss[1][lt,i,j]=-ss[1][lt,i,j]
                        #####################################################################
                        # print('-')
                        # if lt==div(model.Nt,2)+1
                        #     Gt1_,G01_,G0t1_,Gt01_=G4(model,ss[1],lt-1,div(model.Nt,2))
                        # else
                        #     Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt-1,div(model.Nt,2))
                        # end
                        # Gt1_=model.eK*Gt1_*model.eKinv
                        # Gt01_=model.eK*Gt01_
                        # G0t1_=G0t1_*model.eKinv
                        # GM_=GroverMatrix(G01_[index[:],index[:]],G02[index[:],index[:]])
                        # gmInv_=inv(GM_)
                        # detg_=det(GM_)

                        # for jj in size(ss[1])[3]:-1:j
                        #     E=zeros(model.Ns)
                        #     for ii in 1:size(ss[1])[2]
                        #         x,y=model.nnidx[ii,jj]
                        #         E[x]=ss[1][lt,ii,jj]
                        #         E[y]=-ss[1][lt,ii,jj]
                        #     end
                        #     Gt1_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *Gt1_* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        #     Gt01_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:]*Gt01_
                        #     G0t1_=G0t1_*model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        # end
    
                        # if norm(Gt1-Gt1_)+norm(G01-G01_)+norm(Gt01-Gt01_)+norm(G0t1-G0t1_)+
                        #    norm(gmInv-gmInv_)+abs(detg-detg_)>1e-3
                        #     println('\n',norm(Gt1-Gt1_),'\n',norm(G01-G01_),'\n',norm(Gt01-Gt01_),'\n',norm(G0t1-G0t1_))
                        #     println(norm(gmInv_-gmInv)," ",abs(detg-detg_))
                        #     error("s1:  $lt  $j:,,,asdasdasd")
                        # end
                        #####################################################################
                    end

                    # 更新s2
                    E=[-2*ss[2][lt,i,j] , 2*ss[2][lt,i,j]]
                    Δ=uv*diagm(exp.(model.α.*E))*uv'-I(2)
                    r=I(2)+Δ*(I(2)-Gt2[subidx,subidx])

                    b=Gt02[subidx,index[:]]*gmInv
                    a=(2*G01[index[:],index[:]]-I2)*G0t2[index[:],subidx]
                    Tau=b*a

                    p=det(r)^(1-λ)*det(r+Δ*Tau)^λ

                    if p<0
                        println("Negative Sign: $(p)")
                    end

                    if rand(rng)<p
                        rho=inv(r+Δ*Tau)*Δ
                        gmInv-=gmInv* a * rho * b
                        detg*=det(I(2)+inv(r)*Δ*Tau) 
    
                        G02+= (G0t2[:,subidx] /r*Δ *(Gt02[subidx,:]))
                        Gt02+=(Gt2[:,subidx] /r*Δ *(Gt02[subidx,:]))
                        G0t2-=(G0t2[:,subidx] /r*Δ  * (I1-Gt2)[subidx,:]  )
                        Gt2-= (Gt2[:,subidx] /r*Δ *((I1-Gt2)[subidx,:]) )         
                        ss[2][lt,i,j]=-ss[2][lt,i,j]
                        #####################################################################
                        # print('*')
                        # if lt==div(model.Nt,2)+1
                        #     Gt2_,G02_,G0t2_,Gt02_=G4(model,ss[2],lt-1,div(model.Nt,2))
                        # else
                        #     Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt-1,div(model.Nt,2))
                        # end
                        # Gt2_=model.eK*Gt2_*model.eKinv
                        # Gt02_=model.eK*Gt02_
                        # G0t2_=G0t2_*model.eKinv
                        # GM_=GroverMatrix(G01[index[:],index[:]],G02_[index[:],index[:]])
                        # gmInv_=inv(GM_)
                        # detg_=det(GM_)

                        # for jj in size(ss[1])[3]:-1:j
                        #     E=zeros(model.Ns)
                        #     for ii in 1:size(ss[1])[2]
                        #         x,y=model.nnidx[ii,jj]
                        #         E[x]=ss[2][lt,ii,jj]
                        #         E[y]=-ss[2][lt,ii,jj]
                        #     end
                        #     Gt2_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:] *Gt2_* model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        #     Gt02_=model.UV[jj,:,:]'*diagm(exp.(model.α*E))*model.UV[jj,:,:]*Gt02_
                        #     G0t2_=G0t2_*model.UV[jj,:,:]'*diagm(exp.(-model.α*E))*model.UV[jj,:,:]
                        # end
    
                        # if norm(Gt2-Gt2_)+norm(G02-G02_)+norm(Gt02-Gt02_)+norm(G0t2-G0t2_)+
                        #    norm(gmInv-gmInv_)+abs(detg-detg_)>1e-3
                        #     println('\n',norm(Gt2-Gt2_),'\n',norm(G02-G02_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t2-G0t2_))
                        #     println(norm(gmInv_-gmInv)," ",abs(detg-detg_))
                        #     error("s2:  $lt  $x:,,,asdasdasd")
                        # end
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
            tmpO+=(detg)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
        end

        O[loop+1]=tmpO/counter
        tmpO=counter=0
    end
    return ss
end
