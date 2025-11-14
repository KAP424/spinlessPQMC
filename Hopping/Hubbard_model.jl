# using density channel ±1,±2 HS transformation

# mutable
struct _Hubbard_Para
    Lattice::String
    t::Float64
    U::Float64
    site::Vector{Int64}
    Θ::Float64
    Ns::Int64
    Nt::Int64
    K::Array{Float64,2}
    BatchSize::Int64
    WrapTime::Int64
    Δt::Float64
    α::Float64
    γ::Vector{Float64}
    η::Vector{Float64}
    Pt::Array{Float64,2}
    HalfeK::Array{Float64,2}
    eK::Array{Float64,2}
    HalfeKinv::Array{Float64,2}
    eKinv::Array{Float64,2}
    nnidx::Matrix{Tuple{Int64, Int64}}
    UV::Array{Float64, 3}
    nodes::Vector{Int64}
end


function Hubbard_Para(t,U,Lattice::String,site,Δt,Θ,BatchSize,Initial::String)
    Nt::Int64=2*cld(Θ,Δt)
    WrapTime::Int64=div(BatchSize,2)
    
    α = sqrt(Δt * U / 2)
    γ = [1 + sqrt(6) / 3, 1 + sqrt(6) / 3, 1 - sqrt(6) / 3, 1 - sqrt(6) / 3]
    η = [sqrt(2 * (3 - sqrt(6))), -sqrt(2 * (3 - sqrt(6))), sqrt(2 * (3 + sqrt(6))), -sqrt(2 * (3 + sqrt(6)))]
    
    K=K_Matrix(Lattice,site)
    Ns=size(K)[1]
    if Lattice=="SQUARE"
        Ns=prod(site)
        if length(site)==1
            nnidx=fill((0, 0), div(Ns,2), 2)
            count=1
            for i in 1:2:Ns
                nn=nn2idx(Lattice,site,i)
                for j in eachindex(nn)
                    nnidx[count,j]=(i,nn[j])
                end
                count+=1
            end
        elseif length(site)==2
            nnidx=fill((0, 0), div(Ns,2), 4)
            count=1
            for x in 1:site[1]
                for y in 1:site[2]
                    if (x+y)%2==1
                        i=x+(y-1)*site[1]
                        nn=nn2idx(Lattice,site,i)
                        for j in eachindex(nn)
                            nnidx[count,j]=(i,nn[j])
                        end
                        count+=1
                    end
                end
            end
        end
    elseif  occursin("HoneyComb", Lattice)
        Ns=prod(site)*2
        nnidx=fill((0, 0), div(Ns,2), 3)
        count=1
        for i in 1:2:Ns
            nn=nn2idx(Lattice,site,i)
            for j in eachindex(nn)
                nnidx[count,j]=(i,nn[j])
            end
            count+=1
        end
    end

    E,V=LAPACK.syevd!('V', 'L',-t.*K[:,:])
    HalfeK=V*Diagonal(exp.(-Δt.*E./2))*V'
    eK=V*Diagonal(exp.(-Δt.*E))*V'
    HalfeKinv=V*Diagonal(exp.(Δt.*E./2))*V'
    eKinv=V*Diagonal(exp.(Δt.*E))*V'


    Pt=zeros(Float64,Ns,Int(Ns/2))
    if Initial=="H0"
        KK=K[:,:]
        # 交错化学势，打开gap，去兼并
        μ=1e-3
        if occursin("HoneyComb", Lattice)
            KK+=μ*Diagonal(repeat([-1, 1], div(Ns, 2)))
        elseif Lattice=="SQUARE"
            for i in 1:Ns
                x,y=i_xy(Lattice,site,i)
                KK[i,i]+=μ*(-1)^(x+y)
            end
        end

        # hopping 扰动，避免能级简并
        # KK[KK .!= 0] .+=( rand(size(KK)...) * 1e-3)[KK.!= 0]
        # KK=(KK+KK')./2
        
        E,V=LAPACK.syevd!('V', 'L',KK[:,:])
        # Pt.=V[:,1:div(Ns,2)]
        Pt.=V[:,div(Ns,2)+1:end]
    elseif Initial=="V" 
        if occursin("HoneyComb", Lattice)
            for i in 1:div(Ns,2)
                Pt[i*2-1,i]=1
            end
        else
            count=1
            for i in 1:Ns
                x,y=i_xy(Lattice,site,i)
                if (x+y)%2==1
                    Pt[i,count]=1
                    count+=1
                end
            end
        end
    end

    # Pt=HalfeKinv*Pt

    a,b=size(nnidx)
    s=ones(Int8,Nt,a,b)
    UV=zeros(Float64,Ns,Ns,b)
    
    for j in 1:b
        for i in 1:size(s)[2]
            x,y=nnidx[i,j]
            UV[x,x,j]=UV[x,y,j]=UV[y,x,j]=-2^0.5/2
            UV[y,y,j]=2^0.5/2
        end
    end

    if div(Nt, 2) % BatchSize == 0
        nodes = collect(0:BatchSize:Nt)
    else
        nodes = vcat(0, reverse(collect(div(Nt, 2) - BatchSize:-BatchSize:1)), collect(div(Nt, 2):BatchSize:Nt), Nt)
    end

    return _Hubbard_Para(Lattice,t,U,site,Θ,Ns,Nt,K,BatchSize,WrapTime,Δt,α,γ,η,Pt,HalfeK,eK,HalfeKinv,eKinv,nnidx,UV,nodes)

end



if abspath(PROGRAM_FILE) == @__FILE__
    push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/spinlessPQMC/Hopping/")
    using KAPDQMC_tV
    using LinearAlgebra,Random

    model=KAPDQMC_tV.Hubbard_Para(1.0,4.0,"HoneyComb120",[3,3],0.1,2.0,10,"H0")
    println("Hubbard model initialized.")
    s=Initial_s(model,MersenneTwister(1234))
    
    # TEST For diag transformation of nn interaction UV*Diagonal(s)*UV' = V
    lt=1
    tmpVV=zeros(Float64,3,model.Ns,model.Ns)
    for j in 3:-1:1
        tmpN=zeros(Float64,model.Ns)
        tmpV=zeros(Float64,model.Ns,model.Ns)
        println("---------------")
        for i in 1:div(model.Ns,2)
            # println(i," ",j,": ",model.nnidx[i,j])
            x,y=model.nnidx[i,j]
            tmpN[x]=model.η[s[lt,i,j]]
            tmpN[y]=-model.η[s[lt,i,j]]
            
            tmpV[x,y]=model.η[s[lt,i,j]]
            tmpV[y,x]=model.η[s[lt,i,j]]
        end
        tmpVV[j,:,:]=tmpV[:,:]
        @assert norm(model.UV[:,:,j]*Diagonal(tmpN)*model.UV[:,:,j]'-tmpV)<1e-5
        @assert norm(model.UV[:,:,j]'*model.UV[:,:,j]-I(model.Ns))<1e-5
    end

    # TEST for not comute for V1,V2,V3
    @assert norm(tmpVV[1,:,:]*tmpVV[2,:,:]-tmpVV[2,:,:]*tmpVV[1,:,:])>1e-5
    @assert norm(tmpVV[1,:,:]*tmpVV[3,:,:]-tmpVV[3,:,:]*tmpVV[1,:,:])>1e-5
    @assert norm(tmpVV[2,:,:]*tmpVV[3,:,:]-tmpVV[3,:,:]*tmpVV[2,:,:])>1e-5

    # TEST for comute for V_i,V_i^′
    for j in 3:-1:1
        tmp=tmpVV[j,:,:]
        for i in 1:div(model.Ns,2)
            x,y=model.nnidx[i,j]
            tmp[x,y]=-tmp[x,y]
            tmp[y,x]=-tmp[y,x]
            @assert norm(tmp * tmpVV[j,:,:] - tmpVV[j,:,:]*tmp) <1e-5
        end
    end

    # TEST for diag transformation of update UV*Diagonal(s)*UV' = V
    uv=[-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]
    for sx_f in 1:4
        for sx_i in 1:4
            tmp22= (model.η[sx_f]-model.η[sx_i])*[0 1 ; 1.0 0]
            tmp2=(model.η[sx_f]-model.η[sx_i])*[1,-1]
            println(norm(uv*Diagonal(tmp2)*uv' .- tmp22))
        end
    end
end



