# using density channel ±1 HS transformation

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
    Pt::Array{ComplexF64,2}
    HalfeK::Array{ComplexF64,2}
    eK::Array{ComplexF64,2}
    HalfeKinv::Array{ComplexF64,2}
    eKinv::Array{ComplexF64,2}
    nnidx::Matrix{Tuple{Int64, Int64}}
    UV::Array{ComplexF64, 3}
end



function Hubbard_Para(t,U,Lattice::String,site,Δt,Θ,BatchSize,Initial::String)
    Nt::Int64=2*cld(Θ,Δt)
    WrapTime::Int64=div(BatchSize,2)
    
    α::Float64=acosh(exp(Δt*U/2)) 
    
    K=K_Matrix(Lattice,site)
    H0=t*1im*UpperTriangular(K)
    H0=(H0+H0')/4
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
    elseif  Lattice=="HoneyComb"
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
    # 交错化学势，打开gap，去兼并
    # μ=0.0
    # if Lattice=="HoneyComb"
    #     K+=μ*diagm(repeat([-1, 1], div(Ns, 2)))
    # elseif Lattice=="SQUARE"
    #     for i in 1:Ns
    #         x,y=i_xy(Lattice,site,i)
    #         K[i,i]+=μ*(-1)^(x+y)
    #     end
    # end
    # K[K .!= 0] .+=( rand(size(K)...) * 0.1)[K.!= 0]
    # K=(K+K')./2

    E,V=eigen(H0)
    HalfeK=V*diagm(exp.(-Δt.*E./2))*V'
    eK=V*diagm(exp.(-Δt.*E))*V'
    HalfeKinv=V*diagm(exp.(Δt.*E./2))*V'
    eKinv=V*diagm(exp.(Δt.*E))*V'

    if Initial=="H0"
        Pt=V[:,1:div(Ns,2)]
    elseif Initial=="V" 
        Pt=zeros(Float64,Ns,Int(Ns/2))
        for i in 1:Int(Ns/2)
            Pt[i*2,i]=1
        end
    end

    a,b=size(nnidx)
    s=ones(Int8,Nt,a,b)
    UV=zeros(ComplexF64,3,Ns,Ns)
    
    lt=1
    for j in 1:size(s)[3]
        V=zeros(ComplexF64,Ns,Ns)
        for i in 1:size(s)[2]
            x,y=nnidx[i,j]
            V[x,y]=s[lt,i,j]*1im/4
            V[y,x]=-s[lt,i,j]*1im/4
        end
        _,UV[j,:,:]=eigen(V)
        UV[j,:,:]=UV[j,:,:]'
    end

    return _Hubbard_Para(Lattice,t,U,site,Θ,Ns,Nt,K,BatchSize,WrapTime,Δt,α,Pt,HalfeK,eK,HalfeKinv,eKinv,nnidx,UV)

end


# function setμ(model::_Hubbard_Para,μ)
#     # fig1:1d
#     # km=abs(acos(μ/2))
#     # N_particle=Int(round(km/π*model.Ns))
    
#     # fig2:2d-circle Fermi surface
#     N_particle=Int(round( μ^2/4/π *model.Ns ))
    
#     E,V=eigen(model.K)
#     model.Pt=V[:,1:N_particle]
# end

