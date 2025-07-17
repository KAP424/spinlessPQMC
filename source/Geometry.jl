

function xy_i(Lattice::String,site::Vector{Int64},x::Int64,y::Int64)::Int64
    if Lattice=="SQUARE"
        if length(site)==2
            if x>site[1] || y>site[2]
                println()
                error("Error :($x,$y) Out of Lattice Range$(site)!")
            end
            return x+(y-1)*site[1]
        else
            error("Error : Dimension wrong")
        end


    elseif  Lattice=="HoneyComb"
        if x>site[1] ||y>site[2]
            error("Error : Out of Lattice Range!")
        end
        return 2*(x+(y-1)*site[1])
    else
        error("Not support $Lattice")
    end
end

function i_xy(Lattice::String,site::Vector{Int64},i::Int64)
    if  Lattice=="SQUARE"
        if i>prod(site)
            error("Out of Lattice Range")
        else
            x=mod1(i,site[1])
            y=cld(i,site[1])
            return x::Int64,y::Int64
        end
    elseif  Lattice=="HoneyComb"
        j=div(i,2)
        return mod1(j,site[1]),div(j,site[1])+1
    else
        error("Not support $Lattice")
    end
end


function K_Matrix(Lattice::String,site::Vector{Int64})
    
    if Lattice=="SQUARE"
        Ns=prod(site)
        K=zeros(Float64,Ns,Ns)
        if length(site)==1
            for i in 1:Ns
                K[i,mod1(i+1,Ns)]+=1
                K[i,mod1(i-1,Ns)]+=1   
            end

        elseif length(site)==2                
            for i in 1:Ns
                x,y=i_xy(Lattice,site,i)
                K[i,xy_i(Lattice,site,mod1(x+1,site[1]),y)]+=1
                K[i,xy_i(Lattice,site,mod1(x-1,site[1]),y)]+=1
                K[i,xy_i(Lattice,site,x,mod1(y+1,site[2]) )]+=1
                K[i,xy_i(Lattice,site,x,mod1(y-1,site[2]))]+=1

            end
        end

        if norm(K-K')>1e-5
            error("K Matrix ERROR, Not symmetric matrix")
        end
        return K

            
    elseif Lattice=="HoneyComb"
        Ns=prod(site)*2
        K=zeros(Float64,Ns,Ns)

        # x,y--period
        for i in 1:site[1]
            for j in 1:site[2]
            # -------------------------------------------------------------

                k=xy_i(Lattice,site,i,j)
                k1=xy_i(Lattice,site,i,mod1(j+1,site[2]))-1
                k2=xy_i(Lattice,site,mod1(i-1,site[1]), j)-1
                K[k,k-1]=K[k,k1]=K[k,k2]=1

                k=k-1
                k1=xy_i(Lattice,site,mod1(i+1,site[1]),j)
                k2=xy_i(Lattice,site,i,mod1(j-1,site[1]))
                K[k,k+1]=K[k,k1]=K[k,k2]=1

            # -------------------------------------------------------------

                # k=xy_i(Lattice,site,i,j)
                # k1=xy_i(Lattice,site,i                 ,mod1(j+1,site[2]))-1
                # k2=xy_i(Lattice,site,mod1(i-1,site[1]), mod1(j+1,site[2]))-1
                # K[k,k-1]=K[k,k1]=K[k,k2]=1

                # k=k-1
                # k1=xy_i(Lattice,site,mod1(i+1,site[1]),mod1(j-1,site[1]))
                # k2=xy_i(Lattice,site,i,                mod1(j-1,site[1]))
                # K[k,k+1]=K[k,k1]=K[k,k2]=1
            # -------------------------------------------------------------

            end
        end

        # x,y--open
        # for i in 1:site[1]
        #     for j in 1:site[2]
        #         k=xy_i(Lattice,site,i,j)
        #         K[k,k-1]=1
        #         if j+1<=site[2]
        #             k1=xy_i(Lattice,site,i                 ,mod1(j+1,site[2]))-1
        #             K[k,k1]=1
        #             if i-1>0
        #                 k2=xy_i(Lattice,site,mod1(i-1,site[1]), mod1(j+1,site[2]))-1
        #                 K[k,k2]=1
        #             end
        #         end

        #         k=k-1
        #         K[k,k+1]=1
        #         if j-1>0
        #             k2=xy_i(Lattice,site,i,                mod1(j-1,site[1]))
        #             K[k,k2]=1
        #             if i+1<=site[1]
        #                 k1=xy_i(Lattice,site,mod1(i+1,site[1]),mod1(j-1,site[1]))
        #                 K[k,k1]=1
        #             end
        #         end
        #     end
        # end

        # x--open ,y--period
        # for i in 1:site[1]
        #     for j in 1:site[2]
        #         k=xy_i(Lattice,site,i,j)
        #         K[k,k-1]=1
        #         k1=xy_i(Lattice,site,i                 ,mod1(j+1,site[2]))-1
        #         K[k,k1]=1
        #         if i-1>0
        #             k2=xy_i(Lattice,site,mod1(i-1,site[1]), mod1(j+1,site[2]))-1
        #             K[k,k2]=1
        #         end

        #         k=k-1
        #         K[k,k+1]=1
        #         k2=xy_i(Lattice,site,i,                mod1(j-1,site[1]))
        #         K[k,k2]=1
        #         if i+1<=site[1]
        #             k1=xy_i(Lattice,site,mod1(i+1,site[1]),mod1(j-1,site[1]))
        #             K[k,k1]=1
        #         end
        #     end
        # end

        if norm(K-K')>1e-5
            println(K)
            error("K Matrix ERROR, Not symmetric matrix")
        end
        return K
    else
        error("Not support $Lattice")
    end

end


function area_index(Lattice::String,site::Vector{Int64},area::Tuple{Vector{Int64}, Vector{Int64}})::Vector{Int64}
    if Lattice=="SQUARE"
        if length(site)==1
            index=[x for x in area[1][1]:area[2][1]]
            return index

        elseif length(site)==2
            counter=1
            index=zeros(Int64,prod(area[2]-area[1]+[1,1]))
            for lx in area[1][1]:area[2][1]
                for ly in area[1][2]:area[2][2]
                    index[counter]=xy_i(Lattice,site,lx,ly)
                    counter+=1
                end
            end
            return index
        end
    elseif Lattice=="HoneyComb"

        L=site[1]
        if area[1][1]==-1
            println("zigzag")
            index=collect(4:2:xy_i(Lattice,site,L-1,1))
            
            for i in 2:div(2*L,3)
                index=vcat(collect(xy_i(Lattice,site,2,i)-1   :1:  xy_i(Lattice,site,L-i,i) ),index)
            end
            return index        
        
        elseif area[1][1]==-2
            # Bearded triangle(A) 
            index=Vector{Int64}()
            
            for i in 2:div(2*L,3)
                index=vcat(xy_i(Lattice,site,2,i)-1,index)
                index=vcat(collect(xy_i(Lattice,site,3,i)-1   :1:  xy_i(Lattice,site,L-i+1,i)-1 ),index)
            end
            index=vcat(xy_i(Lattice,site,2,div(2*L,3)+1)-1,index)
            return index 
        
        else
            counter=1
            index=zeros(Int64,2*prod(area[2]-area[1]+[1,1]))
            for lx in area[1][1]:area[2][1]
                for ly in area[1][2]:area[2][2]
                    index[counter]=xy_i(Lattice,site,lx,ly)-1
                    index[counter+1]=index[counter]+1
                    counter+=2
                end
            end
            return index
        end
    end

end
