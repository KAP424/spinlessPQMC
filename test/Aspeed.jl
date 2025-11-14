using LinearAlgebra
using BenchmarkTools

mutable struct testStruct
    a::Int64
    A::Array{Float64,3}
    uv
end

uv=[-2^0.5/2 -2^0.5/2;-2^0.5/2 2^0.5/2]
As=testStruct(3,rand(1000,1000,3),uv)


RR=rand(1000,1000)


@btimed copyto!(view(As.A,:, :,1),RR)

function sss()
    return 1,2
end


A=rand(10,2)

for i in axes(A,1)
    println(A[i,1], " ", A[i,2])
end