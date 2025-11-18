using LinearAlgebra,LinearAlgebra.LAPACK
using BenchmarkTools
using Random

rng=MersenneTwister(1234)

elements = (1, 2, 3, 4)
samplers_dict = Dict{UInt8, Random.Sampler}()
for excluded in elements
    allowed = [i for i in elements if i != excluded]
    samplers_dict[excluded] = Random.Sampler(rng, allowed)
end

samplers_dict


rng2=MersenneTwister(5678)

sx = rand(rng2,  samplers_dict[2])


A=rand(5,3)

tau = Vector{Float64}(undef, 3)

@btime LAPACK.geqrf!($A,$tau)

@btime LAPACK.orgqr!($A,$tau)

A=Diagonal([1.0,2.0,3.0,4.0])

A[2,2]

A=rand(4,4)

 reverse(axes(A, 2))

for i in axes(A,1)
    A[i,i] += 1e-3 * rand()
end

