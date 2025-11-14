using LinearAlgebra
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
