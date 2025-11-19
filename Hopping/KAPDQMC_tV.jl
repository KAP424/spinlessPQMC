# H=∑_{<ij>}c_i c_j^† +∑_{<ij>} n_i n_j
module KAPDQMC_tV
    using Base.Filesystem
    using LinearAlgebra,LinearAlgebra.BLAS,LinearAlgebra.LAPACK
    using DelimitedFiles
    using Random
    using Statistics

    include("Geometry.jl")
    export K_Matrix,xy_i,i_xy,area_index,nn2idx,nnidx_F

    include("Hubbard_model.jl")
    export Hubbard_Para,UpdateBuffer

    include("GreenMatrix.jl")
    export Gτ,G4,Initial_s,G12FF,GroverMatrix

    include("phy_update.jl")
    export phy_update,Poss,phy_measure

    # include("EE_update.jl")
    # export ctrl_EEicr

    include("SCEE.jl")
    export ctrl_SCEEicr

    # include("disorder_operate.jl")
    # export DO_icr
    
end
