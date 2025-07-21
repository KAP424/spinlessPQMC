# H=∑_{<ij>}c_i c_j^† +∑_{<ij>} n_i n_j
module KAPDQMC_spinless_CDW
    using Base.Filesystem
    using LinearAlgebra
    using DelimitedFiles
    using Random
    using Statistics

    include("../source/Geometry.jl")
    export K_Matrix,xy_i,i_xy,area_index,nn2idx

    include("../source/Hubbard_model.jl")
    export Hubbard_Para,_Hubbard_Para,setμ

    include("GreenMatrix.jl")
    export Gτ,G4,Initial_s,G12FF,GroverMatrix

    # include("phy_measure.jl")
    # export EK,NN,Magnetism,CzzofSpin

    include("phy_update.jl")
    export phy_update,Poss

    # include("EE_update.jl")
    # export ctrl_EEicr,EE_dir,EEICR

    # include("SCEE.jl")
    # export ctrl_SCEEicr

    # include("disorder_operate.jl")
    # export DO_icr
    
end
