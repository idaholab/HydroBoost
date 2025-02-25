#=
HydroBoost Model

Current version: 1.0
Last update: 09.18.2024

Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov

=#

module HydroBoost_Sim

using JuMP
using MathOptInterface
using LinearAlgebra
const _MOI = MathOptInterface
using CSV
using FileIO

# Solvers
#using CPLEX
using HiGHS

# Create and register module level logger
using Memento
const _LOGGER = Memento.getlogger(@__MODULE__)
__init__() = Memento.register(_LOGGER)

include("core/model_type.jl")
include("core/base_functions.jl")
include("core/model_handler.jl")

include("model_structure/HydroBoost.jl")

end
