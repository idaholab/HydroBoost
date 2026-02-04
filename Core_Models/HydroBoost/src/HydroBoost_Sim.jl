#=

HydroBoost Model

Main Authors:
Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov
Carlos Josue Lopez; Argonne National Laboratory; clopezsalgado@anl.gov
Alberto Grimaldi; Argonne National Laboratory; agrimaldi@anl.gov

Current version: 2.0
Last update: 07.31.2025

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

# Select the model structure by uncommenting the corresponding include line below

# include("model_structure/HydroBoost_Model_A.jl")
include("model_structure/HydroBoost_Model_B.jl")
# include("model_structure/HydroBoost_Model_C.jl")

end

# Use the following instructions on the terminal to run the selected model:

### Model A ###

# import Pkg; Pkg.activate("."); Pkg.instantiate()
# include("Core_Models/HydroBoost/src/HydroBoost_Sim.jl")
# using .HydroBoost_Sim
# result = HydroBoost_Sim.execute_HydroBoost_model("Model_A")

### Model B ###

# import Pkg; Pkg.activate("."); Pkg.instantiate()
# include("Core_Models/HydroBoost/src/HydroBoost_Sim.jl")
# using .HydroBoost_Sim
# result = HydroBoost_Sim.execute_HydroBoost_model("Model_B")

### Model C ###

# import Pkg; Pkg.activate("."); Pkg.instantiate()
# include("Core_Models/HydroBoost/src/HydroBoost_Sim.jl")
# using .HydroBoost_Sim
# result = HydroBoost_Sim.execute_HydroBoost_model("Model_C")