#=
HydroBoost Model

Current version: 1.0
Last update: 09.18.2024

Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov

=#

##### Top Level Abstract Types #####

"Root of the ALEAF model formulation type hierarchy"
abstract type Abstract_ALEAF_Model end

"individual model type"
abstract type Abstract_HydroBoost_Model <: Abstract_ALEAF_Model end

"constructor for ALEAF_Model_Structure_for_HydroBoost"
mutable struct ALEAF_Model_Structure_HydroBoost <: Abstract_ALEAF_Model
    model::Dict{Symbol,<:Any}
    model_type::String

    setting::Dict{String,<:Any}
    solution::Dict{String,<:Any}

    ref::Dict{Symbol,<:Any}
    var::Dict{Symbol,<:Any}
    con::Dict{Symbol,<:Any}

    sol::Dict{Symbol,<:Any}

    cnw::Int
end
