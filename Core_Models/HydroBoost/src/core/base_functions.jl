#=
HydroBoost Model

Current version: 1.0
Last update: 09.18.2024

Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov

=#

function check_and_create_path(Report_path::String)
    if !ispath(Report_path)
        mkpath(Report_path)
    end
end

function read_xlsx_return_dict_string_any(file_location::String, tap_name::String; first_row_value=1)

    data = DataFrame(XLSX.readtable(file_location, tap_name, first_row=first_row_value))

    li = Dict{String, Any}()
    for row in eachrow(data)
        di = Dict{String, Any}()
        for name in names(row)
            di[string(name)] = row[name]
        end
        li[string(DataFrames.row(row))] = di
    end

    return li

end

function get_index(am::Abstract_ALEAF_Model, key1::Symbol, nw::Int)
    return keys(am.ref[:nw][nw][key1])
end

function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, key2::String, idx)
    return am.ref[:nw][nw][key1][idx][key2]
end

function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, idx::Int, key2::String)
    return am.ref[:nw][nw][key1][idx][key2]
end

function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, key2::String, idx1, idx2)
    return am.ref[:nw][nw][key1][idx1][key2][idx2]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, idx1, idx2, key2::String)
    return am.ref[:nw][nw][key1][idx1][idx2][key2]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, key2::String, key3::String, d::Int, h::Int, t::Int)
    return am.ref[:nw][nw][key1][d][key2][string(h)][string(t)][key3]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, key2::String, key3::String, key4::String, d::Int, h::Int, t::Int, y::Int)
    return am.ref[:nw][nw][key1][y][key2][string(d)][key3][string(h)][string(t)][key4]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, bus::Int, line::Int)
    return am.ref[:nw][nw][key1][bus][line]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, key2::String)
    return am.ref[:nw][nw][key1][key2]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, idx::Int)
    return am.ref[:nw][nw][key1][idx]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol)
    return am.ref[:nw][nw][key1]
end


function parameter(am::Abstract_ALEAF_Model, key1::String)
    return am.setting[key1]
end


function parameter(am::Abstract_ALEAF_Model, key1::String, key2::String)
    return am.setting[key1][key2]
end


function add_sol_component(aim::Abstract_ALEAF_Model, nw::Int, comp_name::Symbol, field_name::Symbol, comp_ids, variables)
    for i in comp_ids
        @assert !haskey(dict_pointer(aim.sol[:nw][nw], comp_name, i), field_name)
        dict_pointer(aim.sol[:nw][nw], comp_name, i)[field_name] = variables[i]
    end
end


""
function dict_pointer(dict::Dict, args...)
    for arg in args
        if haskey(dict, arg)
            dict = dict[arg]
        else
            dict = dict[arg] = Dict()
        end
    end
    return dict
end
