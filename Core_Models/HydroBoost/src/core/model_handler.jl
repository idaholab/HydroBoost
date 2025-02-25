#=
HydroBoost Model

Current version: 1.0
Last update: 09.18.2024

Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov

=#


function collect_result_distributed(JuMP_model::JuMP.AbstractModel, solution_list::Dict{Symbol,<:Any}, solve_time)

    result_count = _MOI.get(JuMP_model, _MOI.ResultCount())
    solution = Dict{String,Any}()
    if result_count > 0
        solution = collect_solution_distributed(solution_list)

        gap = 0.0
        try
            gap = JuMP.relative_gap(JuMP_model)
        catch
            gap = 0.0
        end

        result = Dict{String,Any}(
            "result_count" => result_count,
            "optimizer" => JuMP.solver_name(JuMP_model),
            "termination_status" => JuMP.termination_status(JuMP_model),
            "primal_status" => JuMP.primal_status(JuMP_model),
            "dual_status" => JuMP.dual_status(JuMP_model),
            "objective" => JuMP.objective_value(JuMP_model),
            "objective_lb" => JuMP.objective_bound(JuMP_model),
            "solve_time" => solve_time,
            "solution" => solution,
            "relative_gap" => gap
            )        
    else
        Memento.warn(_LOGGER, "model has no results, solution dict will be empty")

        result = Dict{String,Any}(
            "result_count" => result_count,
            "optimizer" => JuMP.solver_name(JuMP_model),
            "termination_status" => JuMP.termination_status(JuMP_model),
            "primal_status" => JuMP.primal_status(JuMP_model),
            "dual_status" => JuMP.dual_status(JuMP_model),
            "objective" => 0,
            "objective_lb" => 0,
            "solve_time" => solve_time,
            "solution" => 0,
            "relative_gap" => 0.0
            )     
    end

    return result
end


"this will collect solution for multiple networks"
function collect_solution_distributed(solution_list::Dict{Symbol,<:Any})

    sol = Dict{String, Any}()
    for solution_category in keys(solution_list)
        sol[string(solution_category)] = Dict{String, Any}()
        for var_idx in keys(solution_list[solution_category])
            sol[string(solution_category)][string(var_idx)] = Dict{String, Any}()
        end
    end
    
    for solution_category in keys(solution_list)
        # for var_idx in keys(solution_list[solution_category])
        Threads.@threads for var_idx in collect(keys(solution_list[solution_category]))
            for var in keys(solution_list[solution_category][var_idx])
                sol[string(solution_category)][string(var_idx)][string(var)] = JuMP.value(solution_list[solution_category][var_idx][var])
            end
        end
    end
    
    # for (k,v) in sol
    #     sol[k] = v
    # end

    return sol
end



function initialize_model_instance_HydroBoost(model_type::Type, data::Dict{String,<:Any}, day_id)
    @assert model_type <: Abstract_ALEAF_Model

    setting = Dict{String,Any}()
    # jump_model::JuMP.AbstractModel=JuMP.Model()

    ref = network_ref_initialize_HydroBoost(model_type, data, day_id) # reference data

    var = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    con = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    sol = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    jump_model = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    for (nw_id, nw) in ref[:nw]
        var[:nw][nw_id] = Dict{Symbol,Any}()
        con[:nw][nw_id] = Dict{Symbol,Any}()
        sol[:nw][nw_id] = Dict{Symbol,Any}()
        jump_model[:nw][nw_id] = JuMP.Model()
    end

    solution = Dict{String,Any}()
    cnw = 1
    
    imo = ALEAF_Model_Structure_HydroBoost(
        jump_model,
        string(model_type),
        setting,
        solution, 
        ref,
        var,
        con,
        sol,
        cnw
    )

    return imo
end


function network_ref_initialize_HydroBoost(model_type::Type, data::Dict{String,<:Any}, day_id)

    refs = Dict{Symbol,Any}()

    nws_data = Dict("0" => data)
    
    nws = refs[:nw] = Dict{Int,Any}()

    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        ref = nws[nw_id] = Dict{Symbol,Any}()

        for (key, item) in nw_data
            if key == "repdays"
                ref[Symbol(key)] = item[day_id]
            else
                if isa(item, Dict{String,Any})
                    item_lookup = Dict{Int,Any}([(parse(Int, k), v) for (k,v) in item])
                    ref[Symbol(key)] = item_lookup
                elseif isa(item, Dict{Int64,Any})
                    item_lookup = Dict{Int,Any}([(k, v) for (k,v) in item])
                    ref[Symbol(key)] = item_lookup
                else
                    ref[Symbol(key)] = item
                end
            end
        end
    end

    return refs
end

