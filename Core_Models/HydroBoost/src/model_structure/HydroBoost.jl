#=
HydroBoost Model

Current version: 1.0
Last update: 09.18.2024

Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov

=#

using XLSX
using DataFrames
using Infiltrator
using Dates
using JSON
using Statistics
using DelimitedFiles


""
function execute_HydroBoost_model(project_id::String)

    # Read setting
    ALEAF_setting = read_ALEAF_HydroBoost_setting(project_id)
    ALEAF_setting["ALEAF_model_type"] = Abstract_HydroBoost_Model

    # Check cases to run
    Memento.info(_LOGGER, "-- Current case: $project_id")

    # Generate Network Data
    start_time = time()
    network_data = generate_networkdata_HydroBoost(ALEAF_setting, project_id)
    Memento.info(_LOGGER, "[HydroBoost Model]:\tGenerate network data. Time(sec): $(round(time() - start_time, digits=2))")

    # Build and Run HydroBoost model instances 
    solutions = build_and_run_daily_HydroBoost(ALEAF_setting, network_data)

    # Report Solutions
    export_HydroBoost_results(ALEAF_setting, network_data, solutions, project_id)

end


function build_and_run_daily_HydroBoost(ALEAF_setting, network_data)

    daily_solutions = Dict{String, Any}()

    # run daily HydroBoost sequantially
    # run day_id 1 first (without prior day solution)
    day_idx = 1
    daily_solutions[string(day_idx)] = build_and_run_HydroBoost_for_each_day(day_idx, ALEAF_setting, network_data)        
    
    # run remaining days using solution from day-1 as initial status
    # for day_id in [2, 3, 4]
    for day_id in eachindex([i for i in 1:365])
        if day_id != 1  # ignore day_id 1
            daily_solutions[string(day_id)] = build_and_run_HydroBoost_for_each_day(day_id, ALEAF_setting, network_data; prior_day_solution=daily_solutions[string(day_id-1)]["solution"])        
        end
    end

    return daily_solutions

end


function build_and_run_HydroBoost_for_each_day(day_id, ALEAF_setting, network_data; prior_day_solution::Dict{String,<:Any} = Dict{String,Any}())

    start_time = time()

    # build common ALEAF model instance structure
    ALEAF_model_instance = initialize_model_instance_HydroBoost(ALEAF_setting["ALEAF_model_type"], network_data, day_id)

    # Add reference data 
    add_ref_HydroBoost_model!(ALEAF_model_instance, ALEAF_setting)

    # build optimization model instance 
    build_HydroBoost_optimization_model_instance!(ALEAF_model_instance, day_id, prior_day_solution)

    # Define JuMP_model, solution_list, and solver setting
    JuMP_model = ALEAF_model_instance.model[:nw][0]
    solution_list = ALEAF_model_instance.sol[:nw][0]
    solver_setting = ALEAF_model_instance.setting["Solver Setting"]
    
    # Export ALEAF_model_instance model
    if ALEAF_setting["Simulation Setting"]["export_model_lp_flag"] == true
        output_path = parameter(ALEAF_model_instance, 0, :output_path)
        file_name = ALEAF_setting["Simulation Setting"]["model_lp_file_name_value"]
        file_name = string(file_name, "_day_", day_id, ".lp")
        JuMP.write_to_file(JuMP_model, string(output_path, file_name))
        Memento.info(_LOGGER, "[HydroBoost Model]:\tExport HydroBoost_model_instance model (lp format)")
    end
    preparation_time = round(time() - start_time, digits=2)

    start_time = time()
    result_HydroBoost = solve_model_HydroBoost!(JuMP_model, solution_list, solver_setting)
    solution_time = round(time() - start_time, digits=2)

    total_time = round(preparation_time + solution_time, digits=2)

    Memento.info(_LOGGER, "[HydroBoost Model]:\tSolved Day ID: $day_id, Prep (sec): $preparation_time, Solution (sec): $solution_time, Total (sec): $total_time")

    ALEAF_model_instance = 0.0 
    
    return result_HydroBoost
    
end



####################################################
#------ DEFINIE HydroBoost Optimization Model
####################################################


function build_HydroBoost_optimization_model_instance!(am::Abstract_ALEAF_Model, day_id, prior_day_solution)

    # Reset a JuMP model
    JuMP_model = JuMP.Model()

    const_name_flag = am.setting["Simulation Setting"]["const_name_flag"]
    if const_name_flag == false
        JuMP.set_string_names_on_creation(JuMP_model, false)
    end

    ###################################
    #------ Pre-processing
    ###################################
    ids_i = [(i) for (i) in get_index(am, :gen_index, 0)]
    ids_j_hydro = [(j) for (j) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO"]
    ids_i_sto = [(i) for (i) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"]

    ids_l = [(l) for (l) in 1:am.setting["Simulation Setting"]["num_hydropower_performance_segment_value"]]
    ids_h = [(h) for (h) in am.setting["run_H"]]
    ids_t = [(t) for (t) in am.setting["run_T"]]
            
    ###################################
    #------ DEFINE DECISION VARIABLES
    ###################################
    
    # Storage Variables
    HydroBoost_variable_iht_binary(JuMP_model, am, :storage, "u_B_iht", ids_i_sto, ids_h, ids_t)                                             # Binary variable driven to 1 when BESS is set to charging mode, and 0 otherwise      
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "p_B_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)          # Energy discharged from BESS and accounted for at point of delivery in period t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "p_B_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)          # Energy contributing to charge BESS in period t before accounting for losses [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "e_B_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)            # Storage device state of charge in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)         # Reserve for regulation up in discharging mode [MW]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)         # Reserve for regulation down in discharging mode [MW]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)         # Reserve for regulation up in charging mode [MW]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)         # Reserve for regulation down in charging mode [MW]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)           # BESS reserve sold to market for regulation up [MW]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)           # BESS reserve sold to market for regulation down [MW]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)         # Reserve for spinning reserve in discharging mode [MW]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)         # Reserve for spinning reserve in charging mode [MW]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0)           # BESS reserve sold to market for spinning reserve [MW]
    
    HydroBoost_variable_ht_real(JuMP_model, am, :storage, "p_B_DT_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0)                      # Total power from hydroelectric plant in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :storage, "p_B_CT_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0)                      # Total charging from hydroelectric plant in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :storage, "p_GB_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0)                        # Energy from grid allocated to charge BESS (before considering any round-trip losses) [MWh]
    
    # Hydropower Dispatch Variables    
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "u_H_jht", ids_j_hydro, ids_h, ids_t)                                           # Binary variable which is equal to 1 if hydro generator j is on-line in hour t, and 0 otherwise
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "a_H_jht", ids_j_hydro, ids_h, ids_t)                                           # Binary variable which is equal to 1 if hydro generator j is started at beginning of hour t, and 0 otherwise
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "z_H_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0,  bounded_upper=true, upper_bound=1.0)                                           
                                                                                                                                             # Binary variable which is equal to 1 if hydro generator j is shut down at beginning of hour t, and 0 otherwise
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "p_H_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0)          # Power output of hydroelectric generator j in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_HT_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0)                        # Total power from hydroelectric plant in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_HB_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0)                        # Hydroelectric power allocated to charge the BESS in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_HG_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0)                        # Hydroelectric power delivered to the grid in hour t [MWh]

    # Hydropower Water Uses Variables
    HydroBoost_variable_ilht_real(JuMP_model, am, :water, "u_l_jht", ids_j_hydro, ids_l, ids_h, ids_t; bounded_lower=true, lower_bound=0)  # Water discharge of block ℓ of hydro generator j in hour t 
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "u_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0)            # Water discharge of hydro generator j in hour t
    HydroBoost_variable_ilht_binary(JuMP_model, am, :water, "w_l_jht", ids_j_hydro, ids_l, ids_h, ids_t)                                   # Binary variable which is equal to 1 if water discharged by hydro generator j has exceeded block ℓ in hour t
    HydroBoost_variable_ht_real(JuMP_model, am, :water, "s_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0)                           # Spillage of reservoir in hour t
    HydroBoost_variable_ht_real(JuMP_model, am, :water, "e_H_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0)                         # Water volume of reservoir in period t

    # Hydropower Rough Zone Variables 
    if am.setting["Simulation Setting"]["hydropower_rough_zone_flag"] == true
        HydroBoost_variable_iht_real(JuMP_model, am, :rough_zone, "y_l_plus_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0)    # Slack variables which take non-zero value if operation outside rough zone ℓ of unit j cannot be honored in hour t
        HydroBoost_variable_iht_real(JuMP_model, am, :rough_zone, "y_l_minus_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0)   # Slack variables which take non-zero value if operation outside rough zone ℓ of unit j cannot be honored in hour t
        HydroBoost_variable_iht_binary(JuMP_model, am, :rough_zone, "phi_l_jht", ids_j_hydro, ids_h, ids_t)                                       # Auxiliary binary variable required for representation of rough zone constraint
    end

    ###################################
    #------ DEFINE Constraints
    ###################################

    # BESS constraints 
    
    for h in ids_h 
        for t in ids_t

            # total discharge and charge 
            HydroBoost_constraint_total_ES_power_discharge_ht(JuMP_model, am, "HydroBoost_constraint_total_ES_power_discharge_ht", ids_i_sto, h, t; const_name_flag)
            HydroBoost_constraint_total_ES_power_charge_ht(JuMP_model, am, "HydroBoost_constraint_total_ES_power_charge_ht", ids_i_sto, h, t; const_name_flag)

            for i in ids_i_sto
            
                # SOC balance 
                HydroBoost_constraint_ES_SOC_Balance_Inter_Hour_iht(JuMP_model, am, "HydroBoost_constraint_ES_SOC_Balance_Inter_Hour_iht", i, h, t, day_id, prior_day_solution; const_name_flag)

                # SOC bounds 
                HydroBoost_constraint_ES_SOC_Bounds_DN_iht(JuMP_model, am, "HydroBoost_constraint_ES_SOC_Bounds_DN_iht", i, h, t; const_name_flag)
                HydroBoost_constraint_ES_SOC_Bounds_UP_iht(JuMP_model, am, "HydroBoost_constraint_ES_SOC_Bounds_UP_iht", i, h, t; const_name_flag)

                # Operational Limits in discharging mode
                HydroBoost_constraint_ES_power_Bounds_UP_iht(JuMP_model, am, "HydroBoost_constraint_ES_power_Bounds_UP_iht", i, h, t; const_name_flag)
                HydroBoost_constraint_ES_power_Bounds_DN_iht(JuMP_model, am, "HydroBoost_constraint_ES_power_Bounds_DN_iht", i, h, t; const_name_flag)
                HydroBoost_constraint_ES_RU_D_cap_iht(JuMP_model, am, "HydroBoost_constraint_ES_RU_D_cap_iht", i, h, t; const_name_flag)
                HydroBoost_constraint_ES_RD_D_cap_iht(JuMP_model, am, "HydroBoost_constraint_ES_RD_D_cap_iht", i, h, t; const_name_flag)
                HydroBoost_constraint_ES_SR_D_cap_iht(JuMP_model, am, "HydroBoost_constraint_ES_SR_D_cap_iht", i, h, t; const_name_flag)

                # Operational Limits in charging mode
                HydroBoost_constraint_ES_charge_Bounds_UP_iht(JuMP_model, am, "HydroBoost_constraint_ES_charge_Bounds_UP_iht", i, h, t; const_name_flag)
                HydroBoost_constraint_ES_charge_Bounds_DN_iht(JuMP_model, am, "HydroBoost_constraint_ES_charge_Bounds_DN_iht", i, h, t; const_name_flag)
                HydroBoost_constraint_ES_RD_C_cap_iht(JuMP_model, am, "HydroBoost_constraint_ES_RD_C_cap_iht", i, h, t; const_name_flag)
                HydroBoost_constraint_ES_RU_C_cap_iht(JuMP_model, am, "HydroBoost_constraint_ES_RU_C_cap_iht", i, h, t; const_name_flag)
                HydroBoost_constraint_ES_SR_C_cap_iht(JuMP_model, am, "HydroBoost_constraint_ES_SR_C_cap_iht", i, h, t; const_name_flag)

                # Ancillary services offered to market 
                HydroBoost_constraint_ES_reg_up_sale_iht(JuMP_model, am, "HydroBoost_constraint_ES_reg_up_sale_iht", i, h, t; const_name_flag)
                HydroBoost_constraint_ES_reg_dn_sale_iht(JuMP_model, am, "HydroBoost_constraint_ES_reg_dn_sale_iht", i, h, t; const_name_flag)
                HydroBoost_constraint_ES_spin_sale_iht(JuMP_model, am, "HydroBoost_constraint_ES_spin_sale_iht", i, h, t; const_name_flag)

            end
        end
    end

    # Hydropower system constraints 
    for h in ids_h
        for t in ids_t

            # Generation of hydropower generators constraints
            HydroBoost_constraint_total_hydro_power_generation_jht(JuMP_model, am, "HydroBoost_constraint_total_hydro_power_generation_jht", ids_j_hydro, h, t; const_name_flag)

            for j in ids_j_hydro
                # Commitment of hydropower generators constraints
                HydroBoost_constraint_hydro_commitment_status_jht(JuMP_model, am, "HydroBoost_constraint_hydro_commitment_status_jht", j, h, t, day_id, prior_day_solution; const_name_flag)
                HydroBoost_constraint_hydro_start_up_shut_down_bound_jht(JuMP_model, am, "HydroBoost_constraint_hydro_start_up_shut_down_bound_jht", j, h, t; const_name_flag)

                # Generation of hydropower generators constraints
                HydroBoost_constraint_hydro_power_generation_jht(JuMP_model, am, "HydroBoost_constraint_hydro_power_generation_jht", ids_l, j, h, t; const_name_flag)
            end

            # Water discharge constraints 
            for j in ids_j_hydro

                HydroBoost_constraint_hydro_total_water_use_jht(JuMP_model, am, "HydroBoost_constraint_hydro_total_water_use_jht", ids_l, j, h, t; const_name_flag)

                for l in ids_l
                    HydroBoost_constraint_hydro_power_water_discharge_bounds_up_ljht(JuMP_model, am, "HydroBoost_constraint_hydro_power_water_discharge_bounds_up_ljht", j, l, h, t; const_name_flag)
                    HydroBoost_constraint_hydro_power_water_discharge_bounds_down_ljht(JuMP_model, am, "HydroBoost_constraint_hydro_power_water_discharge_bounds_down_ljht", j, l, h, t; const_name_flag)
                end
            end

            # Ramping constraints
            for j in ids_j_hydro
                HydroBoost_constraint_hydro_power_ramping_bound_up_jht(JuMP_model, am, "HydroBoost_constraint_hydro_power_ramping_bound_up_jht", j, h, t, day_id, prior_day_solution; const_name_flag)
                HydroBoost_constraint_hydro_power_ramping_bound_dn_jht(JuMP_model, am, "HydroBoost_constraint_hydro_power_ramping_bound_dn_jht", j, h, t, day_id, prior_day_solution; const_name_flag)
            end

        end
    end


    # Water Discharge Constraints
    for h in ids_h
        for t in ids_t

            # Spillage of reservoir constraints
            HydroBoost_constraint_hydro_power_water_spillage_ht(JuMP_model, am, "HydroBoost_constraint_hydro_power_water_spillage_ht", ids_j_hydro, h, t; const_name_flag)

            # Water release lower bounds constraints
            HydroBoost_constraint_hydro_power_water_release_lower_bound_ht(JuMP_model, am, "HydroBoost_constraint_hydro_power_water_release_lower_bound_ht", ids_j_hydro, h, t; const_name_flag)

            # Water balance constraints
            HydroBoost_constraint_hydro_water_balance_ht(JuMP_model, am, "HydroBoost_constraint_hydro_water_balance_ht", ids_j_hydro, h, t, day_id, prior_day_solution; const_name_flag)

            # Resorvior storage constraints
            HydroBoost_constraint_hydro_resorvoir_volume_lower_bound_ht(JuMP_model, am, "HydroBoost_constraint_hydro_resorvoir_volume_lower_bound_ht", h, t; const_name_flag)
            HydroBoost_constraint_hydro_resorvoir_volume_upper_bound_ht(JuMP_model, am, "HydroBoost_constraint_hydro_resorvoir_volume_upper_bound_ht", h, t; const_name_flag)

        end
    end

    # Initial and last-period resorvior volumne constraints
    HydroBoost_constraint_hydro_resorvoir_uses_limit(JuMP_model, am, "HydroBoost_constraint_hydro_resorvoir_uses_limit", ids_h, ids_t; const_name_flag)

    # Coupling constraints
    for h in ids_h
        for t in ids_t
            HydroBoost_constraint_hydropower_power_generation_coupling_ht(JuMP_model, am, "HydroBoost_constraint_hydropower_power_generation_coupling_ht", h, t; const_name_flag)
            HydroBoost_constraint_ES_power_generation_coupling_ht(JuMP_model, am, "HydroBoost_constraint_ES_power_generation_coupling_ht", h, t; const_name_flag)
        end
    end

    # Interconnection constraints
    for h in ids_h
        for t in ids_t
            HydroBoost_constraint_outflow_interconnection_limit_ht(JuMP_model, am, "HydroBoost_constraint_outflow_interconnection_limit_ht", h, t; const_name_flag)
            HydroBoost_constraint_inflow_interconnection_limit_ht(JuMP_model, am, "HydroBoost_constraint_inflow_interconnection_limit_ht", h, t; const_name_flag)
        end
    end

    # Rough zone constraints
    if am.setting["Simulation Setting"]["hydropower_rough_zone_flag"] == true

        rough_zone_segment_number = am.setting["Simulation Setting"]["hydropower_rough_zone_segment_number_value"]

        for h in ids_h
            for t in ids_t
                for j in ids_j_hydro
                    HydroBoost_constraint_hydro_rough_zone_y_l_minus_bound_jht(JuMP_model, am, "HydroBoost_constraint_hydro_rough_zone_y_l_minus_bound_jht", rough_zone_segment_number, j, h, t; const_name_flag)
                    HydroBoost_constraint_hydro_rough_zone_y_l_plus_bound_jht(JuMP_model, am, "HydroBoost_constraint_hydro_rough_zone_y_l_plus_bound_jht", rough_zone_segment_number, j, h, t; const_name_flag)
                    HydroBoost_constraint_hydro_rough_zone_y_l_minus_big_M_jht(JuMP_model, am, "HydroBoost_constraint_hydro_rough_zone_y_l_minus_big_M_jht", rough_zone_segment_number, j, h, t; const_name_flag)
                    HydroBoost_constraint_hydro_rough_zone_y_l_plus_big_M_jht(JuMP_model, am, "HydroBoost_constraint_hydro_rough_zone_y_l_plus_big_M_jht", rough_zone_segment_number, j, h, t; const_name_flag)
                    HydroBoost_constraint_hydro_rough_zone_y_l_plus_upper_bound_jht(JuMP_model, am, "HydroBoost_constraint_hydro_rough_zone_y_l_plus_upper_bound_jht", rough_zone_segment_number, j, h, t; const_name_flag)
                    HydroBoost_constraint_hydro_rough_zone_y_l_minus_upper_bound_jht(JuMP_model, am, "HydroBoost_constraint_hydro_rough_zone_y_l_minus_upper_bound_jht", rough_zone_segment_number, j, h, t; const_name_flag)
                end
            end
        end
    end
    
    ###################################
    #------ DEFINIE Objective
    ###################################
    HydroBoost_objective_function(JuMP_model, am, ids_j_hydro, ids_i_sto, ids_h, ids_t, ids_l)

    # export_JuMP_model_lp_file(JuMP_model) 

    am.model[:nw][0] = JuMP_model

end


###################################
#------ DEFINIE Objective
###################################


function HydroBoost_objective_function(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, ids_j_hydro, ids_i_sto, ids_h, ids_t, ids_l)

    objective = JuMP.AffExpr(0.0)    

    # Energy revenue 
    for h in ids_h
        for t in ids_t
            JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :p_HG_ht, (h,t)))
            JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :p_B_DT_ht, (h,t)))
            JuMP.add_to_expression!(objective, -1.00001 * am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :p_GB_ht, (h,t)))
        end
    end

    # AS Revenue
    for h in ids_h
        for t in ids_t
            for i in ids_i_sto
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Regulation_up"], get_variable(am, :r_RU_iht, (i,h,t)))
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Regulation_down"], get_variable(am, :r_RD_iht, (i,h,t)))
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Spin"], get_variable(am, :r_SR_iht, (i,h,t)))
            end
        end
    end

    # Adjusted revenue due to regulation deployments
    for h in ids_h
        for t in ids_t
            for i in ids_i_sto
                reg_up_signal = 0.0 # TODO
                reg_down_signal = 0.0   # TODO

                JuMP.add_to_expression!(objective, reg_up_signal * am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :r_RU_iht, (i,h,t)))
                JuMP.add_to_expression!(objective, - reg_down_signal * am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :r_RD_iht, (i,h,t)))                
            end
        end
    end

    # Operational costs
    for j in ids_j_hydro
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
        START_UP_COST = parameter(am, bus_idx, :gen_bus, tech_idx, "START_UP_COST")
        SHUT_DN_COST = parameter(am, bus_idx, :gen_bus, tech_idx, "SHUT_DN_COST")

        for h in ids_h
            for t in ids_t
                JuMP.add_to_expression!(objective, - START_UP_COST, get_variable(am, :a_H_jht, (j,h,t)))
                JuMP.add_to_expression!(objective, - SHUT_DN_COST, get_variable(am, :z_H_jht, (j,h,t)))
            end
        end
    end

    # Rough Zone Operation Penalty
    if am.setting["Simulation Setting"]["hydropower_rough_zone_flag"] == true
        for j in ids_j_hydro
            for h in ids_h
                for t in ids_t
                    JuMP.add_to_expression!(objective, - am.setting["Simulation Setting"]["hyropower_rough_zone_operation_penalty_value"], get_variable(am, :y_l_plus_jht, (j,h,t)))
                    JuMP.add_to_expression!(objective, - am.setting["Simulation Setting"]["hyropower_rough_zone_operation_penalty_value"], get_variable(am, :y_l_minus_jht, (j,h,t)))
                end
            end
        end
    end
    
    #DEFINE OBJECTIVE FUNCTION
    return JuMP.@objective(JuMP_model, Max, objective     
    )

end


##################################################
#------ Define Constraints
##################################################


function HydroBoost_constraint_hydro_rough_zone_y_l_plus_big_M_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    big_M = 100000  # Todo: adjust

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    tag_id = string("Water Flow_", rough_zone_segment_number)
    U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
    
    # variable
    phi_l_jht = get_variable(am, :phi_l_jht, (j,h,t))
    y_l_plus_jht = get_variable(am, :y_l_plus_jht, (j,h,t))
    u_l_jht = get_variable(am, :u_l_jht, (j,rough_zone_segment_number,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        (U_bar_l_j - u_l_jht - y_l_plus_jht) - big_M * (1-phi_l_jht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end


function HydroBoost_constraint_hydro_rough_zone_y_l_minus_big_M_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    big_M = 100000  # Todo: adjust
    
    # variable
    phi_l_jht = get_variable(am, :phi_l_jht, (j,h,t))
    y_l_minus_jht = get_variable(am, :y_l_minus_jht, (j,h,t))
    u_l_jht = get_variable(am, :u_l_jht, (j,rough_zone_segment_number,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        (u_l_jht - y_l_minus_jht) - big_M * phi_l_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end


function HydroBoost_constraint_hydro_rough_zone_y_l_plus_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    tag_id = string("Water Flow_", rough_zone_segment_number)
    U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
    
    # variable
    y_l_plus_jht = get_variable(am, :y_l_plus_jht, (j,h,t))
    u_l_jht = get_variable(am, :u_l_jht, (j,rough_zone_segment_number,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        U_bar_l_j - u_l_jht - y_l_plus_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end


function HydroBoost_constraint_hydro_rough_zone_y_l_minus_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    y_l_minus_jht = get_variable(am, :y_l_minus_jht, (j,h,t))
    u_l_jht = get_variable(am, :u_l_jht, (j,rough_zone_segment_number,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        u_l_jht - y_l_minus_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end


function HydroBoost_constraint_hydro_rough_zone_y_l_plus_upper_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    tag_id = string("Water Flow_", rough_zone_segment_number)
    U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
    
    # variable
    y_l_plus_jht = get_variable(am, :y_l_plus_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        y_l_plus_jht - U_bar_l_j
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end


function HydroBoost_constraint_hydro_rough_zone_y_l_minus_upper_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    tag_id = string("Water Flow_", rough_zone_segment_number)
    U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
        
    # variable
    y_l_minus_jht = get_variable(am, :y_l_minus_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        y_l_minus_jht - U_bar_l_j
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

function HydroBoost_constraint_inflow_interconnection_limit_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    inflow_limits = am.setting["Simulation Setting"]["Interconnection Limits Inflow"]

    # variable
    p_GB_ht = get_variable(am, :p_GB_ht, (h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_GB_ht - inflow_limits
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end


function HydroBoost_constraint_outflow_interconnection_limit_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    outflow_limits = am.setting["Simulation Setting"]["Interconnection Limits Outflow"]

    # variable
    p_B_DT_ht = get_variable(am, :p_B_DT_ht, (h,t))
    p_HG_ht = get_variable(am, :p_HG_ht, (h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        (p_B_DT_ht + p_HG_ht) - outflow_limits
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end


function HydroBoost_constraint_ES_power_generation_coupling_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    p_B_CT_ht = get_variable(am, :p_B_CT_ht, (h,t))
    p_HB_ht = get_variable(am, :p_HB_ht, (h,t))
    p_GB_ht = get_variable(am, :p_GB_ht, (h,t))
            
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_B_CT_ht - (p_HB_ht + p_GB_ht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end


function HydroBoost_constraint_hydropower_power_generation_coupling_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    p_HT_ht = get_variable(am, :p_HT_ht, (h,t))
    p_HG_ht = get_variable(am, :p_HG_ht, (h,t))
    p_HB_ht = get_variable(am, :p_HB_ht, (h,t))
        
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_HT_ht - (p_HG_ht + p_HB_ht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end


function HydroBoost_constraint_hydro_power_ramping_bound_dn_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int, day_id, prior_day_solution; const_name_flag::Bool=false)
    
    if day_id == 1
        if h >= 2
            
            # set
            
            # parameter  
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
            MAX_RD = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RD")
            MAX_RD_SHUT_DN = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RD_SHUT_DN")
            
            # variable
            p_H_jht = get_variable(am, :p_H_jht, (j,h,t))
            prior_p_H_jht = get_variable(am, :p_H_jht, (j,h-1,t))
            u_H_jht = get_variable(am, :u_H_jht, (j,h,t))
            z_H_jht = get_variable(am, :z_H_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                - (p_H_jht- prior_p_H_jht) - MAX_RD * u_H_jht - MAX_RD_SHUT_DN * z_H_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

        end
    else
        if h == 1

            # set
            
            # parameter  
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
            MAX_RD = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RD")
            MAX_RD_SHUT_DN = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RD_SHUT_DN")

            idx_string = string("(", j, ", ", 24,  ", ", t,")")
            
            # variable
            p_H_jht = get_variable(am, :p_H_jht, (j,h,t))
            prior_p_H_jht = prior_day_solution["hydro"][idx_string]["p_H_jht"]   # get prior day solution for hour 1
            u_H_jht = get_variable(am, :u_H_jht, (j,h,t))
            z_H_jht = get_variable(am, :z_H_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                - (p_H_jht- prior_p_H_jht) - MAX_RD * u_H_jht - MAX_RD_SHUT_DN * z_H_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

            
        else #if h >= 2
            
            # set
            
            # parameter  
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
            MAX_RD = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RD")
            MAX_RD_SHUT_DN = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RD_SHUT_DN")
            
            # variable
            p_H_jht = get_variable(am, :p_H_jht, (j,h,t))
            prior_p_H_jht = get_variable(am, :p_H_jht, (j,h-1,t))
            u_H_jht = get_variable(am, :u_H_jht, (j,h,t))
            z_H_jht = get_variable(am, :z_H_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                - (p_H_jht- prior_p_H_jht) - MAX_RD * u_H_jht - MAX_RD_SHUT_DN * z_H_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    
            
    
        end
    end


end


function HydroBoost_constraint_hydro_power_ramping_bound_up_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int, day_id, prior_day_solution; const_name_flag::Bool=false)
    
    if day_id == 1
        if h >= 2
            
            # set
            
            # parameter  
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
            MAX_RU = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RU")
            MAX_RU_START_UP = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RU_START_UP")
            
            # variable
            p_H_jht = get_variable(am, :p_H_jht, (j,h,t))
            prior_p_H_jht = get_variable(am, :p_H_jht, (j,h-1,t))
            prior_u_H_jht = get_variable(am, :u_H_jht, (j,h-1,t))
            a_H_jht = get_variable(am, :a_H_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                (p_H_jht- prior_p_H_jht) - MAX_RU * prior_u_H_jht - MAX_RU_START_UP * a_H_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

        end
    else
        if h == 1
            # set
            
            # parameter  
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
            MAX_RU = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RU")
            MAX_RU_START_UP = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RU_START_UP")

            idx_string = string("(", j, ", ", 24,  ", ", t,")")
            
            # variable
            p_H_jht = get_variable(am, :p_H_jht, (j,h,t))
            prior_p_H_jht = prior_day_solution["hydro"][idx_string]["p_H_jht"]   # get prior day solution for hour 1
            prior_u_H_jht = prior_day_solution["hydro"][idx_string]["u_H_jht"]   # get prior day solution for hour 1
            a_H_jht = get_variable(am, :a_H_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                (p_H_jht- prior_p_H_jht) - MAX_RU * prior_u_H_jht - MAX_RU_START_UP * a_H_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end 



        else # if h >= 2
            # set
            
            # parameter  
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
            MAX_RU = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RU")
            MAX_RU_START_UP = parameter(am, bus_idx, :gen_bus, tech_idx, "MAX_RU_START_UP")
            
            # variable
            p_H_jht = get_variable(am, :p_H_jht, (j,h,t))
            prior_p_H_jht = get_variable(am, :p_H_jht, (j,h-1,t))
            prior_u_H_jht = get_variable(am, :u_H_jht, (j,h-1,t))
            a_H_jht = get_variable(am, :a_H_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                (p_H_jht- prior_p_H_jht) - MAX_RU * prior_u_H_jht - MAX_RU_START_UP * a_H_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    
        end
    end

end


function HydroBoost_constraint_hydro_resorvoir_uses_limit(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_h, ids_t; const_name_flag::Bool=false)
    
    # Set/Index
    last_hour = last(ids_h)
    last_min = last(ids_t)
    
    # parameter
    VMAX = am.ref[:nw][0][:hydro_resorvior][1]["VMAX"]
    MAX_Use_Per_Simulation_Percent = am.ref[:nw][0][:hydro_resorvior][1]["MAX_Use_Per_Simulation_Percent"]

    # variable
    e_H_ht = get_variable(am, :e_H_ht, (last_hour,last_min))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_H_ht - VMAX * (1-MAX_Use_Per_Simulation_Percent)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($last_hour,$last_min)")) end    

end


function HydroBoost_constraint_hydro_resorvoir_volume_upper_bound_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter
    VMAX = am.ref[:nw][0][:hydro_resorvior][1]["VMAX"]

    # variable
    e_H_ht = get_variable(am, :e_H_ht, (h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_H_ht - VMAX
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end


function HydroBoost_constraint_hydro_resorvoir_volume_lower_bound_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter
    VMIN = am.ref[:nw][0][:hydro_resorvior][1]["VMIN"]

    # variable
    e_H_ht = get_variable(am, :e_H_ht, (h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_H_ht - VMIN
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end


function HydroBoost_constraint_hydro_water_balance_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_j_hydro, h::Int, t::Int, day_id, prior_day_solution; const_name_flag::Bool=false)
    
    # parameter
    Water_Inflow = am.ref[:nw][0][:repdays]["data"][h][t]["Inflow"]
    CF = 3600 * 0.0000229569 
    
    # variable
    s_ht = get_variable(am, :s_ht, (h,t))
    e_H_ht = get_variable(am, :e_H_ht, (h,t))
    VMAX = am.ref[:nw][0][:hydro_resorvior][1]["VMAX"]
    
    global prior_e_H_ht = am.ref[:nw][0][:hydro_resorvior][1]["VINI_Percent"] * VMAX    # this value will be used for day id 1 at hour 1
    if h == 1
        if day_id > 1
            idx_string = string("(", 24,  ", ", t,")")
            global prior_e_H_ht = prior_day_solution["water"][idx_string]["e_H_ht"]   # get prior day solution for hour 1
        end
    else    # if h >= 2
        global prior_e_H_ht = get_variable(am, :e_H_ht, (h-1,t))

    end

    sum_water_use = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :u_jht, (j,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        (e_H_ht - prior_e_H_ht) - CF * (Water_Inflow - sum_water_use - s_ht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end


function HydroBoost_constraint_hydro_power_water_release_lower_bound_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_j_hydro, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    Water_Release_Requirement = am.ref[:nw][0][:repdays]["data"][h][t]["Water_Release_Requirement"]
    
    # variable
    s_ht = get_variable(am, :s_ht, (h,t))

    sum_water_use = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :u_jht, (j,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        sum_water_use + s_ht - Water_Release_Requirement
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end


function HydroBoost_constraint_hydro_power_water_spillage_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_j_hydro, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    s_ht = get_variable(am, :s_ht, (h,t))

    sum_water_use = JuMP.AffExpr(0.0)
    for j in ids_j_hydro

        bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
        W_j = parameter(am, bus_idx, :gen_bus, tech_idx, "START_UP_Water_Use")

        JuMP.add_to_expression!(sum_water_use, W_j, get_variable(am, :a_H_jht, (j,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        s_ht - sum_water_use
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end


function HydroBoost_constraint_hydro_power_water_discharge_bounds_down_ljht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, l::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
     # parameter  
     bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
     tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
     tag_id = string("Water Flow_", l)
     U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
     
     # variable
     u_l_jht = get_variable(am, :u_l_jht, (j,l,h,t))
     w_l_jht = get_variable(am, :w_l_jht, (j,l,h,t))
     
     # constraint
     expr = JuMP.@expression(JuMP_model,  
         u_l_jht - U_bar_l_j * w_l_jht
         )   
     constraint = JuMP.@constraint(JuMP_model, 
         expr >= 0
         )    
     if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$l,$h,$t)")) end    

end





function HydroBoost_constraint_hydro_power_water_discharge_bounds_up_ljht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, l::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    if l == 1   # first water block

        # parameter  
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
        tag_id = string("Water Flow_", l)
        U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
        
        # variable
        u_l_jht = get_variable(am, :u_l_jht, (j,l,h,t))
        u_H_jht = get_variable(am, :u_H_jht, (j,h,t))
        
        # constraint
        expr = JuMP.@expression(JuMP_model,  
            u_l_jht - U_bar_l_j * u_H_jht
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$l,$h,$t)")) end    


    else

        # parameter  
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
        tag_id = string("Water Flow_", l)
        U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
        
        # variable
        u_l_jht = get_variable(am, :u_l_jht, (j,l,h,t))
        prior_w_l_jht = get_variable(am, :w_l_jht, (j,l-1,h,t))
        
        # constraint
        expr = JuMP.@expression(JuMP_model,  
            u_l_jht - U_bar_l_j * prior_w_l_jht
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$l,$h,$t)")) end    

    end
end


function HydroBoost_constraint_total_ES_power_charge_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_i_sto, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    p_B_CT_ht = get_variable(am, :p_B_CT_ht, (h,t))
        
    sum_power_output = JuMP.AffExpr(0.0)
    for i in ids_i_sto
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :p_B_C_iht, (i,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_B_CT_ht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

function HydroBoost_constraint_total_ES_power_discharge_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_i_sto, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    p_B_DT_ht = get_variable(am, :p_B_DT_ht, (h,t))
        
    sum_power_output = JuMP.AffExpr(0.0)
    for i in ids_i_sto
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :p_B_D_iht, (i,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_B_DT_ht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

function HydroBoost_constraint_total_hydro_power_generation_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_j_hydro, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    p_HT_ht = get_variable(am, :p_HT_ht, (h,t))
        
    sum_power_output = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :p_H_jht, (j,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_HT_ht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end



function HydroBoost_constraint_hydro_total_water_use_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_l, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    Min_Water_Release = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Water_Release")
    
    # variable
    u_jht = get_variable(am, :u_jht, (j,h,t))
    u_H_jht = get_variable(am, :u_H_jht, (j,h,t))
    
    sum_water_use = JuMP.AffExpr(0.0)
    for l in ids_l
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :u_l_jht, (j,l,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        u_jht - sum_water_use - Min_Water_Release * u_H_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end


function HydroBoost_constraint_hydro_power_generation_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_l, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    P0 = parameter(am, bus_idx, :gen_bus, tech_idx, "PMIN (MW)")
    
    # variable
    p_H_jht = get_variable(am, :p_H_jht, (j,h,t))
    u_H_jht = get_variable(am, :u_H_jht, (j,h,t))
    
    sum_power_output = JuMP.AffExpr(0.0)
    for l in ids_l
        tag_id = string("Water_Power_Conversion_", l)
        rho_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)

        JuMP.add_to_expression!(sum_power_output, rho_l_j, get_variable(am, :u_l_jht, (j,l,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_H_jht - P0 * u_H_jht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end


function HydroBoost_constraint_hydro_start_up_shut_down_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    
    # variable
    a_H_jht = get_variable(am, :a_H_jht, (j,h,t))
    z_H_jht = get_variable(am, :z_H_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        a_H_jht + z_H_jht - 1
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end


function HydroBoost_constraint_hydro_commitment_status_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int, day_id, prior_day_solution; const_name_flag::Bool=false)

    
    if day_id == 1
        if h >= 2   # u_H_jht of h=1 for day_id 1 will be optimized
        
            # parameter  
            
            # variable
            u_H_jht = get_variable(am, :u_H_jht, (j,h,t))
            previous_u_H_jht = get_variable(am, :u_H_jht, (j,h-1,t))
            a_H_jht = get_variable(am, :a_H_jht, (j,h,t))
            z_H_jht = get_variable(am, :z_H_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                u_H_jht - previous_u_H_jht - a_H_jht + z_H_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

        end
    else
        if h == 1

            # parameter  
            idx_string = string("(", j, ", ", 24,  ", ", t,")")

            # variable
            u_H_jht = get_variable(am, :u_H_jht, (j,h,t))
            previous_u_H_jht = prior_day_solution["hydro"][idx_string]["u_H_jht"]   # get prior day solution for hour 1
            a_H_jht = get_variable(am, :a_H_jht, (j,h,t))
            z_H_jht = get_variable(am, :z_H_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                u_H_jht - previous_u_H_jht - a_H_jht + z_H_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    



        else

            # parameter  
            
            # variable
            u_H_jht = get_variable(am, :u_H_jht, (j,h,t))
            previous_u_H_jht = get_variable(am, :u_H_jht, (j,h-1,t))
            a_H_jht = get_variable(am, :a_H_jht, (j,h,t))
            z_H_jht = get_variable(am, :z_H_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                u_H_jht - previous_u_H_jht - a_H_jht + z_H_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    



        end



    end

end


function HydroBoost_constraint_ES_spin_sale_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    
    # variable
    r_SR_iht = get_variable(am, :r_SR_iht, (i,h,t))
    r_SR_D_iht = get_variable(am, :r_SR_D_iht, (i,h,t))
    r_SR_C_iht = get_variable(am, :r_SR_C_iht, (i,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_SR_iht - r_SR_D_iht - r_SR_C_iht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_reg_dn_sale_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    
    # variable
    r_RD_iht = get_variable(am, :r_RD_iht, (i,h,t))
    r_RD_D_iht = get_variable(am, :r_RD_D_iht, (i,h,t))
    r_RD_C_iht = get_variable(am, :r_RD_C_iht, (i,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RD_iht - r_RD_D_iht - r_RD_C_iht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_reg_up_sale_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    
    # variable
    r_RU_iht = get_variable(am, :r_RU_iht, (i,h,t))
    r_RU_D_iht = get_variable(am, :r_RU_D_iht, (i,h,t))
    r_RU_C_iht = get_variable(am, :r_RU_C_iht, (i,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RU_iht - r_RU_D_iht - r_RU_C_iht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_SR_C_cap_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    ramp_cap = parameter(am, bus_idx, :gen_bus, tech_idx, "Maximum SP")

    # variable
    r_SR_C_iht = get_variable(am, :r_SR_C_iht, (i,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_SR_C_iht - ramp_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_RU_C_cap_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    ramp_cap = parameter(am, bus_idx, :gen_bus, tech_idx, "Maximum RU")

    # variable
    r_RU_C_iht = get_variable(am, :r_RU_C_iht, (i,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RU_C_iht - ramp_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_RD_C_cap_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    ramp_cap = parameter(am, bus_idx, :gen_bus, tech_idx, "Maximum RD")

    # variable
    r_RD_C_iht = get_variable(am, :r_RD_C_iht, (i,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
    r_RD_C_iht - ramp_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_charge_Bounds_DN_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
        
    # variable
    p_B_C_iht = get_variable(am, :p_B_C_iht, (i,h,t))
    r_RU_C_iht = get_variable(am, :r_RU_C_iht, (i,h,t))
    r_SR_C_iht = get_variable(am, :r_SR_C_iht, (i,h,t))

    u_B_iht = get_variable(am, :u_B_iht, (i,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_B_C_iht - r_RU_C_iht - r_SR_C_iht 
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_charge_Bounds_UP_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    P_C_max = parameter(am, bus_idx, :gen_bus, tech_idx, "Max_Charge_MW")

    # variable
    p_B_C_iht = get_variable(am, :p_B_C_iht, (i,h,t))
    r_RD_C_iht = get_variable(am, :r_RD_C_iht, (i,h,t))
    
    u_B_iht = get_variable(am, :u_B_iht, (i,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_B_C_iht + r_RD_C_iht - P_C_max * (1 - u_B_iht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_SR_D_cap_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    ramp_cap = parameter(am, bus_idx, :gen_bus, tech_idx, "Maximum SP")

    # variable
    r_SR_D_iht = get_variable(am, :r_SR_D_iht, (i,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_SR_D_iht - ramp_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_RD_D_cap_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    ramp_cap = parameter(am, bus_idx, :gen_bus, tech_idx, "Maximum RD")

    # variable
    r_RD_D_iht = get_variable(am, :r_RD_D_iht, (i,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RD_D_iht - ramp_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_RU_D_cap_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    ramp_cap = parameter(am, bus_idx, :gen_bus, tech_idx, "Maximum RU")

    # variable
    r_RU_D_iht = get_variable(am, :r_RU_D_iht, (i,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RU_D_iht - ramp_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_power_Bounds_DN_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  

    # variable
    p_B_D_iht = get_variable(am, :p_B_D_iht, (i,h,t))
    r_RD_D_iht = get_variable(am, :r_RD_D_iht, (i,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_B_D_iht - r_RD_D_iht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_power_Bounds_UP_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    P_B_max = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")

    # variable
    p_B_D_iht = get_variable(am, :p_B_D_iht, (i,h,t))
    r_RU_D_iht = get_variable(am, :r_RU_D_iht, (i,h,t))
    r_SR_D_iht = get_variable(am, :r_SR_D_iht, (i,h,t))

    u_B_iht = get_variable(am, :u_B_iht, (i,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_B_D_iht + r_RU_D_iht + r_SR_D_iht - P_B_max * u_B_iht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_SOC_Bounds_UP_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    BATEFF_C = parameter(am, bus_idx, :gen_bus, tech_idx, "Charging Efficiency") * 0.01 # convert to percent
    Max_SOC_MWh = parameter(am, bus_idx, :gen_bus, tech_idx, "Max_SOC_MWh")
    
    # variable
    e_B_iht = get_variable(am, :e_B_iht, (i,h,t))
    p_B_C_iht = get_variable(am, :p_B_C_iht, (i,h,t))
    r_RD_iht = get_variable(am, :r_RD_iht, (i,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_B_iht - Max_SOC_MWh + BATEFF_C * (p_B_C_iht + r_RD_iht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_SOC_Bounds_DN_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    BATEFF_D = parameter(am, bus_idx, :gen_bus, tech_idx, "Discharging Efficiency") * 0.01 # convert to percent
    Min_SOC_MWh = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_SOC_MWh")

    # variable
    e_B_iht = get_variable(am, :e_B_iht, (i,h,t))
    p_B_D_iht = get_variable(am, :p_B_D_iht, (i,h,t))
    r_RU_iht = get_variable(am, :r_RU_iht, (i,h,t))
    r_SR_iht = get_variable(am, :r_SR_iht, (i,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_B_iht - Min_SOC_MWh - (1/BATEFF_D) * (p_B_D_iht + r_RU_iht + r_SR_iht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


function HydroBoost_constraint_ES_SOC_Balance_Inter_Hour_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int, day_id, prior_day_solution; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    BATEFF_C = parameter(am, bus_idx, :gen_bus, tech_idx, "Charging Efficiency") * 0.01 # convert to percent
    BATEFF_D = parameter(am, bus_idx, :gen_bus, tech_idx, "Discharging Efficiency") * 0.01 # convert to percent
    Max_SOC_MWh = parameter(am, bus_idx, :gen_bus, tech_idx, "Max_SOC_MWh")
    Min_SOC_MWh = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_SOC_MWh")

    reg_up_signal = 0.0     # TODO
    reg_down_signal = 0.0   # TODO

    global prior_e_B_iht = 0
    if (h==1)
        if day_id == 1
            if am.setting["Simulation Setting"]["storage initialization option"] == "Minimum"
                global prior_e_B_iht = Min_SOC_MWh
            elseif am.setting["Simulation Setting"]["storage initialization option"] == "Middle"
                global prior_e_B_iht = Max_SOC_MWh / 2
            elseif am.setting["Simulation Setting"]["storage initialization option"] == "Maximum"
                global prior_e_B_iht = Max_SOC_MWh
            end
        else    # get prior day solution
            idx_string = string("(", i, ", ", 24,  ", ", t,")")
            global prior_e_B_iht = prior_day_solution["storage"][idx_string]["e_B_iht"]
        end
    else
        global prior_e_B_iht = get_variable(am, :e_B_iht, (i,h-1,t)) 
    end
    
    # variable
    e_B_iht = get_variable(am, :e_B_iht, (i,h,t))
    p_B_C_iht = get_variable(am, :p_B_C_iht, (i,h,t))
    p_B_D_iht = get_variable(am, :p_B_D_iht, (i,h,t))
    r_RU_D_iht = get_variable(am, :r_RU_D_iht, (i,h,t))
    r_RD_D_iht = get_variable(am, :r_RD_D_iht, (i,h,t))
    r_RU_C_iht = get_variable(am, :r_RU_C_iht, (i,h,t))
    r_RD_C_iht = get_variable(am, :r_RD_C_iht, (i,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_B_iht - prior_e_B_iht - BATEFF_C * (p_B_C_iht + reg_down_signal * r_RD_C_iht - reg_up_signal * r_RU_C_iht) + (1/BATEFF_D) * (p_B_D_iht + reg_up_signal * r_RU_D_iht - reg_down_signal * r_RD_D_iht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end


##################################################
#------ Define Variables 
##################################################

function HydroBoost_variable_ilht_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, ids_i, ids_l, ids_h, ids_t; bounded_lower::Bool=false, lower_bound=0, bounded_upper::Bool=false, upper_bound=0, report::Bool=true)
    
    ids = [(i, l, h, t) for (i) in ids_i for l in ids_l for h in ids_h for t in ids_t]

    var = am.var[:nw][0][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    # lower bound
    if bounded_lower
        for idx in ids
            JuMP.set_lower_bound(var[idx], lower_bound)
        end
    end

    # upper bound
    if bounded_upper
        for idx in ids
            JuMP.set_upper_bound(var[idx], upper_bound)
        end
    end

    report && add_sol_component(am, 0, category, Symbol(variable_name), ids, var)    
end



function HydroBoost_variable_ilht_binary(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, ids_i, ids_l, ids_h, ids_t; bounded_lower::Bool=false, lower_bound=0, bounded_upper::Bool=false, upper_bound=0, report::Bool=true)
    
    ids = [(i, l, h, t) for (i) in ids_i for l in ids_l for h in ids_h for t in ids_t]

    var = am.var[:nw][0][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=true
    )

    # lower bound
    if bounded_lower
        for idx in ids
            JuMP.set_lower_bound(var[idx], lower_bound)
        end
    end

    # upper bound
    if bounded_upper
        for idx in ids
            JuMP.set_upper_bound(var[idx], upper_bound)
        end
    end

    report && add_sol_component(am, 0, category, Symbol(variable_name), ids, var)    
end


function HydroBoost_variable_ht_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, ids_h, ids_t; bounded_lower::Bool=false, lower_bound=0, bounded_upper::Bool=false, upper_bound=0, report::Bool=true)
    
    ids = [(h, t) for h in ids_h for t in ids_t]

    var = am.var[:nw][0][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    # lower bound
    if bounded_lower
        for idx in ids
            JuMP.set_lower_bound(var[idx], lower_bound)
        end
    end

    # upper bound
    if bounded_upper
        for idx in ids
            JuMP.set_upper_bound(var[idx], upper_bound)
        end
    end

    report && add_sol_component(am, 0, category, Symbol(variable_name), ids, var)    
end

function HydroBoost_variable_iht_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, ids_i, ids_h, ids_t; bounded_lower::Bool=false, lower_bound=0, bounded_upper::Bool=false, upper_bound=0, report::Bool=true)
    
    ids = [(i, h, t) for (i) in ids_i for h in ids_h for t in ids_t]

    var = am.var[:nw][0][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    # lower bound
    if bounded_lower
        for idx in ids
            JuMP.set_lower_bound(var[idx], lower_bound)
        end
    end

    # upper bound
    if bounded_upper
        for idx in ids
            JuMP.set_upper_bound(var[idx], upper_bound)
        end
    end

    report && add_sol_component(am, 0, category, Symbol(variable_name), ids, var)    
end


function HydroBoost_variable_iht_binary(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, ids_i, ids_h, ids_t; bounded_lower::Bool=false, lower_bound=0, bounded_upper::Bool=false, upper_bound=0, report::Bool=true)
    
    ids = [(i, h, t) for (i) in ids_i for h in ids_h for t in ids_t]

    var = am.var[:nw][0][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=true
    )

    # lower bound
    if bounded_lower
        for idx in ids
            JuMP.set_lower_bound(var[idx], lower_bound)
        end
    end

    # upper bound
    if bounded_upper
        for idx in ids
            JuMP.set_upper_bound(var[idx], upper_bound)
        end
    end

    report && add_sol_component(am, 0, category, Symbol(variable_name), ids, var)    
end


##################################################
#------ Utility Functions 
##################################################


function export_HydroBoost_results(ALEAF_setting::Dict{String,<:Any}, network_data, daily_solutions, project_id)
    
    # build common ALEAF model instance structure to get reference data needed for reporting
    ALEAF_model_instance = initialize_model_instance_HydroBoost(ALEAF_setting["ALEAF_model_type"], network_data, 1)

    # Add reference data 
    add_ref_HydroBoost_model!(ALEAF_model_instance, ALEAF_setting)
    
    # define filename
    case_name = project_id
    times = Dates.format(now(), "yyyy_mm_dd_HH_MM") 
    output_path = network_data["output_path"]
    filename = string("ALEAF_HydroBoost_",case_name,"_",times,".json")

    # output dict
    ALEAF_solution = Dict{String, Any}()
    ALEAF_solution["model result"] = daily_solutions
    ALEAF_solution["network data"] = network_data
    ALEAF_solution["setting"] = ALEAF_setting

    # create output file (JSON format)
    # stringdata = JSON.json(ALEAF_solution)
    # open(string(output_path, filename), "w") do f
    #     write(f, stringdata)
    # end
    # Memento.info(_LOGGER, string("-- ", string(output_path, filename), ": file saved"))

    # report output
    HydroBoost_report_result_storage_dispatch(ALEAF_model_instance, ALEAF_setting, daily_solutions, project_id, network_data)
    HydroBoost_report_result_hydro_dispatch(ALEAF_model_instance, ALEAF_setting, daily_solutions, project_id, network_data)
    HydroBoost_report_result_plant_dispatch(ALEAF_model_instance, ALEAF_setting, daily_solutions, project_id, network_data)
    
end


function HydroBoost_report_result_plant_dispatch(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data)
    
    start_time = time()

    # Create and open a file
    case_name = project_id
    output_path = network_data["output_path"]
    dispatch_file_name = string("ALEAF_HydroBoost_",case_name,"__plant_dispatch.csv")
    
    # Write Label 
    dispatch_label_list = ["day", "hour", "time", 
    "p_B_DT_ht", "p_B_CT_ht", "p_GB_ht", "p_HT_ht", "p_HB_ht", "p_HG_ht", "s_ht", "u_ht", "e_H_ht", "I_ht", "LMP", "Reg Up Price", "Reg Dn Price", "Spin Price"]

    # Write Outputs
    ids_j_hydro = [(j) for (j) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO"]
    ids_i_sto = [(i) for (i) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"]

    ids_l = [(l) for (l) in 1:am.setting["Simulation Setting"]["num_hydropower_performance_segment_value"]]
    ids_h = [(h) for (h) in 1:am.setting["Simulation Setting"]["num_hours_per_day_value"]]
    ids_t = [(t) for (t) in 1:am.setting["Simulation Setting"]["num_sub_period_value"]]
    
    num_days = 365
    # num_days = 4
    
    dispatch_output_list = Array{Any}(undef, num_days*length(ids_h)*length(ids_t)+1, length(dispatch_label_list))
    dispatch_output_list[1,:] = dispatch_label_list

    # for day_id in [1, 2, 3, 4]
    for day_id in eachindex([i for i in 1:365])    
        for hour_id in ids_h
            for time_id in ids_t
                
                idx_ht = string("(", hour_id, ", ", time_id, ")")

                # get solutions
                p_B_DT_ht = daily_solutions[string(day_id)]["solution"]["storage"][idx_ht]["p_B_DT_ht"]
                p_B_CT_ht = daily_solutions[string(day_id)]["solution"]["storage"][idx_ht]["p_B_CT_ht"]
                p_GB_ht = daily_solutions[string(day_id)]["solution"]["storage"][idx_ht]["p_GB_ht"]
                p_HT_ht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_HT_ht"]
                p_HB_ht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_HB_ht"]
                p_HG_ht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_HG_ht"]
                s_ht = daily_solutions[string(day_id)]["solution"]["water"][idx_ht]["s_ht"]
                e_H_ht = daily_solutions[string(day_id)]["solution"]["water"][idx_ht]["e_H_ht"]

                u_ht = 0.0
                for j in ids_j_hydro
                    u_ht += daily_solutions[string(day_id)]["solution"]["water"][string("(", j, ", ", hour_id, ", ", time_id, ")")]["u_jht"]
                end

                I_ht = network_data["repdays"][day_id]["data"][hour_id][time_id]["Inflow"]

                LMP = network_data["repdays"][day_id]["data"][hour_id][time_id]["DA_LMP"]
                Reg_up_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Regulation_up"]
                Reg_dn_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Regulation_down"]
                Spin_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Spin"]
                
                row_id = (day_id-1)*length(ids_h)*length(ids_t) + (hour_id-1)*length(ids_t) + (time_id) + 1
                
                dispatch_output_list[row_id, :] = [day_id, hour_id, time_id, 
                    p_B_DT_ht, p_B_CT_ht, p_GB_ht, p_HT_ht, p_HB_ht, p_HG_ht, s_ht, u_ht, e_H_ht, I_ht, LMP, Reg_up_price, Reg_dn_price, Spin_price]
            end
        end
    end
        
    CSV.write(string(output_path, dispatch_file_name), Tables.table(dispatch_output_list), writeheader=false)
        
    Memento.info(_LOGGER, "[ALEAF HydroBoost Model]]:\t - Plant dispatch solution reporting, Total Time(sec): $(round(time() - start_time, digits=2))")

end


function HydroBoost_report_result_hydro_dispatch(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data)
    
    start_time = time()

    # Create and open a file
    case_name = project_id
    output_path = network_data["output_path"]
    dispatch_file_name = string("ALEAF_HydroBoost_",case_name,"__hydro_dispatch.csv")
    
    # Write Label 
    dispatch_label_list = ["day", "hour", "time", "unit_id", "UnitGroup", "Unit_Category", "Unit_Type", 
    "u_H_jht", "a_H_jht", "z_H_jht", "p_H_jht", "u_jht"]
    
    # Write Outputs
    ids_j_hydro = [(j) for (j) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO"]
    ids_i_sto = [(i) for (i) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"]

    ids_l = [(l) for (l) in 1:am.setting["Simulation Setting"]["num_hydropower_performance_segment_value"]]
    ids_h = [(h) for (h) in 1:am.setting["Simulation Setting"]["num_hours_per_day_value"]]
    ids_t = [(t) for (t) in 1:am.setting["Simulation Setting"]["num_sub_period_value"]]

    # expand the label
    push!(dispatch_label_list, string("u_", 0, "_jht"))
    for l in ids_l
        push!(dispatch_label_list, string("u_", l, "_jht"))
    end
    
    num_days = 365
    # num_days = 4
    
    dispatch_output_list = Array{Any}(undef, num_days*length(ids_i_sto)*length(ids_h)*length(ids_t)+1, length(dispatch_label_list))
    dispatch_output_list[1,:] = dispatch_label_list

    # for day_id in [1, 2, 3, 4]
    for day_id in eachindex([i for i in 1:365])    
        for hour_id in ids_h
            for time_id in ids_t
                for unit_idx in eachindex(ids_j_hydro)

                    i = ids_j_hydro[unit_idx]

                    idx_idhty = string("(", i, ", ", hour_id, ", ", time_id, ")")

                    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
                    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
                    
                    # unit info
                    UNITGROUP = parameter(am, bus_idx, :gen_bus, tech_idx, "UNITGROUP")
                    Unit_Category = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_CATEGORY")
                    Unit_Type = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE")
                    
                    # get solutions
                    u_H_jht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_idhty]["u_H_jht"]
                    a_H_jht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_idhty]["a_H_jht"]
                    z_H_jht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_idhty]["z_H_jht"]
                    p_H_jht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_idhty]["p_H_jht"]
                    u_jht = daily_solutions[string(day_id)]["solution"]["water"][idx_idhty]["u_jht"]
                    
                    row_id = (day_id-1)*length(ids_j_hydro)*length(ids_h)*length(ids_t) + (hour_id-1)*length(ids_j_hydro)*length(ids_t) + (time_id-1)*length(ids_j_hydro) + (unit_idx) + 1
                    
                    output_list = [day_id, hour_id, time_id, i, UNITGROUP, Unit_Category, Unit_Type, 
                    u_H_jht, a_H_jht, z_H_jht, p_H_jht, u_jht]

                    # minimum water release 
                    Min_Water_Release = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Water_Release")
                    u_0_jht = u_H_jht * Min_Water_Release
                    push!(output_list, u_0_jht)

                    # water release blocks
                    for l in ids_l
                        push!(output_list, daily_solutions[string(day_id)]["solution"]["water"][string("(", i, ", ", l, ", ", hour_id, ", ", time_id, ")")]["u_l_jht"])
                    end

                    dispatch_output_list[row_id, :] = output_list

                end
            end
        end
    end
        
    CSV.write(string(output_path, dispatch_file_name), Tables.table(dispatch_output_list), writeheader=false)
        
    Memento.info(_LOGGER, "[ALEAF HydroBoost Model]]:\t - Hydro dispatch solution reporting, Total Time(sec): $(round(time() - start_time, digits=2))")

end


function HydroBoost_report_result_storage_dispatch(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data)
    
    start_time = time()

    # Create and open a file
    case_name = project_id
    output_path = network_data["output_path"]
    dispatch_file_name = string("ALEAF_HydroBoost_",case_name,"__storage_dispatch.csv")
    
    # Write Label 
    dispatch_label_list = ["day", "hour", "time", "unit_id", "UnitGroup", "Unit_Category", "Unit_Type", 
    "u_B_iht", "p_B_D_iht", "p_B_C_iht", "e_B_iht", "r_RU_D_iht", "r_RD_D_iht", "r_RU_C_iht", "r_RD_C_iht", "r_RU_iht", "r_RD_iht", "r_SR_D_iht", "r_SR_C_iht", "r_SR_iht"]
    
    # Write Outputs
    ids_j_hydro = [(j) for (j) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO"]
    ids_i_sto = [(i) for (i) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"]

    ids_l = [(l) for (l) in 1:am.setting["Simulation Setting"]["num_hydropower_performance_segment_value"]]
    ids_h = [(h) for (h) in 1:am.setting["Simulation Setting"]["num_hours_per_day_value"]]
    ids_t = [(t) for (t) in 1:am.setting["Simulation Setting"]["num_sub_period_value"]]
    
    num_days = 365
    # num_days = 4
    
    dispatch_output_list = Array{Any}(undef, num_days*length(ids_i_sto)*length(ids_h)*length(ids_t)+1, length(dispatch_label_list))
    dispatch_output_list[1,:] = dispatch_label_list

    # for day_id in [1, 2, 3, 4]
    for day_id in eachindex([i for i in 1:365])    
        for hour_id in ids_h
            for time_id in ids_t
                for unit_idx in eachindex(ids_i_sto)

                    i = ids_i_sto[unit_idx]

                    idx_idhty = string("(", i, ", ", hour_id, ", ", time_id, ")")

                    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
                    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
                    
                    # unit info
                    UNITGROUP = parameter(am, bus_idx, :gen_bus, tech_idx, "UNITGROUP")
                    Unit_Category = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_CATEGORY")
                    Unit_Type = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE")
                    
                    # get solutions
                    u_B_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["u_B_iht"]
                    p_B_D_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["p_B_D_iht"]
                    p_B_C_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["p_B_C_iht"]
                    e_B_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["e_B_iht"]
                    r_RU_D_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["r_RU_D_iht"]
                    r_RD_D_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["r_RD_D_iht"]
                    r_RU_C_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["r_RU_C_iht"]
                    r_RD_C_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["r_RD_C_iht"]
                    r_RU_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["r_RU_iht"]
                    r_RD_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["r_RD_iht"]
                    r_SR_D_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["r_SR_D_iht"]
                    r_SR_C_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["r_SR_C_iht"]
                    r_SR_iht = daily_solutions[string(day_id)]["solution"]["storage"][idx_idhty]["r_SR_iht"]

                    row_id = (day_id-1)*length(ids_i_sto)*length(ids_h)*length(ids_t) + (hour_id-1)*length(ids_i_sto)*length(ids_t) + (time_id-1)*length(ids_i_sto) + (unit_idx) + 1
                    
                    dispatch_output_list[row_id, :] = [day_id, hour_id, time_id, i, UNITGROUP, Unit_Category, Unit_Type, 
                        u_B_iht, p_B_D_iht, p_B_C_iht, e_B_iht, r_RU_D_iht, r_RD_D_iht, r_RU_C_iht, r_RD_C_iht, r_RU_iht, r_RD_iht, r_SR_D_iht, r_SR_C_iht, r_SR_iht]

                end
            end
        end
    end
        
    CSV.write(string(output_path, dispatch_file_name), Tables.table(dispatch_output_list), writeheader=false)
        
    Memento.info(_LOGGER, "[ALEAF HydroBoost Model]]:\t - Storage dispatch solution reporting, Total Time(sec): $(round(time() - start_time, digits=2))")

end


function solve_model_HydroBoost!(JuMP_model::JuMP.AbstractModel, solution_list::Dict{Symbol,<:Any}, solver_setting)
    
    solver_name = solver_setting["solver_name"]
    sol_iteration = 1
    if solver_name  == "CPLEX"
        optimizer = CPLEX.Optimizer
        if solver_setting["1"]["Value"] == false  # CPLEX needs direct mode
            JuMP.set_optimizer(JuMP_model, optimizer)
        end
    elseif solver_name == "HiGHS"
        optimizer = HiGHS.Optimizer
        JuMP.set_optimizer(JuMP_model, optimizer)     
    end

    "set solver setting"
    for setting_id in keys(solver_setting)
        if setting_id ∉ ["optimizer", "solver_name", "1"]
            if solver_setting[setting_id]["Flag"] == true
                JuMP.set_optimizer_attribute(JuMP_model, solver_setting[setting_id]["Parameter"], solver_setting[setting_id]["Value"])
            end
        end
    end

    "solve"
    start_time = time()
    try
        _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(JuMP_model, )
    catch
        println("[ALEAF HydroBoost Model]:\tERROR found in the optimizer")
    end
    solve_time = time() - start_time

    "check solution"
    if (JuMP.termination_status(JuMP_model) == _MOI.OPTIMAL) || (JuMP.termination_status(JuMP_model) == _MOI.LOCALLY_SOLVED)
        result = collect_result_distributed(JuMP_model, solution_list, solve_time)
        
        return result 
        
    elseif JuMP.termination_status(JuMP_model) == _MOI.TIME_LIMIT
        if JuMP.has_values(JuMP_model)
           
            result = collect_result_distributed(JuMP_model, solution_list, solve_time)
            try 
                gap = round(result["relative_gap"], digits=5)
            catch
                gap = 100.0
            end
            
            println("[ALEAF HydroBoost Model]:\tTime Limit. Sub-optimal solution exists; solution_time: $solve_time, gap: $gap")
            return result
        else
            println("[ALEAF HydroBoost Model]:\tTime Limit Reached without solution; solution_time: $solve_time")
            result = Dict()
            return result
        end
    else
        if JuMP.has_values(JuMP_model)
            result = collect_result_distributed(JuMP_model, solution_list, solve_time)
            gap = 0.0
            try 
                gap = round(result["relative_gap"], digits=5)
            catch
            end
            println("[ALEAF HydroBoost Model]:\tSub-optimal solution exists; solution_time: $solve_time, gap: $gap")
            return result
        else
            println("[ALEAF HydroBoost Model]:\tFailed to find a solution; solution_time: $solve_time")
            result = Dict()
            return result
        end
    end
end


function export_JuMP_model_lp_file(JuMP_model)

    file_name = "HydroBoost_test.lp"
    JuMP.write_to_file(JuMP_model, file_name)
    Memento.info(_LOGGER, "[ALEAF HydroBoost Model]:\tExport HydroBoost_model_instance model (lp format)")

end


function get_variable(am::Abstract_ALEAF_Model, key1::Symbol, idx)
    return am.var[:nw][0][key1][idx]
end


function add_ref_HydroBoost_model!(am::Abstract_ALEAF_Model, ALEAF_setting::Dict{String, Any})
    
    # model setting
    add_model_setting_HydroBoost!(am, ALEAF_setting)
    
    # run period
    add_run_period_HydroBoost!(am)

    # plant data update 
    update_individual_plant_info_HydroBoost!(am) 
    
    # add gen index
    add_gen_index_ref_HydroBoost!(am)

end


function add_gen_index_ref_HydroBoost!(am::Abstract_ALEAF_Model)

    # add gen_index
    am.ref[:nw][0][:gen_index] = Dict{Int64, Any}()
    gen_idx = 1
    num_sub_area = 1 # Assuming single bus system (might be extended in the future)
    
    for sub_area_idx in 1:num_sub_area
        for genco_tech_id in sort!(collect(keys(am.ref[:nw][sub_area_idx][:gen_bus])))
            UNIT_TYPE = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["UNIT_TYPE"]
            UNIT_GROUP = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["UNITGROUP"]
            UNIT_CATEGORY = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["UNIT_CATEGORY"]
            CAP = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["PMAX"]
            am.ref[:nw][0][:gen_index][gen_idx] = Dict{String, Any}("bus_idx" => sub_area_idx, 
                                                                    "genco_tech_id" => genco_tech_id, 
                                                                    "genco_tech_UNIT_Type" => UNIT_TYPE, 
                                                                    "UNIT_GROUP" => UNIT_GROUP,
                                                                    "UNIT_CATEGORY" => UNIT_CATEGORY,
                                                                    "CAP" => CAP
                                                                    )
            am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["gen_idx"] = gen_idx
            gen_idx += 1
        end
    end

    am.ref[:nw][0][:gen_index] = sort(am.ref[:nw][0][:gen_index])

    # remove unnecessary dics
    for idx in keys(am.ref[:nw])
        if idx != 0
            delete!(am.ref[:nw][idx], :gen_technology_hydro)
            delete!(am.ref[:nw][idx], :gen_technology_BESS)
        end
    end

end


function update_individual_plant_info_HydroBoost!(am::Abstract_ALEAF_Model)

    # Assuming single bus system (might be extended in the future)
    sub_area_id = 1
    
    # add sub_area dict first
    am.ref[:nw][sub_area_id] = Dict{Symbol, Any}()
    am.ref[:nw][sub_area_id][:gen_bus] = Dict{Int64, Any}()

    # allocate plant data
    global plant_number = 1
    # 1) hydro
    for plant_id in keys(am.ref[:nw][0][:gen_technology_hydro])
        am.ref[:nw][sub_area_id][:gen_bus][plant_number] = deepcopy(am.ref[:nw][0][:gen_technology_hydro][plant_id])
        global plant_number += 1
    end
    # 2) BESS
    for plant_id in keys(am.ref[:nw][0][:gen_technology_BESS])
        am.ref[:nw][sub_area_id][:gen_bus][plant_number] = deepcopy(am.ref[:nw][0][:gen_technology_BESS][plant_id])
        global plant_number += 1
    end
end


function add_run_period_HydroBoost!(am::Abstract_ALEAF_Model)

    "Hour"
    num_hours_per_day = am.setting["Simulation Setting"]["num_hours_per_day_value"] * am.setting["Simulation Setting"]["look_ahead_days_value"]
    index_list = [] 
    for h = 1:num_hours_per_day
        append!(index_list, h)
    end
    am.setting["run_H"] = index_list

    "Sub-Hour"
    index_list = []
    append!(index_list, 1)
    am.setting["run_T"] = index_list
        
end


function add_model_setting_HydroBoost!(am::Abstract_ALEAF_Model, ALEAF_setting::Dict{String, Any}; nw::Int=am.cnw)

    # Add LC_GEP setting to the model instance
    am.setting["Simulation Setting"] = deepcopy(ALEAF_setting["Simulation Setting"])

    # Solver setting
    if ALEAF_setting["Simulation Setting"]["solver_name"] == "CPLEX"
        am.setting["Solver Setting"] = ALEAF_setting["CPLEX Setting"]
        am.setting["Solver Setting"]["solver_name"] = "CPLEX"
        am.setting["Solver Setting"]["optimizer"] = CPLEX.Optimizer
    elseif ALEAF_setting["Simulation Setting"]["solver_name"] == "HiGHS"
        am.setting["Solver Setting"] = ALEAF_setting["HiGHS Setting"]
        am.setting["Solver Setting"]["solver_name"] = "HiGHS"
        am.setting["Solver Setting"]["optimizer"] = HiGHS.Optimizer    
    end
end


function allocate_hydro_sub_hourly_timeseries_data!(timeSeries, network_data, ALEAF_setting, data_label_list)
    
    # rename water release requirment 
    rename!(timeSeries, "Water Release Requirment" => "Water_Release_Requirement")
    
    # convert date-time string to date time format
    timeSeries[!, "Date-time"] = DateTime.(timeSeries[!, "Date-Time"], "yyyy-mm-dd HH:MM:SS")

    # add month and day info 
    timeSeries[!, :Month] = Dates.month.(timeSeries[!, "Date-time"])
    timeSeries[!, :Day] = Dates.day.(timeSeries[!, "Date-time"])
    timeSeries[!, :Hour] = Dates.hour.(timeSeries[!, "Date-time"])
    
    num_hours_per_day = ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    num_sub_period = ALEAF_setting["Simulation Setting"]["num_sub_period_value"]
    
    for day_id in keys(network_data["repdays"])

        global hour_idx = 1
        if !haskey(network_data["repdays"][day_id], "data") network_data["repdays"][day_id]["data"] = Dict{Any, Any}() end
        if !haskey(network_data["repdays"][day_id]["data"], hour_idx) network_data["repdays"][day_id]["data"][hour_idx] = Dict{Any, Any}() end
        
        for sim_day_id in network_data["repdays"][day_id]["look_ahead_simulation_day_idx_list"]

            month = network_data["repdays"][sim_day_id]["Month"]
            day = network_data["repdays"][sim_day_id]["Day"]

            # Select timeseries data
            sub_dataframe = timeSeries[(timeSeries.Month.==month).&(timeSeries.Day.==day),:]

            # allocate/update timeseries data (hour/sub-hour)
            daily_data_list = Dict()
            for data_label_idx in data_label_list
                daily_data_list[data_label_idx] = sub_dataframe[:, data_label_idx]
            end

            for h = 1:num_hours_per_day
                for t = 1:num_sub_period
                    if !haskey(network_data["repdays"][day_id]["data"][h] , t) network_data["repdays"][day_id]["data"][h][t] = Dict{Any, Any}() end
                    for data_label_idx in data_label_list
                        network_data["repdays"][day_id]["data"][hour_idx][t][data_label_idx] = daily_data_list[data_label_idx][h]
                    end
                    global hour_idx += 1
                end
            end
            daily_data_list = 0
            
        end
    end

end


function allocate_hydro_daily_timeseries_data!(timeSeries, network_data, data_label_list)

    # add month and day info 
    timeSeries[!, :Month] = Dates.month.(timeSeries[!, "Date-time"])
    timeSeries[!, :Day] = Dates.day.(timeSeries[!, "Date-time"])
    
    for day_id in keys(network_data["repdays"])

        network_data["repdays"][day_id]["daily_data"] = Dict{Any, Any}()

        for sim_day_id in network_data["repdays"][day_id]["look_ahead_simulation_day_idx_list"]

            month = network_data["repdays"][sim_day_id]["Month"]
            day = network_data["repdays"][sim_day_id]["Day"]

            # Select timeseries data
            sub_dataframe = timeSeries[(timeSeries.Month.==month).&(timeSeries.Day.==day),:]

            # allocate/update timeseries data (daily)
            for data_label_idx in data_label_list
                if !haskey(network_data["repdays"][day_id]["daily_data"], data_label_idx) network_data["repdays"][day_id]["daily_data"][data_label_idx] = [] end
                push!(network_data["repdays"][day_id]["daily_data"][data_label_idx],  sub_dataframe[:, data_label_idx][1])
            end
        end
    end

end


function allocate_market_price_sub_hourly_timeseries_data!(timeSeries, network_data, ALEAF_setting, data_label)

    num_hours_per_day = ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    num_sim_hours_per_day = ALEAF_setting["Simulation Setting"]["look_ahead_days_value"] * num_hours_per_day
    num_sub_period = ALEAF_setting["Simulation Setting"]["num_sub_period_value"]

    for day_id in keys(network_data["repdays"])
        # Select timeseries data
        sub_dataframe = timeSeries[:,day_id+1]  # add 1 due to the index

        # allocate/update timeseries data (hour/sub-hour)
        id_dic = 1
        if !haskey(network_data["repdays"][day_id], "data") network_data["repdays"][day_id]["data"] = Dict{Any, Any}() end
        for h = 1:num_sim_hours_per_day
            if !haskey(network_data["repdays"][day_id]["data"], h) network_data["repdays"][day_id]["data"][h] = Dict{Any, Any}() end
            for t = 1:num_sub_period
                if !haskey(network_data["repdays"][day_id]["data"][h] , t) network_data["repdays"][day_id]["data"][h][t] = Dict{Any, Any}() end
                network_data["repdays"][day_id]["data"][h][t][data_label] = sub_dataframe[id_dic]
                id_dic += 1
            end
        end
    end

end


function read_ALEAF_HydroBoost_setting(project_id::String)

    ALEAF_setting = Dict{String, Any}()

    setting_file_location = string("Core_Models/HydroBoost/generated_data/", project_id, "/")
    setting_file_name = string(project_id, ".xlsx")

    ALEAF_setting["Simulation Setting"] = Dict{String, Any}()
    data = DataFrame(XLSX.readtable(string(setting_file_location, setting_file_name), "Simulation Setting"))
    for row in eachrow(data)
        ALEAF_setting["Simulation Setting"][values(row)[1]] = values(row)[2]
    end

    if ALEAF_setting["Simulation Setting"]["solver_name"] == "HiGHS"
        ALEAF_setting["HiGHS Setting"] = read_xlsx_return_dict_string_any(string(setting_file_location, setting_file_name), "HiGHS Setting")
    elseif ALEAF_setting["Simulation Setting"]["solver_name"] == "CPLEX"
        ALEAF_setting["CPLEX Setting"] = read_xlsx_return_dict_string_any(string(setting_file_location, setting_file_name), "CPLEX Setting")
    end

    ALEAF_setting["PTC per Tech"] = read_xlsx_return_dict_string_any(string(setting_file_location, setting_file_name), "PTC per Tech")

    tab_names_with_category_list = ["Gen Technology - BESS", "Gen Technology - Hydro", "Hydro Resorvior"]
    for data_category in tab_names_with_category_list
        ALEAF_setting[data_category] = read_xlsx_return_dict_string_any(string(setting_file_location, setting_file_name), data_category; first_row_value = 2)
    end
    
    return ALEAF_setting
end


function generate_networkdata_HydroBoost(ALEAF_setting::Dict{String,<:Any}, project_id; year_id::Int64=1)

    # Generate network data
    network_data = Dict{String, Any}()
    data_location = string(pwd(), "/Core_Models/HydroBoost/generated_data/", project_id, "/")
    
    # Generate network_data for a given project_id
    network_data["case_name"] = project_id
    network_data["model_type"] = "MILP"

    # Gen Technology
    network_data["gen_technology_hydro"] = deepcopy(ALEAF_setting["Gen Technology - Hydro"])
    network_data["gen_technology_BESS"] = deepcopy(ALEAF_setting["Gen Technology - BESS"])
    network_data["hydro_resorvior"] = deepcopy(ALEAF_setting["Hydro Resorvior"])

    # day group setting for each month 
    num_days_each_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    total_num_days = 365
    if ALEAF_setting["Simulation Setting"]["simulation_year_value"] % 4 == 0 # check leap year
        num_days_each_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        total_num_days = 366
    end

    # day group definition
    network_data["repdays"] = Dict()
    num_days_per_simulation = ALEAF_setting["Simulation Setting"]["look_ahead_days_value"] 
    global start_day_id = 1
    global month = 1
    for day_id in 1:total_num_days
        network_data["repdays"][day_id] = Dict()
        network_data["repdays"][day_id]["Month"] = month
        network_data["repdays"][day_id]["Day"] = start_day_id
        
        network_data["repdays"][day_id]["look_ahead_simulation_day_idx_list"] = []
        for sim_day_id in 1:num_days_per_simulation
            if day_id + sim_day_id - 1 <= total_num_days
                push!(network_data["repdays"][day_id]["look_ahead_simulation_day_idx_list"], day_id + sim_day_id - 1)
            else
                push!(network_data["repdays"][day_id]["look_ahead_simulation_day_idx_list"], day_id + sim_day_id - 1 - total_num_days)
            end
        end
        
        if start_day_id < num_days_each_month[month]
            global start_day_id += 1
        else
            global start_day_id = 1
            global month += 1
        end
    end

    # Collect market price time-series data ["DA_LMP", "Regulation_down", "Regulation_up", "Spin"]
    
    forecasting_method = ALEAF_setting["Simulation Setting"]["Market_price_forecasting_method"]     # Perfect_foresight, Mean_persistence, Additive_model_with_regressors, Additive_model_no_regressors, Autoregressive_with_regressors, Autoregressive_no_regressors, Manual_Forecast

    if forecasting_method in ["Perfect_foresight", "Mean_persistence", "Manual_Forecast"]
        
        market_price_time_series_data_list = ["DA_LMP", "Regulation_down", "Regulation_up", "Spin"]
        market_price_data_path = string(data_location, "Market/", forecasting_method, "/")

        for data_label in market_price_time_series_data_list
            timeSeries_df = CSV.read(string(market_price_data_path, data_label, ".csv"), DataFrame)
            allocate_market_price_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, data_label)
        end

    else

        # read LMP
        timeSeries_df = CSV.read(string(data_location, "Market/", forecasting_method, ".csv"), DataFrame)
        allocate_market_price_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, "DA_LMP")
        
        # read AS prices from Mean_persistence folder
        market_price_time_series_data_list = ["Regulation_down", "Regulation_up", "Spin"]
        market_price_data_path = string(data_location, "Market/Mean_persistence/")

        for data_label in market_price_time_series_data_list
            timeSeries_df = CSV.read(string(market_price_data_path, data_label, ".csv"), DataFrame)
            allocate_market_price_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, data_label)
        end

    end

    
    
    # collect hydro time-series data ["Hydro - Daily", "Hydro - Hourly"]
    hydro_data_path = string(data_location, "Hydro/")
    
    # --- "Hydro - Daily"
    timeSeries_df = CSV.read(string(hydro_data_path, "Daily_flow.csv"), DataFrame)
    allocate_hydro_daily_timeseries_data!(timeSeries_df, network_data, ["Vmin", "Vmax"])

    # --- "Hydro - Hourly"
    timeSeries_df = CSV.read(string(hydro_data_path, "Hourly_flow.csv"), DataFrame)

    allocate_hydro_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, ["Inflow", "Diversion", "Water_Release_Requirement"])
    
    #------------------------------------------

    # Define output path
    output_path = string(pwd(), "/Simulation_Results/")
    check_and_create_path(output_path)

    output_path = string(output_path, project_id, "/")
    check_and_create_path(output_path)

    network_data["output_path"] = output_path
    
    return network_data

end


