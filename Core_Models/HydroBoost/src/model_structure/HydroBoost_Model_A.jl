#=

HydroBoost Model A: Hydroelectric System Without Reservoir

Main Authors:
Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov
Carlos Josue Lopez; Argonne National Laboratory; clopezsalgado@anl.gov
Alberto Grimaldi; Argonne National Laboratory; agrimaldi@anl.gov

Current version: 2.0
Last update: 07.31.2025

=#

using XLSX
using DataFrames
using Infiltrator
using Dates
using JSON
using Statistics
using DelimitedFiles
using Base.Iterators: collect

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
#------ Define HydroBoost optimization model
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

    ###################################
    #------ Define sets & indices
    ###################################

    ids_i = [(i) for (i) in get_index(am, :gen_index, 0)]
    ids_j_hydro = [(j) for (j) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) === "HYDRO"]
    ids_i_sto = [(i) for (i) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) === "STORAGE"]
    
    ids_k_ren_idx = [k for k in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", k) == "RENEWABLE"]
    id_to_name = Dict{Int,String}()
    for k in ids_k_ren_idx
        bus_idx  = parameter(am, 0, :gen_index, k, "bus_idx")
        tech_idx = parameter(am, 0, :gen_index, k, "genco_tech_id")
        id_to_name[k] = parameter(am, bus_idx, :gen_bus, tech_idx, "Tech_ID")
    end
    ids_k_ren = [id_to_name[k] for k in ids_k_ren_idx ]

    ids_l = [(l) for (l) in 1:am.setting["Simulation Setting"]["num_hydropower_performance_segment_value"]]
    ids_h = [(h) for (h) in am.setting["run_H"]]
    ids_t = [(t) for (t) in am.setting["run_T"]]

    ###################################
    #------ Define decision variables (TOT = 39 variables)
    ###################################
    
    # Storage plant (BESS): dispatch variables
    HydroBoost_variable_iht_binary(JuMP_model, am, :storage, "u_B_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0,  bounded_upper=true, upper_bound=1.0)               # Binary variable driven to 1 when BESS is set to charging mode, and 0 otherwise      
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "e_B_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                     # Storage device state of charge in hour t [MWh]

    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "p_B_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Power discharged from BESS and accounted for at point of delivery in period t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "p_B_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Power contributing to charge BESS in period t before accounting for losses [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :storage, "p_B_DT_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                               # Power discharged from BESS and accounted for at point of delivery in period t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :storage, "p_B_CT_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                               # Power contributing to charge BESS in period t before accounting for losses [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :storage, "p_Gr_St_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                              # Power from grid allocated to charge storage (BESS) (before considering any round-trip losses) [MWh]
    
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for regulation up in discharging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for regulation up in charging mode provided by BESS unit i in hour t  [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                    # Reserve for regulation up provided by BESS unit i in hour t [MWh]
    
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for regulation down in discharging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for regulation down in charging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                    # Reserve for regulation down provided by BESS unit i in hour t [MWh]
    
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for spinning reserve in discharging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for spinning reserve in charging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                    # Reserve for spinning reserve provided by BESS unit i in hour t [MWh]
    
    # r_RU_ht: Reserve for regulation up provided by all BESS units in hour t [MWh] is already included in the objective function
    # r_RD_ht: Reserve for regulation down provided by all BESS units in hour t [MWh] is already included in the objective function
    # r_SR_ht: Reserve for spinning reserve provided by all BESS units in hour t [MWh] is already included in the objective function

    # Hydro plant: power dispatch variables    
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "u_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)             # Binary variable which is equal to 1 if hydro generator j is on-line in hour t, and 0 otherwise
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "a_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)             # Binary variable which is equal to 1 if hydro generator j is started at beginning of hour t, and 0 otherwise
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "z_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)             # Binary variable which is equal to 1 if hydro generator j is shut down at beginning of hour t, and 0 otherwise
    
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "p_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                     # Power output of hydro generator j in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_G_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                    # Day-ahead total power setpoint of dispatchable hydro plant in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_GB_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                   # Hydro power allocated to charge the BESS in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_G_Gr_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                 # Dispatchable hydro power sold to the grid in the day-ahead market in hour t [MWh] 

    # Hydro plant: water flow variables
    HydroBoost_variable_ilht_real(JuMP_model, am, :water, "q_Gl_jht", ids_j_hydro, ids_l, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                            # Water discharge of block ℓ of hydro generator j in hour t [ft^3/s]
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "q_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                     # Water discharge of hydro generator j in hour t [ft^3/s]
    HydroBoost_variable_ilht_binary(JuMP_model, am, :water, "w_Gl_jht", ids_j_hydro, ids_l, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)    # Binary variable which is equal to 1 if water discharged by hydro generator j has exceeded block ℓ in hour t
    HydroBoost_variable_ht_real(JuMP_model, am, :water, "s_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                      # Spillage of reservoir in hour t [ft^3/s]
    HydroBoost_variable_ht_real(JuMP_model, am, :water, "e_H_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                    # Water volume of reservoir in period t [A-F]

    # Hydro plant: rough zone variables
    if am.setting["Simulation Setting"]["hydropower_rough_zone_flag"] == true
        HydroBoost_variable_iht_real(JuMP_model, am, :rough_zone, "y_Gl_plus_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                      # Slack variables (upper bound) which take non-zero value if operation outside rough zone ℓ of unit j cannot be honored in hour t [ft^3/s]
        HydroBoost_variable_iht_real(JuMP_model, am, :rough_zone, "y_Gl_minus_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                     # Slack variables (lower bound) which take non-zero value if operation outside rough zone ℓ of unit j cannot be honored in hour t [ft^3/s]
        HydroBoost_variable_iht_binary(JuMP_model, am, :rough_zone, "phi_Gl_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0) # Auxiliary binary variable required for representation of rough zone constraint
    end

    # Non-dispatchable generators (RES plant): variables
    HydroBoost_variable_iht_real(JuMP_model, am, :renewable, "p_R_kht", ids_k_ren, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Power output of non-dispatchable generator 𝑘 in hour 𝑡 [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :renewable, "p_R_Gr_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                             # Power of non-dispatchable generators sent to grid in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :renewable, "p_Rspill_kht", ids_k_ren, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                              # Curtailed power of non-dispatchable generator 𝑘 in hour 𝑡 [MWh]
    HydroBoost_variable_ht_real(JuMP_model,am, :renewable, "p_R_St_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                              # Power of non−dispatchable generators allocated to storage (BESS) in hour 𝑡 [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :renewable, "p_R_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                # Power of all non-dispatchable generators in hour t [MWh]
    
    ###################################
    #------ Call of constraints in the build-loop
    ###################################

    ### BESS constraints: from (6) to (28) ###
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

    ### Hydro power system constraints: from (29) to (48) ###
    for h in ids_h
        for t in ids_t

            # Generation of hydro power generators constraints
            HydroBoost_constraint_total_hydro_power_generation_jht(JuMP_model, am, "HydroBoost_constraint_total_hydro_power_generation_jht", ids_j_hydro, h, t; const_name_flag)

            for j in ids_j_hydro
                # Commitment of hydro power generators constraints
                HydroBoost_constraint_hydro_commitment_status_jht(JuMP_model, am, "HydroBoost_constraint_hydro_commitment_status_jht", j, h, t, day_id, prior_day_solution; const_name_flag)
                HydroBoost_constraint_hydro_start_up_shut_down_bound_jht(JuMP_model, am, "HydroBoost_constraint_hydro_start_up_shut_down_bound_jht", j, h, t; const_name_flag)

                # Generation of hydro power generators constraints
                HydroBoost_constraint_hydro_power_generation_jht(JuMP_model, am, "HydroBoost_constraint_hydro_power_generation_jht", ids_l, j, h, t; const_name_flag)
            end

            # Water discharge of plant constraints
            for j in ids_j_hydro

                HydroBoost_constraint_hydro_total_water_use_jht(JuMP_model, am, "HydroBoost_constraint_hydro_total_water_use_jht", ids_l, j, h, t; const_name_flag)

                for l in ids_l
                    HydroBoost_constraint_hydro_power_water_discharge_bounds_up_ljht(JuMP_model, am, "HydroBoost_constraint_hydro_power_water_discharge_bounds_up_ljht", j, l, h, t; const_name_flag)
                    HydroBoost_constraint_hydro_power_water_discharge_bounds_down_ljht(JuMP_model, am, "HydroBoost_constraint_hydro_power_water_discharge_bounds_down_ljht", j, l, h, t; const_name_flag)
                end
            end
        end
    end

    # Water discharge constraints
    for h in ids_h
        for t in ids_t

            # Spillage of reservoir constraints
            HydroBoost_constraint_hydro_power_water_spillage_ht(JuMP_model, am, "HydroBoost_constraint_hydro_power_water_spillage_ht", ids_j_hydro, h, t; const_name_flag)

            # Water balance constraints
            HydroBoost_constraint_hydro_water_balance_ht(JuMP_model, am, "HydroBoost_constraint_hydro_water_balance_ht", ids_j_hydro, h, t, day_id, prior_day_solution; const_name_flag)

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

    ### RES generators constraints: from (49) to (51) ###
    for h in ids_h
        for t in ids_t
            for k in ids_k_ren
                HydroBoost_constraint_RES_balance_kht(JuMP_model, am,"HydroBoost_constraint_RES_balance_kht", k, h, t; const_name_flag)
            end
        end
    end

    for h in ids_h
        for t in ids_t
            HydroBoost_constraint_RES_total_ht(JuMP_model, am, "HydroBoost_constraint_RES_total_ht", ids_k_ren, h, t; const_name_flag)
        end
    end

    ### Coupling constraints: from (52) to (56)###
    for h in ids_h
        for t in ids_t
            HydroBoost_constraint_Hydro_power_generation_coupling_ht(JuMP_model, am, "HydroBoost_constraint_Hydro_power_generation_coupling_ht", h, t; const_name_flag)
            HydroBoost_constraint_RES_power_generation_coupling_ht(JuMP_model, am, "HydroBoost_constraint_RES_power_generation_coupling_ht", h, t; const_name_flag)
            HydroBoost_constraint_ES_power_generation_coupling_ht(JuMP_model, am, "HydroBoost_constraint_ES_power_generation_coupling_ht", h, t; const_name_flag)
        end
    end

    for h in ids_h
        for t in ids_t
            HydroBoost_constraint_outflow_interconnection_limit_ht(JuMP_model, am, "HydroBoost_constraint_outflow_interconnection_limit_ht", h, t; const_name_flag)
            HydroBoost_constraint_inflow_interconnection_limit_ht(JuMP_model, am, "HydroBoost_constraint_inflow_interconnection_limit_ht", h, t; const_name_flag)
        end
    end

    ###################################
    #------ Call of the objective function in the build-loop
    ###################################
    HydroBoost_objective_function(JuMP_model, am, ids_j_hydro, ids_i_sto, ids_h, ids_t, ids_l)

    # export_JuMP_model_lp_file(JuMP_model) 

    am.model[:nw][0] = JuMP_model

end

###################################
#------ Define objective function
###################################

### Objective Function ###

# Constraints (1-5)

function HydroBoost_objective_function(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, ids_j_hydro, ids_i_sto, ids_h, ids_t, ids_l)

    objective = JuMP.AffExpr(0.0)    

    # Energy revenue 
    for h in ids_h
        for t in ids_t
            JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :p_G_Gr_ht, (h,t)))
            JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :p_B_DT_ht, (h,t)))
            JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"],get_variable(am, :p_R_Gr_ht, (h,t)))
            JuMP.add_to_expression!(objective, -1.00001 * am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :p_Gr_St_ht, (h,t)))
        end
    end

    # AS revenue
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
                JuMP.add_to_expression!(objective, - START_UP_COST, get_variable(am, :a_G_jht, (j,h,t)))
                JuMP.add_to_expression!(objective, - SHUT_DN_COST, get_variable(am, :z_G_jht, (j,h,t)))
            end
        end
    end

    # Rough zone operation penalty
    if am.setting["Simulation Setting"]["hydropower_rough_zone_flag"] == true
        for j in ids_j_hydro
            for h in ids_h
                for t in ids_t
                    JuMP.add_to_expression!(objective, - am.setting["Simulation Setting"]["hyropower_rough_zone_operation_penalty_value"], get_variable(am, :y_Gl_plus_jht, (j,h,t)))
                    JuMP.add_to_expression!(objective, - am.setting["Simulation Setting"]["hyropower_rough_zone_operation_penalty_value"], get_variable(am, :y_Gl_minus_jht, (j,h,t)))
                end
            end
        end
    end
    
    # Define objective function
    return JuMP.@objective(JuMP_model, Max, objective     
    )

end

##################################################
#------ Define constraints
##################################################

### Battery Energy Storage System (BESS) ###

# (a) BESS energy balance equation, considering losses and regulation signal (deployment of up and down regulation):

# Constraint (6)
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
        e_B_iht - prior_e_B_iht - (BATEFF_C * (p_B_C_iht + reg_down_signal * r_RD_C_iht - reg_up_signal * r_RU_C_iht)) + ((1/BATEFF_D) * (p_B_D_iht + reg_up_signal * r_RU_D_iht - reg_down_signal * r_RD_D_iht))
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end

# (b) Limits on storage device state of charge, considering the provision of ancillary services:

# Constraint (7)
function HydroBoost_constraint_ES_SOC_Bounds_DN_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    BATEFF_D = parameter(am, bus_idx, :gen_bus, tech_idx, "Discharging Efficiency") * 0.01 # convert to percent
    Min_SOC_MWh = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_SOC_MWh")

    # variable
    e_B_iht = get_variable(am, :e_B_iht, (i,h,t))
    r_RU_D_iht = get_variable(am, :r_RU_D_iht, (i,h,t))
    r_SR_D_iht = get_variable(am, :r_SR_D_iht, (i,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_B_iht - Min_SOC_MWh - (1/BATEFF_D) * (r_RU_D_iht + r_SR_D_iht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end

# Constraint (8)
function HydroBoost_constraint_ES_SOC_Bounds_UP_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    BATEFF_C = parameter(am, bus_idx, :gen_bus, tech_idx, "Charging Efficiency") * 0.01 # convert to percent
    Max_SOC_MWh = parameter(am, bus_idx, :gen_bus, tech_idx, "Max_SOC_MWh")
    
    # variable
    e_B_iht = get_variable(am, :e_B_iht, (i,h,t))
    r_RD_C_iht = get_variable(am, :r_RD_C_iht, (i,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_B_iht - Max_SOC_MWh + BATEFF_C * r_RD_C_iht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$h,$t)")) end    

end

# (c) Operational limits in discharging mode:

# Constraint (9)
function HydroBoost_constraint_ES_power_Bounds_UP_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    P_B_max = parameter(am, bus_idx, :gen_bus, tech_idx, "Max_Power")

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

# Constraint (10)
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

# Constraint (11)
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

# Constraint (12)
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

# Constraint (13)
function HydroBoost_constraint_ES_SR_D_cap_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    ramp_cap = parameter(am, bus_idx, :gen_bus, tech_idx, "Maximum SR")

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

# (d) Operational limits in charging mode:

# Constraint (14)
function HydroBoost_constraint_ES_charge_Bounds_UP_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    P_C_max = parameter(am, bus_idx, :gen_bus, tech_idx, "Max_Charge")

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

# Constraint (15)
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

# Constraint (16)
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

# Constraint (17)
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

# Constraint (18)
function HydroBoost_constraint_ES_SR_C_cap_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    ramp_cap = parameter(am, bus_idx, :gen_bus, tech_idx, "Maximum SR")

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

# (e) Discharging and charging power of BESS:

# Constraint (19)
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

# Constraint (20)
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

#(f) Intial and final storage of BESS
# Constraint (21)
# Constraint (22)
# Missing constraints present in the mathematical formulation. Are they already included in constrint (6)?

# (g) Allocation of ancillary services provided by the BESS:

# Constraint (23)
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

# Constraint (24)
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

# Constraint (25)
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

# Constraint (26), (27), and (28) are already included in the objective function formulation.

### Hydroelectric System without Reservoir ###

# h Commitment of hydro generators:

# Constraint (29)
function HydroBoost_constraint_hydro_commitment_status_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int, day_id, prior_day_solution; const_name_flag::Bool=false)

    if day_id == 1
        if h >= 2   # u_G_jht of h=1 for day_id 1 will be optimized
        
            # parameter  
            
            # variable
            u_G_jht = get_variable(am, :u_G_jht, (j,h,t))
            previous_u_G_jht = get_variable(am, :u_G_jht, (j,h-1,t))
            a_G_jht = get_variable(am, :a_G_jht, (j,h,t))
            z_G_jht = get_variable(am, :z_G_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                u_G_jht - previous_u_G_jht - a_G_jht + z_G_jht
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
            u_G_jht = get_variable(am, :u_G_jht, (j,h,t))
            previous_u_G_jht = prior_day_solution["hydro"][idx_string]["u_G_jht"]   # get prior day solution for hour 1
            a_G_jht = get_variable(am, :a_G_jht, (j,h,t))
            z_G_jht = get_variable(am, :z_G_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                u_G_jht - previous_u_G_jht - a_G_jht + z_G_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

        else

            # parameter  
            
            # variable
            u_G_jht = get_variable(am, :u_G_jht, (j,h,t))
            previous_u_G_jht = get_variable(am, :u_G_jht, (j,h-1,t))
            a_G_jht = get_variable(am, :a_G_jht, (j,h,t))
            z_G_jht = get_variable(am, :z_G_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                u_G_jht - previous_u_G_jht - a_G_jht + z_G_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

        end

    end

end

# Constraint (30)
function HydroBoost_constraint_hydro_start_up_shut_down_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    
    # variable
    a_G_jht = get_variable(am, :a_G_jht, (j,h,t))
    z_G_jht = get_variable(am, :z_G_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        a_G_jht + z_G_jht - 1
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraints (31) and (32) already included in section "Define decision variables".

# (i) Water discharge of plant:

# Constraints (33) and (35)
function HydroBoost_constraint_hydro_power_water_discharge_bounds_up_ljht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, l::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    if l == 1   # first water block

        # parameter  
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
        tag_id = string("Water Flow_", l)
        U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
        
        # variable
        q_Gl_jht = get_variable(am, :q_Gl_jht, (j,l,h,t))
        u_G_jht = get_variable(am, :u_G_jht, (j,h,t))
        
        # constraint
        expr = JuMP.@expression(JuMP_model,  
            q_Gl_jht - U_bar_l_j * u_G_jht
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
        q_Gl_jht = get_variable(am, :q_Gl_jht, (j,l,h,t))
        prior_w_Gl_jht = get_variable(am, :w_Gl_jht, (j,l-1,h,t))
        
        # constraint
        expr = JuMP.@expression(JuMP_model,  
            q_Gl_jht - U_bar_l_j * prior_w_Gl_jht
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$l,$h,$t)")) end    

    end
end

# Constraints (34) and (36)
function HydroBoost_constraint_hydro_power_water_discharge_bounds_down_ljht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, l::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
     # parameter  
     bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
     tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
     tag_id = string("Water Flow_", l)
     U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
     
     # variable
     q_Gl_jht = get_variable(am, :q_Gl_jht, (j,l,h,t))
     w_Gl_jht = get_variable(am, :w_Gl_jht, (j,l,h,t))
     
     # constraint
     expr = JuMP.@expression(JuMP_model,  
         q_Gl_jht - U_bar_l_j * w_Gl_jht
         )   
     constraint = JuMP.@constraint(JuMP_model, 
         expr >= 0
         )    
     if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$l,$h,$t)")) end    

end

# Constraint (37)
function HydroBoost_constraint_hydro_total_water_use_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_l, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    Min_Water_Release = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Water_Release")
    
    # variable
    q_G_jht = get_variable(am, :q_G_jht, (j,h,t))
    u_G_jht = get_variable(am, :u_G_jht, (j,h,t))
    
    sum_water_use = JuMP.AffExpr(0.0)
    for l in ids_l
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :q_Gl_jht, (j,l,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        q_G_jht - sum_water_use - Min_Water_Release * u_G_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# (h) Generation of hydro power:

# Constraint (38)
function HydroBoost_constraint_hydro_power_generation_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_l, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    P0 = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Power")
    
    # variable
    p_G_jht = get_variable(am, :p_G_jht, (j,h,t))
    u_G_jht = get_variable(am, :u_G_jht, (j,h,t))
    
    sum_power_output = JuMP.AffExpr(0.0)
    for l in ids_l
        tag_id = string("Water_Power_Conversion_", l)
        rho_l_j  = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)

        JuMP.add_to_expression!(sum_power_output, rho_l_j, get_variable(am, :q_Gl_jht, (j,l,h,t)))
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_G_jht - P0 * u_G_jht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (39)
function HydroBoost_constraint_total_hydro_power_generation_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_j_hydro, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    p_G_ht = get_variable(am, :p_G_ht, (h,t))
        
    sum_power_output = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :p_G_jht, (j,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_G_ht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

# (k) Water balance equation:

# Constraint (40)
function HydroBoost_constraint_hydro_water_balance_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_j_hydro, h::Int, t::Int, day_id, prior_day_solution; const_name_flag::Bool=false)
    
    # parameter
    Water_Inflow = am.ref[:nw][0][:repdays]["data"][h][t]["Inflow"]
    CF = 3600 * 0.0000229569
    
    # variable
    s_ht = get_variable(am, :s_ht, (h,t))
    e_H_ht = get_variable(am, :e_H_ht, (h,t))
    VMAX = am.ref[:nw][0][:hydro_reservoir][1]["VMAX"]
    
    global prior_e_H_ht = am.ref[:nw][0][:hydro_reservoir][1]["phi_ini"] * VMAX    # this value will be used for day id 1 at hour 1
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
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :q_G_jht, (j,h,t)))        
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

# (l) Spillage bounds:

# Constraint (41)
function HydroBoost_constraint_hydro_power_water_spillage_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_j_hydro, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    s_ht = get_variable(am, :s_ht, (h,t))

    sum_water_use = JuMP.AffExpr(0.0)
    for j in ids_j_hydro

        bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
        W_j = parameter(am, bus_idx, :gen_bus, tech_idx, "START_UP_Water_Use")
        tau_j = parameter(am, bus_idx, :gen_bus, tech_idx, "START_UP_Time")

        coeff = W_j * (tau_j / 60)

        JuMP.add_to_expression!(sum_water_use, coeff, get_variable(am, :a_G_jht, (j,h,t)))        
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

# (m) Rough zone constraints:

# Constraint (42)
function HydroBoost_constraint_hydro_rough_zone_y_l_minus_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    y_Gl_minus_jht = get_variable(am, :y_Gl_minus_jht, (j,h,t))
    q_Gl_jht = get_variable(am, :q_Gl_jht, (j,rough_zone_segment_number,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        q_Gl_jht - y_Gl_minus_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (43)
function HydroBoost_constraint_hydro_rough_zone_y_l_plus_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    tag_id = string("Water Flow_", rough_zone_segment_number)
    U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
    
    # variable
    y_Gl_plus_jht = get_variable(am, :y_Gl_plus_jht, (j,h,t))
    q_Gl_jht = get_variable(am, :q_Gl_jht, (j,rough_zone_segment_number,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        U_bar_l_j - q_Gl_jht - y_Gl_plus_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (44)
function HydroBoost_constraint_hydro_rough_zone_y_l_minus_big_M_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    big_M = 100000  # Todo: adjust
    
    # variable
    phi_Gl_jht = get_variable(am, :phi_Gl_jht, (j,h,t))
    y_Gl_minus_jht = get_variable(am, :y_Gl_minus_jht, (j,h,t))
    q_Gl_jht = get_variable(am, :q_Gl_jht, (j,rough_zone_segment_number,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        (q_Gl_jht - y_Gl_minus_jht) - big_M * phi_Gl_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (45)
function HydroBoost_constraint_hydro_rough_zone_y_l_plus_big_M_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    big_M = 100000  # Todo: adjust

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    tag_id = string("Water Flow_", rough_zone_segment_number)
    U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
    
    # variable
    phi_Gl_jht = get_variable(am, :phi_Gl_jht, (j,h,t))
    y_Gl_plus_jht = get_variable(am, :y_Gl_plus_jht, (j,h,t))
    q_Gl_jht = get_variable(am, :q_Gl_jht, (j,rough_zone_segment_number,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        (U_bar_l_j - q_Gl_jht - y_Gl_plus_jht) - big_M * (1-phi_Gl_jht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraints (46) and (47) already included in section "Define decision variables".

# Constraint (48-A)
function HydroBoost_constraint_hydro_rough_zone_y_l_plus_upper_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    tag_id = string("Water Flow_", rough_zone_segment_number)
    U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
    
    # variable
    y_Gl_plus_jht = get_variable(am, :y_Gl_plus_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        y_Gl_plus_jht - U_bar_l_j
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (48-B)
function HydroBoost_constraint_hydro_rough_zone_y_l_minus_upper_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, rough_zone_segment_number, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    tag_id = string("Water Flow_", rough_zone_segment_number)
    U_bar_l_j = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
        
    # variable
    y_Gl_minus_jht = get_variable(am, :y_Gl_minus_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        y_Gl_minus_jht - U_bar_l_j
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

### RES generators constraints ###

# (n) Non-dispatchable renewable generation constraints:

# Constraint (49)
function HydroBoost_constraint_RES_balance_kht(JuMP_model::JuMP.AbstractModel,am::Abstract_ALEAF_Model, const_name::String, k::String, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter
    P_Rfor = am.ref[:nw][0][:repdays]["data"][h][t][("P_Rfor", k)]

    # variable
    p_R   = get_variable(am, :p_R_kht,    (k,h,t))
    p_spill = get_variable(am, :p_Rspill_kht, (k,h,t))

    # constraint p_R + p_spill == P_Rfor
    expr = JuMP.@expression(JuMP_model, p_R + p_spill)
    c = @constraint(JuMP_model, expr == P_Rfor)
    
    if const_name_flag
        set_name(c, string(const_name, "_($k,$h,$t)"))
    end
end

# Constraint (50)
function HydroBoost_constraint_RES_total_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, k, h::Int, t::Int; const_name_flag::Bool=false,)

    # variable
    expr = JuMP.@expression(JuMP_model, sum(get_variable(am, :p_R_kht, (k_i,h,t)) for k_i in k))
    pTot = get_variable(am, :p_R_ht,(h,t))

    # constraint pTot == sum(p_R)
    c = JuMP.@constraint(JuMP_model, pTot == expr)

    if const_name_flag
        set_name(c, string(const_name, "_($h,$t)"))
    end
end

# Constraint (51) already included in section "Define decision variables".

### Coupling constraints ###

# (o) Hydro + BESS + RES power balances and interconnection limits:

# Constraint (52)
function HydroBoost_constraint_Hydro_power_generation_coupling_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    p_G_ht = get_variable(am, :p_G_ht, (h,t))
    p_G_Gr_ht = get_variable(am, :p_G_Gr_ht, (h,t))
    p_GB_ht = get_variable(am, :p_GB_ht, (h,t))
        
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_G_ht - (p_G_Gr_ht + p_GB_ht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

# Constraint (53)
function HydroBoost_constraint_RES_power_generation_coupling_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    p_R_ht = get_variable(am, :p_R_ht, (h,t))
    p_R_Gr_ht = get_variable(am, :p_R_Gr_ht, (h,t))
    p_R_St_ht = get_variable(am, :p_R_St_ht, (h,t))
        
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_R_ht - (p_R_Gr_ht + p_R_St_ht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

# Constraint (54)
function HydroBoost_constraint_ES_power_generation_coupling_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    p_B_CT_ht = get_variable(am, :p_B_CT_ht, (h,t))
    p_GB_ht = get_variable(am, :p_GB_ht, (h,t))
    p_R_St_ht = get_variable(am, :p_R_St_ht, (h,t))
    p_Gr_St_ht = get_variable(am, :p_Gr_St_ht, (h,t))
            
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_B_CT_ht - (p_GB_ht + p_R_St_ht + p_Gr_St_ht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

# Constraint (55)
function HydroBoost_constraint_outflow_interconnection_limit_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    outflow_limits = am.setting["Simulation Setting"]["Interconnection Limits Outflow"]

    # variable
    p_B_DT_ht = get_variable(am, :p_B_DT_ht, (h,t))
    p_G_Gr_ht = get_variable(am, :p_G_Gr_ht, (h,t))
    p_R_Gr_ht = get_variable(am, :p_R_Gr_ht, (h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        (p_B_DT_ht + p_G_Gr_ht + p_R_Gr_ht) - outflow_limits
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

# Constraint (56)
function HydroBoost_constraint_inflow_interconnection_limit_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    inflow_limits = am.setting["Simulation Setting"]["Interconnection Limits Inflow"]

    # variable
    p_Gr_St_ht = get_variable(am, :p_Gr_St_ht, (h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_Gr_St_ht - inflow_limits
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

##################################################
#------ Define variables bounds
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
    HydroBoost_report_result_RES_dispatch(ALEAF_model_instance, ALEAF_setting, daily_solutions, project_id, network_data)
    HydroBoost_report_result_objective(ALEAF_model_instance, ALEAF_setting, daily_solutions, project_id, network_data)

end

function HydroBoost_report_result_plant_dispatch(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data)
    
    start_time = time()

    # Create and open a file
    case_name = project_id
    output_path = network_data["output_path"]
    dispatch_file_name = string("ALEAF_HydroBoost_",case_name,"__plant_dispatch.csv")
    
    # Write Label 
    dispatch_label_list = ["day", "hour", "time", "p_B_DT_ht", "p_B_CT_ht", "p_Gr_St_ht", "p_G_ht", "p_GB_ht", "p_G_Gr_ht", "p_R_ht", "p_R_Gr_ht", "p_R_St_ht", "s_ht", "q_G_jht", "e_H_ht", "I_ht", "LMP", "Reg Up Price", "Reg Dn Price", "Spin Price"]

    # Write Outputs
    ids_j_hydro = [(j) for (j) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO"]
    ids_i_sto = [(i) for (i) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"]
    ids_k_ren = [(k) for k in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", k) == "RENEWABLE"]

    ids_l = [(l) for (l) in 1:am.setting["Simulation Setting"]["num_hydropower_performance_segment_value"]]
    ids_h = [(h) for (h) in 1:am.setting["Simulation Setting"]["num_hours_per_day_value"]]
    ids_t = [(t) for (t) in 1:am.setting["Simulation Setting"]["num_sub_period_value"]]
    
    
    num_days = 365
    # num_days = 4
    
    dispatch_output_list = Array{Any}(undef, num_days*length(ids_h)*length(ids_t)+1, length(dispatch_label_list))
    dispatch_output_list[1,:] = dispatch_label_list

    # for day_id in [1, 2, 3, 4]
    for day_id in eachindex([i for i in 1:num_days])    
        for hour_id in ids_h
            for time_id in ids_t
                
                idx_ht = string("(", hour_id, ", ", time_id, ")")

                # get solutions

                # storage
                p_B_DT_ht = daily_solutions[string(day_id)]["solution"]["storage"][idx_ht]["p_B_DT_ht"]
                p_B_CT_ht = daily_solutions[string(day_id)]["solution"]["storage"][idx_ht]["p_B_CT_ht"]
                p_Gr_St_ht = daily_solutions[string(day_id)]["solution"]["storage"][idx_ht]["p_Gr_St_ht"]

                # hydro
                p_G_ht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_G_ht"]
                p_GB_ht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_GB_ht"]
                p_G_Gr_ht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_G_Gr_ht"]

                # RES
                p_R_ht = daily_solutions[string(day_id)]["solution"]["renewable"][idx_ht]["p_R_ht"]
                p_R_Gr_ht = daily_solutions[string(day_id)]["solution"]["renewable"][idx_ht]["p_R_Gr_ht"]
                p_R_St_ht = daily_solutions[string(day_id)]["solution"]["renewable"][idx_ht]["p_R_St_ht"]

                # water flow
                s_ht = daily_solutions[string(day_id)]["solution"]["water"][idx_ht]["s_ht"]
                e_H_ht = daily_solutions[string(day_id)]["solution"]["water"][idx_ht]["e_H_ht"]

                q_G_jht = 0.0
                for j in ids_j_hydro
                    q_G_jht += daily_solutions[string(day_id)]["solution"]["water"][string("(", j, ", ", hour_id, ", ", time_id, ")")]["q_G_jht"]
                end

                I_ht = network_data["repdays"][day_id]["data"][hour_id][time_id]["Inflow"]

                # prices
                LMP = network_data["repdays"][day_id]["data"][hour_id][time_id]["DA_LMP"]
                Reg_up_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Regulation_up"]
                Reg_dn_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Regulation_down"]
                Spin_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Spin"]
                
                row_id = (day_id-1)*length(ids_h)*length(ids_t) + (hour_id-1)*length(ids_t) + (time_id) + 1
                
                dispatch_output_list[row_id, :] = [day_id, hour_id, time_id, 
                    p_B_DT_ht, p_B_CT_ht, p_Gr_St_ht, p_G_ht, p_GB_ht, p_G_Gr_ht, p_R_ht, p_R_Gr_ht, p_R_St_ht, s_ht, q_G_jht, e_H_ht, I_ht, LMP, Reg_up_price, Reg_dn_price, Spin_price]
            end
        end
    end
        
    CSV.write(string(output_path, dispatch_file_name), Tables.table(dispatch_output_list), writeheader=false)
        
    Memento.info(_LOGGER, "[ALEAF HydroBoost Model]:\t - Plant dispatch solution reporting, Total Time(sec): $(round(time() - start_time, digits=2))")

end

function HydroBoost_report_result_hydro_dispatch(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data)
    start_time = time()

    # Build output filename
    case_name          = project_id
    output_path        = network_data["output_path"]
    dispatch_file_name = string(output_path, "ALEAF_HydroBoost_", case_name, "__hydro_dispatch.csv")

    # Define CSV header
    dispatch_label_list = [
        "day", "hour", "time", "unit_id",
        "UnitGroup", "Unit_Category", "Unit_Type",
        "u_G_jht", "a_G_jht", "z_G_jht", "p_G_jht", "q_G_jht",
        "q_0_jht"
    ]
    num_segments = ALEAF_setting["Simulation Setting"]["num_hydropower_performance_segment_value"]
    for l in 1:num_segments
        push!(dispatch_label_list, "q_G$(l)_jht")
    end

    # Collect hydro unit indices
    ids_j_hydro = [j for j in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO"]

    # Map each hydro index to its Tech_ID
    id_to_name_hydro = Dict{Int,String}()
    for j in ids_j_hydro
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx",      j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
        id_to_name_hydro[j] = parameter(am, bus_idx, :gen_bus, tech_idx, "Tech_ID")
    end

    # Time indices
    ids_h    = 1:ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    ids_t    = 1:ALEAF_setting["Simulation Setting"]["num_sub_period_value"]
    num_days = 365

    # Preallocate output array (including header)
    dispatch_output_list = Array{Any}(undef, num_days * length(ids_j_hydro) * length(ids_h) * length(ids_t) + 1, length(dispatch_label_list))
    dispatch_output_list[1, :] = dispatch_label_list

    # Fill rows with solution data
    row = 2
    for day_id in 1:num_days, h in ids_h, t in ids_t, unit_idx in eachindex(ids_j_hydro)
        j   = ids_j_hydro[unit_idx]
        key = string("(", j, ", ", h, ", ", t, ")")

        # Textual unit ID
        unit_id = id_to_name_hydro[j]

        # Unit metadata
        bus_idx       = parameter(am, 0, :gen_index, "bus_idx",      j)
        tech_idx      = parameter(am, 0, :gen_index, "genco_tech_id", j)
        UNITGROUP     = parameter(am, bus_idx, :gen_bus, tech_idx, "UNITGROUP")
        Unit_Category = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_CATEGORY")
        Unit_Type     = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE")

        # Get decision and flow variables
        u_G_jht = daily_solutions[string(day_id)]["solution"]["hydro"][key]["u_G_jht"]
        a_G_jht = daily_solutions[string(day_id)]["solution"]["hydro"][key]["a_G_jht"]
        z_G_jht = daily_solutions[string(day_id)]["solution"]["hydro"][key]["z_G_jht"]
        p_G_jht = daily_solutions[string(day_id)]["solution"]["hydro"][key]["p_G_jht"]
        q_G_jht = daily_solutions[string(day_id)]["solution"]["water"][key]["q_G_jht"]

        # Minimum water release
        Min_Water_Release = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Water_Release")
        q_0_jht = u_G_jht * Min_Water_Release

        # Build row
        row_data = Any[day_id, h, t, unit_id, UNITGROUP, Unit_Category, Unit_Type, u_G_jht, a_G_jht, z_G_jht, p_G_jht, q_G_jht, q_0_jht]

        # Append each performance segment flow
        for l in 1:num_segments
            key_l = string("(", j, ", ", l, ", ", h, ", ", t, ")")
            push!(row_data,
                daily_solutions[string(day_id)]["solution"]["water"][key_l]["q_Gl_jht"]
            )
        end

        dispatch_output_list[row, :] = row_data
        row += 1
    end

    # Write CSV without extra header
    CSV.write(dispatch_file_name, Tables.table(dispatch_output_list); writeheader=false)
    Memento.info(_LOGGER, "[ALEAF HydroBoost Model]: - Hydro dispatch solution reporting, Total Time(sec): $(round(time() - start_time, digits=2))")
end

function HydroBoost_report_result_storage_dispatch(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data)
    start_time = time()

    # Build output filename
    case_name          = project_id
    output_path        = network_data["output_path"]
    dispatch_file_name = string(output_path, "ALEAF_HydroBoost_", case_name, "__storage_dispatch.csv")

    # Define CSV header
    dispatch_label_list = [
        "day", "hour", "time", "unit_id",
        "UnitGroup", "Unit_Category", "Unit_Type",
        "u_B_iht", "p_B_D_iht", "p_B_C_iht", "e_B_iht",
        "r_RU_D_iht", "r_RD_D_iht", "r_RU_C_iht", "r_RD_C_iht",
        "r_RU_iht", "r_RD_iht", "r_SR_D_iht", "r_SR_C_iht", "r_SR_iht"
    ]

    # Collect storage unit indices
    ids_i_sto = [
        i for i in get_index(am, :gen_index, 0)
        if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
    ]

    # Map each storage index to its Tech_ID from the Excel sheet
    id_to_name_sto = Dict{Int,String}()
    for i in ids_i_sto
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx",      i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        id_to_name_sto[i] = parameter(am, bus_idx, :gen_bus, tech_idx, "Tech_ID")
    end

    # Time indices
    ids_h    = 1:ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    ids_t    = 1:ALEAF_setting["Simulation Setting"]["num_sub_period_value"]
    num_days = 365

    # Preallocate output array (including header)
    dispatch_output_list = Array{Any}(
        undef,
        num_days * length(ids_i_sto) * length(ids_h) * length(ids_t) + 1,
        length(dispatch_label_list)
    )
    dispatch_output_list[1, :] = dispatch_label_list

    # Fill rows with solution data
    row = 2
    for day_id in 1:num_days, h in ids_h, t in ids_t, unit_idx in eachindex(ids_i_sto)
        i   = ids_i_sto[unit_idx]
        key = string("(", i, ", ", h, ", ", t, ")")

        # Textual unit ID
        unit_id = id_to_name_sto[i]

        # Unit metadata
        bus_idx       = parameter(am, 0, :gen_index, "bus_idx",      i)
        tech_idx      = parameter(am, 0, :gen_index, "genco_tech_id", i)
        UNITGROUP     = parameter(am, bus_idx, :gen_bus, tech_idx, "UNITGROUP")
        Unit_Category = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_CATEGORY")
        Unit_Type     = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE")

        # Get solution values
        sol     = daily_solutions[string(day_id)]["solution"]["storage"][key]
        u_B_iht = sol["u_B_iht"]
        p_B_D_iht = sol["p_B_D_iht"]
        p_B_C_iht = sol["p_B_C_iht"]
        e_B_iht = sol["e_B_iht"]
        r_RU_D_iht = sol["r_RU_D_iht"]
        r_RD_D_iht = sol["r_RD_D_iht"]
        r_RU_C_iht = sol["r_RU_C_iht"]
        r_RD_C_iht = sol["r_RD_C_iht"]
        r_RU_iht = sol["r_RU_iht"]
        r_RD_iht = sol["r_RD_iht"]
        r_SR_D_iht = sol["r_SR_D_iht"]
        r_SR_C_iht = sol["r_SR_C_iht"]
        r_SR_iht = sol["r_SR_iht"]

        # Populate the row
        dispatch_output_list[row, :] = [
            day_id, h, t, unit_id,
            UNITGROUP, Unit_Category, Unit_Type,
            u_B_iht, p_B_D_iht, p_B_C_iht, e_B_iht,
            r_RU_D_iht, r_RD_D_iht, r_RU_C_iht, r_RD_C_iht,
            r_RU_iht, r_RD_iht, r_SR_D_iht, r_SR_C_iht, r_SR_iht
        ]
        row += 1
    end

    # Write CSV without extra header
    CSV.write(dispatch_file_name, Tables.table(dispatch_output_list); writeheader=false)
    Memento.info(_LOGGER, "[ALEAF HydroBoost Model]: - Storage dispatch solution reporting, Total Time(sec): $(round(time() - start_time, digits=2))")
end

function HydroBoost_report_result_RES_dispatch(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data)
    start_time = time()
    case_name       = project_id
    output_path     = network_data["output_path"]
    dispatch_file_name = string("ALEAF_HydroBoost_", case_name, "__RES_dispatch.csv")
    dispatch_label_list = ["day", "hour", "time", "res_id", "UnitGroup", "Unit_Category", "Unit_Type", "p_R_kht", "p_Rspill_kht"]

    ids_k_ren = [ k for k in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", k) == "RENEWABLE" ]

    id_to_name = Dict{Int,String}()
    for k in ids_k_ren
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx",      k)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", k)
        id_to_name[k] = parameter(am, bus_idx, :gen_bus, tech_idx, "Tech_ID")
    end

    ids_h = 1:ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    ids_t = 1:ALEAF_setting["Simulation Setting"]["num_sub_period_value"]
    num_days = length(keys(network_data["repdays"]))

    total_rows = num_days * length(ids_h) * length(ids_t) * length(ids_k_ren) + 1
    dispatch_output_list = Array{Any}(undef, total_rows, length(dispatch_label_list))
    dispatch_output_list[1, :] = dispatch_label_list

    for (day_idx, day_id) in enumerate(sort(collect(keys(network_data["repdays"]))))
        for (h_idx, h) in enumerate(ids_h)
            for (t_idx, t) in enumerate(ids_t)
                for (unit_idx, k) in enumerate(ids_k_ren)
                    kname = id_to_name[k]
                    idx_kht = string("(\"", kname, "\", ", h, ", ", t, ")")

                    bus_idx  = parameter(am, 0, :gen_index,      "bus_idx",      k)
                    tech_idx = parameter(am, 0, :gen_index,      "genco_tech_id", k)
                    UNITGROUP     = parameter(am, bus_idx, :gen_bus, tech_idx, "UNITGROUP")
                    Unit_Category = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_CATEGORY")
                    Unit_Type     = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE")

                    sol_res      = daily_solutions[string(day_id)]["solution"]["renewable"][idx_kht]
                    p_R_kht      = sol_res["p_R_kht"]
                    p_Rspill_kht = sol_res["p_Rspill_kht"]

                    row_id = (day_idx-1)*length(ids_h)*length(ids_t)*length(ids_k_ren) +
                             (h_idx-1)*length(ids_t)*length(ids_k_ren) +
                             (t_idx-1)*length(ids_k_ren) +
                             unit_idx + 1

                    dispatch_output_list[row_id, :] = [day_id, h, t, kname, UNITGROUP, Unit_Category, Unit_Type, p_R_kht, p_Rspill_kht]
                end
            end
        end
    end

    CSV.write(string(output_path, dispatch_file_name), Tables.table(dispatch_output_list); writeheader=false)

    Memento.info(_LOGGER, "[ALEAF HydroBoost Model]: - RES dispatch solution reporting, Total Time(sec): $(round(time() - start_time, digits=2))")
end

function HydroBoost_report_result_objective(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data)
    start_time = time()

    # ---------- File path ----------
    case_name   = project_id
    output_path = network_data["output_path"]
    csv_file    = string(output_path, "ALEAF_HydroBoost_", case_name, "__objective_function.csv")

    # ---------- Indices ----------
    ids_h    = 1:ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    ids_t    = 1:ALEAF_setting["Simulation Setting"]["num_sub_period_value"]
    num_days = length(keys(network_data["repdays"]))

    # ---------- Unit sets ----------
    ids_i_sto = [i for i in get_index(am, :gen_index, 0)
                 if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"]

    ids_j_hydro = [j for j in get_index(am, :gen_index, 0)
                   if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO" &&
                   begin
                       bus_idx  = parameter(am, 0, :gen_index, "bus_idx",      j)
                       tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
                       parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE") == "GENERATOR"
                   end]

    # ---------- Rough-zone config ----------
    rz_flag    = get(am.setting["Simulation Setting"], "hydropower_rough_zone_flag", false)
    rz_penalty = get(am.setting["Simulation Setting"], "hyropower_rough_zone_operation_penalty_value", 0.0) # key name as in your code

    # ---------- Safe getter ----------
    get_or0(dict, key) = haskey(dict, key) ? dict[key] : 0.0

    # ---------- Header (add *_cum columns) ----------
    header = ["day","hour","time",
              "obj_ht",
              "energy_market","as_bess","as_hydro","startup_shutdown_cost","rough_zone_penalty","reg_adj_bess","reg_adj_hydro",
              "energy_market_cum","as_bess_cum","as_hydro_cum","startup_shutdown_cost_cum","rough_zone_penalty_cum",
              "reg_adj_bess_cum","reg_adj_hydro_cum","obj_cumulative"]

    n_rows = num_days * length(ids_h) * length(ids_t) + 1
    out = Array{Any}(undef, n_rows, length(header))
    out[1,:] = header
    row = 2

    # ---------- Cumulative trackers (run across the whole year) ----------
    cum_obj       = 0.0
    cum_energy    = 0.0
    cum_as_bess   = 0.0
    cum_as_hydro  = 0.0
    cum_susd      = 0.0
    cum_rz        = 0.0
    cum_adj_bess  = 0.0
    cum_adj_hydro = 0.0

    # ---------- Main loops ----------
    for day_id in 1:num_days
        day_key = string(day_id)
        for h in ids_h, t in ids_t
            idx_ht = string("(", h, ", ", t, ")")

            # Prices
            LMP   = network_data["repdays"][day_id]["data"][h][t]["DA_LMP"]
            Pru   = network_data["repdays"][day_id]["data"][h][t]["Regulation_up"]
            Prd   = network_data["repdays"][day_id]["data"][h][t]["Regulation_down"]
            Pspin = network_data["repdays"][day_id]["data"][h][t]["Spin"]

            # Energy market term
            p_G_Gr_ht  = daily_solutions[day_key]["solution"]["hydro"][idx_ht]["p_G_Gr_ht"]
            p_B_DT_ht  = daily_solutions[day_key]["solution"]["storage"][idx_ht]["p_B_DT_ht"]
            p_R_Gr_ht  = daily_solutions[day_key]["solution"]["renewable"][idx_ht]["p_R_Gr_ht"]
            p_Gr_St_ht = daily_solutions[day_key]["solution"]["storage"][idx_ht]["p_Gr_St_ht"]
            energy_market = LMP*(p_G_Gr_ht + p_B_DT_ht + p_R_Gr_ht - 1.00001*p_Gr_St_ht)

            # AS: BESS
            as_bess = 0.0
            reg_adj_bess = 0.0
            for i in ids_i_sto
                key_iht = string("(", i, ", ", h, ", ", t, ")")
                sol_s = daily_solutions[day_key]["solution"]["storage"][key_iht]
                rRU = sol_s["r_RU_iht"]; rRD = sol_s["r_RD_iht"]; rSR = sol_s["r_SR_iht"]
                as_bess += Pru*rRU + Prd*rRD + Pspin*rSR
                reg_up_signal   = 0.0
                reg_down_signal = 0.0
                reg_adj_bess += reg_up_signal*LMP*rRU - reg_down_signal*LMP*rRD
            end

            # AS: Hydro
            as_hydro = 0.0
            reg_adj_hydro = 0.0
            for j in ids_j_hydro
                key_jht = string("(", j, ", ", h, ", ", t, ")")
                sol_h = daily_solutions[day_key]["solution"]["hydro"][key_jht]
                rRUg = sol_h["r_RU_G_jht"]; rRDg = sol_h["r_RD_G_jht"]; rSRg = sol_h["r_SR_G_jht"]
                as_hydro += Pru*rRUg + Prd*rRDg + Pspin*rSRg
                reg_up_signal   = 0.0
                reg_down_signal = 0.0
                reg_adj_hydro += reg_up_signal*LMP*rRUg - reg_down_signal*LMP*rRDg
            end

            # Startup/Shutdown (Hydro)
            startup_shutdown_cost = 0.0
            for j in ids_j_hydro
                key_jht = string("(", j, ", ", h, ", ", t, ")")
                sol_h = daily_solutions[day_key]["solution"]["hydro"][key_jht]
                a_G = sol_h["a_G_jht"]; z_G = sol_h["z_G_jht"]
                bus_idx  = parameter(am, 0, :gen_index, "bus_idx",      j)
                tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
                c_su = parameter(am, bus_idx, :gen_bus, tech_idx, "START_UP_COST")
                c_sd = parameter(am, bus_idx, :gen_bus, tech_idx, "SHUT_DN_COST")
                startup_shutdown_cost += -(c_su*a_G + c_sd*z_G)
            end

            # Rough-zone (if enabled)
            rough_zone_penalty = 0.0
            if rz_flag && haskey(daily_solutions[day_key]["solution"], "rough_zone")
                for j in ids_j_hydro
                    key_jht = string("(", j, ", ", h, ", ", t, ")")
                    if haskey(daily_solutions[day_key]["solution"]["rough_zone"], key_jht)
                        sol_rz = daily_solutions[day_key]["solution"]["rough_zone"][key_jht]
                        y_plus  = get_or0(sol_rz, "y_Gl_plus_jht")
                        y_minus = get_or0(sol_rz, "y_Gl_minus_jht")
                        rough_zone_penalty += -rz_penalty*(y_plus + y_minus)
                    end
                end
            end

            # Objective at (h,t)
            obj_ht = energy_market + as_bess + as_hydro + startup_shutdown_cost +
                     rough_zone_penalty + reg_adj_bess + reg_adj_hydro

            # ---------- Update cumulatives ----------
            cum_obj       += obj_ht
            cum_energy    += energy_market
            cum_as_bess   += as_bess
            cum_as_hydro  += as_hydro
            cum_susd      += startup_shutdown_cost
            cum_rz        += rough_zone_penalty
            cum_adj_bess  += reg_adj_bess
            cum_adj_hydro += reg_adj_hydro

            # ---------- Write row ----------
            out[row, :] = Any[
                day_id, h, t,
                obj_ht,
                energy_market, as_bess, as_hydro, startup_shutdown_cost, rough_zone_penalty, reg_adj_bess, reg_adj_hydro,
                cum_energy, cum_as_bess, cum_as_hydro, cum_susd, cum_rz,
                cum_adj_bess, cum_adj_hydro, cum_obj
            ]
            row += 1
        end
    end

    CSV.write(csv_file, Tables.table(out); writeheader=false)

    Memento.info(_LOGGER,
        "[ALEAF HydroBoost Model]: - Objective function CSV (with cumulatives per component) written, " *
        "Total Time(sec): $(round(time() - start_time, digits=2))")
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
            CAP = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["Max_Power"]
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
    # 3) RES
    for plant_id in keys(am.ref[:nw][0][:gen_technology_RES])
        am.ref[:nw][sub_area_id][:gen_bus][plant_number] = deepcopy(am.ref[:nw][0][:gen_technology_RES][plant_id])
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

function allocate_res_sub_hourly_timeseries_data!(timeSeries, network_data, ALEAF_setting, data_label_list)
    # convert date-time string to DateTime
    timeSeries[!, "Date-time"] = DateTime.(timeSeries[!, "Date-Time"], "yyyy-mm-dd HH:MM:SS")

    # add month, day, hour info
    timeSeries[!, :Month] = Dates.month.(timeSeries[!, "Date-time"])
    timeSeries[!, :Day]   = Dates.day.(timeSeries[!, "Date-time"])
    timeSeries[!, :Hour]  = Dates.hour.(timeSeries[!, "Date-time"])

    num_hours_per_day = ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    num_sub_period    = ALEAF_setting["Simulation Setting"]["num_sub_period_value"]

    for day_id in keys(network_data["repdays"])
        
        global hour_idx = 1
        if !haskey(network_data["repdays"][day_id], "data") network_data["repdays"][day_id]["data"] = Dict{Any, Any}() end
        if !haskey(network_data["repdays"][day_id]["data"], hour_idx) network_data["repdays"][day_id]["data"][hour_idx] = Dict{Any, Any}() end

        for sim_day_id in network_data["repdays"][day_id]["look_ahead_simulation_day_idx_list"]
            
            month = network_data["repdays"][sim_day_id]["Month"]
            day   = network_data["repdays"][sim_day_id]["Day"]

            # select time-series data for this day
            sub_dataframe = timeSeries[(timeSeries.Month .== month).&(timeSeries.Day   .== day), :]

            # prepare daily data list
            daily_data_list = Dict{String, Vector{Float64}}()
            for data_label in data_label_list
                col = sub_dataframe[:, data_label]
                daily_data_list[data_label] = Float64.(coalesce.(col, 0.0))
            end

            for h in 1:num_hours_per_day
                for t in 1:num_sub_period
                    if !haskey(network_data["repdays"][day_id]["data"][h], t)
                        network_data["repdays"][day_id]["data"][h][t] = Dict{Any, Any}()
                    end
                    for data_label in data_label_list
                        network_data["repdays"][day_id]["data"][hour_idx][t][("P_Rfor", data_label)] =daily_data_list[data_label][h]
                    end
                    global hour_idx += 1
                end
            end
            daily_data_list = 0

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

    tab_names_with_category_list = ["Gen Technology - BESS", "Gen Technology - Hydro", "Gen Technology - RES", "Hydro Reservoir"]
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
    network_data["gen_technology_RES"] = deepcopy(ALEAF_setting["Gen Technology - RES"])
    network_data["hydro_reservoir"] = deepcopy(ALEAF_setting["Hydro Reservoir"])
    
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
    
    # collect RES profiles time-series data ["RES - Hourly"]
    res_data_path = string(data_location, "RES/")

    # --- "RES - Hourly"
    timeSeries_df = CSV.read(string(res_data_path, "RES_profiles.csv"), DataFrame)
    allocate_res_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, ["RES_1_PV", "RES_2_PV", "RES_3_WIND", "RES_4_WIND"])

    #------------------------------------------

    # Define output path
    output_path = string(pwd(), "/Simulation_Results/")
    check_and_create_path(output_path)

    output_path = string(output_path, project_id, "/")
    check_and_create_path(output_path)

    network_data["output_path"] = output_path
    
    return network_data

end