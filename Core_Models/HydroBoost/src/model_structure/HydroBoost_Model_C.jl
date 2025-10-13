#=

HydroBoost Model C: Pumped Storage Hydro (PSH) Salina Plant

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

    Memento.info(_LOGGER, "-- Current case: $project_id")

    case_set = sort(collect(keys(filter(x->(x.second["Run_Flag"] == true),
                                        ALEAF_setting["Simulation Configuration"]))))

    for case_id in case_set
        # REMOVE this line (it breaks the String keys "1","3"):
        # case_id = parse(Int64, case_id)

        Memento.info(_LOGGER, "-- Current case ID: $case_id")

        # update ALEAF_setting dictionary 
        ALEAF_setting_of_case_id = update_simulation_params(ALEAF_setting, case_id)

        # NEW: derive the human label for folder name (e.g., "Test 1")
        case_label = ALEAF_setting["Simulation Configuration"][case_id]["Case_ID"]

        # Generate Network Data  --> pass the 3rd positional arg (case_label)
        start_time = time()
        network_data = generate_networkdata_HydroBoost(ALEAF_setting_of_case_id, project_id, case_label)
        Memento.info(_LOGGER, "[HydroBoost Model]:\tGenerate network data. Time(sec): $(round(time() - start_time, digits=2))")

        solutions = build_and_run_daily_HydroBoost(ALEAF_setting_of_case_id, network_data)
        export_HydroBoost_results(ALEAF_setting_of_case_id, network_data, solutions, project_id)
    end
end

function update_simulation_params(ALEAF_setting, case_id)
    # copy ALEAF_setting dictionary first
    ALEAF_setting_of_case_id = deepcopy(ALEAF_setting)

    # keys in "Simulation Configuration" are "1","2","3","4" -> normalize to String
    key = string(case_id)

    # update parameters for the given case_id
    ALEAF_setting_of_case_id["Simulation Setting"]["look_ahead_days_value"]         = ALEAF_setting["Simulation Configuration"][key]["Look_Ahead_Days"]
    ALEAF_setting_of_case_id["Simulation Setting"]["Interconnection Limits Inflow"] = ALEAF_setting["Simulation Configuration"][key]["Interconnection Limits Inflow"]
    ALEAF_setting_of_case_id["Simulation Setting"]["Interconnection Limits Outflow"]= ALEAF_setting["Simulation Configuration"][key]["Interconnection Limits Outflow"]
    ALEAF_setting_of_case_id["Simulation Setting"]["Market_price_forecasting_method"]= ALEAF_setting["Simulation Configuration"][key]["Market_Price_File_ID"]

    for (uid, unit_any) in ALEAF_setting_of_case_id["Gen Technology - BESS"]
        unit = unit_any
        unit["Max_Power"]             = ALEAF_setting["Simulation Configuration"][key]["Max_Power"]
        unit["Max_SOC_MWh"]           = ALEAF_setting["Simulation Configuration"][key]["Max_SOC_MWh"]
        unit["Min_SOC_MWh"]           = ALEAF_setting["Simulation Configuration"][key]["Min_SOC_MWh"]
        unit["Max_Charge"]            = ALEAF_setting["Simulation Configuration"][key]["Max_Charge"]  
        unit["Roundtrip Efficiency"]  = ALEAF_setting["Simulation Configuration"][key]["Roundtrip Efficiency"]
        unit["Charging Efficiency"]   = ALEAF_setting["Simulation Configuration"][key]["Charging Efficiency"]
        unit["Discharging Efficiency"]= ALEAF_setting["Simulation Configuration"][key]["Discharging Efficiency"]
    end

    return ALEAF_setting_of_case_id
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
    ids_j_hydro = [j for j in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO" &&
    begin
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
        parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE") == "GENERATOR"
    end]

    ids_j_pump = [j for j in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO" &&
    begin
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
        parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE") == "PUMP"
    end]                                                                                                                                  

    ids_i_sto = [(i) for (i) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) === "STORAGE"]

    # Mapping unit -> UNITGROUP (read from :gen_bus)
    group_of = Dict{Int,String}()
    for j in vcat(ids_j_hydro, ids_j_pump)
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
        group_of[j] = parameter(am, bus_idx, :gen_bus, tech_idx, "UNITGROUP")
    end

    ids_g = unique(values(group_of))
    gens_of_g  = Dict{String,Vector{Int}}( g => [j for j in ids_j_hydro if group_of[j] == g] for g in ids_g )
    pumps_of_g = Dict{String,Vector{Int}}( g => [j for j in ids_j_pump  if group_of[j] == g] for g in ids_g )

    ids_k_ren_idx = [k for k in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", k) == "RENEWABLE"]
    id_to_name = Dict{Int,String}()
    for k in ids_k_ren_idx
        bus_idx  = parameter(am, 0, :gen_index, k, "bus_idx")
        tech_idx = parameter(am, 0, :gen_index, k, "genco_tech_id")
        id_to_name[k] = parameter(am, bus_idx, :gen_bus, tech_idx, "Tech_ID")
    end
    ids_k_ren = [id_to_name[k] for k in ids_k_ren_idx ]                                                                                                                                  # RES generator index

    ids_l = [(l) for (l) in 1:am.setting["Simulation Setting"]["num_hydropower_performance_segment_value"]]
    ids_h = [(h) for (h) in am.setting["run_H"]]
    ids_t = [(t) for (t) in am.setting["run_T"]]

    ###################################
    #------ Define decision variables (TOT = 57 variables)
    ###################################
    
    # Storage plant (BESS): dispatch variables
    HydroBoost_variable_iht_binary(JuMP_model, am, :storage, "u_B_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0,  bounded_upper=true, upper_bound=1.0)                # Binary variable driven to 1 when BESS is set to charging mode, and 0 otherwise      
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "e_B_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                      # State of charge of BESS unit i in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :storage, "e_B_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                   # State of charge of all BESS units in hour t [MWh]

    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "p_B_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                    # Power discharged from BESS and accounted for at point of delivery in period t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "p_B_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                    # Power contributing to charge BESS in period t before accounting for losses [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :storage, "p_B_DT_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                # Power discharged from BESS and accounted for at point of delivery in period t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :storage, "p_B_CT_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                # Power contributing to charge BESS in period t before accounting for losses [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :storage, "p_Gr_St_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                               # Power from grid allocated to charge storage (BESS) (before considering any round-trip losses) [MWh]
    
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Reserve for regulation up in discharging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Reserve for regulation up in charging mode provided by BESS unit i in hour t  [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                     # Reserve for regulation up provided by BESS unit i in hour t [MWh]
    
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Reserve for regulation down in discharging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Reserve for regulation down in charging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                     # Reserve for regulation down provided by BESS unit i in hour t [MWh]
    
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Reserve for spinning reserve in discharging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Reserve for spinning reserve in charging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                     # Reserve for spinning reserve provided by BESS unit i in hour t [MWh]
    
    # r_RU_ht: Reserve for regulation up provided by all BESS units in hour t [MWh] is already included in the objective function
    # r_RD_ht: Reserve for regulation down provided by all BESS units in hour t [MWh] is already included in the objective function
    # r_SR_ht: Reserve for spinning reserve provided by all BESS units in hour t [MWh] is already included in the objective function

    # Hydro plant: power dispatch variables    
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "u_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)              # Binary variable which is equal to 1 if hydro generator j is on-line in hour t, and 0 otherwise
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "a_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)              # Binary variable which is equal to 1 if hydro generator j is started at beginning of hour t, and 0 otherwise
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "z_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)              # Binary variable which is equal to 1 if hydro generator j is shut down at beginning of hour t, and 0 otherwise

    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "u_P_jht", ids_j_pump, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)               # Binary variable which is equal to 1 if hydro pump j is on-line in hour t, and 0 otherwise
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "a_P_jht", ids_j_pump, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)               # Binary variable which is equal to 1 if hydro pump j is started at beginning of hour t, and 0 otherwise
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "z_P_jht", ids_j_pump, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)               # Binary variable which is equal to 1 if hydro pump j is shut down at beginning of hour t, and 0 otherwise

    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "y_PUMP_ght", ids_g, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)                 # Binary variable which is equal to 1 if hydro group g is in pump mode at (h,t); 0 = generator mode

    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "p_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                      # Power output of hydro generator j in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_G_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                     # Day-ahead total power setpoint of dispatchable hydro plant in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_GB_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                    # Hydro power allocated to charge the BESS in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_G_Gr_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                  # Dispatchable hydro power sold to the grid in the day-ahead market in hour t [MWh] 

    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "p_P_jht", ids_j_pump, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                       # Power consumption by hydro pump j in pumping mode in hour t
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_P_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                     # Power consumption of all hydro pumps in hour t [MWh]
    
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_RU_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Hydro reserve for regulation up provided by hydro generator j in hour t, when in generation mode [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_RD_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Hydro reserve for regulation down provided by hydro generator j in hour t, when in generation mode [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_SR_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Hydro reserve for spinning reserve provided by hydro generator j in hour t, when in generation mode [MWh]

    # r_RU_G_ht: Hydro reserve for regulation up provided by all hydro generators in hour t, when in generation mode [MWh] is already included in the objective function
    # r_RD_G_ht: Hydro reserve for regulation down provided by all hydro generators in hour t, when in generation mode [MWh] is already included in the objective function
    # r_SR_G_ht: Hydro reserve for spinning reserve provided by all hydro generators in hour t, when in generation mode [MWh] is already included in the objective function
    
    # r_RU_GB_ht: Hydro + BESS reserve for regulation up provided by all hydro generators in hour t [MWh] is already included in the objective function
    # r_RD_GB_ht: Hydro + BESS reserve for regulation down provided by all hydro generators in hour t [MWh] is already included in the objective function
    # r_SR_GB_ht: Hydro + BESS reserve for spinning reserve provided by all hydro generators in hour t [MWh] is already included in the objective function

    # Hydro plant: water flow variables
    HydroBoost_variable_ilht_real(JuMP_model, am, :water, "q_Gl_jht", ids_j_hydro, ids_l, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                             # Water discharge of block ℓ of the piecewise linear function of hydro generator j in hour t when in generation mode [ft^3/s]
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "q_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                      # Water discharge of hydro generator j in hour t [ft^3/s]
    HydroBoost_variable_ilht_binary(JuMP_model, am, :water, "w_Gl_jht", ids_j_hydro, ids_l, ids_h, ids_t; bounded_lower=true, lower_bound=0.0, bounded_upper=true, upper_bound=1.0)      # Binary variable which is equal to 1 if water discharged by hydro generator j has exceeded block ℓ in hour t when in generation mode
    HydroBoost_variable_ht_real(JuMP_model, am, :water, "e_H_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                     # Water volume of reservoir in period t [A-F]
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "q_P_jht", ids_j_pump, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                       # Water flow rate pumped through hydro pump j in pumping mode at hour t [ft^3/s]
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "q_P_TOT_jht", ids_j_pump, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Total water flow rate pumped through hydro pump j in pumping mode at hour t [ft^3/s]

    HydroBoost_variable_iht_real(JuMP_model, am, :water, "q_RU_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Flow equivalent of regulation up allocated to hydro generator j in hour t, when in generation mode [ft^3/s]
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "q_RD_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Flow equivalent of regulation down allocated to hydro generator j in hour t, when in generation mode [ft^3/s]
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "q_SR_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Flow equivalent of spinning reserve allocated to hydro generator j in hour t, when in generation mode [ft^3/s]
    
    # Hydro plant: rough zone variables
    if am.setting["Simulation Setting"]["hydropower_rough_zone_flag"] == true
        HydroBoost_variable_iht_real(JuMP_model, am, :rough_zone, "y_Gl_plus_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                       # Slack variables (upper bound) which take non-zero value if operation outside rough zone ℓ of unit j cannot be honored in hour t, when operating in generation mode [ft^3/s]
        HydroBoost_variable_iht_real(JuMP_model, am, :rough_zone, "y_Gl_minus_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                      # Slack variables (lower bound) which take non-zero value if operation outside rough zone ℓ of unit j cannot be honored in hour t, when operating in generation mode [ft^3/s]
        HydroBoost_variable_iht_binary(JuMP_model, am, :rough_zone, "phi_Gl_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0, bounded_upper=true, upper_bound=1.0)   # Auxiliary binary variable required for representation of rough zone constraints in generation mode
    end                                                            
    
    
    # Non-dispatchable generators (RES plant): variables
    HydroBoost_variable_iht_real(JuMP_model, am, :renewable, "p_R_kht", ids_k_ren, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                    # Power output of non-dispatchable generator 𝑘 in hour 𝑡 [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :renewable, "p_R_Gr_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                              # Power of non-dispatchable generators sent to grid in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :renewable, "p_Rspill_kht", ids_k_ren, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                               # Curtailed power of non-dispatchable generator 𝑘 in hour 𝑡 [MWh]
    HydroBoost_variable_ht_real(JuMP_model,am, :renewable, "p_R_St_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                               # Power of non−dispatchable generators allocated to storage (BESS) in hour 𝑡 [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :renewable, "p_R_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                 # Power of all non-dispatchable generators in hour t [MWh]
    
    ###################################
    #------ Call of constraints in the build-loop
    ###################################


    ### BESS constraints: from (6) to (29) ###
    for h in ids_h
        for t in ids_t

            # total discharge and charge 
            HydroBoost_constraint_total_ES_power_discharge_ht(JuMP_model, am, "HydroBoost_constraint_total_ES_power_discharge_ht", ids_i_sto, h, t; const_name_flag)
            HydroBoost_constraint_total_ES_power_charge_ht(JuMP_model, am, "HydroBoost_constraint_total_ES_power_charge_ht", ids_i_sto, h, t; const_name_flag)
            HydroBoost_constraint_total_ES_soc_ht(JuMP_model, am, "HydroBoost_constraint_total_ES_soc_ht", ids_i_sto, h, t; const_name_flag)

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

    ### Hydro power system constraints: from (30) to (76) ###
    for h in ids_h
        for t in ids_t

            # Generation of hydro power generators constraints
            HydroBoost_constraint_total_hydro_power_generation_jht(JuMP_model, am, "HydroBoost_constraint_total_hydro_power_generation_jht", ids_j_hydro, h, t; const_name_flag)

            for j in ids_j_hydro
                # Commitment of hydro power generators constraints
                HydroBoost_constraint_hydro_generators_commitment_status_jht(JuMP_model, am, "HydroBoost_constraint_hydro_generators_commitment_status_jht", j, h, t, day_id, prior_day_solution; const_name_flag)
                HydroBoost_constraint_hydro_generators_start_up_shut_down_bound_jht(JuMP_model, am, "HydroBoost_constraint_hydro_generators_start_up_shut_down_bound_jht", j, h, t; const_name_flag)

                # Operational limits of hydro generators considering generation and ancillary services:
                HydroBoost_constraint_hydro_power_Bounds_UP_jht(JuMP_model, am, "HydroBoost_constraint_hydro_power_Bounds_UP_jht", j, h, t; const_name_flag)
                HydroBoost_constraint_hydro_power_Bounds_DN_jht(JuMP_model, am, "HydroBoost_constraint_hydro_power_Bounds_DN_jht", j, h, t; const_name_flag)
                HydroBoost_constraint_hydro_power_RU_cap_iht(JuMP_model, am, "HydroBoost_constraint_hydro_power_RU_cap_iht", j, h, t; const_name_flag)
                HydroBoost_constraint_hydro_power_RD_cap_iht(JuMP_model, am, "HydroBoost_constraint_hydro_power_RD_cap_iht", j, h, t; const_name_flag)
                HydroBoost_constraint_hydro_power_SR_cap_iht(JuMP_model, am, "HydroBoost_constraint_hydro_power_SR_cap_iht", j, h, t; const_name_flag)

                # Generation of hydro power and allocation of ancillary services constraints
                HydroBoost_constraint_hydro_power_generation_jht(JuMP_model, am, "HydroBoost_constraint_hydro_power_generation_jht", ids_l, j, h, t; const_name_flag)
                HydroBoost_constraint_hydro_power_generation_reg_up_jht(JuMP_model, am, "HydroBoost_constraint_hydro_power_generation_reg_up_jht", j, h, t; const_name_flag)
                HydroBoost_constraint_hydro_power_generation_reg_down_jht(JuMP_model, am, "HydroBoost_constraint_hydro_power_generation_reg_down_jht", j, h, t; const_name_flag)
                HydroBoost_constraint_hydro_power_generation_spin_res_jht(JuMP_model, am, "HydroBoost_constraint_hydro_power_generation_spin_res_jht", j, h, t; const_name_flag)

            end

            # Consumption of hydro power pumps constraints
            HydroBoost_constraint_total_hydro_pump_power_consumption_jht(JuMP_model, am, "HydroBoost_constraint_total_hydro_pump_power_consumption_jht", ids_j_pump, h, t; const_name_flag)
            
            for j in ids_j_pump
                # Commitment of hydro power generators constraints
                HydroBoost_constraint_hydro_pumps_commitment_status_jht(JuMP_model, am, "HydroBoost_constraint_hydro_pumps_commitment_status_jht", j, h, t, day_id, prior_day_solution; const_name_flag)
                HydroBoost_constraint_hydro_pumps_start_up_shut_down_bound_jht(JuMP_model, am, "HydroBoost_constraint_hydro_pumps_start_up_shut_down_bound_jht", j, h, t; const_name_flag)
                # Consumption of hydro power, when in pumping mode:
                HydroBoost_constraint_hydro_pump_power_consumption_jht(JuMP_model, am, "HydroBoost_constraint_hydro_pump_power_consumption_jht", j, h, t; const_name_flag)
                HydroBoost_constraint_hydro_pump_power_consumption_upper_limit_jht(JuMP_model, am, "HydroBoost_constraint_hydro_pump_power_consumption_upper_limit_jht", j, h, t; const_name_flag)
                HydroBoost_constraint_total_hydro_pump_flow_jht(JuMP_model, am, "HydroBoost_constraint_total_hydro_pump_flow_jht", j, h, t; const_name_flag)
            end
            
            for g in ids_g
                # Mutual exclusivity constraint
                HydroBoost_constraint_hydro_gen_pump_mutex_ght(JuMP_model, am, "HydroBoost_constraint_hydro_gen_pump_mutex_ght",g, h, t, gens_of_g, pumps_of_g; const_name_flag)
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

            # Flow equivalent of ancillary services provided by hydro generators constraints
            for j in ids_j_hydro
                HydroBoost_constraint_equivalent_hydro_flow_RU_SR_jht(JuMP_model, am, "HydroBoost_constraint_equivalent_hydro_flow_RU_SR_jht", ids_l, j, h, t; const_name_flag)
                HydroBoost_constraint_equivalent_hydro_flow_RD_jht(JuMP_model, am, "HydroBoost_constraint_equivalent_hydro_flow_RD_jht", j, h, t; const_name_flag)
            end

        end
    end

    # Water discharge constraints
    for h in ids_h
        for t in ids_t

            # Volume limits
            HydroBoost_constraint_hydro_reservoir_volume_lower_bound_ht(JuMP_model, am, "HydroBoost_constraint_hydro_reservoir_volume_lower_bound_ht", ids_j_hydro, h, t; const_name_flag)
            HydroBoost_constraint_hydro_reservoir_volume_upper_bound_ht(JuMP_model, am, "HydroBoost_constraint_hydro_reservoir_volume_upper_bound_ht", ids_j_hydro, h, t; const_name_flag)

            # Water balance constraints
            HydroBoost_constraint_hydro_water_balance_ht(JuMP_model, am, "HydroBoost_constraint_hydro_water_balance_ht", ids_j_hydro, ids_j_pump, h, t, day_id, prior_day_solution; const_name_flag)

        end
    end

    # Initial and End-of-Period water volume constraints
    # HydroBoost_constraint_hydro_reservoir_uses_ini_limit(JuMP_model, am, "HydroBoost_constraint_hydro_reservoir_uses_ini_limit", ids_h, ids_t; const_name_flag) --> Already included in the water balance equation - constraint (64)
    HydroBoost_constraint_hydro_reservoir_uses_end_limit(JuMP_model, am, "HydroBoost_constraint_hydro_reservoir_uses_end_limit", ids_h, ids_t; const_name_flag)

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

    ### RES generators constraints: from (77) to (79) ###
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

    ### Coupling constraints: from (80) to (87) ###
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

    # AS revenue provided by BESS
    for h in ids_h
        for t in ids_t
            for i in ids_i_sto
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Regulation_up"], get_variable(am, :r_RU_iht, (i,h,t)))
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Regulation_down"], get_variable(am, :r_RD_iht, (i,h,t)))
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Spin"], get_variable(am, :r_SR_iht, (i,h,t)))
            end
        end
    end

    # AS revenue provided by hydro generators
    for h in ids_h
        for t in ids_t
            for j in ids_j_hydro
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Regulation_up"], get_variable(am, :r_RU_G_jht, (j,h,t)))
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Regulation_down"], get_variable(am, :r_RD_G_jht, (j,h,t)))
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Spin"], get_variable(am, :r_SR_G_jht, (j,h,t)))
            end
        end
    end

    # Adjusted revenue due to regulation deployments provided by BESS
    for h in ids_h
        for t in ids_t
            for i in ids_i_sto
                reg_up_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RU"]
                reg_down_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RD"]

                JuMP.add_to_expression!(objective, reg_up_signal * am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :r_RU_iht, (i,h,t)))
                JuMP.add_to_expression!(objective, - reg_down_signal * am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :r_RD_iht, (i,h,t)))                
            end
        end
    end

    # Adjusted revenue due to regulation deployments provided by hydro generators
    for h in ids_h
        for t in ids_t
            for j in ids_j_hydro
                reg_up_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RU"]
                reg_down_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RD"]

                JuMP.add_to_expression!(objective, reg_up_signal * am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :r_RU_G_jht, (j,h,t)))
                JuMP.add_to_expression!(objective, - reg_down_signal * am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :r_RD_G_jht, (j,h,t)))                
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

    reg_up_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RU"]
    reg_down_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RD"]

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

# Constraint (21)
function HydroBoost_constraint_total_ES_soc_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_i_sto, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter

    # variable
    e_B_ht = get_variable(am, :e_B_ht, (h,t))

    sum_soc = JuMP.AffExpr(0.0)
    for i in ids_i_sto
        JuMP.add_to_expression!(sum_soc, 1.0, get_variable(am, :e_B_iht, (i,h,t)))
    end

    # constraint
    expr = JuMP.@expression(JuMP_model,
        e_B_ht - sum_soc
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )

    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end 
end

# (f) Intial and final storage of BESS
# Constraint (22) and (23) are already included in constraint (6).

# (g) Allocation of ancillary services provided by the BESS:

# Constraint (24)
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

# Constraint (25)
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

# Constraint (26)
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

# Constraint (27), (28), and (29) are already included in the objective function formulation.

### Pumped Storage Hydro (PSH) System ###

# (h) Commitment of hydro generators and pumps:

# Constraint (30)
function HydroBoost_constraint_hydro_generators_commitment_status_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int, day_id, prior_day_solution; const_name_flag::Bool=false)

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

# Constraint (31)
function HydroBoost_constraint_hydro_generators_start_up_shut_down_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)

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

# Constraints (32) and (33) are already included in section "Define decision variables".

# Constraint (34)
function HydroBoost_constraint_hydro_pumps_commitment_status_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int, day_id, prior_day_solution; const_name_flag::Bool=false)

    if day_id == 1
        if h >= 2   # u_P_jht of h=1 for day_id 1 will be optimized
        
            # parameter  
            
            # variable
            u_P_jht = get_variable(am, :u_P_jht, (j,h,t))
            previous_u_P_jht = get_variable(am, :u_P_jht, (j,h-1,t))
            a_P_jht = get_variable(am, :a_P_jht, (j,h,t))
            z_P_jht = get_variable(am, :z_P_jht, (j,h,t))

            # constraint
            expr = JuMP.@expression(JuMP_model,  
                u_P_jht - previous_u_P_jht - a_P_jht + z_P_jht
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
            u_P_jht = get_variable(am, :u_P_jht, (j,h,t))
            previous_u_P_jht = prior_day_solution["hydro"][idx_string]["u_P_jht"]   # get prior day solution for hour 1
            a_P_jht = get_variable(am, :a_P_jht, (j,h,t))
            z_P_jht = get_variable(am, :z_P_jht, (j,h,t))

            # constraint
            expr = JuMP.@expression(JuMP_model,  
                u_P_jht - previous_u_P_jht - a_P_jht + z_P_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

        else

            # parameter  
            
            # variable
            u_P_jht = get_variable(am, :u_P_jht, (j,h,t))
            previous_u_P_jht = get_variable(am, :u_P_jht, (j,h-1,t))
            a_P_jht = get_variable(am, :a_P_jht, (j,h,t))
            z_P_jht = get_variable(am, :z_P_jht, (j,h,t))

            # constraint
            expr = JuMP.@expression(JuMP_model,  
                u_P_jht - previous_u_P_jht - a_P_jht + z_P_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

        end

    end

end

# Constraint (35)
function HydroBoost_constraint_hydro_pumps_start_up_shut_down_bound_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    
    # variable
    a_P_jht = get_variable(am, :a_P_jht, (j,h,t))
    z_P_jht = get_variable(am, :z_P_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        a_P_jht + z_P_jht - 1
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraints (36) and (37) already included in section "Define decision variables".

# Constraint (38): group-level operating mode exclusivity (either pumps OR generators)
function HydroBoost_constraint_hydro_gen_pump_mutex_ght(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, g::String, h::Int, t::Int, gens_of_g::Dict{String,Vector{Int}}, pumps_of_g::Dict{String,Vector{Int}}; const_name_flag::Bool=false)

    # Group-level operating-mode exclusivity:
    # y_PUMP_ght = 1 → pump mode allowed, generators off
    # y_PUMP_ght = 0 → generator mode allowed, pumps off

    # g: the group id (string) for a hydro plant  (UNITGROUP).
    # gens_of_g[g]: list of generator unit indices in group g
    # pumps_of_g[g]: list of pump unit indices in group g

        # Mode binary for this group and time (1 = pump mode, 0 = generator mode)
    yP = get_variable(am, :y_PUMP_ght, (g,h,t))

    # Collect unit lists (may be empty)
    gens  = get(gens_of_g,  g, Int[])
    pumps = get(pumps_of_g, g, Int[])

    # If the group has no units at all, skip
    if isempty(gens) && isempty(pumps)
        return
    end

    # Build sum(u_G) as an affine expression (robust even if gens is empty)
    sum_uG = JuMP.AffExpr(0.0)
    for j in gens
        JuMP.add_to_expression!(sum_uG, 1.0, get_variable(am, :u_G_jht, (j,h,t)))
    end

    # Build sum(u_P) as an affine expression (robust even if pumps is empty)
    sum_uP = JuMP.AffExpr(0.0)
    for j in pumps
        JuMP.add_to_expression!(sum_uP, 1.0, get_variable(am, :u_P_jht, (j,h,t)))
    end

    # Natural big-M values = number of units in the group
    MG = length(gens)
    MP = length(pumps)

    # If yP = 0 (generation mode) → generators allowed, pumps forbidden
    c1 = @constraint(JuMP_model, sum_uG <= MG * (1 - yP))
    # If yP = 1 (pump mode)    → pumps allowed, generators forbidden
    c2 = @constraint(JuMP_model, sum_uP <= MP * yP)

    if const_name_flag
        JuMP.set_name(c1, string(const_name, "_G(", g, ",", h, ",", t, ")"))
        JuMP.set_name(c2, string(const_name, "_P(", g, ",", h, ",", t, ")"))
    end
end

# (i) Operational limits of hydro generators considering generation and ancillary services:

# Constraint (39)
function HydroBoost_constraint_hydro_power_Bounds_UP_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    P_G_max = parameter(am, bus_idx, :gen_bus, tech_idx, "Max_Power")

    # variable
    p_G_jht = get_variable(am, :p_G_jht, (j,h,t))
    r_RU_G_jht = get_variable(am, :r_RU_G_jht, (j,h,t))
    r_SR_G_jht = get_variable(am, :r_SR_G_jht, (j,h,t))
    
    u_G_jht = get_variable(am, :u_G_jht, (j,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_G_jht + r_RU_G_jht + r_SR_G_jht - P_G_max * u_G_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (40)
function HydroBoost_constraint_hydro_power_Bounds_DN_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    P_G_min = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Power")

    # variable
    p_G_jht = get_variable(am, :p_G_jht, (j,h,t))
    r_RD_G_jht = get_variable(am, :r_RD_G_jht, (j,h,t))

    u_G_jht = get_variable(am, :u_G_jht, (j,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_G_jht - r_RD_G_jht - P_G_min * u_G_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (41)
function HydroBoost_constraint_hydro_power_RU_cap_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    ramp_cap = parameter(am, bus_idx, :gen_bus, tech_idx, "Maximum RU")

    # variable
    r_RU_G_jht = get_variable(am, :r_RU_G_jht, (j,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RU_G_jht - ramp_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (42)
function HydroBoost_constraint_hydro_power_RD_cap_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    ramp_cap = parameter(am, bus_idx, :gen_bus, tech_idx, "Maximum RD")

    # variable
    r_RD_G_jht = get_variable(am, :r_RD_G_jht, (j,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RD_G_jht - ramp_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (43)
function HydroBoost_constraint_hydro_power_SR_cap_iht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    ramp_cap = parameter(am, bus_idx, :gen_bus, tech_idx, "Maximum SR")

    # variable
    r_SR_G_jht = get_variable(am, :r_SR_G_jht, (j,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_SR_G_jht - ramp_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# (j) Ramping constraints when in generation mode:

# Constraint (44)
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
            p_G_jht = get_variable(am, :p_G_jht, (j,h,t))
            prior_p_G_jht = get_variable(am, :p_G_jht, (j,h-1,t))
            prior_u_G_jht = get_variable(am, :u_G_jht, (j,h-1,t))
            a_G_jht = get_variable(am, :a_G_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                (p_G_jht- prior_p_G_jht) - MAX_RU * prior_u_G_jht - MAX_RU_START_UP * a_G_jht
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
            p_G_jht = get_variable(am, :p_G_jht, (j,h,t))
            prior_p_G_jht = prior_day_solution["hydro"][idx_string]["p_G_jht"]   # get prior day solution for hour 1
            prior_u_G_jht = prior_day_solution["hydro"][idx_string]["u_G_jht"]   # get prior day solution for hour 1
            a_G_jht = get_variable(am, :a_G_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                (p_G_jht- prior_p_G_jht) - MAX_RU * prior_u_G_jht - MAX_RU_START_UP * a_G_jht
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
            p_G_jht = get_variable(am, :p_G_jht, (j,h,t))
            prior_p_G_jht = get_variable(am, :p_G_jht, (j,h-1,t))
            prior_u_G_jht = get_variable(am, :u_G_jht, (j,h-1,t))
            a_G_jht = get_variable(am, :a_G_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                (p_G_jht- prior_p_G_jht) - MAX_RU * prior_u_G_jht - MAX_RU_START_UP * a_G_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    
        end
    end

end

# Constraint (45)
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
            p_G_jht = get_variable(am, :p_G_jht, (j,h,t))
            prior_p_G_jht = get_variable(am, :p_G_jht, (j,h-1,t))
            u_G_jht = get_variable(am, :u_G_jht, (j,h,t))
            z_G_jht = get_variable(am, :z_G_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                - (p_G_jht- prior_p_G_jht) - MAX_RD * u_G_jht - MAX_RD_SHUT_DN * z_G_jht
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
            p_G_jht = get_variable(am, :p_G_jht, (j,h,t))
            prior_p_G_jht = prior_day_solution["hydro"][idx_string]["p_G_jht"]   # get prior day solution for hour 1
            u_G_jht = get_variable(am, :u_G_jht, (j,h,t))
            z_G_jht = get_variable(am, :z_G_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                - (p_G_jht- prior_p_G_jht) - MAX_RD * u_G_jht - MAX_RD_SHUT_DN * z_G_jht
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
            p_G_jht = get_variable(am, :p_G_jht, (j,h,t))
            prior_p_G_jht = get_variable(am, :p_G_jht, (j,h-1,t))
            u_G_jht = get_variable(am, :u_G_jht, (j,h,t))
            z_G_jht = get_variable(am, :z_G_jht, (j,h,t))
            
            # constraint
            expr = JuMP.@expression(JuMP_model,  
                - (p_G_jht- prior_p_G_jht) - MAX_RD * u_G_jht - MAX_RD_SHUT_DN * z_G_jht
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    
            
    
        end
    end


end

# (k) Water discharge of plant:

# Constraints (46) and (48)
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

# Constraints (47) and (49)
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

# Constraint (50)
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

# (l) Flow equivalent of ancillary services provided by hydro generators:

# Constraint (51)
function HydroBoost_constraint_equivalent_hydro_flow_RU_SR_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_l, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    Min_Water_Release = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Water_Release")
    sum_U_bar = 0.0
    for l in ids_l
        tag_id = string("Water Flow_", l)
        U_bar = parameter(am, bus_idx, :gen_bus, tech_idx, tag_id)
        sum_U_bar += U_bar
    end
    
    # variable
    q_G_jht = get_variable(am, :q_G_jht, (j,h,t))
    q_RU_G_jht = get_variable(am, :q_RU_G_jht, (j,h,t))
    q_SR_G_jht = get_variable(am, :q_SR_G_jht, (j,h,t))
    u_G_jht = get_variable(am, :u_G_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        q_G_jht + q_RU_G_jht + q_SR_G_jht - (Min_Water_Release + sum_U_bar) * u_G_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (52)
function HydroBoost_constraint_equivalent_hydro_flow_RD_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    Min_Water_Release = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Water_Release")

    # variable
    q_G_jht = get_variable(am, :q_G_jht, (j,h,t))
    q_RD_G_jht = get_variable(am, :q_RD_G_jht, (j,h,t))
    u_G_jht = get_variable(am, :u_G_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        q_G_jht - q_RD_G_jht - Min_Water_Release * u_G_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# (m) Generation of hydro power and allocation of ancillary services, when in generation mode::

# Constraint (53)
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

# Constraint (54)
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

# Constraint (55)
function HydroBoost_constraint_hydro_power_generation_reg_up_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    rho_ave_j = parameter(am, bus_idx, :gen_bus, tech_idx, "Average_Water_Power_Conversion_Coeff")
    
    # variable
    r_RU_G_jht = get_variable(am, :r_RU_G_jht, (j,h,t))
    q_RU_G_jht = get_variable(am, :q_RU_G_jht, (j,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RU_G_jht - rho_ave_j * q_RU_G_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (56)
function HydroBoost_constraint_hydro_power_generation_reg_down_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    rho_ave_j = parameter(am, bus_idx, :gen_bus, tech_idx, "Average_Water_Power_Conversion_Coeff")
    
    # variable
    r_RD_G_jht = get_variable(am, :r_RD_G_jht, (j,h,t))
    q_RD_G_jht = get_variable(am, :q_RD_G_jht, (j,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RD_G_jht - rho_ave_j * q_RD_G_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (57)
function HydroBoost_constraint_hydro_power_generation_spin_res_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    rho_ave_j = parameter(am, bus_idx, :gen_bus, tech_idx, "Average_Water_Power_Conversion_Coeff")
    
    # variable
    r_SR_G_jht = get_variable(am, :r_SR_G_jht, (j,h,t))
    q_SR_G_jht = get_variable(am, :q_SR_G_jht, (j,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_SR_G_jht - rho_ave_j * q_SR_G_jht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (58), (59), and (60) are already included in the objective function formulation.

# (n) Consumption of hydro power, when in pumping mode:

# Constraint (61)
function HydroBoost_constraint_hydro_pump_power_consumption_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
    P0_P = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Power")
    rho_ave_j = parameter(am, bus_idx, :gen_bus, tech_idx, "Average_Water_Power_Conversion_Coeff")
    
    # variable
    p_P_jht = get_variable(am, :p_P_jht, (j,h,t))
    q_P_jht = get_variable(am, :q_P_jht, (j,h,t))
    u_P_jht = get_variable(am, :u_P_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_P_jht - rho_ave_j * q_P_jht - u_P_jht * P0_P
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    

end

# Constraint (62)
function HydroBoost_constraint_total_hydro_pump_power_consumption_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_j_pump, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter

    # variable
    p_P_ht = get_variable(am, :p_P_ht, (h,t))

    sum_power_output = JuMP.AffExpr(0.0)
    for j in ids_j_pump
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :p_P_jht, (j,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_P_ht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

# Constraint (63)
function HydroBoost_constraint_hydro_pump_power_consumption_upper_limit_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
    Q_P_MAX = parameter(am, bus_idx, :gen_bus, tech_idx, "Max_Water_Release")
    P0_P = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Power")
    rho_ave_j = parameter(am, bus_idx, :gen_bus, tech_idx, "Average_Water_Power_Conversion_Coeff")

    # variable
    q_P_jht = get_variable(am, :q_P_jht, (j,h,t))
    u_P_jht = get_variable(am, :u_P_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        q_P_jht + (1 / rho_ave_j) * P0_P * u_P_jht - u_P_jht * Q_P_MAX
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    
end

# Constraint (64)
function HydroBoost_constraint_total_hydro_pump_flow_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, h::Int, t::Int; const_name_flag::Bool=false)
    
    # set
    
    # parameter  
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
    P0_P = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Power")
    rho_ave_j = parameter(am, bus_idx, :gen_bus, tech_idx, "Average_Water_Power_Conversion_Coeff")

    # variable
    q_P_TOT_jht = get_variable(am, :q_P_TOT_jht, (j,h,t))
    q_P_jht = get_variable(am, :q_P_jht, (j,h,t))
    u_P_jht = get_variable(am, :u_P_jht, (j,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        q_P_TOT_jht - q_P_jht - (1 / rho_ave_j) * P0_P * u_P_jht
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$h,$t)")) end    
end


# (o) Water balance equation:

# Constraint (65)
function HydroBoost_constraint_hydro_water_balance_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_j_hydro, ids_j_pump, h::Int, t::Int, day_id, prior_day_solution; const_name_flag::Bool=false)

    # constants / parameters
    CF = 3600 * 0.0000229569 # Factor to convert water flow units in [ft^3/s] into net volume units in [A-F/h]         
    reg_up_signal   = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RU"]
    reg_down_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RD"]
    C_0 = am.ref[:nw][0][:repdays]["data"][h][t]["Diversion_C0"]
    C_1 = am.ref[:nw][0][:repdays]["data"][h][t]["Diversion_C1"]

    # variables
    e_H_ht = get_variable(am, :e_H_ht, (h,t))
    VMIN = am.ref[:nw][0][:hydro_reservoir][1]["VMIN"]
    VMAX = am.ref[:nw][0][:hydro_reservoir][1]["VMAX"]
    
    global prior_e_H_ht = (1 - am.ref[:nw][0][:hydro_reservoir][1]["phi_ini"]) * (VMAX - VMIN) + VMIN   # this value will be used for day id 1 at hour 1
    if h == 1
        if day_id > 1
            idx_string = string("(", 24,  ", ", t,")")
            global prior_e_H_ht = prior_day_solution["water"][idx_string]["e_H_ht"]   # get prior day solution for hour 1
        end
    else    # if h >= 2
        global prior_e_H_ht = get_variable(am, :e_H_ht, (h-1,t))

    end

    # sum of pump contributions (positive in balance)
    sum_pump = JuMP.AffExpr(0.0)
    for j in ids_j_pump
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)

        rho_ave_j = parameter(am, bus_idx, :gen_bus, tech_idx, "Average_Water_Power_Conversion_Coeff") 
        P0_P      = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Power")                           

        q_P_jht = get_variable(am, :q_P_jht, (j,h,t))
        u_P_jht = get_variable(am, :u_P_jht, (j,h,t))

        JuMP.add_to_expression!(sum_pump, 1.0, q_P_jht)
        JuMP.add_to_expression!(sum_pump, (1.0 / rho_ave_j) * P0_P, u_P_jht)
    end

    # sum of generator withdrawals (negative in balance)
    sum_gen = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        q_G_jht     = get_variable(am, :q_G_jht,     (j,h,t))
        q_RU_G_jht  = get_variable(am, :q_RU_G_jht,  (j,h,t))
        q_RD_G_jht  = get_variable(am, :q_RD_G_jht,  (j,h,t))

        JuMP.add_to_expression!(sum_gen, 1.0, q_G_jht)
        JuMP.add_to_expression!(sum_gen, reg_up_signal,   q_RU_G_jht)
        JuMP.add_to_expression!(sum_gen, - reg_down_signal, q_RD_G_jht)
    end

    # constraint
    expr = JuMP.@expression(JuMP_model,
        (e_H_ht - prior_e_H_ht) - CF * (sum_pump - sum_gen - C_0 - C_1 * e_H_ht)
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr == 0
        )
    if const_name_flag
        JuMP.set_name(constraint, string(const_name, "_(", h, ",", t, ")"))
    end
end

# (p) Volume limits:

# Constraint (66)
function HydroBoost_constraint_hydro_reservoir_volume_lower_bound_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_j_hydro, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter
    VMIN = am.ref[:nw][0][:hydro_reservoir][1]["VMIN"]
    CF = 3600 * 0.0000229569

    # variable
    e_H_ht = get_variable(am, :e_H_ht, (h,t))

    sum_water_use = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :q_RU_G_jht, (j,h,t)))
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :q_SR_G_jht, (j,h,t)))
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_H_ht - VMIN - CF * sum_water_use
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

# Constraint (67)
function HydroBoost_constraint_hydro_reservoir_volume_upper_bound_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_j_hydro, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter
    VMAX = am.ref[:nw][0][:hydro_reservoir][1]["VMAX"]
    CF = 3600 * 0.0000229569

    # variable
    e_H_ht = get_variable(am, :e_H_ht, (h,t))

    sum_water_use = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :q_RD_G_jht, (j,h,t)))
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_H_ht - VMAX + CF * sum_water_use
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

# (q) Initial and End-of-Period water volume constraints:

# Constraints (68) --> Already included in the water balance equation - constraint (65)
# function HydroBoost_constraint_hydro_reservoir_uses_ini_limit(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_h, ids_t; const_name_flag::Bool=false)
    
#     # Set/Index
#     first_hour = first(ids_h)
#     first_min = first(ids_t)

#     # parameter
#     VMIN = am.ref[:nw][0][:hydro_reservoir][1]["VMIN"]
#     VMAX = am.ref[:nw][0][:hydro_reservoir][1]["VMAX"]
#     phi_ini = am.ref[:nw][0][:hydro_reservoir][1]["phi_ini"]

#     # variable
#     e_H_ht = get_variable(am, :e_H_ht, (first_hour,first_min))

#     # constraint
#     expr = JuMP.@expression(JuMP_model,  
#         e_H_ht - VMIN - (VMAX - VMIN) * (1 - phi_ini)
#         )   
#     constraint = JuMP.@constraint(JuMP_model, 
#         expr == 0
#         )    
#     if const_name_flag JuMP.set_name(constraint, string(const_name, "_($first_hour,$first_min)")) end    

# end

# Constraints (69)
function HydroBoost_constraint_hydro_reservoir_uses_end_limit(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_h, ids_t; const_name_flag::Bool=false)

    # Set/Index
    last_hour = last(ids_h)
    last_min = last(ids_t)

    # parameter
    VMIN = am.ref[:nw][0][:hydro_reservoir][1]["VMIN"]
    VMAX = am.ref[:nw][0][:hydro_reservoir][1]["VMAX"]
    phi_end = am.ref[:nw][0][:hydro_reservoir][1]["phi_end"]

    # variable
    e_H_ht = get_variable(am, :e_H_ht, (last_hour,last_min))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_H_ht - VMIN - (VMAX - VMIN) * (1 - phi_end)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($last_hour,$last_min)")) end    

end

# (r) Rough zone constraints:

# Constraint (70)
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

# Constraint (71)
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

# Constraint (72)
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

# Constraint (73)
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

# Constraints (74) and (75) are already included in section "Define decision variables".

# Constraint (76-A)
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

# Constraint (76-B)
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

# (s) Non-dispatchable renewable generation constraints:

# Constraint (77)
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

# Constraint (78)
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

# Constraint (79) is already included in section "Define decision variables".

### Coupling constraints ###

# (t) Ancillary services offered to the market:

# Constraints (80), (81), and (82) are already included in the objective function formulation.

# (u) Hydro + BESS + RES power balances and interconnection limits:

# Constraint (83)
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

# Constraint (84)
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

# Constraint (85)
function HydroBoost_constraint_ES_power_generation_coupling_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    
    # variable
    p_B_CT_ht = get_variable(am, :p_B_CT_ht, (h,t))
    p_P_ht = get_variable(am, :p_P_ht, (h,t))
    p_GB_ht = get_variable(am, :p_GB_ht, (h,t))
    p_R_St_ht = get_variable(am, :p_R_St_ht, (h,t))
    p_Gr_St_ht = get_variable(am, :p_Gr_St_ht, (h,t))
            
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_B_CT_ht + p_P_ht - (p_GB_ht + p_R_St_ht + p_Gr_St_ht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

# Constraint (86)
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

# Constraint (87)
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
    dispatch_label_list = ["day", "hour", "time", "LMP", "Reg Up Price", "Reg Dn Price", "Spin Price", "Delta RU", "Delta RD", "p_G_ht", "p_P_ht", "p_GB_ht", "p_G_Gr_ht", "q_G_jht", "q_P_TOT_jht", "e_H_ht", "p_B_DT_ht", "p_B_CT_ht", "e_B_ht",  "p_R_ht", "p_R_Gr_ht", "p_R_St_ht", "p_Gr_St_ht"]

    # Write Outputs
    ids_j_hydro = [j for j in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO" &&
    begin
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
        parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE") == "GENERATOR"
    end]

    ids_j_pump = [j for j in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO" &&
    begin
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
        parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE") == "PUMP"
    end]

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
                e_B_ht = daily_solutions[string(day_id)]["solution"]["storage"][idx_ht]["e_B_ht"]
                p_Gr_St_ht = daily_solutions[string(day_id)]["solution"]["storage"][idx_ht]["p_Gr_St_ht"]

                # hydro
                p_G_ht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_G_ht"]
                p_P_ht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_P_ht"]
                p_GB_ht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_GB_ht"]
                p_G_Gr_ht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_G_Gr_ht"]

                # RES
                p_R_ht = daily_solutions[string(day_id)]["solution"]["renewable"][idx_ht]["p_R_ht"]
                p_R_Gr_ht = daily_solutions[string(day_id)]["solution"]["renewable"][idx_ht]["p_R_Gr_ht"]
                p_R_St_ht = daily_solutions[string(day_id)]["solution"]["renewable"][idx_ht]["p_R_St_ht"]

                # water flow
                e_H_ht = daily_solutions[string(day_id)]["solution"]["water"][idx_ht]["e_H_ht"]

                q_G_jht = 0.0
                for j in ids_j_hydro
                    q_G_jht += daily_solutions[string(day_id)]["solution"]["water"][string("(", j, ", ", hour_id, ", ", time_id, ")")]["q_G_jht"]
                end

                q_P_TOT_jht = 0.0
                for j in ids_j_pump
                    q_P_TOT_jht += daily_solutions[string(day_id)]["solution"]["water"][string("(", j, ", ", hour_id, ", ", time_id, ")")]["q_P_TOT_jht"]
                end

                # prices
                LMP = network_data["repdays"][day_id]["data"][hour_id][time_id]["DA_LMP"]
                Reg_up_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Regulation_up"]
                Reg_dn_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Regulation_down"]
                Spin_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Spin"]
                Delta_RU = network_data["repdays"][day_id]["data"][hour_id][time_id]["Delta_RU"]
                Delta_RD = network_data["repdays"][day_id]["data"][hour_id][time_id]["Delta_RD"]
                
                row_id = (day_id-1)*length(ids_h)*length(ids_t) + (hour_id-1)*length(ids_t) + (time_id) + 1

                dispatch_output_list[row_id, :] = [day_id, hour_id, time_id, LMP, Reg_up_price, Reg_dn_price, Spin_price, Delta_RU, Delta_RD,
                p_G_ht, p_P_ht, p_GB_ht, p_G_Gr_ht, q_G_jht, q_P_TOT_jht, e_H_ht, p_B_DT_ht, p_B_CT_ht, e_B_ht, p_R_ht, p_R_Gr_ht, p_R_St_ht, p_Gr_St_ht]
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
    dispatch_label_list = ["day", "hour", "time", "unit_id", "UnitGroup", "Unit_Category", "Unit_Type", "u_G_jht", "u_P_jht", "a_G_jht", "a_P_jht", "z_G_jht", "z_P_jht",
    "p_G_jht", "p_P_jht", "q_G_jht", "q_P_jht", "q_P_TOT_jht", "q_RU_G_jht", "q_RD_G_jht", "q_SR_G_jht", "r_RU_G_jht", "r_RD_G_jht", "r_SR_G_jht", "q_0_jht"]

    num_segments = ALEAF_setting["Simulation Setting"]["num_hydropower_performance_segment_value"]
    for l in 1:num_segments
        push!(dispatch_label_list, "q_G$(l)_jht")
    end

    # Collect hydro unit indices
    ids_j_hydro = [j for j in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO" &&
    begin
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
        parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE") == "GENERATOR"
    end]

    ids_j_pump = [j for j in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO" &&
    begin
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
        parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE") == "PUMP"
    end]

    ids_j_hydro_all = vcat(ids_j_hydro, ids_j_pump)

    # Map each hydro index to its Tech_ID
    id_to_name_hydro = Dict{Int,String}()
    for j in ids_j_hydro_all
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
        id_to_name_hydro[j] = parameter(am, bus_idx, :gen_bus, tech_idx, "Tech_ID")
    end

    # Time indices
    ids_h    = 1:ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    ids_t    = 1:ALEAF_setting["Simulation Setting"]["num_sub_period_value"]
    num_days = 365

    # Preallocate output array (including header)
    n_units = length(ids_j_hydro_all)
    n_rows  = num_days * n_units * length(ids_h) * length(ids_t) + 1
    n_cols  = length(dispatch_label_list)
    dispatch_output_list = Array{Any}(undef, n_rows, n_cols)
    dispatch_output_list[1, :] = dispatch_label_list

    # Fill rows with solution data
    row = 2

    # Hydro generators
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

        # Hydro and water solution slices
        sol_h      = daily_solutions[string(day_id)]["solution"]["hydro"][key]
        sol_w      = daily_solutions[string(day_id)]["solution"]["water"][key]
        u_G_jht    = sol_h["u_G_jht"]
        a_G_jht    = sol_h["a_G_jht"]
        z_G_jht    = sol_h["z_G_jht"]
        p_G_jht    = sol_h["p_G_jht"]
        q_G_jht    = sol_w["q_G_jht"]
        q_RU_G_jht = sol_w["q_RU_G_jht"]
        q_RD_G_jht = sol_w["q_RD_G_jht"]
        q_SR_G_jht = sol_w["q_SR_G_jht"]
        r_RU_G_jht = sol_h["r_RU_G_jht"]
        r_RD_G_jht = sol_h["r_RD_G_jht"]
        r_SR_G_jht = sol_h["r_SR_G_jht"]

        # pump fields not applicable → 0.0
        u_P_jht = 0.0; a_P_jht = 0.0; z_P_jht = 0.0; p_P_jht = 0.0; q_P_jht = 0.0; q_P_TOT_jht = 0.0

        # Minimum water release calculation
        Min_Water_Release = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Water_Release")
        q_0_jht = u_G_jht * Min_Water_Release

        # Assemble the row
        row_data = Any[day_id, h, t, unit_id, UNITGROUP, Unit_Category, Unit_Type, u_G_jht, u_P_jht, a_G_jht, a_P_jht, z_G_jht, z_P_jht, p_G_jht, p_P_jht, q_G_jht, q_P_jht, q_P_TOT_jht, q_RU_G_jht, q_RD_G_jht, q_SR_G_jht, r_RU_G_jht, r_RD_G_jht, r_SR_G_jht, q_0_jht]

        # Append each flow segment
        for l in 1:num_segments
            key_l = string("(", j, ", ", l, ", ", h, ", ", t, ")")
            push!(row_data, daily_solutions[string(day_id)]["solution"]["water"][key_l]["q_Gl_jht"])
        end

        dispatch_output_list[row, :] = row_data
        row += 1
    end

    # Hydro pumps
    for day_id in 1:num_days, h in ids_h, t in ids_t, unit_idx in eachindex(ids_j_pump)
        j   = ids_j_pump[unit_idx]
        key = string("(", j, ", ", h, ", ", t, ")")

        # Textual unit ID
        unit_id = id_to_name_hydro[j]

        # Unit metadata
        # Unit metadata
        bus_idx       = parameter(am, 0, :gen_index, "bus_idx",      j)
        tech_idx      = parameter(am, 0, :gen_index, "genco_tech_id", j)
        UNITGROUP     = parameter(am, bus_idx, :gen_bus, tech_idx, "UNITGROUP")
        Unit_Category = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_CATEGORY")
        Unit_Type     = parameter(am, bus_idx, :gen_bus, tech_idx, "UNIT_TYPE")


        # Hydro and water solution slices
        sol_h      = daily_solutions[string(day_id)]["solution"]["hydro"][key]
        sol_w      = daily_solutions[string(day_id)]["solution"]["water"][key]
        u_P_jht    = sol_h["u_P_jht"]
        a_P_jht    = sol_h["a_P_jht"]
        z_P_jht    = sol_h["z_P_jht"]
        p_P_jht    = sol_h["p_P_jht"]
        q_P_jht    = sol_w["q_P_jht"]
        q_P_TOT_jht = sol_w["q_P_TOT_jht"]

        # pump fields not applicable → 0.0
        u_G_jht = 0.0; a_G_jht = 0.0; z_G_jht = 0.0; p_G_jht = 0.0; q_G_jht = 0.0; q_RU_G_jht = 0.0; q_RD_G_jht = 0.0; q_SR_G_jht = 0.0; r_RU_G_jht = 0.0; r_RD_G_jht = 0.0; r_SR_G_jht = 0.0; q_0_jht = 0.0

        # Assemble the row
        row_data = Any[day_id, h, t, unit_id, UNITGROUP, Unit_Category, Unit_Type, u_G_jht, u_P_jht, a_G_jht, a_P_jht, z_G_jht, z_P_jht, p_G_jht, p_P_jht, q_G_jht, q_P_jht, q_P_TOT_jht, q_RU_G_jht, q_RD_G_jht, q_SR_G_jht, r_RU_G_jht, r_RD_G_jht, r_SR_G_jht, q_0_jht]

        # Append each flow segment
        for l in 1:num_segments
            push!(row_data, 0.0)
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
            reg_up_signal   = network_data["repdays"][day_id]["data"][h][t]["Delta_RU"]
            reg_down_signal = network_data["repdays"][day_id]["data"][h][t]["Delta_RD"]

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
                reg_up_signal   = network_data["repdays"][day_id]["data"][h][t]["Delta_RU"]
                reg_down_signal = network_data["repdays"][day_id]["data"][h][t]["Delta_RD"]
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
                reg_up_signal   = network_data["repdays"][day_id]["data"][h][t]["Delta_RU"]
                reg_down_signal = network_data["repdays"][day_id]["data"][h][t]["Delta_RD"]
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

    ALEAF_setting["Simulation Configuration"] = read_xlsx_return_dict_string_any(string(setting_file_location, setting_file_name), "Simulation Configuration"; first_row_value = 2)

    return ALEAF_setting
end

function generate_networkdata_HydroBoost(ALEAF_setting::Dict{String,<:Any}, project_id, case_id; year_id::Int64=1)

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

    # Collect market price time-series data ["DA_LMP", "Regulation_down", "Regulation_up", "Spin", "Delta_RU", "Delta_RD"]

    forecasting_method = ALEAF_setting["Simulation Setting"]["Market_price_forecasting_method"]  # Perfect_foresight, Mean_persistence, Additive_model_with_regressors, Additive_model_no_regressors, Autoregressive_with_regressors, Autoregressive_no_regressors, Manual_Forecast

    if forecasting_method in ["Perfect_foresight", "Mean_persistence", "Manual_Forecast"]

        market_price_time_series_data_list = ["DA_LMP", "Regulation_down", "Regulation_up", "Spin", "Delta_RU", "Delta_RD"]
        market_price_data_path = string(data_location, "Market/", forecasting_method, "/")

        for data_label in market_price_time_series_data_list
            timeSeries_df = CSV.read(string(market_price_data_path, data_label, ".csv"), DataFrame)
            allocate_market_price_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, data_label)
        end

    else
        # read LMP from selected method
        timeSeries_df = CSV.read(string(data_location, "Market/", forecasting_method, ".csv"), DataFrame)
        allocate_market_price_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, "DA_LMP")

        # read AS prices + signals from Mean_persistence folder
        market_price_time_series_data_list = ["Regulation_down", "Regulation_up", "Spin", "Delta_RU", "Delta_RD"]
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

    allocate_hydro_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, ["Diversion_C0", "Diversion_C1", "Water_Release_Requirement"])
    
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

    output_path = string(output_path, case_id, "/")
    check_and_create_path(output_path)

    network_data["output_path"] = output_path
    
    return network_data

end