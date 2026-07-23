#=

HydroBoost Model B: Hydroelectric System With Reservoir (Cascading Modeling) 

Main Authors:
Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov
Carlos Josue Lopez; Argonne National Laboratory; clopezsalgado@anl.gov
Alberto Grimaldi; Argonne National Laboratory; agrimaldi@anl.gov

Current version: 3.0
Last update from Alberto Grimaldi: 03.10.2026

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
    ALEAF_setting_of_case_id["Simulation Setting"]["look_ahead_days_value"]           = ALEAF_setting["Simulation Configuration"][key]["Look_Ahead_Days"]
    ALEAF_setting_of_case_id["Simulation Setting"]["Interconnection Limits Inflow"]   = ALEAF_setting["Simulation Configuration"][key]["Interconnection Limits Inflow"]
    ALEAF_setting_of_case_id["Simulation Setting"]["Interconnection Limits Outflow"]  = ALEAF_setting["Simulation Configuration"][key]["Interconnection Limits Outflow"]
    ALEAF_setting_of_case_id["Simulation Setting"]["Market_price_forecasting_method"] = ALEAF_setting["Simulation Configuration"][key]["Market_Price_File_ID"]

    for (uid, unit_any) in ALEAF_setting_of_case_id["Gen Technology - BESS"]
        unit = unit_any
        unit["Max_Power"]              = ALEAF_setting["Simulation Configuration"][key]["Max_Power"]
        unit["Max_SOC_MWh"]            = ALEAF_setting["Simulation Configuration"][key]["Max_SOC_MWh"]
        unit["Min_SOC_MWh"]            = ALEAF_setting["Simulation Configuration"][key]["Min_SOC_MWh"]
        unit["Max_Charge"]             = ALEAF_setting["Simulation Configuration"][key]["Max_Charge"]  
        unit["Roundtrip Efficiency"]   = ALEAF_setting["Simulation Configuration"][key]["Roundtrip Efficiency"]
        unit["Charging Efficiency"]    = ALEAF_setting["Simulation Configuration"][key]["Charging Efficiency"]
        unit["Discharging Efficiency"] = ALEAF_setting["Simulation Configuration"][key]["Discharging Efficiency"]
        unit["Maximum RU"]             = ALEAF_setting["Simulation Configuration"][key]["Maximum RU"]
        unit["Maximum RD"]             = ALEAF_setting["Simulation Configuration"][key]["Maximum RD"]
        unit["Maximum SR"]             = ALEAF_setting["Simulation Configuration"][key]["Maximum SR"]
        unit["AET"]                    = ALEAF_setting["Simulation Configuration"][key]["AET"]
    end

    return ALEAF_setting_of_case_id
end

function build_and_run_daily_HydroBoost(ALEAF_setting, network_data)

    daily_solutions = Dict{String, Any}()
    num_days = length(keys(network_data["repdays"]))

    # run daily HydroBoost sequantially
    # run day_id 1 first (without prior day solution)
    day_idx = 1
    daily_solutions[string(day_idx)] = build_and_run_HydroBoost_for_each_day(day_idx, ALEAF_setting, network_data)        
    
    # run remaining days using solution from day-1 as initial status
    # for day_id in [2, 3, 4]
    for day_id in 2:num_days
        daily_solutions[string(day_id)] = build_and_run_HydroBoost_for_each_day(day_id, ALEAF_setting, network_data; prior_day_solution=daily_solutions[string(day_id-1)]["solution"])
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
    ids_p = [(p) for (p) in get_index(am, :plant_index, 0)]                                                                                                                             # plant index, ### NEW for cascading modeling
	
    ids_k_ren_idx = [k for k in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", k) == "RENEWABLE"]
    id_to_name = Dict{Int,String}()
	#println("id_to_name:")
    for k in ids_k_ren_idx
        bus_idx  = parameter(am, 0, :gen_index, k, "bus_idx")
        tech_idx = parameter(am, 0, :gen_index, k, "genco_tech_id")
        id_to_name[k] = parameter(am, bus_idx, :gen_bus, tech_idx, "Tech_ID")
		#println(id_to_name[k])
    end
    ids_k_ren = [id_to_name[k] for k in ids_k_ren_idx ]# name of renewable
	#println("ids_k_ren=$ids_k_ren")

    ids_l = [(l) for (l) in 1:am.setting["Simulation Setting"]["num_hydropower_performance_segment_value"]]
    ids_h = [(h) for (h) in am.setting["run_H"]]
    ids_t = [(t) for (t) in am.setting["run_T"]]

    ###################################
    #------ Define decision variables (TOT = 51 variables)
    ###################################
    
    # Storage plant (BESS): dispatch variables
    HydroBoost_variable_iht_binary(JuMP_model, am, :storage, "u_B_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0,  bounded_upper=true, upper_bound=1.0)               # Binary variable driven to 1 when BESS is set to charging mode, and 0 otherwise      
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "e_B_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                     # Storage device state of charge in hour t [MWh]

    ### REMOVED ### # HydroBoost_variable_ht_real(JuMP_model, am, :storage, "e_B_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                # State of charge of all BESS units in hour t [MWh]
    ### REMOVED ### HydroBoost_variable_iht_real(JuMP_model, am, :storage, "e_B_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                         # Plant state of charge of all BESS units in hour t [MWh]

    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "p_B_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Power discharged from BESS and accounted for at point of delivery in period t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "p_B_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Power contributing to charge BESS in period t before accounting for losses [MWh]
    
    # HydroBoost_variable_ht_real(JuMP_model, am, :storage, "p_B_DT_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                             # Power discharged from BESS and accounted for at point of delivery in period t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "p_B_DT_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0)                                                        # Plant power discharged from BESS and accounted for at point of delivery in period t [MWh]	
	
    # HydroBoost_variable_ht_real(JuMP_model, am, :storage, "p_B_CT_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                             # Power contributing to charge BESS in period t before accounting for losses [MWh]
	HydroBoost_variable_iht_real(JuMP_model, am, :storage, "p_B_CT_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0)                                                        # Plant power from grid allocated to charge storage (BESS) (before considering any round-trip losses) [MWh]
	
    # HydroBoost_variable_iht_real(JuMP_model, am, :storage, "p_Gr_St_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0)                                                     # Plant power purchased from grid allocated to charge storage (BESS) [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :storage,  "p_Gr_St_ht",ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                              # NEW: System-wide bus bar variables (h,t indexed only -- no plant index) - Grid → Storage (system)

    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for regulation up in discharging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for regulation up in charging mode provided by BESS unit i in hour t  [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                    # Reserve for regulation up provided by BESS unit i in hour t [MWh]
    
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for regulation down in discharging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for regulation down in charging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                    # Reserve for regulation down provided by BESS unit i in hour t [MWh]
    
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_D_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for spinning reserve in discharging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_C_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Reserve for spinning reserve in charging mode provided by BESS unit i in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_iht", ids_i_sto, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                    # Reserve for spinning reserve provided by BESS unit i in hour t [MWh]
    
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RU_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                        # Plant reserve for regulation up provided by all BESS units in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_RD_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                        # Plant reserve for regulation down provided by all BESS units in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :storage, "r_SR_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                        # Plant reserve for spinning reserve provided by all BESS units in hour t [MWh]

    # Hydro plant: power dispatch variables    
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "u_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)             # Binary variable which is equal to 1 if hydro generator j is on-line in hour t, and 0 otherwise
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "a_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)             # Binary variable which is equal to 1 if hydro generator j is started at beginning of hour t, and 0 otherwise
    HydroBoost_variable_iht_binary(JuMP_model, am, :hydro, "z_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)             # Binary variable which is equal to 1 if hydro generator j is shut down at beginning of hour t, and 0 otherwise
    
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "p_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                     # Power output of hydro generator j in hour t [MWh]
    
    # HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_G_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                  # Day-ahead total power setpoint of dispatchable hydro plant in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "p_G_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                           # Plant day-ahead total power setpoint of dispatchable hydro plant in hour t [MWh]
    
    # HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_GB_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                 # Hydro power allocated to charge the BESS in hour t [MWh]
    # HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "p_GB_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                        # Plant hydro power allocated to charge the BESS in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro,   "p_GB_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                 # NEW: System-wide bus bar variables (h,t indexed only -- no plant index) - Hydro → Storage (system)
    
    # HydroBoost_variable_ht_real(JuMP_model, am, :hydro, "p_G_Gr_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                               # Dispatchable hydro power sold to the grid in the day-ahead market in hour t [MWh] 
    # HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "p_G_Gr_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                      # Plant dispatchable hydro power sold to the grid in the day-ahead market in hour t [MWh] 
    HydroBoost_variable_ht_real(JuMP_model, am, :hydro,   "p_G_Gr_ht",  ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                              # NEW: System-wide bus bar variables (h,t indexed only -- no plant index) - Hydro → Grid  (system)

    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_RU_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Hydro reserve for regulation up provided by hydro generator j in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_RD_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Hydro reserve for regulation down provided by hydro generator j in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_SR_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Hydro reserve for spinning reserve provided by hydro generator j in hour t [MWh]

    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_RU_G_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                        # Plant hydro reserve for regulation up provided by all hydro generators in hour t, when in generation mode [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_RD_G_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                        # Plant hydro reserve for regulation down provided by all hydro generators in hour t, when in generation mode [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_SR_G_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                        # Plant hydro reserve for spinning reserve provided by all hydro generators in hour t, when in generation mode [MWh]
    
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_RU_GB_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                       # Plant hydro + BESS reserve for regulation up provided by all hydro generators in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_RD_GB_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                       # Plant hydro + BESS reserve for regulation down provided by all hydro generators in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :hydro, "r_SR_GB_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                       # Plant hydro + BESS reserve for spinning reserve provided by all hydro generators in hour t [MWh]

    # Hydro plant: water flow variables
    HydroBoost_variable_ilht_real(JuMP_model, am, :water, "q_Gl_jht", ids_j_hydro, ids_l, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                            # Water discharge of block ℓ of hydro generator j in hour t [ft^3/s]
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "q_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                     # Water discharge of hydro generator j in hour t [ft^3/s]
    HydroBoost_variable_ilht_binary(JuMP_model, am, :water, "w_Gl_jht", ids_j_hydro, ids_l, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0)    # Binary variable which is equal to 1 if water discharged by hydro generator j has exceeded block ℓ in hour t
    
    # HydroBoost_variable_ht_real(JuMP_model, am, :water, "s_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                    # Spillage of reservoir in hour t [ft^3/s]
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "s_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                             # Plant spillage of reservoir in hour t [ft^3/s]
    
    # HydroBoost_variable_ht_real(JuMP_model, am, :water, "e_H_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                                  # Water volume of reservoir in period t [A-F]
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "e_H_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                           # Plant water volume of reservoir in period t [A-F]
    
	HydroBoost_variable_iht_real(JuMP_model, am, :water, "q_RU_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Flow equivalent of regulation up allocated to hydro generator j in hour t [ft^3/s]
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "q_RD_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Flow equivalent of regulation down allocated to hydro generator j in hour t [ft^3/s]
    HydroBoost_variable_iht_real(JuMP_model, am, :water, "q_SR_G_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Flow equivalent of spinning reserve allocated to hydro generator j in hour t [ft^3/s]
    
    # Hydro plant: rough zone variables
    if am.setting["Simulation Setting"]["hydropower_rough_zone_flag"] == true
        HydroBoost_variable_iht_real(JuMP_model, am, :rough_zone, "y_Gl_plus_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                      # Slack variables (upper bound) which take non-zero value if operation outside rough zone ℓ of unit j cannot be honored in hour t [ft^3/s]
        HydroBoost_variable_iht_real(JuMP_model, am, :rough_zone, "y_Gl_minus_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                     # Slack variables (lower bound) which take non-zero value if operation outside rough zone ℓ of unit j cannot be honored in hour t [ft^3/s]
        HydroBoost_variable_iht_binary(JuMP_model, am, :rough_zone, "phi_Gl_jht", ids_j_hydro, ids_h, ids_t; bounded_lower=true, lower_bound=0.0,  bounded_upper=true, upper_bound=1.0) # Auxiliary binary variable required for representation of rough zone constraint
    end

    # Non-dispatchable generators (RES plant): variables
	println("renewable_ids_k_ren:$ids_k_ren")
    HydroBoost_variable_iht_real(JuMP_model, am, :renewable, "p_R_kht", ids_k_ren, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                   # Power output of non-dispatchable generator 𝑘 in hour 𝑡 [MWh]
    
    # HydroBoost_variable_ht_real(JuMP_model, am, :renewable, "p_R_Gr_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                           # Power of non-dispatchable generators sent to grid in hour t [MWh]
    # HydroBoost_variable_iht_real(JuMP_model, am, :renewable, "p_R_Gr_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Plant power of non-dispatchable generators sent to grid in hour t [MWh]
    HydroBoost_variable_ht_real(JuMP_model, am, :renewable,"p_R_Gr_ht",  ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                             # NEW: System-wide bus bar variables (h,t indexed only -- no plant index) - RES → Grid  (system)

    HydroBoost_variable_iht_real(JuMP_model, am, :renewable, "p_Rspill_kht", ids_k_ren, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                              # Curtailed power of non-dispatchable generator 𝑘 in hour 𝑡 [MWh]
    
    # HydroBoost_variable_ht_real(JuMP_model,am, :renewable, "p_R_St_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                            # Power of non−dispatchable generators allocated to storage (BESS) in hour 𝑡 [MWh]
    # HydroBoost_variable_iht_real(JuMP_model, am, :renewable, "p_R_St_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                  # Plant power of non−dispatchable generators allocated to storage (BESS) in hour 𝑡 [MWh]
    HydroBoost_variable_ht_real(JuMP_model,am, :renewable, "p_R_St_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                              # NEW: System-wide bus bar variables (h,t indexed only -- no plant index) - RES → Storage (system)

    # HydroBoost_variable_ht_real(JuMP_model, am, :renewable, "p_R_ht", ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                              # Power of all non-dispatchable generators in hour t [MWh]
    HydroBoost_variable_iht_real(JuMP_model, am, :renewable, "p_R_pht", ids_p, ids_h, ids_t; bounded_lower=true, lower_bound=0.0)                                                       # Plant power output of non-dispatchable generator 𝑘 in hour 𝑡 [MWh]    
    
    ###################################
    #------ Call of constraints in the build-loop
    ###################################

    ### BESS constraints: from (6) to (29) ###
    for h in ids_h
        for t in ids_t

            # Total discharge and charge power of BESS units
			for p in ids_p # added
                HydroBoost_constraint_total_ES_power_discharge_pht(JuMP_model, am, "HydroBoost_constraint_total_ES_power_discharge_pht", p, h, t; const_name_flag)  # --- MOD
                HydroBoost_constraint_total_ES_power_charge_pht(JuMP_model, am, "HydroBoost_constraint_total_ES_power_charge_pht", p, h, t; const_name_flag)        # --- MOD
                ### REMOVED HydroBoost_constraint_total_ES_soc_pht(JuMP_model, am, "HydroBoost_constraint_total_ES_soc_pht", ids_p_sto, p, h, t; const_name_flag)		
            end

            # Total reg up, reg down, and spinning reserve provided by BESS units
            for p in ids_p # added
	            HydroBoost_constraint_total_ES_reg_up_reserve_pht(JuMP_model, am, "HydroBoost_constraint_total_ES_reg_up_reserve_pht", p, h, t; const_name_flag)    # --- MOD
	            HydroBoost_constraint_total_ES_reg_dn_reserve_pht(JuMP_model, am, "HydroBoost_constraint_total_ES_reg_dn_reserve_pht", p, h, t; const_name_flag)    # --- MOD
                HydroBoost_constraint_total_ES_spin_reserve_pht(JuMP_model, am, "HydroBoost_constraint_total_ES_spin_reserve_pht", p, h, t; const_name_flag)        # --- MOD
            end
			
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

    # (29) BESS energy throughput over the look-ahead window
    if am.setting["Simulation Setting"]["storage_AET_limit_flag"] == true
        for i in ids_i_sto
            HydroBoost_constraint_ES_AET_i(JuMP_model, am, "HydroBoost_constraint_ES_AET_i", i; const_name_flag)
        end
    end

    ### Hydro power system constraints: from (30) to (70) ###
    for h in ids_h
        for t in ids_t
            for p in ids_p
				# Generation of hydro power generators constraints
				HydroBoost_constraint_total_hydro_power_generation_jht(JuMP_model, am, "HydroBoost_constraint_total_hydro_power_generation_jht", p, h, t; const_name_flag)

                # Total reg up, reg down, and spinning reserve provided by hydro generators
                HydroBoost_constraint_total_reg_up_reserve_jht(JuMP_model, am, "HydroBoost_constraint_total_reg_up_reserve_jht", p, h, t; const_name_flag)
                HydroBoost_constraint_total_reg_dn_reserve_jht(JuMP_model, am, "HydroBoost_constraint_total_reg_dn_reserve_jht", p, h, t; const_name_flag)
                HydroBoost_constraint_total_spin_reserve_jht(JuMP_model, am, "HydroBoost_constraint_total_spin_reserve_jht", p, h, t; const_name_flag)
            end
            for j in ids_j_hydro
                # Commitment of hydro power generators constraints
                HydroBoost_constraint_hydro_commitment_status_jht(JuMP_model, am, "HydroBoost_constraint_hydro_commitment_status_jht", j, h, t, day_id, prior_day_solution; const_name_flag)
                HydroBoost_constraint_hydro_start_up_shut_down_bound_jht(JuMP_model, am, "HydroBoost_constraint_hydro_start_up_shut_down_bound_jht", j, h, t; const_name_flag)

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
            for p in ids_p# added
			    #ids_p_hydro = parameter(am, 0, :hydro_index, p)# gen_index for all the hydro units in plant p, added
				
				# Water storage and spillage bounds constraints
				HydroBoost_constraint_hydro_reservoir_volume_lower_bound_pht(JuMP_model, am, "HydroBoost_constraint_hydro_reservoir_volume_lower_bound_pht", p, h, t; const_name_flag)
				HydroBoost_constraint_hydro_reservoir_volume_upper_bound_pht(JuMP_model, am, "HydroBoost_constraint_hydro_reservoir_volume_upper_bound_pht", p, h, t; const_name_flag)
				HydroBoost_constraint_hydro_power_water_spillage_pht(JuMP_model, am, "HydroBoost_constraint_hydro_power_water_spillage_pht", p, h, t; const_name_flag)
				HydroBoost_constraint_hydro_power_water_release_lower_bound_pht(JuMP_model, am, "HydroBoost_constraint_hydro_power_water_release_lower_bound_pht", p, h, t; const_name_flag)

				# Water balance constraints
				HydroBoost_constraint_hydro_water_balance_pht(JuMP_model, am, "HydroBoost_constraint_hydro_water_balance_pht", p, h, t, day_id, prior_day_solution; const_name_flag)
            end
        end
    end

    # Initial and End-of-Period water volume constraints
    # HydroBoost_constraint_hydro_reservoir_uses_ini_limit(JuMP_model, am, "HydroBoost_constraint_hydro_reservoir_uses_ini_limit", ids_h, ids_t; const_name_flag) --> Already included in the water balance equation - constraint (55)
    for p in ids_p# added
	    HydroBoost_constraint_hydro_reservoir_uses_end_limit(JuMP_model, am, "HydroBoost_constraint_hydro_reservoir_uses_end_limit", p, ids_h, ids_t; const_name_flag)
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

    # Startup and shutdown daily limit of hydro generators
    num_days = am.setting["Simulation Setting"]["look_ahead_days_value"]

    for d in 1:num_days
        for j in ids_j_hydro
            HydroBoost_constraint_hydro_generators_daily_caps_jd(JuMP_model, am, "HydroBoost_constraint_hydro_generators_daily_caps_jd", j, d, ids_t; const_name_flag)
        end
    end

    ### RES generators constraints: from (71) to (73) ###
    for h in ids_h
        for t in ids_t
            for k in ids_k_ren
                HydroBoost_constraint_RES_balance_kht(JuMP_model, am,"HydroBoost_constraint_RES_balance_kht", k, h, t; const_name_flag)
            end
        end
    end


    for h in ids_h
        for t in ids_t
		    for p in ids_p
                HydroBoost_constraint_RES_total_pht(JuMP_model, am, "HydroBoost_constraint_RES_total_pht", p, h, t; const_name_flag)
			end
        end
    end

    ### Coupling constraints: from (74) to (81) --- 2ND MOD (June 2026) ###
    for h in ids_h
        for t in ids_t

            # (74)(75)(76): per-plant AS sale constraints
            for p in ids_p
                HydroBoost_constraint_reg_up_sale_pht(JuMP_model, am, "HydroBoost_constraint_reg_up_sale_pht", p, h, t; const_name_flag)
                HydroBoost_constraint_reg_dn_sale_pht(JuMP_model, am, "HydroBoost_constraint_reg_dn_sale_pht", p, h, t; const_name_flag)
                HydroBoost_constraint_spin_sale_pht(JuMP_model, am, "HydroBoost_constraint_spin_sale_pht", p, h, t; const_name_flag)
            end

            # (77)(78)(79): system-level power coupling at common bus bar (called once per h,t)
            HydroBoost_constraint_Hydro_power_generation_coupling_ht(JuMP_model, am, "HydroBoost_constraint_Hydro_power_generation_coupling_ht", ids_p, h, t; const_name_flag)
            HydroBoost_constraint_RES_power_generation_coupling_ht(JuMP_model, am, "HydroBoost_constraint_RES_power_generation_coupling_ht", ids_p, h, t; const_name_flag)
            HydroBoost_constraint_ES_power_generation_coupling_ht(JuMP_model, am, "HydroBoost_constraint_ES_power_generation_coupling_ht", ids_p, h, t; const_name_flag)

        end
    end

    # (80)(81): system-wide interconnection limits at common bus bar (called once per h,t)
    for h in ids_h
        for t in ids_t
            HydroBoost_constraint_outflow_interconnection_limit_ht(JuMP_model, am, "HydroBoost_constraint_outflow_interconnection_limit_ht", ids_p, h, t; const_name_flag)
            HydroBoost_constraint_inflow_interconnection_limit_ht(JuMP_model, am, "HydroBoost_constraint_inflow_interconnection_limit_ht", ids_p, h, t; const_name_flag)
        end
    end

    ###################################
    #------ Call of the objective function in the build-loop
    ###################################
    HydroBoost_objective_function(JuMP_model, am, ids_p, ids_j_hydro, ids_i_sto, ids_h, ids_t, ids_l)

    # export_JuMP_model_lp_file(JuMP_model) 

    am.model[:nw][0] = JuMP_model

end

###################################
#------ Define objective function
###################################

### Objective Function ###

# Constraints (1-5)

function HydroBoost_objective_function(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, ids_p, ids_j_hydro, ids_i_sto, ids_h, ids_t, ids_l)

    objective = JuMP.AffExpr(0.0)    

    # Energy revenue (plant-level, consistent with base Model_B formulation) --- 2ND MOD (June 2026)
    for h in ids_h
        for t in ids_t

            # system-level variables: called once per (h,t)
            JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :p_G_Gr_ht, (h,t)))
            JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :p_R_Gr_ht, (h,t)))
            JuMP.add_to_expression!(objective, -1.00001 * am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :p_Gr_St_ht, (h,t)))

            # p_B_DT_pht: still per-plant, summed over all plants
            for p in ids_p
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["DA_LMP"], get_variable(am, :p_B_DT_pht, (p,h,t)))
            end

        end
    end

    # AS revenue provided by BESS
    for h in ids_h
        for t in ids_t
            for p in ids_p
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Regulation_up"], get_variable(am, :r_RU_pht, (p,h,t)))
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Regulation_down"], get_variable(am, :r_RD_pht, (p,h,t)))
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Spin"], get_variable(am, :r_SR_pht, (p,h,t)))
            end
        end
    end

    # AS revenue provided by hydro generators
    for h in ids_h
        for t in ids_t
            for p in ids_p
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Regulation_up"], get_variable(am, :r_RU_G_pht, (p,h,t)))
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Regulation_down"], get_variable(am, :r_RD_G_pht, (p,h,t)))
                JuMP.add_to_expression!(objective, am.ref[:nw][0][:repdays]["data"][h][t]["Spin"], get_variable(am, :r_SR_G_pht, (p,h,t)))
            end
        end
    end

    # Adjusted revenue due to regulation deployments provided by BESS
    for h in ids_h
        for t in ids_t
            for p in ids_p
                reg_up_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RU"]
                reg_down_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RD"]

                JuMP.add_to_expression!(objective, reg_up_signal * am.ref[:nw][0][:repdays]["data"][h][t]["RT_LMP"], get_variable(am, :r_RU_pht, (p,h,t)))
                JuMP.add_to_expression!(objective, - reg_down_signal * am.ref[:nw][0][:repdays]["data"][h][t]["RT_LMP"], get_variable(am, :r_RD_pht, (p,h,t)))
            end
        end
    end

    # Adjusted revenue due to regulation deployments provided by hydro generators
    for h in ids_h
        for t in ids_t
            for p in ids_p
                reg_up_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RU"]
                reg_down_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RD"]

                JuMP.add_to_expression!(objective, reg_up_signal * am.ref[:nw][0][:repdays]["data"][h][t]["RT_LMP"], get_variable(am, :r_RU_G_pht, (p,h,t)))
                JuMP.add_to_expression!(objective, - reg_down_signal * am.ref[:nw][0][:repdays]["data"][h][t]["RT_LMP"], get_variable(am, :r_RD_G_pht, (p,h,t)))
            end
        end
    end

    # Variable O&M cost of hydro generators
    for j in ids_j_hydro
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)  
        VAR_OM_COST = parameter(am, bus_idx, :gen_bus, tech_idx, "VAR_OM_COST")

        for h in ids_h
            for t in ids_t
                reg_up_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RU"]
                reg_down_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RD"]

                JuMP.add_to_expression!(objective, - VAR_OM_COST, get_variable(am, :p_G_jht, (j,h,t)))
                JuMP.add_to_expression!(objective, - VAR_OM_COST * reg_up_signal, get_variable(am, :r_RU_G_jht, (j,h,t)))
                JuMP.add_to_expression!(objective, VAR_OM_COST * reg_down_signal, get_variable(am, :r_RD_G_jht, (j,h,t)))  
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

# Constraint (19) --- MOD
function HydroBoost_constraint_total_ES_power_discharge_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)# added
    
    # parameter
    ids_i_sto = parameter(am, 0, :sto_index, p)
    
    # variable
    p_B_DT_pht = get_variable(am, :p_B_DT_pht, (p,h,t))
        
    sum_power_output = JuMP.AffExpr(0.0)
    for i in ids_i_sto
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :p_B_D_iht, (i,h,t)))   
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_B_DT_pht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )  
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end  

end

# Constraint (20) --- MOD
function HydroBoost_constraint_total_ES_power_charge_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter
    ids_i_sto = parameter(am, 0, :sto_index, p)
    
    # variable
    p_B_CT_pht = get_variable(am, :p_B_CT_pht, (p,h,t))
        
    sum_power_output = JuMP.AffExpr(0.0)
    for i in ids_i_sto
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :p_B_C_iht, (i,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_B_CT_pht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

###REMOVED ### # # Constraint (21) --- MOD
# function HydroBoost_constraint_total_ES_soc_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_i_sto, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)

#     # parameter

#     # variable
#     e_B_pht = get_variable(am, :e_B_pht, (p,h,t))

#     sum_soc = JuMP.AffExpr(0.0)
#     for i in ids_i_sto
#         JuMP.add_to_expression!(sum_soc, 1.0, get_variable(am, :e_B_iht, (i,h,t)))
#     end

#     # constraint
#     expr = JuMP.@expression(JuMP_model,
#         e_B_pht - sum_soc
#         )
#     constraint = JuMP.@constraint(JuMP_model, 
#         expr == 0
#         )

#     if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end 
# end

# (f) Intial and final storage of BESS
# Constraint (21) and (22) are already included in constraint (6).

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

# Constraint (26)
function HydroBoost_constraint_total_ES_reg_up_reserve_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter
    ids_i_sto = parameter(am, 0, :sto_index, p)
    
    # variable
    r_RU_pht = get_variable(am, :r_RU_pht, (p,h,t))
        
    sum_power_output = JuMP.AffExpr(0.0)
    for i in ids_i_sto
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :r_RU_iht, (i,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RU_pht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# Constraint (27)
function HydroBoost_constraint_total_ES_reg_dn_reserve_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter
    ids_i_sto = parameter(am, 0, :sto_index, p)
    
    # variable
    r_RD_pht = get_variable(am, :r_RD_pht, (p,h,t))
        
    sum_power_output = JuMP.AffExpr(0.0)
    for i in ids_i_sto
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :r_RD_iht, (i,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RD_pht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# Constraint (28)
function HydroBoost_constraint_total_ES_spin_reserve_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter
    ids_i_sto = parameter(am, 0, :sto_index, p)
    
    # variable
    r_SR_pht = get_variable(am, :r_SR_pht, (p,h,t))
        
    sum_power_output = JuMP.AffExpr(0.0)
    for i in ids_i_sto
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :r_SR_iht, (i,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_SR_pht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# (h) BESS annual energy throughput limit over the optimization window:

# Constraint (29)
function HydroBoost_constraint_ES_AET_i(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, i::Int; const_name_flag::Bool=false)

    # sets of the optimization window
    ids_h = am.setting["run_H"]
    ids_t = am.setting["run_T"]

    # parameters
    bus_idx  = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)

    AET = parameter(am, bus_idx, :gen_bus, tech_idx, "AET")
    # BATEFF_C = parameter(am, bus_idx, :gen_bus, tech_idx, "Charging Efficiency")   * 0.01
    BATEFF_D = parameter(am, bus_idx, :gen_bus, tech_idx, "Discharging Efficiency") * 0.01

    # If no annual energy throughput limit is specified, skip the constraint
    if AET <= 0
        return
    end

    # Scale annual AET to the optimization window based on look_ahead_days_value
    look_ahead_days = am.setting["Simulation Setting"]["look_ahead_days_value"]

    year_days = 365.0
    if haskey(am.setting, "simulation_year_value") && (am.setting["simulation_year_value"] % 4 == 0)
        year_days = 366.0
    end

    AET_scale = look_ahead_days / year_days

    # Build total energy throughput over the optimization window
    sum_AET = JuMP.AffExpr(0.0)

    for h in ids_h
        for t in ids_t

            # dispatch variables
            # p_B_C_iht = get_variable(am, :p_B_C_iht, (i,h,t))
            p_B_D_iht = get_variable(am, :p_B_D_iht, (i,h,t))

            # reserve variables (aggregated RU / RD)
            # r_RD_iht = get_variable(am, :r_RD_iht, (i,h,t))
            r_RU_iht = get_variable(am, :r_RU_iht, (i,h,t))
            
            # regulation signals
            # reg_down_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RD"]
            reg_up_signal   = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RU"]

            # Charging contribution (energy going into the BESS)
            # JuMP.add_to_expression!(sum_AET, BATEFF_C * p_B_C_iht)

            # Discharging contribution (energy going out of the BESS)
            JuMP.add_to_expression!(sum_AET, (1.0 / BATEFF_D) * p_B_D_iht)

            # Regulation contribution
            # JuMP.add_to_expression!(sum_AET, reg_down_signal * r_RD_iht)
            JuMP.add_to_expression!(sum_AET, (1.0 / BATEFF_D) * reg_up_signal * r_RU_iht)
        
        end
    end

    # Annual throughput limit scaled to the optimization window
    expr = JuMP.@expression(JuMP_model,
        sum_AET - AET * AET_scale         # sum_AET - 2.0 * AET * AET_scale
    )

    constraint = JuMP.@constraint(JuMP_model, expr <= 0)

    if const_name_flag
        JuMP.set_name(constraint, string(const_name, "_($i)"))
    end

end

### Hydroelectric System with Reservoir ###

# (i) Commitment of hydro generators:

# Constraint (30)
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

# Constraint (31)
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

# Constraints (32) and (33) already included in section "Define decision variables".

# (j) Operational limits of hydro generators considering generation and ancillary services:

# Constraint (34)
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

# Constraint (35)
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

# Constraint (36)
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

# Constraint (37)
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

# Constraint (38)
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

# (k) Ramping constraints:

# Constraint (39)
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

# Constraint (40)
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

# (l) Water discharge of plant:

# Constraints (41) and (43)
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

# Constraints (42) and (44)
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

# Constraint (45)
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

# (m) Flow equivalent of ancillary services provided by hydro generators:

# Constraint (46)
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

# Constraint (47)
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

# (n) Generation of hydro power and allocation of ancillary services:

# Constraint (48)
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

# Constraint (49) --- MOD
function HydroBoost_constraint_total_hydro_power_generation_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
    ids_j_hydro = parameter(am, 0, :hydro_index, p)# gen_index for all the hydro units in plant p, added
	
    # variable
    p_G_pht = get_variable(am, :p_G_pht, (p,h,t))
        
    sum_power_output = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :p_G_jht, (j,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        p_G_pht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# Constraint (50)
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

# Constraint (51)
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

# Constraint (52)
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

# Constraint (53)
function HydroBoost_constraint_total_reg_up_reserve_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter
    ids_j_hydro = parameter(am, 0, :hydro_index, p)

    # variable
    r_RU_G_pht = get_variable(am, :r_RU_G_pht, (p,h,t))

    sum_power_output = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :r_RU_G_jht, (j,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RU_G_pht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# Constraint (54)
function HydroBoost_constraint_total_reg_dn_reserve_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter
    ids_j_hydro = parameter(am, 0, :hydro_index, p)

    # variable
    r_RD_G_pht = get_variable(am, :r_RD_G_pht, (p,h,t))

    sum_power_output = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :r_RD_G_jht, (j,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RD_G_pht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# Constraint (55)
function HydroBoost_constraint_total_spin_reserve_jht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter
    ids_j_hydro = parameter(am, 0, :hydro_index, p)

    # variable
    r_SR_G_pht = get_variable(am, :r_SR_G_pht, (p,h,t))

    sum_power_output = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_power_output, 1.0, get_variable(am, :r_SR_G_jht, (j,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_SR_G_pht - sum_power_output
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# (o) Water balance equation:

# Constraint (56) --- MOD
function HydroBoost_constraint_hydro_water_balance_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int, day_id, prior_day_solution; const_name_flag::Bool=false)
    
    # parameter
	#println("p=$p")
    Water_Inflow = am.ref[:nw][0][:repdays]["data"][h][t]["Inflow_"*string(p)]# revised
    CF = 3600 * 0.0000229569 # Factor to convert water flow units in [ft^3/s] into net volume units in [A-F/h]
    C_0 = am.ref[:nw][0][:repdays]["data"][h][t]["Diversion_C0_"*string(p)]
    C_1 = am.ref[:nw][0][:repdays]["data"][h][t]["Diversion_C1_"*string(p)]
	ids_j_hydro = parameter(am, 0, :hydro_index, p)# gen_index for all the hydro units in plant p, added
    #println("plant p:$p, hydro_id:$ids_j_hydro")
    reg_up_signal   = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RU"]
    reg_down_signal = am.ref[:nw][0][:repdays]["data"][h][t]["Delta_RD"]
    
    # variable
    s_pht = get_variable(am, :s_pht, (p,h,t))
	e_H_pht = get_variable(am, :e_H_pht, (p,h,t))# revised
    VMIN = am.ref[:nw][0][:hydro_reservoir][p]["VMIN"]
    VMAX = am.ref[:nw][0][:hydro_reservoir][p]["VMAX"]
    
    global prior_e_H_pht = (1 - am.ref[:nw][0][:hydro_reservoir][p]["phi_ini"]) * (VMAX - VMIN) + VMIN    # this value will be used for day id 1 at hour 1
    if h == 1
        if day_id > 1
			idx_string = string("(", p, ", ",24,  ", ", t,")")# revised
			global prior_e_H_pht = prior_day_solution["water"][idx_string]["e_H_pht"]   # get prior day solution for hour 1, revised
        end
    else    # if h >= 2
        global prior_e_H_pht = get_variable(am, :e_H_pht, (p,h-1,t))
    end

    sum_water_use = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :q_G_jht, (j,h,t)))
        JuMP.add_to_expression!(sum_water_use, reg_up_signal, get_variable(am, :q_RU_G_jht, (j,h,t)))
        JuMP.add_to_expression!(sum_water_use, -reg_down_signal, get_variable(am, :q_RD_G_jht, (j,h,t)))
    end
    
	########################
	# upstream spillage+discharge
	upstream_water_release = JuMP.AffExpr(0.0)# add	
	
	# upstream reservoir id
    id_upstream_raw = am.ref[:nw][0][:hydro_reservoir][p]["upstream_reservoir"]# add
	
    # Parse upstream reservoir id.
    # Allowed inputs: missing/""/0 (meaning no upstream), numeric, or strings like "Res_3".
    id_upstream = 0
    if !(id_upstream_raw isa Missing)
        if id_upstream_raw isa Integer
            id_upstream = Int(id_upstream_raw)
        elseif id_upstream_raw isa AbstractFloat
            id_upstream = Int(round(id_upstream_raw))
        else
            id_upstream_str = strip(string(id_upstream_raw))
            if isempty(id_upstream_str) || id_upstream_str == "0"
                id_upstream = 0
            else
                match_result = match(r"(?i)res_(\d+)", id_upstream_str)
                if match_result !== nothing
                    id_upstream = parse(Int, match_result.captures[1])
                else
                    # Fall back to numeric parsing if it's a plain number in a string
                    try
                        id_upstream = parse(Int, id_upstream_str)
                    catch
                        error("Invalid upstream_reservoir value for plant $(p): $(id_upstream_str). Expected missing/0 or Res_<n>.")
                    end
                end
            end
        end
    end
	
	h_lag = am.ref[:nw][0][:hydro_reservoir][p]["water_delay_time"]# add	
	#println("Current plant id:$p, upsteam id:$id_upstream, Water delay time:$h_lag")
	
	if id_upstream !=0 # the reservior has upstream reserviors, add	
	   	# upstream hydro unit ids:
    	ids_upstream_hydro = parameter(am, 0, :hydro_index, id_upstream)# gen_index for all the storages in plant p, add	
		
		# distance
		h_upstream = h - h_lag# add	
        
		#println("Current plant id:$p, upsteam id:$id_upstream, Water delay time:$h_lag, upstream time:$h_upstream")
		#println("Upstream plants:$ids_upstream_hydro")

		if h_upstream>=1# if the desired upstream spillage and discharge are within the same day, add			
			for j in ids_upstream_hydro# add	
				#JuMP.add_to_expression!(upstream_water_release,1.0, get_variable(am, :u_jht, (j,h_upstream,t)))# add	
				JuMP.add_to_expression!(upstream_water_release, 1.0, get_variable(am, :q_G_jht, (j,h_upstream,t)))
                JuMP.add_to_expression!(upstream_water_release, reg_up_signal, get_variable(am, :q_RU_G_jht, (j,h_upstream,t)))
                JuMP.add_to_expression!(upstream_water_release, -reg_down_signal, get_variable(am, :q_RD_G_jht, (j,h_upstream,t)))
			end# add	
			
			#spillage
			JuMP.add_to_expression!(upstream_water_release,1.0, get_variable(am, :s_pht, (id_upstream,h_upstream,t)))# add	
		else# if the desired upstream spillage and discharge are in the previous day, add		
		    if day_id>1# starting from second day, add	
                # discharge
			    for j in ids_upstream_hydro# add	
                    idx_string = string("(", j, ", ",24+h_upstream,  ", ", t,")")# add	
					prior_q_G_jht = prior_day_solution["water"][idx_string]["q_G_jht"]# add	
					#println("prior_q_G_jht:$prior_q_G_jht (j:$j h:$(24+h_upstream)")# add	
					
					prior_q_RU_G_jht = prior_day_solution["water"][idx_string]["q_RU_G_jht"]# add	
					prior_q_RD_G_jht = prior_day_solution["water"][idx_string]["q_RD_G_jht"]# add
					
					JuMP.add_to_expression!(upstream_water_release, 1.0, prior_q_G_jht)
					JuMP.add_to_expression!(upstream_water_release, reg_up_signal, prior_q_RU_G_jht)
					JuMP.add_to_expression!(upstream_water_release, -reg_down_signal, prior_q_RD_G_jht)					
			    end# add	
				
			    # spillage
			    idx_string = string("(", id_upstream, ", ",24+h_upstream,  ", ", t,")")# add	
			    prior_s_pht = prior_day_solution["water"][idx_string]["s_pht"]# add	
				
				#println("prior_s_pht:$prior_s_pht (p:$id_upstream h:$(24+h_upstream) prior_e_H_pht:$prior_e_H_pht Water_Inflow:$Water_Inflow at hour:$h")# add	
				JuMP.add_to_expression!(upstream_water_release,1.0, prior_s_pht)# add		
            #else# day 1
			    #h_upstream = 1

			    #for j in ids_upstream_hydro
				    #JuMP.add_to_expression!(upstream_water_release,1.0, get_variable(am, :u_jht, (j,h_upstream,t)))
			    #end
			
    			##spillage
	     		#JuMP.add_to_expression!(upstream_water_release,1.0, get_variable(am, :s_pht, (id_upstream,h_upstream,t)))				
			end# add	  			
		end# add
        #if h_upstream == 0 and day_id == 1
		    #println("prior_s_pht:$prior_s_pht (p:$id_upstream h:$(24+h_upstream) upstream_water_releaser:$upstream_water_release")# add
		#end 		
	end# add	
	########################
	
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        (e_H_pht - prior_e_H_pht) - CF * (Water_Inflow - sum_water_use - s_pht - C_0 - C_1 * e_H_pht+upstream_water_release)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
	if id_upstream == 1
	   if day_id == 1 
	      if h_upstream == 0
	         #println("prior_e_H_pht:$prior_e_H_pht CF:$CF Water_Inflow:$Water_Inflow sum_water_use:$sum_water_use s_pht:$s_pht C_0:$C_0 C_1:$C_1 upstream_water_release$upstream_water_release")# add
	      end
	   end
	end 
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end    

end

# (p) Water storage and spillage bounds:

# Constraint (57) --- MOD
function HydroBoost_constraint_hydro_reservoir_volume_lower_bound_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter
    VMIN = am.ref[:nw][0][:hydro_reservoir][p]["VMIN"]
    CF = 3600 * 0.0000229569 # Factor to convert water flow units in [ft^3/s] into net volume units in [A-F/s]
	
	ids_j_hydro = parameter(am, 0, :hydro_index, p)# gen_index for all the storages in plant p, added, DONE

    # variable
    e_H_pht = get_variable(am, :e_H_pht, (p,h,t))

    sum_water_use = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :q_RU_G_jht, (j,h,t)))
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :q_SR_G_jht, (j,h,t)))
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_H_pht - VMIN - CF * sum_water_use
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# Constraint (58) --- MOD
function HydroBoost_constraint_hydro_reservoir_volume_upper_bound_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)

    # parameter
    VMAX = am.ref[:nw][0][:hydro_reservoir][p]["VMAX"]
    CF = 3600 * 0.0000229569 # Factor to convert water flow units in [ft^3/s] into net volume units in [A-F/s]
	ids_j_hydro = parameter(am, 0, :hydro_index, p)# gen_index for all the hydro units in plant p, added, DONE

    # variable
    e_H_pht = get_variable(am, :e_H_pht, (p,h,t))

    sum_water_use = JuMP.AffExpr(0.0)
	
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :q_RD_G_jht, (j,h,t)))
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_H_pht - VMAX + CF * sum_water_use
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# Constraint (59) --- MOD
function HydroBoost_constraint_hydro_power_water_spillage_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter
    ids_j_hydro = parameter(am, 0, :hydro_index, p)# gen_index for all the hydro units in plant p, added, DONE
	
    # variable
    s_pht = get_variable(am, :s_pht, (p,h,t))

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
        s_pht - sum_water_use
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# Constraint (60) --- MOD
function HydroBoost_constraint_hydro_power_water_release_lower_bound_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter  
	keys_in_data = keys(am.ref[:nw][0][:repdays]["data"][h][t])
    #println("Keys available in am.ref[:nw][0][:repdays]['data'][$h][$t]: ", keys_in_data)
    Water_Release_Requirement = am.ref[:nw][0][:repdays]["data"][h][t]["Water_Release_Requirement_"*string(p)]
	#Water_Release_Requirement = am.ref[:nw][0][:repdays]["data"][h][t]["Water_Release_Requirement"]
	ids_j_hydro = parameter(am, 0, :hydro_index, p)# gen_index for all the hydro units in plant p, added, DONE
    
    # variable
    s_pht = get_variable(am, :s_pht, (p,h,t))

    sum_water_use = JuMP.AffExpr(0.0)
    for j in ids_j_hydro
        JuMP.add_to_expression!(sum_water_use, 1.0, get_variable(am, :q_G_jht, (j,h,t)))        
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        sum_water_use + s_pht - Water_Release_Requirement
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# (q) Initial and End-of-Period water volume constraints:

# Constraints (61) --> Already included in the water balance equation - constraint (57)
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

# Constraints (62) --- MOD
function HydroBoost_constraint_hydro_reservoir_uses_end_limit(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, ids_h, ids_t; const_name_flag::Bool=false)
    
    # Set/Index
    last_hour = last(ids_h)
    last_min = last(ids_t)

    # parameter
    VMIN = am.ref[:nw][0][:hydro_reservoir][p]["VMIN"]
    VMAX = am.ref[:nw][0][:hydro_reservoir][p]["VMAX"]
    phi_end = am.ref[:nw][0][:hydro_reservoir][p]["phi_end"]

    # variable
    e_H_pht = get_variable(am, :e_H_pht, (p,last_hour,last_min))

    # constraint
    expr = JuMP.@expression(JuMP_model,  
        e_H_pht - VMIN - (VMAX - VMIN) * (1 - phi_end)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($last_hour,$last_min)")) end    

end

# (r) Rough zone constraints:

# Constraint (63)
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

# Constraint (64)
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

# Constraint (65)
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

# Constraint (66)
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

# Constraints (67) and (68) already included in section "Define decision variables".

# Constraint (69-A)
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

# Constraint (69-B)
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

# (s) Startup and shutdown daily limit of hydro generators and hydro pumps

# Constraint (70)
function HydroBoost_constraint_hydro_generators_daily_caps_jd(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, j::Int, d::Int, ids_t; const_name_flag::Bool=false)

    # parameter
    bus_idx  = parameter(am, 0, :gen_index, "bus_idx", j)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
    SUcap_d = parameter(am, bus_idx, :gen_bus, tech_idx, "Max_Daily_Startup")

    num_hours_per_day = am.setting["Simulation Setting"]["num_hours_per_day_value"]
    h_start = (d-1) * num_hours_per_day + 1
    h_end   = d * num_hours_per_day

    # variables (sum on the day window)
    sum_a = JuMP.AffExpr(0.0)
    for h in h_start:h_end
        for t in ids_t
            JuMP.add_to_expression!(sum_a, 1.0, get_variable(am, :a_G_jht, (j,h,t)))
        end
    end

    # constraint
    expr1 = JuMP.@expression(JuMP_model, sum_a - SUcap_d)
    c1 = JuMP.@constraint(JuMP_model, expr1 <= 0)
    if const_name_flag JuMP.set_name(c1, string(const_name, "_SU($j,$d)")) end

end

### RES generators constraints ###

# (t) Non-dispatchable renewable generation constraints:

# Constraint (71)
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

# Constraint (72) --- MOD
function HydroBoost_constraint_RES_total_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false,)

    k = parameter(am, 0, :renw_index, p)# gen_index for all the hydro units in plant p, added, DONE
    # variable
    expr = JuMP.@expression(JuMP_model, sum(get_variable(am, :p_R_kht, (k_i,h,t)) for k_i in k))
    pTot = get_variable(am, :p_R_pht,(p,h,t))

    # constraint pTot == sum(p_R)
    c = JuMP.@constraint(JuMP_model, pTot == expr)

    if const_name_flag
        set_name(c, string(const_name, "_($p,$h,$t)"))
    end
end

# Constraint (73) already included in section "Define decision variables".

### Coupling constraints ###

# (u) Ancillary services offered to the market:

# Constraint (74)
function HydroBoost_constraint_reg_up_sale_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)

    # variable
    r_RU_GB_pht = get_variable(am, :r_RU_GB_pht, (p,h,t))
    r_RU_pht = get_variable(am, :r_RU_pht, (p,h,t))
    r_RU_G_pht = get_variable(am, :r_RU_G_pht, (p,h,t))
   
    # constraint
    expr = JuMP.@expression(JuMP_model,  
        r_RU_GB_pht - r_RU_pht - r_RU_G_pht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end    

end

# Constraint (75)
function HydroBoost_constraint_reg_dn_sale_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)
    
    # variable
    r_RD_GB_pht = get_variable(am, :r_RD_GB_pht, (p,h,t))
    r_RD_pht = get_variable(am, :r_RD_pht, (p,h,t))
    r_RD_G_pht = get_variable(am, :r_RD_G_pht, (p,h,t))
    
    # constraint
    expr = JuMP.@expression(JuMP_model,
        r_RD_GB_pht - r_RD_pht - r_RD_G_pht
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr == 0
        )
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end

end


# Constraint (76)
function HydroBoost_constraint_spin_sale_pht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, p::Int64, h::Int, t::Int; const_name_flag::Bool=false)

    # variable
    r_SR_GB_pht = get_variable(am, :r_SR_GB_pht, (p,h,t))
    r_SR_pht = get_variable(am, :r_SR_pht, (p,h,t))
    r_SR_G_pht = get_variable(am, :r_SR_G_pht, (p,h,t))

    # constraint
    expr = JuMP.@expression(JuMP_model,
        r_SR_GB_pht - r_SR_pht - r_SR_G_pht
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr == 0
        )
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($p,$h,$t)")) end

end

# (v) Hydro + BESS + RES power balances and interconnection limits:

# Constraint (77) --- MOD --- 2ND MOD (June 2026)
function HydroBoost_constraint_Hydro_power_generation_coupling_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_p, h::Int, t::Int; const_name_flag::Bool=false)

    # variable
    p_G_Gr_ht = get_variable(am, :p_G_Gr_ht, (h,t))
    p_GB_ht   = get_variable(am, :p_GB_ht,   (h,t))

    sum_G = JuMP.AffExpr(0.0)
    for p in ids_p
        JuMP.add_to_expression!(sum_G, 1.0, get_variable(am, :p_G_pht, (p,h,t)))
    end

    # constraint
    expr = JuMP.@expression(JuMP_model, sum_G - p_G_Gr_ht - p_GB_ht)
    constraint = JuMP.@constraint(JuMP_model, expr == 0)
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end
end

# Constraint (78) --- MOD --- 2ND MOD (June 2026)
function HydroBoost_constraint_RES_power_generation_coupling_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_p, h::Int, t::Int; const_name_flag::Bool=false)
    
    # variable
    p_R_Gr_ht = get_variable(am, :p_R_Gr_ht, (h,t))
    p_R_St_ht = get_variable(am, :p_R_St_ht, (h,t))

    sum_R = JuMP.AffExpr(0.0)
    for p in ids_p
        JuMP.add_to_expression!(sum_R, 1.0, get_variable(am, :p_R_pht, (p,h,t)))
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model, sum_R - p_R_Gr_ht - p_R_St_ht)
    constraint = JuMP.@constraint(JuMP_model, expr == 0)
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end
end

# Constraint (79) --- MOD --- 2ND MOD (June 2026)
function HydroBoost_constraint_ES_power_generation_coupling_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_p, h::Int, t::Int; const_name_flag::Bool=false)
    
    # variable
    p_GB_ht    = get_variable(am, :p_GB_ht,    (h,t))
    p_R_St_ht  = get_variable(am, :p_R_St_ht,  (h,t))
    p_Gr_St_ht = get_variable(am, :p_Gr_St_ht, (h,t))

    sum_BC = JuMP.AffExpr(0.0)
    for p in ids_p
        JuMP.add_to_expression!(sum_BC, 1.0, get_variable(am, :p_B_CT_pht, (p,h,t)))
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model, sum_BC - p_GB_ht - p_R_St_ht - p_Gr_St_ht)
    constraint = JuMP.@constraint(JuMP_model, expr == 0)
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end
end

# Constraint (80) --- MOD --- 2ND MOD (June 2026)
function HydroBoost_constraint_outflow_interconnection_limit_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_p, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter 
    outflow_limits = am.setting["Simulation Setting"]["Interconnection Limits Outflow"]

    # variable
    p_G_Gr_ht = get_variable(am, :p_G_Gr_ht, (h,t))
    p_R_Gr_ht = get_variable(am, :p_R_Gr_ht, (h,t))

    sum_BD_RU_SR = JuMP.AffExpr(0.0)
    for p in ids_p
        JuMP.add_to_expression!(sum_BD_RU_SR, 1.0, get_variable(am, :p_B_DT_pht,  (p,h,t)))
        JuMP.add_to_expression!(sum_BD_RU_SR, 1.0, get_variable(am, :r_RU_GB_pht, (p,h,t)))
        JuMP.add_to_expression!(sum_BD_RU_SR, 1.0, get_variable(am, :r_SR_GB_pht, (p,h,t)))
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model,
        sum_BD_RU_SR + p_G_Gr_ht + p_R_Gr_ht - outflow_limits)
    constraint = JuMP.@constraint(JuMP_model, expr <= 0)
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($h,$t)")) end
end

# Constraint (81) --- MOD --- 2ND MOD (June 2026)
function HydroBoost_constraint_inflow_interconnection_limit_ht(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, ids_p, h::Int, t::Int; const_name_flag::Bool=false)
    
    # parameter
    inflow_limits = am.setting["Simulation Setting"]["Interconnection Limits Inflow"]
    
    # variable
    p_Gr_St_ht = get_variable(am, :p_Gr_St_ht, (h,t))

    sum_RD = JuMP.AffExpr(0.0)
    for p in ids_p
        JuMP.add_to_expression!(sum_RD, 1.0, get_variable(am, :r_RD_GB_pht, (p,h,t)))
    end
    
    # constraint
    expr = JuMP.@expression(JuMP_model, p_Gr_St_ht + sum_RD - inflow_limits)
    constraint = JuMP.@constraint(JuMP_model, expr <= 0)
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
    HydroBoost_report_reservoir_volume_validation(ALEAF_model_instance, ALEAF_setting, daily_solutions, project_id, network_data)
    HydroBoost_report_result_RES_dispatch(ALEAF_model_instance, ALEAF_setting, daily_solutions, project_id, network_data)
    HydroBoost_report_result_objective(ALEAF_model_instance, ALEAF_setting, daily_solutions, project_id, network_data)

end

function HydroBoost_report_reservoir_volume_validation(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data)

    start_time = time()

    case_name = project_id
    output_path = network_data["output_path"]
    out_csv_name = string("ALEAF_HydroBoost_", case_name, "__reservoir_volume_validation.csv")
    out_summary_name = string("ALEAF_HydroBoost_", case_name, "__reservoir_volume_validation_summary.txt")

    ids_h = 1:ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    ids_t = 1:ALEAF_setting["Simulation Setting"]["num_sub_period_value"]
    ids_p = [(p) for (p) in get_index(am, :plant_index, 0)]

    # Determine days from available solutions (keys are Strings like "1", "2", ...)
    day_ids = sort(parse.(Int, collect(keys(daily_solutions))))
    num_days = length(day_ids)

    header = [
        "day", "hour", "time", "plant",
        "e_H_solver", "e_H_calc", "error_calc_minus_solver",
        "constraint_residual",
        "prior_e_H", "Water_Inflow", "sum_water_use", "s_pht", "C0", "C1", "upstream_water_release"
    ]

    total_rows = num_days * length(ids_h) * length(ids_t) * length(ids_p) + 1
    out = Array{Any}(undef, total_rows, length(header))
    out[1, :] = header

    CF = 3600 * 0.0000229569 # [ft^3/s] -> [A-F/h]

    row = 2
    # Track simple stats per plant
    stats = Dict{Int, Dict{String, Float64}}()
    for p in ids_p
        stats[p] = Dict(
            "max_abs_error" => 0.0,
            "max_abs_residual" => 0.0,
            "sum_abs_error" => 0.0,
            "sum_abs_residual" => 0.0,
            "count" => 0.0
        )
    end

    for day_id in day_ids
        day_key = string(day_id)
        for h in ids_h
            for t in ids_t
                reg_up_signal   = network_data["repdays"][day_id]["data"][h][t]["Delta_RU"]
                reg_down_signal = network_data["repdays"][day_id]["data"][h][t]["Delta_RD"]

                for p in ids_p
                    idx_pht = string("(", p, ", ", h, ", ", t, ")")

                    # solved variables
                    sol_w = daily_solutions[day_key]["solution"]["water"][idx_pht]
                    e_H_solver = sol_w["e_H_pht"]
                    s_pht = sol_w["s_pht"]

                    # parameters & inputs
                    Water_Inflow = network_data["repdays"][day_id]["data"][h][t]["Inflow_" * string(p)]
                    C0 = network_data["repdays"][day_id]["data"][h][t]["Diversion_C0_" * string(p)]
                    C1 = network_data["repdays"][day_id]["data"][h][t]["Diversion_C1_" * string(p)]

                    # prior volume (matches constraint logic)
                    VMIN = am.ref[:nw][0][:hydro_reservoir][p]["VMIN"]
                    VMAX = am.ref[:nw][0][:hydro_reservoir][p]["VMAX"]
                    prior_e_H = (1 - am.ref[:nw][0][:hydro_reservoir][p]["phi_ini"]) * (VMAX - VMIN) + VMIN
                    if h == 1
                        if day_id > 1
                            prior_key = string("(", p, ", ", 24, ", ", t, ")")
                            prior_e_H = daily_solutions[string(day_id - 1)]["solution"]["water"][prior_key]["e_H_pht"]
                        end
                    else
                        prior_key = string("(", p, ", ", h - 1, ", ", t, ")")
                        prior_e_H = daily_solutions[day_key]["solution"]["water"][prior_key]["e_H_pht"]
                    end

                    # sum water use by hydro units in plant p
                    sum_water_use = 0.0
                    ids_j_hydro = parameter(am, 0, :hydro_index, p)
                    for j in ids_j_hydro
                        key_jht = string("(", j, ", ", h, ", ", t, ")")
                        sol_j = daily_solutions[day_key]["solution"]["water"][key_jht]
                        sum_water_use += sol_j["q_G_jht"]
                        sum_water_use += reg_up_signal * sol_j["q_RU_G_jht"]
                        sum_water_use += -reg_down_signal * sol_j["q_RD_G_jht"]
                    end

                    # upstream water release (matches constraint logic)
                    upstream_water_release = 0.0
                    id_upstream_raw = am.ref[:nw][0][:hydro_reservoir][p]["upstream_reservoir"]
                    id_upstream = 0
                    if !(id_upstream_raw isa Missing)
                        if id_upstream_raw isa Integer
                            id_upstream = Int(id_upstream_raw)
                        elseif id_upstream_raw isa AbstractFloat
                            id_upstream = Int(round(id_upstream_raw))
                        else
                            id_upstream_str = strip(string(id_upstream_raw))
                            if !(isempty(id_upstream_str) || id_upstream_str == "0")
                                match_result = match(r"(?i)res_(\d+)", id_upstream_str)
                                if match_result !== nothing
                                    id_upstream = parse(Int, match_result.captures[1])
                                else
                                    try
                                        id_upstream = parse(Int, id_upstream_str)
                                    catch
                                        id_upstream = 0
                                    end
                                end
                            end
                        end
                    end

                    if id_upstream != 0
                        h_lag_raw = am.ref[:nw][0][:hydro_reservoir][p]["water_delay_time"]
                        h_lag = Int(round(h_lag_raw))
                        h_upstream = h - h_lag

                        ids_upstream_hydro = parameter(am, 0, :hydro_index, id_upstream)

                        if h_upstream >= 1
                            for j in ids_upstream_hydro
                                key_up = string("(", j, ", ", h_upstream, ", ", t, ")")
                                sol_up = daily_solutions[day_key]["solution"]["water"][key_up]
                                upstream_water_release += sol_up["q_G_jht"]
                                upstream_water_release += reg_up_signal * sol_up["q_RU_G_jht"]
                                upstream_water_release += -reg_down_signal * sol_up["q_RD_G_jht"]
                            end
                            spill_key = string("(", id_upstream, ", ", h_upstream, ", ", t, ")")
                            upstream_water_release += daily_solutions[day_key]["solution"]["water"][spill_key]["s_pht"]
                        else
                            if day_id > 1
                                hour_prev = 24 + h_upstream
                                for j in ids_upstream_hydro
                                    key_prev = string("(", j, ", ", hour_prev, ", ", t, ")")
                                    sol_prev = daily_solutions[string(day_id - 1)]["solution"]["water"][key_prev]
                                    upstream_water_release += sol_prev["q_G_jht"]
                                    upstream_water_release += reg_up_signal * sol_prev["q_RU_G_jht"]
                                    upstream_water_release += -reg_down_signal * sol_prev["q_RD_G_jht"]
                                end
                                spill_prev_key = string("(", id_upstream, ", ", hour_prev, ", ", t, ")")
                                upstream_water_release += daily_solutions[string(day_id - 1)]["solution"]["water"][spill_prev_key]["s_pht"]
                            end
                        end
                    end

                    # Derived values
                    # constraint residual using solved e_H (should be ~0)
                    constraint_residual = (e_H_solver - prior_e_H) - CF * (Water_Inflow - sum_water_use - s_pht - C0 - C1 * e_H_solver + upstream_water_release)
                    # calculated e_H using rearranged equation (handles implicit C1*e term)
                    denom = 1 + CF * C1
                    e_H_calc = (prior_e_H + CF * (Water_Inflow - sum_water_use - s_pht - C0 + upstream_water_release)) / denom
                    err = e_H_calc - e_H_solver

                    out[row, :] = [day_id, h, t, p, e_H_solver, e_H_calc, err, constraint_residual,
                                   prior_e_H, Water_Inflow, sum_water_use, s_pht, C0, C1, upstream_water_release]
                    row += 1

                    st = stats[p]
                    abs_err = abs(err)
                    abs_res = abs(constraint_residual)
                    st["max_abs_error"] = max(st["max_abs_error"], abs_err)
                    st["max_abs_residual"] = max(st["max_abs_residual"], abs_res)
                    st["sum_abs_error"] += abs_err
                    st["sum_abs_residual"] += abs_res
                    st["count"] += 1
                end
            end
        end
    end

    CSV.write(string(output_path, out_csv_name), Tables.table(out), writeheader=false)

    # Write a tiny summary text file
    open(string(output_path, out_summary_name), "w") do io
        write(io, "Reservoir volume validation summary for case: $(case_name)\n")
        write(io, "Days: $(num_days), Hours/day: $(length(ids_h)), Sub-periods: $(length(ids_t)), Plants: $(length(ids_p))\n")
        write(io, "\nPer-plant stats (based on e_H_calc vs e_H_solver and constraint residual):\n")
        for p in ids_p
            st = stats[p]
            cnt = max(st["count"], 1.0)
            mae = st["sum_abs_error"] / cnt
            mar = st["sum_abs_residual"] / cnt
            write(io, "Plant $(p): max|error|=$(st["max_abs_error"]), mean|error|=$(mae), max|residual|=$(st["max_abs_residual"]), mean|residual|=$(mar)\n")
        end
    end

    Memento.info(_LOGGER, "[ALEAF HydroBoost Model]:\t - Reservoir volume validation exported, Total Time(sec): $(round(time() - start_time, digits=2))")

end

function HydroBoost_report_result_plant_dispatch(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data) # --- MOD
    
    start_time = time()

    # Create and open a file
    case_name = project_id
    output_path = network_data["output_path"]
    dispatch_file_name = string("ALEAF_HydroBoost_",case_name,"__plant_dispatch.csv")
    
    # Write Label --- 2ND MOD (June 2026): _pht renamed to _ht for 5 system-level variables
    dispatch_label_list = ["day", "hour", "time", "plant", "LMP", "RT_LMP", "Reg Up Price", "Reg Dn Price", "Spin Price", "Delta RU", "Delta RD",
    "p_G_pht", "p_GB_ht", "p_G_Gr_ht", "r_RU_G_pht", "r_RD_G_pht", "r_SR_G_pht", "q_G_jht", "e_H_pht", "Inflow", "s_pht",
    "p_B_DT_pht", "p_B_CT_pht", "p_R_pht", "p_R_Gr_ht", "p_R_St_ht", "p_Gr_St_ht", "r_RU_pht", "r_RD_pht", "r_SR_pht",
    "r_RU_GB_pht", "r_RD_GB_pht", "r_SR_GB_pht"]

    # Write Outputs
    ids_j_hydro = [(j) for (j) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO"]
    ids_i_sto = [(i) for (i) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"]
    ids_k_ren = [(k) for k in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", k) == "RENEWABLE"]

    ids_l = [(l) for (l) in 1:am.setting["Simulation Setting"]["num_hydropower_performance_segment_value"]]
    ids_h = [(h) for (h) in 1:am.setting["Simulation Setting"]["num_hours_per_day_value"]]
    ids_t = [(t) for (t) in 1:am.setting["Simulation Setting"]["num_sub_period_value"]]
    
	ids_p = [(p) for (p) in get_index(am, :plant_index, 0)]    # NEW: plant index

    num_days = length(keys(network_data["repdays"]))

    #dispatch_output_list = Array{Any}(undef, num_days*length(ids_h)*length(ids_t)+1, length(dispatch_label_list))
	dispatch_output_list = Array{Any}(undef, num_days*length(ids_h)*length(ids_t)*length(ids_p)+1, length(dispatch_label_list))
    dispatch_output_list[1,:] = dispatch_label_list

    # for day_id in [1, 2, 3, 4]
    for day_id in 1:num_days
        for hour_id in ids_h
            for time_id in ids_t
				#I_ht = network_data["repdays"][day_id]["data"][hour_id][time_id]["Inflow"]

				# prices
                LMP = network_data["repdays"][day_id]["data"][hour_id][time_id]["DA_LMP"]
                RT_LMP = network_data["repdays"][day_id]["data"][hour_id][time_id]["RT_LMP"]
                Reg_up_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Regulation_up"]
                Reg_dn_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Regulation_down"]
                Spin_price = network_data["repdays"][day_id]["data"][hour_id][time_id]["Spin"]
                Delta_RU = network_data["repdays"][day_id]["data"][hour_id][time_id]["Delta_RU"]
                Delta_RD = network_data["repdays"][day_id]["data"][hour_id][time_id]["Delta_RD"]
 
                for plant_i in eachindex(ids_p)
                    p = ids_p[plant_i]
				    row_id = (day_id-1)*length(ids_h)*length(ids_t)*length(ids_p) + (hour_id-1)*length(ids_t)*length(ids_p) + (time_id-1)*length(ids_p) + (plant_i)+1
				    idx_pht = string("(", p, ", ", hour_id, ", ", time_id, ")") 
					#idx_ht = string("(", hour_id, ", ", time_id, ")")

					# get solutions --- 2ND MOD (June 2026): 5 system-level variables use idx_ht

					# system-level index (h, t) for common bus bar variables
					idx_ht = string("(", hour_id, ", ", time_id, ")")

					# storage
					p_B_DT_pht  = daily_solutions[string(day_id)]["solution"]["storage"][idx_pht]["p_B_DT_pht"]
					p_B_CT_pht  = daily_solutions[string(day_id)]["solution"]["storage"][idx_pht]["p_B_CT_pht"]
					p_Gr_St_ht  = daily_solutions[string(day_id)]["solution"]["storage"][idx_ht]["p_Gr_St_ht"]
                    r_RU_pht    = daily_solutions[string(day_id)]["solution"]["storage"][idx_pht]["r_RU_pht"]
                    r_RD_pht    = daily_solutions[string(day_id)]["solution"]["storage"][idx_pht]["r_RD_pht"]
                    r_SR_pht    = daily_solutions[string(day_id)]["solution"]["storage"][idx_pht]["r_SR_pht"]

					# hydro
					p_G_pht    = daily_solutions[string(day_id)]["solution"]["hydro"][idx_pht]["p_G_pht"]
					p_GB_ht    = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_GB_ht"]
					p_G_Gr_ht  = daily_solutions[string(day_id)]["solution"]["hydro"][idx_ht]["p_G_Gr_ht"]
                    r_RU_G_pht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_pht]["r_RU_G_pht"]
                    r_RD_G_pht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_pht]["r_RD_G_pht"]
                    r_SR_G_pht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_pht]["r_SR_G_pht"]

                    # Hydro + BESS
                    r_RU_GB_pht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_pht]["r_RU_GB_pht"]
                    r_RD_GB_pht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_pht]["r_RD_GB_pht"]
                    r_SR_GB_pht = daily_solutions[string(day_id)]["solution"]["hydro"][idx_pht]["r_SR_GB_pht"]

					# RES
					p_R_pht   = daily_solutions[string(day_id)]["solution"]["renewable"][idx_pht]["p_R_pht"]
					p_R_Gr_ht = daily_solutions[string(day_id)]["solution"]["renewable"][idx_ht]["p_R_Gr_ht"]
					p_R_St_ht = daily_solutions[string(day_id)]["solution"]["renewable"][idx_ht]["p_R_St_ht"]

					# water flow
					s_pht = daily_solutions[string(day_id)]["solution"]["water"][idx_pht]["s_pht"]
					e_H_pht = daily_solutions[string(day_id)]["solution"]["water"][idx_pht]["e_H_pht"]
					Water_Inflow = network_data["repdays"][day_id]["data"][hour_id][time_id]["Inflow_"*string(p)] 

					q_G_jht = 0.0
					ids_j_hydro = parameter(am, 0, :hydro_index, p)
					
					for j in ids_j_hydro
						q_G_jht += daily_solutions[string(day_id)]["solution"]["water"][string("(", j, ", ", hour_id, ", ", time_id, ")")]["q_G_jht"]
					end
					
					#row_id = (day_id-1)*length(ids_h)*length(ids_t) + (hour_id-1)*length(ids_t) + (time_id) + 1
					
					dispatch_output_list[row_id, :] = [day_id, hour_id, time_id, p, LMP, RT_LMP, Reg_up_price, Reg_dn_price, Spin_price, Delta_RU, Delta_RD,
						p_G_pht, p_GB_ht, p_G_Gr_ht, r_RU_G_pht, r_RD_G_pht, r_SR_G_pht, q_G_jht, e_H_pht, Water_Inflow, s_pht,
                        p_B_DT_pht, p_B_CT_pht, p_R_pht, p_R_Gr_ht, p_R_St_ht, p_Gr_St_ht, r_RU_pht, r_RD_pht, r_SR_pht,
                        r_RU_GB_pht, r_RD_GB_pht, r_SR_GB_pht]
				end	
            end
        end
    end
        
    CSV.write(string(output_path, dispatch_file_name), Tables.table(dispatch_output_list), writeheader=false)

    Memento.info(_LOGGER, "[ALEAF HydroBoost Model]:\t - Plant dispatch solution reporting, Total Time(sec): $(round(time() - start_time, digits=2))")

end###stop here

function HydroBoost_report_result_hydro_dispatch(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data)
    start_time = time()

    # Build output filename
    case_name          = project_id
    output_path        = network_data["output_path"]
    dispatch_file_name = string(output_path, "ALEAF_HydroBoost_", case_name, "__hydro_dispatch.csv")

    # Define CSV header
    dispatch_label_list = ["day", "hour", "time", "unit_id", "UnitGroup", "Unit_Category", "Unit_Type", "u_G_jht", "a_G_jht", "z_G_jht", 
    "p_G_jht", "q_G_jht", "q_RU_G_jht", "q_RD_G_jht", "q_SR_G_jht", "r_RU_G_jht", "r_RD_G_jht", "r_SR_G_jht", "q_0_jht"]
    num_segments = ALEAF_setting["Simulation Setting"]["num_hydropower_performance_segment_value"]
    for l in 1:num_segments
        push!(dispatch_label_list, "q_G$(l)_jht")
    end

    # Collect hydro unit indices
    ids_j_hydro = [j for j in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", j) == "HYDRO"]

    # Map each hydro index to its Tech_ID
    id_to_name_hydro = Dict{Int,String}()
    for j in ids_j_hydro
        bus_idx  = parameter(am, 0, :gen_index, "bus_idx", j)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
        id_to_name_hydro[j] = parameter(am, bus_idx, :gen_bus, tech_idx, "Tech_ID")
    end

    # Time indices
    ids_h      = 1:ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    ids_t      = 1:ALEAF_setting["Simulation Setting"]["num_sub_period_value"]
    repday_ids = sort(collect(keys(network_data["repdays"])))
    num_days   = length(repday_ids)

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

        # Minimum water release calculation
        Min_Water_Release = parameter(am, bus_idx, :gen_bus, tech_idx, "Min_Water_Release")
        q_0_jht = u_G_jht * Min_Water_Release

        # Assemble the row
        row_data = Any[day_id, h, t, unit_id, UNITGROUP, Unit_Category, Unit_Type, u_G_jht, a_G_jht, z_G_jht, p_G_jht, q_G_jht, q_RU_G_jht, q_RD_G_jht, q_SR_G_jht, r_RU_G_jht, r_RD_G_jht, r_SR_G_jht, q_0_jht]

        # Append each flow segment
        for l in 1:num_segments
            key_l = string("(", j, ", ", l, ", ", h, ", ", t, ")")
            push!(row_data, daily_solutions[string(day_id)]["solution"]["water"][key_l]["q_Gl_jht"])
        end

        dispatch_output_list[row, :] = row_data
        row += 1
    end

    # Write CSV without extra header
    CSV.write(dispatch_file_name, Tables.table(dispatch_output_list); writeheader=false)
    Memento.info(
        _LOGGER,
        "[ALEAF HydroBoost Model]: - Hydro dispatch solution reporting, Total Time(sec): $(round(time() - start_time, digits=2))"
    )
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
    ids_h      = 1:ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    ids_t      = 1:ALEAF_setting["Simulation Setting"]["num_sub_period_value"]
    repday_ids = sort(collect(keys(network_data["repdays"])))
    num_days   = length(repday_ids)

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

function HydroBoost_report_result_objective(am::Abstract_ALEAF_Model, ALEAF_setting, daily_solutions, project_id, network_data) # --- MOD
    start_time = time()

    # ---------- File path ----------
    case_name   = project_id
    output_path = network_data["output_path"]
    csv_file    = string(output_path, "ALEAF_HydroBoost_", case_name, "__objective_function.csv")

    # ---------- Indices ----------
    ids_h      = 1:ALEAF_setting["Simulation Setting"]["num_hours_per_day_value"]
    ids_t      = 1:ALEAF_setting["Simulation Setting"]["num_sub_period_value"]
    repday_ids = sort(collect(keys(network_data["repdays"])))
    num_days   = length(repday_ids)

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
              "energy_arbitrage", "energy_import_cost",
              "as_bess_reg_up", "as_bess_reg_down", "as_bess_spin", "as_bess",
              "as_hydro_reg_up", "as_hydro_reg_down", "as_hydro_spin", "as_hydro",
              "reg_adj_bess","reg_adj_hydro",
              "startup_shutdown_cost","rough_zone_penalty",
              "vom_cost",                            
              "energy_arbitrage_cum", "cum_energy_import_cost",
              "as_bess_reg_up_cum", "as_bess_reg_down_cum", "as_bess_spin_cum","as_bess_cum",
              "as_hydro_reg_up_cum", "as_hydro_reg_down_cum", "as_hydro_spin_cum", "as_hydro_cum",
              "reg_adj_bess_cum","reg_adj_hydro_cum",
              "startup_shutdown_cost_cum","rough_zone_penalty_cum",
              "vom_cost_cum",                         
              "obj_cumulative"]

    n_rows = num_days * length(ids_h) * length(ids_t) + 1
    out = Array{Any}(undef, n_rows, length(header))
    out[1,:] = header
    row = 2

    # ---------- Cumulative trackers (run across the whole year) ----------
    cum_obj       = 0.0
    cum_energy    = 0.0
    cum_energy_import_cost = 0.0
    cum_as_bess_reg_up = 0.0
    cum_as_bess_reg_down = 0.0
    cum_as_bess_spin = 0.0
    cum_as_bess   = 0.0
    cum_as_hydro_reg_up = 0.0
    cum_as_hydro_reg_down = 0.0
    cum_as_hydro_spin = 0.0
    cum_as_hydro  = 0.0
    cum_adj_bess  = 0.0
    cum_adj_hydro = 0.0
    cum_susd      = 0.0
    cum_rz        = 0.0
    cum_vom       = 0.0  

    # ---------- Main loops ----------
    ids_p = [(p) for (p) in get_index(am, :plant_index, 0)]     # NEW: plant index
	energy_arbitrage = 0.0
    for day_id in repday_ids
        day_key = string(day_id)
        for h in ids_h, t in ids_t
            #idx_ht = string("(", h, ", ", t, ")")
			
            # Prices
            LMP   = network_data["repdays"][day_id]["data"][h][t]["DA_LMP"]
            RT_LMP = network_data["repdays"][day_id]["data"][h][t]["RT_LMP"]
            Pru   = network_data["repdays"][day_id]["data"][h][t]["Regulation_up"]
            Prd   = network_data["repdays"][day_id]["data"][h][t]["Regulation_down"]
            Pspin = network_data["repdays"][day_id]["data"][h][t]["Spin"]
            reg_up_signal   = network_data["repdays"][day_id]["data"][h][t]["Delta_RU"]
            reg_down_signal = network_data["repdays"][day_id]["data"][h][t]["Delta_RD"]

            # Energy arbitrage term --- 2ND MOD (June 2026): system-level variables read once per (h,t)
            energy_arbitrage = 0.0
            energy_import_cost = 0.0

            idx_ht = string("(", h, ", ", t, ")")
            p_G_Gr_ht  = daily_solutions[day_key]["solution"]["hydro"][idx_ht]["p_G_Gr_ht"]
            p_R_Gr_ht  = daily_solutions[day_key]["solution"]["renewable"][idx_ht]["p_R_Gr_ht"]
            p_Gr_St_ht = daily_solutions[day_key]["solution"]["storage"][idx_ht]["p_Gr_St_ht"]
            energy_arbitrage  += LMP * (p_G_Gr_ht + p_R_Gr_ht - 1.00001 * p_Gr_St_ht)
            energy_import_cost = LMP * p_Gr_St_ht

            for p in ids_p
                idx_pht = string("(", p, ", ", h, ", ", t, ")")
                p_B_DT_pht = daily_solutions[day_key]["solution"]["storage"][idx_pht]["p_B_DT_pht"]
                energy_arbitrage += LMP * p_B_DT_pht
            end
			
            # AS: BESS (plant-level, aligned with objective function)
            as_bess_reg_up = 0.0
            as_bess_reg_down = 0.0
            as_bess_spin = 0.0
            as_bess = 0.0
            reg_adj_bess = 0.0

            for p in ids_p
                idx_pht = string("(", p, ", ", h, ", ", t, ")")
                sol_s = daily_solutions[day_key]["solution"]["storage"][idx_pht]
                rRU = sol_s["r_RU_pht"]; rRD = sol_s["r_RD_pht"]; rSR = sol_s["r_SR_pht"]
                as_bess_reg_up += Pru * rRU
                as_bess_reg_down += Prd * rRD
                as_bess_spin += Pspin * rSR
                as_bess += Pru * rRU + Prd * rRD + Pspin * rSR
                reg_adj_bess += reg_up_signal * RT_LMP * rRU - reg_down_signal * RT_LMP * rRD
            end

            # AS: Hydro (plant-level, aligned with objective function)
            as_hydro_reg_up = 0.0
            as_hydro_reg_down = 0.0
            as_hydro_spin = 0.0
            as_hydro = 0.0
            reg_adj_hydro = 0.0
            for p in ids_p
                idx_pht = string("(", p, ", ", h, ", ", t, ")")
                sol_h = daily_solutions[day_key]["solution"]["hydro"][idx_pht]
                rRUg = sol_h["r_RU_G_pht"]; rRDg = sol_h["r_RD_G_pht"]; rSRg = sol_h["r_SR_G_pht"]
                as_hydro_reg_up += Pru * rRUg
                as_hydro_reg_down += Prd * rRDg
                as_hydro_spin += Pspin * rSRg
                as_hydro += Pru * rRUg + Prd * rRDg + Pspin * rSRg
                reg_adj_hydro += reg_up_signal * RT_LMP * rRUg - reg_down_signal * RT_LMP * rRDg
            end

            # Startup/Shutdown cost: hydro generators (G-mode)
            startup_shutdown_cost = 0.0
            for j in ids_j_hydro
                key_jht = string("(", j, ", ", h, ", ", t, ")")
                sol_h = daily_solutions[day_key]["solution"]["hydro"][key_jht]

                # Safe getters for generator start-up/shut-down variables
                a_G = get_or0(sol_h, "a_G_jht")
                z_G = get_or0(sol_h, "z_G_jht")

                bus_idx  = parameter(am, 0, :gen_index, "bus_idx",      j)
                tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
                c_su = parameter(am, bus_idx, :gen_bus, tech_idx, "START_UP_COST")
                c_sd = parameter(am, bus_idx, :gen_bus, tech_idx, "SHUT_DN_COST")

                startup_shutdown_cost += -(c_su * a_G + c_sd * z_G)
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
                        rough_zone_penalty += -rz_penalty * (y_plus + y_minus)
                    end
                end
            end

            # VOM cost
            vom_cost = 0.0
            for j in ids_j_hydro
                key_jht = string("(", j, ", ", h, ", ", t, ")")
                sol_h = daily_solutions[day_key]["solution"]["hydro"][key_jht]

                pG   = sol_h["p_G_jht"]
                rRUg = sol_h["r_RU_G_jht"]
                rRDg = sol_h["r_RD_G_jht"]

                bus_idx  = parameter(am, 0, :gen_index, "bus_idx",      j)
                tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", j)
                λ_vom = parameter(am, bus_idx, :gen_bus, tech_idx, "VAR_OM_COST")

                vom_cost += -(λ_vom * (pG + reg_up_signal * rRUg - reg_down_signal * rRDg))
            end

            # Objective at (h,t)
			#println("energy_market_obj:$energy_market")
            obj_ht = energy_arbitrage + as_bess + as_hydro + reg_adj_bess + reg_adj_hydro + startup_shutdown_cost + rough_zone_penalty + vom_cost

            # ---------- Update cumulatives ----------
            cum_obj               += obj_ht
            cum_energy            += energy_arbitrage
            cum_energy_import_cost+= energy_import_cost
            cum_as_bess_reg_up    += as_bess_reg_up
            cum_as_bess_reg_down  += as_bess_reg_down
            cum_as_bess_spin      += as_bess_spin
            cum_as_bess           += as_bess
            cum_as_hydro_reg_up   += as_hydro_reg_up
            cum_as_hydro_reg_down += as_hydro_reg_down
            cum_as_hydro_spin     += as_hydro_spin
            cum_as_hydro          += as_hydro
            cum_adj_bess          += reg_adj_bess
            cum_adj_hydro         += reg_adj_hydro
            cum_susd              += startup_shutdown_cost
            cum_rz                += rough_zone_penalty
            cum_vom               += vom_cost  

            # ---------- Write row ----------
            out[row, :] = Any[
                day_id, h, t,
                obj_ht,
                energy_arbitrage, energy_import_cost,
                as_bess_reg_up, as_bess_reg_down, as_bess_spin, as_bess,
                as_hydro_reg_up, as_hydro_reg_down, as_hydro_spin, as_hydro,
                reg_adj_bess, reg_adj_hydro,
                startup_shutdown_cost, rough_zone_penalty, 
                vom_cost,                            
                cum_energy, cum_energy_import_cost,
                cum_as_bess_reg_up, cum_as_bess_reg_down, cum_as_bess_spin, cum_as_bess,
                cum_as_hydro_reg_up, cum_as_hydro_reg_down, cum_as_hydro_spin, cum_as_hydro,
                cum_adj_bess, cum_adj_hydro,
                cum_susd, cum_rz,
                cum_vom,                               
                cum_obj
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

    # Helper: fetch a field from Dict{String,Any} allowing for naming variations
    # (different case, underscores, or extra spaces in Excel headers)
    function _get_field(row::Dict{String,Any}, candidates::Vector{String})
        for k in candidates
            if haskey(row, k)
                return row[k]
            end
        end
        norm(s) = lowercase(strip(s))
        cand_norm = Set(norm.(candidates))
        for (k, v) in row
            if norm(k) in cand_norm
                return v
            end
        end
        error("Missing required field. Tried keys: $(candidates). Available keys: $(collect(keys(row)))")
    end

    # add gen_index
    am.ref[:nw][0][:gen_index] = Dict{Int64, Any}()
	
	# add gen_index in each plant
    am.ref[:nw][0][:hydro_index] = Dict{Int64,Vector{Int64}}()
	am.ref[:nw][0][:sto_index] = Dict{Int64,Vector{Int64}}()
	am.ref[:nw][0][:renw_index] = Dict{Int64, Vector{String}}()
	am.ref[:nw][0][:plant_index] = []
	
	
    gen_idx = 1
    num_sub_area = 1 # Assuming single bus system (might be extended in the future)
    
    for sub_area_idx in 1:num_sub_area
        for genco_tech_id in sort!(collect(keys(am.ref[:nw][sub_area_idx][:gen_bus])))
            row = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]
            UNIT_NAME = row["Tech_ID"]
                UNIT_TYPE = row["UNIT_TYPE"]
                UNIT_GROUP = row["UNITGROUP"]
                UNIT_CATEGORY = row["UNIT_CATEGORY"]
                CAP = row["Max_Power"]
			
            unit_id_val = _get_field(row, ["UnitID", "Unit_ID", "UNIT_ID", "unit_id", "UNITID"]) 
            unit_ID = unit_id_val isa Integer ? Int(unit_id_val) : unit_id_val isa AbstractFloat ? Int(round(unit_id_val)) : parse(Int, string(unit_id_val))

            plant_id_val = _get_field(row, ["PlantID", "Plant_ID", "PLANT_ID", "plant_id", "PLANTID"]) 
            plant_id_str = string(plant_id_val)

            # Extract plant index from strings like "Res_1" (case-insensitive),
            # or fall back to numeric parsing if the value is already numeric.
            match_result = match(r"(?i)res_(\d+)", plant_id_str)
            if match_result !== nothing
                plant_ID = parse(Int, match_result.captures[1])
            else
                try
                    plant_ID = parse(Int, strip(plant_id_str))
                catch
                    error("Invalid PlantID format: $(plant_id_str). Expected something like Res_1 or an integer.")
                end
            end

			#println("Tech_ID:$UNIT_NAME, Plant Index: $plant_ID")

			#plant_ID = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["PlantID"]
			
            am.ref[:nw][0][:gen_index][gen_idx] = Dict{String, Any}("bus_idx" => sub_area_idx, 
                                                                    "genco_tech_id" => genco_tech_id, 
                                                                    "genco_tech_UNIT_Type" => UNIT_TYPE, 
                                                                    "UNIT_GROUP" => UNIT_GROUP,
                                                                    "UNIT_CATEGORY" => UNIT_CATEGORY,
                                                                    "CAP" => CAP,
																	"PlantID" => plant_ID,
																	"Unit_ID" => unit_ID
                                                                    )
            am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["gen_idx"] = gen_idx
			
			push!(am.ref[:nw][0][:plant_index],plant_ID)
            
			if UNIT_CATEGORY == "HYDRO"
				if !haskey(am.ref[:nw][0][:hydro_index], plant_ID)
					am.ref[:nw][0][:hydro_index][plant_ID] = Int[]  # Initialize if not exists
				end
				#push!(am.ref[:nw][0][:hydro_index][plant_ID], unit_ID)  # Use push! instead of append
				push!(am.ref[:nw][0][:hydro_index][plant_ID], gen_idx)  # Use push! instead of append
				#println("gen_idx(HYDRO):$gen_idx plant_ID:$plant_ID, unit_ID:$unit_ID")
			elseif UNIT_CATEGORY == "STORAGE"
				if !haskey(am.ref[:nw][0][:sto_index], plant_ID)
					am.ref[:nw][0][:sto_index][plant_ID] = Int[]  # Initialize if not exists
				end
				#push!(am.ref[:nw][0][:sto_index][plant_ID], unit_ID)  # Use push! instead of append
				push!(am.ref[:nw][0][:sto_index][plant_ID], gen_idx)  # Use push! instead of append
				#println("gen_idx(STORAGE):$gen_idx plant_ID:$plant_ID, unit_ID:$unit_ID")
			elseif UNIT_CATEGORY == "RENEWABLE"
				if !haskey(am.ref[:nw][0][:renw_index], plant_ID)
					#println("gen_idx(RENEWABLE):$gen_idx")
				    #println("plant_ID:$plant_ID, unit_ID:$unit_ID")	
                    am.ref[:nw][0][:renw_index][plant_ID] = String[]  # Initialize if not exists
				end		
				#println("gen_idx(RENEWABLE):$gen_idx")
				#println("plant_ID:$plant_ID, unit_ID:$unit_ID")	 				
                #push!(am.ref[:nw][0][:renw_index][plant_ID], gen_idx)  # Use push! instead of append			
                push!(am.ref[:nw][0][:renw_index][plant_ID], UNIT_NAME)  # Use push! instead of append				
			end			
            gen_idx += 1
        end
		# Make sure plant_index contains unique entries
        am.ref[:nw][0][:plant_index] = unique(am.ref[:nw][0][:plant_index])
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
    
    # Normalize Water Release Requirement column names.
    # (The input spreadsheet can have either the legacy typo "Requirment" or the correct "Requirement",
    # and can be either a single column or numbered: "... 1", "... 2", ...)
    if "Water Release Requirment" in names(timeSeries)
        rename!(timeSeries, "Water Release Requirment" => "Water_Release_Requirement")
    elseif "Water Release Requirement" in names(timeSeries)
        rename!(timeSeries, "Water Release Requirement" => "Water_Release_Requirement")
    end

    for col in collect(names(timeSeries))
        m1 = match(r"^Water Release Requirment\s+(\d{1,2})$", col)
        if m1 !== nothing
            new_name = string("Water_Release_Requirement_", m1.captures[1])
            if !(new_name in names(timeSeries))
                rename!(timeSeries, col => new_name)
            end
            continue
        end

        m2 = match(r"^Water Release Requirement\s+(\d{1,2})$", col)
        if m2 !== nothing
            new_name = string("Water_Release_Requirement_", m2.captures[1])
            if !(new_name in names(timeSeries))
                rename!(timeSeries, col => new_name)
            end
        end
    end

    # Normalize data_label_list to match any renames performed above
    normalized_labels = String[]
    for lbl_any in data_label_list
        lbl = string(lbl_any)
        if lbl == "Water Release Requirment" || lbl == "Water Release Requirement"
            push!(normalized_labels, ("Water_Release_Requirement" in names(timeSeries)) ? "Water_Release_Requirement" : lbl)
            continue
        end

        m1 = match(r"^Water Release Requirment\s+(\d{1,2})$", lbl)
        if m1 !== nothing
            new_lbl = string("Water_Release_Requirement_", m1.captures[1])
            push!(normalized_labels, (new_lbl in names(timeSeries)) ? new_lbl : lbl)
            continue
        end

        m2 = match(r"^Water Release Requirement\s+(\d{1,2})$", lbl)
        if m2 !== nothing
            new_lbl = string("Water_Release_Requirement_", m2.captures[1])
            push!(normalized_labels, (new_lbl in names(timeSeries)) ? new_lbl : lbl)
            continue
        end

        push!(normalized_labels, lbl)
    end
    data_label_list = normalized_labels
    
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
                    if !haskey(network_data["repdays"][day_id]["data"][hour_idx], t)
                        network_data["repdays"][day_id]["data"][hour_idx][t] = Dict{Any, Any}()
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

    # Optional test mode: set this in Simulation Setting to run a shorter horizon (e.g., 8 days)
    # ALEAF_setting["Simulation Setting"]["test_num_days_value"] = 8
    test_num_days_value = get(ALEAF_setting["Simulation Setting"], "test_num_days_value", 0)
    test_num_days = try
        Int(round(Float64(test_num_days_value)))
    catch
        0
    end
    if test_num_days > 0
        total_num_days = min(test_num_days, total_num_days)
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

    # Collect market price time-series data ["DA_LMP", "RT_LMP", "Regulation_down", "Regulation_up", "Spin", "Delta_RU", "Delta_RD"]
    
    forecasting_method = ALEAF_setting["Simulation Setting"]["Market_price_forecasting_method"]     # Perfect_foresight, Mean_persistence, Additive_model_with_regressors, Additive_model_no_regressors, Autoregressive_with_regressors, Autoregressive_no_regressors, Manual_Forecast

    if forecasting_method in ["Perfect_foresight", "Mean_persistence", "Manual_Forecast"]
        
        market_price_time_series_data_list = ["DA_LMP", "RT_LMP", "Regulation_down", "Regulation_up", "Spin", "Delta_RU", "Delta_RD"]
        market_price_data_path = string(data_location, "Market/", forecasting_method, "/")

        for data_label in market_price_time_series_data_list
            timeSeries_df = CSV.read(string(market_price_data_path, data_label, ".csv"), DataFrame)
            allocate_market_price_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, data_label)
        end

    else

        # read LMP
        timeSeries_df = CSV.read(string(data_location, "Market/", forecasting_method, ".csv"), DataFrame)
        allocate_market_price_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, "DA_LMP")

        timeSeries_df = CSV.read(string(data_location, "Market/", forecasting_method, ".csv"), DataFrame)
        allocate_market_price_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, "RT_LMP")
        
        # read AS prices from Mean_persistence folder
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

	#column_names = names(timeSeries_df)
    ##println(column_names)
	
	## Define the patterns to match
    ##patterns = r"Diversion|Water Release Requirment|Inflow\(plant\d+\)"
	#patterns = r"Diversion|Water Release Requirment|Inflow"

    ## Use filter to select columns that match the patterns
    #selected_columns = filter(col -> occursin(patterns, col), column_names)
	
	# Replace "Water Release Requirment" with "Water_Release_Requirement"
    #Updated_columns= replace(selected_columns, "Water Release Requirment" => "Water_Release_Requirement")
	
	#println("Updated columns1: ", Updated_columns)
	# Iterate through each column name in timeSeries
    for col in names(timeSeries_df)
	    #println(col)
        # Check if the column name matches the pattern "Water Release Requirment x"
        if occursin(r"^Water Release Requirment \d{1,2}$", col)
            # Extract the number x from the column name
            num = match(r"\d{1,2}", col).match
            # Construct the new column name
            new_name = "Water_Release_Requirement_$num"
            # Rename the column
            rename!(timeSeries_df, col => new_name)
			#println("new_name:$new_name")
        end
    end

    ## Print the updated column names
    #for col in names(timeSeries_df)
	    #println("after change:$col")
    #end

	column_names = names(timeSeries_df)
    #println("column_names:$column_names")
	
	## Define the patterns to match
    ##patterns = r"Diversion|Water Release Requirment|Inflow\(plant\d+\)"
	#patterns = r"Diversion|Water Release Requirment|Inflow"

    ## Use filter to select columns that match the patterns
    #Updated_columns = filter(col -> occursin(patterns, col), column_names)
	#println("Updated_columns2:$Updated_columns")
	
	
    ## Print the selected columns
    ##println("----------------column names---------:$Updated_columns")
    
    allocate_hydro_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, column_names)
	
    #allocate_hydro_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, ["Inflow", "Diversion_C0", "Diversion_C1", "Water_Release_Requirement"])
    
    # collect RES profiles time-series data ["RES - Hourly"]
    res_data_path = string(data_location, "RES/")

    # --- "RES - Hourly"
    timeSeries_df = CSV.read(string(res_data_path, "RES_profiles.csv"), DataFrame)
    # Get the list of column names, excluding 'Date-Time'
    res_list = filter(name -> name != "Date-Time", names(timeSeries_df))
	#println("Column names excluding 'Date-Time': ", res_list)
	
    #allocate_res_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, ["RES_1_PV", "RES_2_PV", "RES_3_WIND", "RES_4_WIND"])
	allocate_res_sub_hourly_timeseries_data!(timeSeries_df, network_data, ALEAF_setting, res_list)

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