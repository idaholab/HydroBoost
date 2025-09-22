"""
Generate O&M Cost Calculator

This script imports Operation & Maintenance cost data from Excel input files
and calculates actual O&M costs based on simulation results.

Supports different model configurations with various asset structures:
- Run-of-River hydro + battery configuration
- Reservoir hydro + battery configuration  
- Reservoir hydro + battery + pumped storage hydro configuration

Features:
- Import O&M cost parameters from Excel
- Calculate hydro O&M costs from dispatch data (startup, shutdown, ramping)
- Generate hourly cost data for spreadsheet integration
- Provide annual cost summaries

Author: HydroBoost Project
"""

import pandas as pd
import numpy as np
import os
import argparse
from typing import Tuple, Dict


def extract_om_cost_rates(hydro_data: pd.DataFrame, battery_data: pd.DataFrame, psh_data: pd.DataFrame) -> Dict[str, float]:
    """
    Extract O&M cost rates from the imported data.
    
    Args:
        hydro_data (pd.DataFrame): Hydro O&M cost data
        battery_data (pd.DataFrame): Battery O&M cost data
        psh_data (pd.DataFrame): PSH O&M cost data
        
    Returns:
        Dict[str, float]: Dictionary of cost rates
    """
    cost_rates = {}
    
    # Extract hydro cost rates
    for _, row in hydro_data.iterrows():
        cost_type = str(row.iloc[0]).strip()
        cost_value = float(row.iloc[1])
        cost_rates[f"hydro_{cost_type}"] = cost_value
    
    # Extract battery cost rates if available
    if not battery_data.empty:
        for _, row in battery_data.iterrows():
            cost_type = str(row.iloc[0]).strip()
            cost_value = float(row.iloc[1])
            cost_rates[f"battery_{cost_type}"] = cost_value
    
    # Extract PSH cost rates if available
    if not psh_data.empty:
        for _, row in psh_data.iterrows():
            cost_type = str(row.iloc[0]).strip()
            cost_value = float(row.iloc[1])
            cost_rates[f"psh_{cost_type}"] = cost_value
    
    return cost_rates


def import_hydro_dispatch_data(filepath: str) -> pd.DataFrame:
    """
    Import hydro dispatch data from CSV file.
    
    Args:
        filepath (str): Path to the hydro dispatch CSV file
        
    Returns:
        pd.DataFrame: Hydro dispatch data
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Hydro dispatch file not found: {filepath}")
    
    # Read CSV with chunking for large files
    try:
        dispatch_data = pd.read_csv(filepath)
        print(f"Successfully imported hydro dispatch data: {len(dispatch_data)} rows")
        return dispatch_data
    except Exception as e:
        raise Exception(f"Error reading hydro dispatch file: {str(e)}")


def import_storage_dispatch_data(filepath: str) -> pd.DataFrame:
    """
    Import battery storage dispatch data from CSV file.
    
    Args:
        filepath (str): Path to the storage dispatch CSV file
        
    Returns:
        pd.DataFrame: Storage dispatch data
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Storage dispatch file not found: {filepath}")
    
    # Read CSV with chunking for large files
    try:
        dispatch_data = pd.read_csv(filepath)
        print(f"Successfully imported storage dispatch data: {len(dispatch_data)} rows")
        return dispatch_data
    except Exception as e:
        raise Exception(f"Error reading storage dispatch file: {str(e)}")


def calculate_startup_shutdown_events(dispatch_data: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate startup and shutdown events based on commitment transitions.
    Only processes hydro generator units (hydro_t*), ignoring pump units.
    
    Args:
        dispatch_data (pd.DataFrame): Hydro dispatch data
        
    Returns:
        pd.DataFrame: Data with startup/shutdown event flags
    """
    result_data = dispatch_data.copy()
    
    # Filter to only hydro generator units (hydro_t*), exclude pumps (hydro_p*)
    generator_units = [unit for unit in dispatch_data['unit_id'].unique() if 'hydro_t' in unit or ('hydro_' in unit and 'hydro_p' not in unit)]
    
    # Group by unit_id to handle multiple hydro generator units
    for unit_id in generator_units:
        unit_mask = dispatch_data['unit_id'] == unit_id
        unit_data = dispatch_data[unit_mask].copy()
        
        # Sort by time to ensure proper order
        unit_data = unit_data.sort_values(['day', 'hour'])
        
        # Use u_G_jht for generator commitment
        commitment = unit_data['u_G_jht'].values
        
        # Initialize event arrays
        startup_events = np.zeros(len(unit_data))
        shutdown_events = np.zeros(len(unit_data))
        
        # Detect transitions
        for i in range(1, len(commitment)):
            # Startup: 0 -> 1 transition
            if commitment[i-1] == 0 and commitment[i] == 1:
                startup_events[i] = 1
            # Shutdown: 1 -> 0 transition  
            elif commitment[i-1] == 1 and commitment[i] == 0:
                shutdown_events[i] = 1
        
        # Update the result dataframe
        result_data.loc[unit_mask, 'startup_event'] = startup_events
        result_data.loc[unit_mask, 'shutdown_event'] = shutdown_events
    
    return result_data


def calculate_ramping_costs(dispatch_data: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate ramping up and down amounts for each time step.
    
    Args:
        dispatch_data (pd.DataFrame): Hydro dispatch data
        
    Returns:
        pd.DataFrame: Data with ramping amounts
    """
    result_data = dispatch_data.copy()
    
    # Group by unit_id to handle multiple hydro units
    for unit_id in dispatch_data['unit_id'].unique():
        unit_mask = dispatch_data['unit_id'] == unit_id
        unit_data = dispatch_data[unit_mask].copy()
        
        # Sort by time to ensure proper order
        unit_data = unit_data.sort_values(['day', 'hour'])
        
        # Calculate power differences
        power = unit_data['p_G_jht'].values
        
        # Initialize ramping arrays
        ramp_up_mw = np.zeros(len(unit_data))
        ramp_down_mw = np.zeros(len(unit_data))
        
        # Calculate ramping for each time step
        for i in range(1, len(power)):
            power_diff = power[i] - power[i-1]
            
            if power_diff > 0:
                ramp_up_mw[i] = power_diff
            elif power_diff < 0:
                ramp_down_mw[i] = abs(power_diff)
        
        # Update the result dataframe
        result_data.loc[unit_mask, 'ramp_up_mw'] = ramp_up_mw
        result_data.loc[unit_mask, 'ramp_down_mw'] = ramp_down_mw
    
    return result_data


def calculate_hourly_om_costs(dispatch_data: pd.DataFrame, cost_rates: Dict[str, float]) -> pd.DataFrame:
    """
    Calculate hourly O&M costs for each time step.
    
    Args:
        dispatch_data (pd.DataFrame): Hydro dispatch data with events and ramping
        cost_rates (Dict[str, float]): O&M cost rates
        
    Returns:
        pd.DataFrame: Data with hourly O&M costs
    """
    result_data = dispatch_data.copy()
    
    # Calculate costs for each row
    startup_costs = result_data['startup_event'] * cost_rates.get('hydro_startup_cost', 0)
    shutdown_costs = result_data['shutdown_event'] * cost_rates.get('hydro_shutdown_cost', 0)
    ramp_up_costs = result_data['ramp_up_mw'] * cost_rates.get('hydro_ramp_up_cost', 0)
    ramp_down_costs = result_data['ramp_down_mw'] * cost_rates.get('hydro_ramp_down_cost', 0)
    
    # Add cost columns
    result_data['startup_cost_$'] = startup_costs
    result_data['shutdown_cost_$'] = shutdown_costs
    result_data['ramp_up_cost_$'] = ramp_up_costs
    result_data['ramp_down_cost_$'] = ramp_down_costs
    result_data['total_hourly_om_cost_$'] = startup_costs + shutdown_costs + ramp_up_costs + ramp_down_costs
    
    return result_data


def generate_om_cost_summary(dispatch_data: pd.DataFrame) -> Dict[str, float]:
    """
    Generate annual O&M cost summary.
    
    Args:
        dispatch_data (pd.DataFrame): Hydro dispatch data with calculated costs
        
    Returns:
        Dict[str, float]: Summary of annual costs
    """
    summary = {
        'total_startup_events': dispatch_data['startup_event'].sum(),
        'total_shutdown_events': dispatch_data['shutdown_event'].sum(),
        'total_ramp_up_mw': dispatch_data['ramp_up_mw'].sum(),
        'total_ramp_down_mw': dispatch_data['ramp_down_mw'].sum(),
        'annual_startup_cost_$': dispatch_data['startup_cost_$'].sum(),
        'annual_shutdown_cost_$': dispatch_data['shutdown_cost_$'].sum(),
        'annual_ramp_up_cost_$': dispatch_data['ramp_up_cost_$'].sum(),
        'annual_ramp_down_cost_$': dispatch_data['ramp_down_cost_$'].sum(),
        'total_annual_om_cost_$': dispatch_data['total_hourly_om_cost_$'].sum()
    }
    
    return summary


def generate_per_unit_cost_summary(dispatch_data: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    """
    Generate per-unit O&M cost summary for each hydro generator.
    
    Args:
        dispatch_data (pd.DataFrame): Hydro dispatch data with calculated costs
        
    Returns:
        Dict[str, Dict[str, float]]: Per-unit summary of annual costs
    """
    per_unit_summary = {}
    
    # Filter to only hydro generator units (hydro_t*), exclude pumps (hydro_p*)
    generator_units = [unit for unit in dispatch_data['unit_id'].unique() if 'hydro_t' in unit or ('hydro_' in unit and 'hydro_p' not in unit)]
    
    for unit_id in generator_units:
        unit_data = dispatch_data[dispatch_data['unit_id'] == unit_id]
        
        per_unit_summary[unit_id] = {
            'startup_events': unit_data['startup_event'].sum(),
            'shutdown_events': unit_data['shutdown_event'].sum(),
            'ramp_up_mw': unit_data['ramp_up_mw'].sum(),
            'ramp_down_mw': unit_data['ramp_down_mw'].sum(),
            'startup_cost_$': unit_data['startup_cost_$'].sum(),
            'shutdown_cost_$': unit_data['shutdown_cost_$'].sum(),
            'ramp_up_cost_$': unit_data['ramp_up_cost_$'].sum(),
            'ramp_down_cost_$': unit_data['ramp_down_cost_$'].sum(),
            'total_om_cost_$': unit_data['total_hourly_om_cost_$'].sum()
        }
    
    return per_unit_summary


def display_om_cost_summary(summary: Dict[str, float], model_type: str):
    """
    Display the O&M cost summary.
    
    Args:
        summary (Dict[str, float]): Cost summary data
        model_type (str): Model type description
    """
    print("\n" + "="*60)
    print(f"ANNUAL HYDRO O&M COST SUMMARY - {model_type}")
    print("="*60)
    
    print(f"\nOPERATIONAL EVENTS:")
    print(f"  Total Startup Events: {summary['total_startup_events']:.0f}")
    print(f"  Total Shutdown Events: {summary['total_shutdown_events']:.0f}")
    print(f"  Total Ramp Up (MW): {summary['total_ramp_up_mw']:.2f}")
    print(f"  Total Ramp Down (MW): {summary['total_ramp_down_mw']:.2f}")
    
    print(f"\nANNUAL COSTS:")
    print(f"  Startup Costs: ${summary['annual_startup_cost_$']:,.2f}")
    print(f"  Shutdown Costs: ${summary['annual_shutdown_cost_$']:,.2f}")
    print(f"  Ramp Up Costs: ${summary['annual_ramp_up_cost_$']:,.2f}")
    print(f"  Ramp Down Costs: ${summary['annual_ramp_down_cost_$']:,.2f}")
    print(f"  TOTAL ANNUAL O&M COST: ${summary['total_annual_om_cost_$']:,.2f}")
    
    print("\n" + "="*60)


def display_per_unit_cost_summary(per_unit_summary: Dict[str, Dict[str, float]], model_type: str):
    """
    Display the per-unit O&M cost breakdown.
    
    Args:
        per_unit_summary (Dict[str, Dict[str, float]]): Per-unit cost summary data
        model_type (str): Model type description
    """
    print("\n" + "="*80)
    print(f"PER-UNIT HYDRO O&M COST BREAKDOWN - {model_type}")
    print("="*80)
    
    for unit_id, unit_costs in per_unit_summary.items():
        print(f"\n{unit_id.upper()}:")
        print("-" * 40)
        print(f"  Startup Events: {unit_costs['startup_events']:.0f}")
        print(f"  Shutdown Events: {unit_costs['shutdown_events']:.0f}")
        print(f"  Ramp Up (MW): {unit_costs['ramp_up_mw']:.2f}")
        print(f"  Ramp Down (MW): {unit_costs['ramp_down_mw']:.2f}")
        print(f"  Startup Costs: ${unit_costs['startup_cost_$']:,.2f}")
        print(f"  Shutdown Costs: ${unit_costs['shutdown_cost_$']:,.2f}")
        print(f"  Ramp Up Costs: ${unit_costs['ramp_up_cost_$']:,.2f}")
        print(f"  Ramp Down Costs: ${unit_costs['ramp_down_cost_$']:,.2f}")
        print(f"  UNIT TOTAL O&M COST: ${unit_costs['total_om_cost_$']:,.2f}")
    
    print("\n" + "="*80)


def import_hydro_data(filepath: str) -> pd.DataFrame:
    """
    Import hydro O&M cost data from cells A3:D6.
    
    Args:
        filepath (str): Path to the Excel file
        
    Returns:
        pd.DataFrame: Hydro O&M cost data
    """
    hydro_data = pd.read_excel(
        filepath,
        sheet_name="O&M Cost",
        skiprows=2,  # Skip rows 1 and 2 to start at A3
        nrows=4,     # Read 4 rows (A3:A6)
        usecols="A:D",  # Read columns A through D
        header=None  # Don't use first row as header
    )
    
    # Clean up the data - remove rows that are all NaN or contain "Battery Unit"
    hydro_data = hydro_data.dropna(how='all')
    hydro_data = hydro_data[~hydro_data.iloc[:, 0].astype(str).str.contains('Battery', na=False)]
    
    return hydro_data


def import_battery_data(filepath: str) -> pd.DataFrame:
    """
    Import battery O&M cost data from cells A8:D8.
    
    Args:
        filepath (str): Path to the Excel file
        
    Returns:
        pd.DataFrame: Battery O&M cost data
    """
    battery_data = pd.read_excel(
        filepath,
        sheet_name="O&M Cost",
        skiprows=7,  # Skip rows 1-7 to start at A8
        nrows=1,     # Read 1 row (A8 only)
        usecols="A:D",  # Read columns A through D
        header=None  # Don't use first row as header
    )
    
    # Clean up the data - remove rows that are all NaN
    battery_data = battery_data.dropna(how='all')
    
    return battery_data


def import_psh_data(filepath: str) -> pd.DataFrame:
    """
    Import Pumped Storage Hydro O&M cost data from cells A10:D13.
    
    Args:
        filepath (str): Path to the Excel file
        
    Returns:
        pd.DataFrame: PSH O&M cost data (empty if not present)
    """
    try:
        psh_data = pd.read_excel(
            filepath,
            sheet_name="O&M Cost",
            skiprows=9,  # Skip rows 1-9 to start at A10
            nrows=4,     # Read 4 rows (A10:A13)
            usecols="A:D",  # Read columns A through D
            header=None  # Don't use first row as header
        )
        
        # Clean up the data - remove rows that are all NaN
        psh_data = psh_data.dropna(how='all')
        
        # Check if this actually contains PSH data (not just empty cells)
        if psh_data.empty or psh_data.iloc[0, 0] != 'startup_cost':
            return pd.DataFrame()  # Return empty DataFrame if no PSH data
        
        return psh_data
        
    except Exception:
        return pd.DataFrame()  # Return empty DataFrame on error


def determine_model_type(filename: str, has_psh: bool) -> str:
    """
    Determine the model type based on filename and PSH presence.
    
    Args:
        filename (str): Input filename
        has_psh (bool): Whether PSH data is present
        
    Returns:
        str: Model type description
    """
    if "Model_A" in filename:
        return "Run-of-River"
    elif "Model_B" in filename:
        return "Reservoir"
    elif "Model_C" in filename or has_psh:
        return "Reservoir with Pumped Storage"
    else:
        # Determine based on PSH presence for custom files
        return "Reservoir with Pumped Storage" if has_psh else "Unknown"


def import_all_om_data(filepath: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, bool, str]:
    """
    Import all O&M cost data from Excel file.
    
    Args:
        filepath (str): Path to the Excel file
        
    Returns:
        Tuple: (hydro_data, battery_data, psh_data, has_PSH, model_type)
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Input file not found: {filepath}")
    
    # Import all data
    hydro_data = import_hydro_data(filepath)
    battery_data = import_battery_data(filepath)
    psh_data = import_psh_data(filepath)
    
    has_PSH = not psh_data.empty
    model_type = determine_model_type(filepath, has_PSH)
    
    print("Successfully imported all O&M cost data")
    
    return hydro_data, battery_data, psh_data, has_PSH, model_type


def display_data_summary(hydro_data: pd.DataFrame, battery_data: pd.DataFrame, 
                        psh_data: pd.DataFrame, has_PSH: bool, model_type: str):
    """
    Display a summary of the imported O&M cost data.
    
    Args:
        hydro_data (pd.DataFrame): Hydro O&M cost data
        battery_data (pd.DataFrame): Battery O&M cost data
        psh_data (pd.DataFrame): PSH O&M cost data
        has_PSH (bool): Flag indicating if PSH asset is present
        model_type (str): Model type description
    """
    print("\n" + "="*60)
    print(f"O&M COST DATA SUMMARY - {model_type}")
    print("="*60)
    
    print("\nHYDRO O&M DATA:")
    print("-" * 30)
    print(hydro_data.to_string(index=False))
    
    print("\nBATTERY O&M DATA:")
    print("-" * 30)
    print(battery_data.to_string(index=False))
    
    if has_PSH and not psh_data.empty:
        print("\nPUMPED STORAGE HYDRO O&M DATA:")
        print("-" * 40)
        print(psh_data.to_string(index=False))
    
    print("\n" + "="*60)


def determine_dispatch_filepath(input_filepath: str) -> str:
    """
    Automatically determine the dispatch file path based on input file.
    
    Args:
        input_filepath (str): Path to the input Excel file
        
    Returns:
        str: Path to the corresponding dispatch CSV file
    """
    # Extract model name from input file
    filename = os.path.basename(input_filepath)
    model_name = filename.replace('.xlsx', '')
    
    # Construct dispatch file path
    dispatch_path = f"Simulation_Results/{model_name}/ALEAF_HydroBoost_{model_name}__hydro_dispatch.csv"
    
    return dispatch_path


def parse_arguments():
    """
    Parse command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Calculate O&M costs from HydroBoost Excel input files and simulation results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python Generate_OM_cost.py Input_Spreadsheets/Model_A.xlsx
  python Generate_OM_cost.py Input_Spreadsheets/Model_C.xlsx
  python Generate_OM_cost.py Input_Spreadsheets/Model_A.xlsx --rates-only
        """
    )
    
    parser.add_argument(
        "filename",
        type=str,
        help="Path to Excel input file to import O&M cost data from"
    )
    
    parser.add_argument(
        "--rates-only",
        action="store_true",
        help="Only import O&M rates, skip dispatch calculation"
    )
    
    return parser.parse_args()


def main():
    """
    Main function for O&M cost calculation.
    """
    args = parse_arguments()
    
    print("HydroBoost O&M Cost Calculator")
    print("=" * 40)
    print(f"Importing O&M rates from: {args.filename}")
    
    try:
        # Import O&M cost parameters
        hydro_data, battery_data, psh_data, has_PSH, model_type = import_all_om_data(args.filename)
        display_data_summary(hydro_data, battery_data, psh_data, has_PSH, model_type)
        
        # Extract cost rates
        cost_rates = extract_om_cost_rates(hydro_data, battery_data, psh_data)
        
        # Skip dispatch calculation if --rates-only flag is used
        if args.rates_only:
            print("\nSkipping dispatch calculation (--rates-only flag used)")
            return
        
        # Automatically determine dispatch file path
        dispatch_path = determine_dispatch_filepath(args.filename)
        
        # Check if dispatch file exists
        if os.path.exists(dispatch_path):
            print(f"\nProcessing hydro dispatch data from: {dispatch_path}")
            
            # Import dispatch data
            dispatch_data = import_hydro_dispatch_data(dispatch_path)
            
            # Calculate startup/shutdown events
            dispatch_data = calculate_startup_shutdown_events(dispatch_data)
            
            # Calculate ramping costs
            dispatch_data = calculate_ramping_costs(dispatch_data)
            
            # Calculate hourly O&M costs
            dispatch_data = calculate_hourly_om_costs(dispatch_data, cost_rates)
            
            # Generate and display summaries
            summary = generate_om_cost_summary(dispatch_data)
            per_unit_summary = generate_per_unit_cost_summary(dispatch_data)
            
            display_om_cost_summary(summary, model_type)
            display_per_unit_cost_summary(per_unit_summary, model_type)
        
        else:
            print(f"\nDispatch file not found: {dispatch_path}")
            print("Only O&M cost rates imported. Run simulation to generate dispatch data.")
            
    except Exception as e:
        print(f"Error: {str(e)}")


if __name__ == "__main__":
    main()
