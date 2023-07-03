# Importing the necessary libraries for the script
# csv - to work with csv files
# logging - to log errors and messages 
# os - provides functions for interacting with the operating system
# math.ceil, math.fabs - provides ceiling function and absolute value function
# time - provides time-related functions
import csv
import logging
import os
from math import ceil, fabs
import time

# Importing Pyomo - a Python-based open-source optimization modeling language
# pyomo.environ - provides a Pythonic interface to define optimization models
# pyomo.gdp - allows to formulate and solve Generalized Disjunctive Programming (GDP) models
# pyomo.util.infeasible - provides utility to log infeasible constraints
import pyomo.environ as pe
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints

# Importing functions from the gdp.dsda.dsda_functions module which seems to be custom functions for solving optimization problems
# external_ref - function to reference external variables
# generate_initialization - function to generate an initial guess for optimization variables
# get_external_information - function to get external information
# initialize_model - function to initialize the optimization model
# solve_complete_external_enumeration - function to solve complete external enumeration
# solve_subproblem - function to solve subproblems
# solve_with_dsda - function to solve with dual sequential decision algorithm (dsda)
# solve_with_gdpopt - function to solve with the GDPopt solver
# solve_with_minlp - function to solve with mixed integer nonlinear programming (MINLP) solver
# visualize_dsda - function to visualize the results of the dsda
from gdp.dsda.dsda_functions import (external_ref, generate_initialization,
                                     get_external_information,
                                     initialize_model,
                                     solve_complete_external_enumeration,
                                     solve_subproblem, solve_with_dsda,
                                     solve_with_gdpopt, solve_with_minlp,
                                     visualize_dsda)

# Importing the function build_small_batch from the module gdp.small_batch.gdp_small_batch
# This function seems to construct the optimization model for your problem
from gdp.small_batch.gdp_small_batch import build_small_batch



def problem_logic_batch(m):
    """
    This function provides the logic of the batch problem.
    [Add reference]

    Parameters
    ----------
    m: Pyomo model
        Original optimization model
    
    Returns
    -------
    logic_expr: List[Tuple[Binary variable, Indicator variable]]
        List of assigned indicator variables with corresponding logical expression over binary variables.
    """
    logic_expr = []  # Initialize an empty list for logical expressions
    for k in m.k:  # Iterate over potential number of parallel units
        for j in m.j:  # Iterate over stages
            # If the disjunction 'k' is selected for stage 'j', append its corresponding binary variable and indicator variable to the logic_expr list
            logic_expr.append([m.Y[k, j], m.Y_exists[k, j].indicator_var])
            # If the disjunction 'k' is not selected for stage 'j', append its corresponding binary variable and indicator variable to the logic_expr list
            logic_expr.append([~m.Y[k, j], m.Y_not_exists[k, j].indicator_var])
    return logic_expr

if __name__ == "__main__":
    # Setting the time limit and the starting point for the model
    timelimit = 900
    starting_point = [3, 3, 3]

    # Enabling global logging to capture all print statements
    globaltee = True
    logging.basicConfig(level=logging.ERROR)

    # Defining the columns of the CSV file for storing results
    csv_columns = ['Method', 'Approach', 'Solver', 'Objective', 'Time', 'Status', 'User_time']
    dict_data = []  # List to store results before writing to CSV
    dir_path = os.path.dirname(os.path.abspath(__file__))  # Getting the directory of the current file
    csv_file = os.path.join(dir_path, "results", "small_batch_results.csv")  # Defining the path of the results file

    # Defining the list of Nonlinear Programming (NLP) and Mixed Integer Nonlinear Programming (MINLP) solvers to be used
    nlps = ['knitro', 'baron']
    minlps = ['antigone', 'baron', 'scip', 'dicopt', 'sbb', 'knitro']

    # Defining options for the NLP and MINLP solvers
    nlp_opts = dict((nlp, {}) for nlp in nlps)
    minlps_opts = dict((minlp, {}) for minlp in minlps)

    # Defining some additional options for certain solvers
    minlps_opts['dicopt']['add_options'] = [ ... ]
    minlps_opts['sbb']['add_options'] = [ ... ]

    # Defining the transformations, external formulations, and strategies for the solvers
    transformations = ['bigm', 'hull']
    ks = ['Infinity', '2']
    strategies = ['LOA', 'LBB']

   # Uncomment this block to initialize the model and solve it using an initialization file
# json_file = os.path.join(
#     dir_path, "gdp/dsda/", "small_batch_initialization.json")
# Check if the initialization file exists
# if os.path.exists(json_file):
#     init_path = json_file  # If the file exists, use it to initialize the model
# else:  # If the file doesn't exist, create a new model and solve it to generate the initialization
#     m = build_small_batch()  # Building the initial model
#     ext_ref = {m.Y: m.k}  # Mapping the binary variables to their external references
#     reformulation_dict, number_of_external_variables, lower_bounds, upper_bounds = get_external_information(
#         m, ext_ref, tee=globaltee)  # Getting the external information
#     m_fixed = external_ref(m=m, x=starting_point, extra_logic_function=problem_logic_batch,
#                            dict_extvar=reformulation_dict, tee=globaltee)  # Fixing the external variables
#     m_solved = solve_subproblem(
#         m=m_fixed, subproblem_solver='baron', timelimit=100, tee=globaltee)  # Solving the fixed model
#     init_path = generate_initialization(
#         m=m_solved, starting_initialization=True, model_name='small_batch')  # Generating the initialization

# Uncomment this block to solve the model using different MINLP solvers and transformations
# for solver in minlps:
#     for transformation in transformations:
#         new_result = {}  # Create an empty dictionary to store the results
#         m = build_small_batch()  # Building the initial model
#         m_init = initialize_model(m, json_path=init_path)  # Initialize the model
#         m_solved = solve_with_minlp(
#             m_init,
#             transformation=transformation,
#             minlp=solver,
#             minlp_options=minlps_opts[solver],
#             timelimit=timelimit,
#             gams_output=False,
#             tee=globaltee,
#         )  # Solve the model
#         # Store the results
#         new_result = {'Method': 'MINLP', 'Approach': transformation, 'Solver': solver, 'Objective': pe.value(
#             m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition, 'User_time': 'NA'}
#         dict_data.append(new_result)  # Add the results to the list
#         print(new_result)  # Print the results

# Uncomment this block to solve the model using the GDPopt method with different NLP solvers and strategies
# for solver in nlps:
#     for strategy in strategies:
#         new_result = {}  # Create an empty dictionary to store the results
#         m = build_small_batch()  # Building the initial model
#         m_init = initialize_model(m, json_path=init_path)  # Initialize the model
#         m_solved = solve_with_gdpopt(
#             m_init,
#             mip='cplex',
#             nlp=solver,
#             nlp_options=nlp_opts[solver],
#             timelimit=timelimit,
#             strategy=strategy,
#             tee=globaltee,
#         )  # Solve the model
#         # Store the results
#         new_result = {'Method': 'GDPopt', 'Approach': strategy, 'Solver': solver, 'Objective': pe.value(
#             m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition, 'User_time': 'NA'}
#         dict_data.append(new_result)  # Add the results to the list
#         print(new_result)  # Print the results


    # D-SDA - MINLP: This is the main block of code that solves the model
    m = build_small_batch()  # Building the initial model
    ext_ref = {m.Y: m.k}  # Mapping the binary variables to their external references
    get_external_information(m, ext_ref, tee=globaltee)  # Getting the external information

    # Iterating over the NLP solvers, external formulations, and transformations
    for solver in nlps:
        for k in ks:
            for transformation in transformations:
                # Solving the model with the D-SDA method
                m_solved, _, _ = solve_with_dsda(
                    model_function=build_small_batch,
                    model_args={},
                    starting_point=starting_point,
                    ext_dict=ext_ref,
                    mip_transformation=True,
                    transformation=transformation,
                    ext_logic=problem_logic_batch,
                    k=k,
                    provide_starting_initialization=True,
                    feasible_model='small_batch',
                    subproblem_solver=solver,
                    subproblem_solver_options=nlp_opts[solver],
                    iter_timelimit=timelimit,
                    timelimit=timelimit,
                    gams_output=False,
                    tee=globaltee,
                    global_tee=globaltee,
                )

                # Storing the results
                new_result = {
                    'Method': str('D-SDA_MIP_' + transformation),
                    'Approach': str('k=' + k),
                    'Solver': solver,
                    'Objective': pe.value(m_solved.obj),
                    'Time': m_solved.dsda_time,
                    'Status': m_solved.dsda_status,
                    'User_time': m_solved.dsda_usertime
                }
                dict_data.append(new_result)
                print(new_result)

    # Writing the results to a CSV file
    try:
        with open(csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in dict_data:
                writer.writerow(data)
    except IOError:
        print("I/O error")

    # Uncomment this block to solve the model using complete enumeration, for each transformation and solver
# for transformation in transformations:  # Iterate over each transformation
#     for solver in ['knitro', 'baron']:  # Iterate over each solver
#         m = build_small_batch()  # Building the initial model
#         ext_ref = {m.Y: m.k}  # Mapping the binary variables to their external references
#         get_external_information(m, ext_ref, tee=False)  # Get external information of the model
#         # Solve the model using complete external enumeration
#         m_solved = solve_complete_external_enumeration(
#             model_function=build_small_batch,  # The function to build the model
#             model_args={},  # Arguments to the model function
#             ext_dict=ext_ref,  # The dictionary mapping the binary variables to their external references
#             ext_logic=problem_logic_batch,  # The logic for the external references
#             feasible_model='small_batch',  # The name of the feasible model
#             mip_transformation=True,  # Whether to transform the problem to a MIP problem
#             transformation=transformation,  # The transformation to use
#             subproblem_solver=solver,  # The solver to use
#             subproblem_solver_options=nlp_opts[solver],  # The options for the solver
#             iter_timelimit=900,  # The time limit for each iteration
#             gams_output=False,  # Whether to output GAMS code
#             tee=globaltee,  # Whether to display solver output
#             global_tee=globaltee,  # Whether to display solver output globally
#             export_csv=True,  # Whether to export the results to a CSV file
#         )

