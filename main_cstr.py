from __future__ import division

# Imports for various utilities
import csv  # To work with csv files
import json  # To work with json files
import os  # To interact with the operating system
from math import ceil, fabs  # Importing ceil and fabs functions from math module

# Imports for visualization and graph operations
import matplotlib.pyplot as plt  # Used for creating static, animated, and interactive visualizations in Python
import networkx as nx  # Allows the creation, manipulation, and study of the structure, dynamics, and functions of complex networks.

# Imports for optimization
import pyomo.environ as pe  # Importing the Pyomo module
import logging  # Logging module to log messages
from pyomo.environ import SolverFactory, Suffix, value  # SolverFactory is a factory class for creating solver plugins. Suffix and value are for accessing solver related information.
from pyomo.gdp import Disjunct, Disjunction  # Imports for GDP model
from pyomo.util.infeasible import log_infeasible_constraints  # Logging infeasible constraints

# Imports for modules in the project
from gdp.cstr.gdp_reactor import build_cstrs  # Module to build the CSTRs
from gdp.dsda.dsda_functions import (  # Importing various functions from dsda_functions module
                                    external_ref,
                                    extvars_gdp_to_mip,
                                    generate_initialization,
                                    get_external_information,
                                    initialize_model,
                                    solve_complete_external_enumeration,
                                    solve_subproblem, solve_with_dsda,
                                    solve_with_gdpopt, solve_with_minlp,
                                    visualize_dsda)



def visualize_cstr_superstructure(m, NT):
    x = list((range(1, NT+1)))  # Creates a list of integers from 1 to NT+1

    # Initialize lists for bypasses (b) and reactors(r) and variable for recycle
    xb = []  # Stores the bypasses
    xr = []  # Stores the reactors
    recycle = 0  # Indicates if recycling is in use

    yp = {}  # Dictionary to store reactor (or not) information
    yr = {}  # Dictionary to store recycle (or not) information

    # Populates yp and yr with binary indicator variable values from the model
    for n in m.N:
        yp[n] = pe.value(pe.value(m.YP[n].get_associated_binary()))
        yr[n] = pe.value(pe.value(m.YR[n].get_associated_binary()))

    # Classifies units as bypasses (b) or reactors(r) and determines recycle status
    for i in x:
        if yp[i] > 0.5:  # If the unit is a reactor
            xr.append(i)
        else:  # If the unit is a bypass
            xb.append(i)
        if yr[i] > 0.5:  # If the unit is a recycle
            recycle = i

    # Creates labels for bypasses (b), reactors(r), input/output (f) and recycle(recy)
    blabels = dict(zip(range(1, len(xb)+1), xb[::-1]))
    rlabels = dict(zip(range(len(xb)+1, len(xb)+1+len(xr)), xr[::-1]))
    flabels = {0: '', NT+1: ''}
    recylabels = {'r1': '', 'r2': '', 'r3': '', 'r4': ''}

    # Creates positions (pos) for bypasses (b), reactors(r), input/output(f) and recycle(recy)
    posb = {i: (i, 0) for i in range(1, len(xb)+1)}
    posr = {i: (i, 0) for i in range(len(xb)+1, len(xb)+1+len(xr))}
    posf = {0: (0.2, 0), NT+1: (NT+1, 0)}
    posrecy = {'r1': (NT+0.5, -0.0009), 'r2': (NT+0.5, 0.008), 'r3': (NT-recycle+0.5, 0.007), 'r4': (NT-recycle+0.5, -0.0009)}

    # Creates flow arrow from input to output
    arcsf = [(0, NT+1)]

    # Declare a directed graph
    graph = nx.DiGraph()

    # Draws labels, nodes and edges for input/output
    nx.draw_networkx_labels(graph, posf, flabels)
    nx.draw_networkx_edges(graph, posf, arcsf, width=1, arrowsize=10)
    nx.draw_networkx(graph, posf, node_size=1, node_color='black', nodelist=flabels, with_labels=True, node_shape='', edgecolors='black')

    # Draws labels and nodes for bypasses
    nx.draw_networkx_labels(graph, posb, blabels)
    nx.draw_networkx(graph, posb, node_size=900, node_color='whitesmoke', width=1.5, nodelist=blabels, with_labels=True, node_shape='s', edgecolors='black', linewidths=0.2)

    # Draws labels and nodes for reactors
    nx.draw_networkx_labels(graph, posr, rlabels)
    nx.draw_networkx(graph, posr, node_size=900, node_color='lightslategray', width=1.5, nodelist=rlabels, with_labels=True, node_shape='s', edgecolors='black', linewidths=1.5)

    # Draws labels, nodes and edges for recycle if it exists
    if recycle != 0:
        arcsrecy = [('r1', 'r2'), ('r3', 'r4')]
        pairs = list(zip(list(arcsrecy), ['R', 'R']))
        edgelabels = dict(pairs)
        nx.draw_networkx_labels(graph, posrecy, recylabels)
        nx.draw_networkx_edges(graph, posrecy, arcsrecy, width=1, arrowsize=10)
        nx.draw_networkx(graph, posrecy, node_size=0, node_color='white', nodelist=recylabels, node_shape='', edgecolors='black')
        nx.draw_networkx_edge_labels(graph, posrecy, edge_labels=edgelabels)

    plt.show()  # Shows the graph



def problem_logic_cstr(m):
    # Initialize the list of logical expressions
    logic_expr = []
    
    # Loop through each unit operation
    for n in m.N:
        
        # If recycling is active at a given stage, then the indicator variable 
        # for the recycle disjunct must be active
        logic_expr.append([m.YR[n], m.YR_is_recycle[n].indicator_var])
        
        # If recycling is not active at a given stage, then the indicator variable 
        # for the "no recycle" disjunct must be active
        logic_expr.append([~m.YR[n], m.YR_is_not_recycle[n].indicator_var])
        
        # The reactor can only be a CSTR if there is no unreacted feed in the 
        # previous stages or if there is unreacted feed in the current stage
        logic_expr.append([pe.lor(pe.land(~m.YF[n2] for n2 in range(1, n)), m.YF[n]), m.YP_is_cstr[n].indicator_var])
        
        # If the reactor is not a CSTR, then it must be a bypass
        logic_expr.append([~pe.lor(pe.land(~m.YF[n2] for n2 in range(1, n)), m.YF[n]), m.YP_is_bypass[n].indicator_var])
        
        # The reactor can only operate if there is no unreacted feed in the previous stages 
        # or if there is unreacted feed in the current stage
        logic_expr.append([pe.lor(pe.land(~m.YF[n2] for n2 in range(1, n)), m.YF[n]), m.YP[n]])
    
    # Return the list of logical expressions
    return logic_expr



if __name__ == "__main__":

    # Number of stages to try for the model
    NTs = range(5, 26, 1)

    # Maximum time allowed for solving each instance
    timelimit = 900

    # Starting point for D-SDA
    starting_point = [1, 1]

    # Enable verbose output
    globaltee = True

    # Setup logging to only show errors
    logging.basicConfig(level=logging.ERROR)

    # Setup columns for csv output
    csv_columns = ['Method', 'Approach', 'Solver', 'Objective', 'Time', 'Status', 'User_time', 'NT']
    dict_data = []

    # Setup file paths
    dir_path = os.path.dirname(os.path.abspath(__file__))
    csv_file = os.path.join(dir_path, "results", "cstr_results.csv")

    # List of NLP solvers to try
    nlps = ['knitro', 'baron']

    # Setup options for each NLP solver
    nlp_opts = dict((nlp, {}) for nlp in nlps)

    # List of MINLP solvers to try
    minlps = ['antigone', 'baron', 'scip', 'dicopt', 'sbb', 'knitro']

    # Setup options for each MINLP solver
    minlps_opts = dict((minlp, {}) for minlp in minlps)
    minlps_opts['dicopt']['add_options'] = [
        'GAMS_MODEL.optfile = 1;'
        '\n'
        '$onecho > dicopt.opt \n'
        'stop 0 \n'
        '*relaxed 2 \n'
        'maxcycles 10000 \n'
        'nlpsolver knitro \n'
        '$offecho \n'
    ]
    minlps_opts['sbb']['add_options'] = [
        'GAMS_MODEL.optfile = 1;'
        '\n'
        '$onecho > sbb.opt \n'
        'rootsolver knitro \n'
        'subsolver knitro \n'
        '$offecho \n'
    ]

    # List of transformations to try
    transformations = ['bigm', 'hull']

    # List of k-values to try for D-SDA
    ks = ['Infinity', '2']

    # List of strategies to try for GDPopt
    strategies = ['LOA', 'GLOA', 'LBB']

    # Loop over each number of stages
    for NT in NTs:
        # Check if initialization exists for this number of stages
        json_file = os.path.join(dir_path, 'gdp/dsda/', 'cstr_' + str(NT) + '_initialization.json')
        if os.path.exists(json_file):
            init_path = json_file
        else:
            # If not, create an initial solution
            m = build_cstrs(NT)
            ext_ref = {m.YF: m.N, m.YR: m.N}
            reformulation_dict, number_of_external_variables, lower_bounds, upper_bounds = get_external_information(m, ext_ref, tee=globaltee)
            m_fixed = external_ref(m=m, x=[1, 1], extra_logic_function=problem_logic_cstr, dict_extvar=reformulation_dict, tee=globaltee)
            m_solved = solve_subproblem(m=m_fixed, subproblem_solver='baron', timelimit=100, tee=True)
            init_path = generate_initialization(m=m_solved, starting_initialization=True, model_name='cstr_'+str(NT))


        # # MINLP
        # for solver in minlps:
        #     for transformation in transformations:
        #         new_result = {}
        #         m = build_cstrs(NT)
        #         m_init = initialize_model(m, json_path=init_path)
        #         m_solved = solve_with_minlp(
        #             m=m_init,
        #             transformation=transformation,
        #             minlp=solver,
        #             minlp_options=minlps_opts[solver],
        #             timelimit=timelimit,
        #             gams_output=False,
        #             tee=globaltee,
        #         )
        #         new_result = {'Method': 'MINLP', 'Approach': transformation, 'Solver': solver, 'Objective': pe.value(
        #             m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition, 'User_time': 'NA', 'NT': NT}
        #         dict_data.append(new_result)
        #         print(new_result)

        # # GDPopt
        # for solver in nlps:
        #     for strategy in strategies:
        #         new_result = {}
        #         m = build_cstrs(NT)
        #         m_init = initialize_model(m, json_path=init_path)
        #         m_solved = solve_with_gdpopt(
        #             m=m_init,
        #             mip='cplex',
        #             nlp=solver,
        #             nlp_options=nlp_opts[solver],
        #             timelimit=timelimit,
        #             strategy=strategy,
        #             tee=globaltee,
        #         )
        #         new_result = {'Method': 'GDPopt', 'Approach': strategy, 'Solver': solver, 'Objective': pe.value(
        #             m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition, 'User_time': 'NA', 'NT': NT}
        #         dict_data.append(new_result)
        #         print(new_result)

         # Setup for D-SDA algorithm
        m = build_cstrs(NT)
        ext_ref = {m.YF: m.N, m.YR: m.N}
        get_external_information(m, ext_ref, tee=False)

        # Loop over each NLP solver
        for solver in nlps:
            # Loop over each k-value
            for k in ks:
                # Loop over each transformation
                for transformation in ['hull','bigm']:
                    new_result = {}
                    # Run the D-SDA algorithm
                    m_solved, _, _ = solve_with_dsda(
                        model_function=build_cstrs,
                        model_args={'NT': NT},
                        starting_point=starting_point,
                        ext_dict=ext_ref,
                        ext_logic=problem_logic_cstr,
                        mip_transformation=True,
                        transformation=transformation,
                        k=k,
                        provide_starting_initialization=True,
                        feasible_model='cstr_' + str(NT),
                        subproblem_solver=solver,
                        subproblem_solver_options=nlp_opts[solver],
                        iter_timelimit=timelimit,
                        timelimit=timelimit,
                        gams_output=False,
                        tee=False,
                        global_tee=False,
                    )
                    # Store the result
                    new_result = {'Method': str('D-SDA_MIP_'+transformation), 'Approach': str('k='+k), 'Solver': solver, 'Objective': pe.value(
                        m_solved.obj), 'Time': m_solved.dsda_time, 'Status': m_solved.dsda_status, 'User_time': m_solved.dsda_usertime, 'NT': NT}
                    dict_data.append(new_result)
                    print(new_result)

        # Try to write the results to a csv file
        try:
            with open(csv_file, 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for data in dict_data:
                    writer.writerow(data)
        except IOError:
            print("I/O error")

    # # Complete enumeration
    # for transformation in ['bigm']:
    #     for solver in ['knitro']:
    #         NT = 25
    #         m = build_cstrs(NT)
    #         ext_ref = {m.YF: m.N, m.YR: m.N}
    #         reformulation_dict, _, _, _ = get_external_information(
    #             m, ext_ref, tee=True)
    #         mip_m, mip_extvars = extvars_gdp_to_mip(m, reformulation_dict)
    #         m_solved = solve_complete_external_enumeration(
    #             model_function=build_cstrs,
    #             model_args={'NT': NT},
    #             ext_dict=ext_ref,
    #             ext_logic=problem_logic_cstr,
    #             mip_transformation=True,
    #             transformation=transformation,
    #             feasible_model='cstr_'+str(NT),
    #             subproblem_solver=solver,
    #             subproblem_solver_options=nlp_opts[solver],
    #             iter_timelimit=900,
    #             gams_output=False,
    #             points=[(j, i) for j in [24,25] for i in range(26)],
    #             tee=False,
    #             global_tee=globaltee,
    #             export_csv=True,
    #         )
