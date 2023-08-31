"""
main_cstr.py

The code imports the CSTR system from the `gdp.cstr.gdp_reactor` module and solves it using different methods (MINLP, GDPopt, DSDA).
The code ...
# TODO: complete description

References:
[1] Linan, David A., et al. "Optimal design of superstructures for placing units and streams with multiple and ordered available locations. Part I: A new mathematical framework." Computers & Chemical Engineering 137, (2020): 106794.
"""
from __future__ import division

import csv
import json
import os
from math import ceil, fabs

import matplotlib.pyplot as plt
import networkx as nx
import pyomo.environ as pyo
import logging
from pyomo.environ import SolverFactory, Suffix, value
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints

from gdp.cstr.gdp_reactor import build_cstrs
from gdp.dsda.dsda_functions import (
    external_ref,
    extvars_gdp_to_mip,
    generate_initialization,
    get_external_information,
    initialize_model,
    solve_complete_external_enumeration,
    solve_subproblem,
    solve_with_dsda,
    solve_with_gdpopt,
    solve_with_minlp,
    visualize_dsda,
)


def visualize_cstr_superstructure(m, NT):
    """
    Visualize the Continuous Stirred Tank Reactor (CSTR) superstructure.

    This function constructs and displays a visual representation of the CSTR system.
    The stages are depicted as nodes, where bypasses are shown in 'whitesmoke' color,
    and reactors are in 'lightslategray'. If recycling exists in the system, it is
    depicted with edges labeled 'R'.

    Args:
        m (pyomo.ConcreteModel) : The Pyomo model instance containing all model variables and parameters
        related to the CSTR system. Particularly, the function expects `m.YP`, `m.YR` variables denoting bypasses and reactors.

        NT (int) : The total number of stages or reactors in the CSTR system.

    Returns:
        None

    Note:
        The function uses the NetworkX library for the visualization, and the
        representation is displayed using `plt.show()`.
    """
    x = list((range(1, NT + 1)))

    # Initialize bypasses (b), reactors(r) and recycle
    xb = []
    xr = []
    recycle = 0

    yp = {}
    yr = {}

    # Use information from solved model
    for n in m.N:
        yp[n] = pyo.value(pyo.value(m.YP[n].get_associated_binary()))
        yr[n] = pyo.value(pyo.value(m.YR[n].get_associated_binary()))

    # Classify units in bypasses (b) or reactors(r) and determine recycle
    for i in x:
        if yp[i] > 0.5:
            xr.append(i)
        else:
            xb.append(i)
        if yr[i] > 0.5:
            recycle = i

    # Create labels for bypasses (b), reactors(r), input/output (f) and recycle(recy)
    blabels = dict(zip(range(1, len(xb) + 1), xb[::-1]))
    rlabels = dict(zip(range(len(xb) + 1, len(xb) + 1 + len(xr)), xr[::-1]))
    flabels = {0: '', NT + 1: ''}
    recylabels = {'r1': '', 'r2': '', 'r3': '', 'r4': ''}

    # Create posicions (pos) for bypasses (b), reactors(r), input/output(f) and recycle(recy)
    posb = {}
    posr = {}
    posf = {0: (0.2, 0), NT + 1: (NT + 1, 0)}
    posrecy = {
        'r1': (NT + 0.5, -0.0009),
        'r2': (NT + 0.5, 0.008),
        'r3': (NT - recycle + 0.5, 0.007),
        'r4': (NT - recycle + 0.5, -0.0009),
    }

    for i in range(1, len(xb) + 1):
        posb[i] = (i, 0)

    for i in range(len(xb) + 1, len(xb) + 1 + len(xr)):
        posr[i] = (i, 0)

    # Create flow arrow from input to output
    arcsf = [(0, NT + 1)]

    # Declare graph
    graph = nx.DiGraph()

    # Graph input/out(f)
    nx.draw_networkx_labels(graph, posf, flabels)
    nx.draw_networkx_edges(graph, posf, arcsf, width=1, arrowsize=10)
    nx.draw_networkx(
        graph,
        posf,
        node_size=1,
        node_color='black',
        nodelist=flabels,
        with_labels=True,
        node_shape='',
        edgecolors='black',
    )

    # Graph bypasses(b)
    nx.draw_networkx_labels(graph, posb, blabels)
    nx.draw_networkx(
        graph,
        posb,
        node_size=900,
        node_color='whitesmoke',
        width=1.5,
        nodelist=blabels,
        with_labels=True,
        node_shape='s',
        edgecolors='black',
        linewidths=0.2,
    )

    # Graph reactors(r)
    nx.draw_networkx_labels(graph, posr, rlabels)
    nx.draw_networkx(
        graph,
        posr,
        node_size=900,
        node_color='lightslategray',
        width=1.5,
        nodelist=rlabels,
        with_labels=True,
        node_shape='s',
        edgecolors='black',
        linewidths=1.5,
    )

    # Graph recycle(recy) if it exists
    if recycle != 0:
        arcsrecy = [('r1', 'r2'), ('r3', 'r4')]
        pairs = list(zip(list(arcsrecy), ['R', 'R']))
        edgelabels = dict(pairs)
        nx.draw_networkx_labels(graph, posrecy, recylabels)
        nx.draw_networkx_edges(graph, posrecy, arcsrecy, width=1, arrowsize=10)
        nx.draw_networkx(
            graph,
            posrecy,
            node_size=0,
            node_color='white',
            nodelist=recylabels,
            node_shape='',
            edgecolors='black',
        )
        nx.draw_networkx_edge_labels(graph, posrecy, edge_labels=edgelabels)

    plt.show()


def problem_logic_cstr(m):
    """
    Generate a set of logical expressions based on model parameters.

    Args:
        m (pyomo.ConcreteModel) : The Pyomo model instance containing all model variables and parameters
        related to the CSTR system. Particularly, the function expects `m.YP`, `m.YR` variables denoting bypasses and reactors.

    Returns:
        logic_expr (list) : A list of logical expressions that are used to define the external variables of the CSTR system.
    """
    logic_expr = []
    for n in m.N:
        logic_expr.append([m.YR[n], m.YR_is_recycle[n].indicator_var])
        logic_expr.append([~m.YR[n], m.YR_is_not_recycle[n].indicator_var])
        logic_expr.append(
            [
                pyo.lor(pyo.land(~m.YF[n2] for n2 in range(1, n)), m.YF[n]),
                m.YP_is_cstr[n].indicator_var,
            ]
        )
        logic_expr.append(
            [
                ~pyo.lor(pyo.land(~m.YF[n2] for n2 in range(1, n)), m.YF[n]),
                m.YP_is_bypass[n].indicator_var,
            ]
        )
        logic_expr.append(
            [pyo.lor(pyo.land(~m.YF[n2] for n2 in range(1, n)), m.YF[n]), m.YP[n]]
        )
    return logic_expr


if __name__ == "__main__":
    # Results
    # NTs = range(5, 26, 1)
    NTs = [5]
    timelimit = 900
    starting_point = [1, 1]

    globaltee = True
    # Setting logging level to ERROR to avoid printing FBBT warning of some constraints not implemented
    logging.basicConfig(level=logging.ERROR)

    csv_columns = [
        'Method',
        'Approach',
        'Solver',
        'Objective',
        'Time',
        'Status',
        'User_time',
        'NT',
    ]
    dict_data = []

    dir_path = os.path.dirname(os.path.abspath(__file__))
    csv_file = os.path.join(dir_path, "results", "cstr_results.csv")

    # nlps = ['knitro', 'baron']# , 'msnlp']
    nlps = ['knitro']

    nlp_opts = dict((nlp, {}) for nlp in nlps)
    # nlp_opts['msnlp']['add_options'] = [
    #     'GAMS_MODEL.optfile = 1;'
    #     '\n'
    #     '$onecho > msnlp.opt \n'
    #     'nlpsolver knitro \n'
    #     '$offecho \n'
    # ]

    minlps = ['antigone', 'baron', 'scip', 'dicopt', 'sbb', 'knitro']

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
    # transformations = ['bigm', 'hull']
    transformations = ['bigm']
    ks = ['Infinity', '2']
    # strategies = ['LOA', 'GLOA', 'LBB']
    strategies = ['LOA']

    for NT in NTs:
        # Create initialization for all methods starting with a single reactor
        json_file = os.path.join(
            dir_path, 'gdp/dsda/', 'cstr_' + str(NT) + '_initialization.json'
        )
        if os.path.exists(json_file):
            init_path = json_file
        else:
            m = build_cstrs(NT)
            ext_ref = {m.YF: m.N, m.YR: m.N}
            (
                reformulation_dict,
                number_of_external_variables,
                lower_bounds,
                upper_bounds,
            ) = get_external_information(m, ext_ref, tee=globaltee)
            m_fixed = external_ref(
                m=m,
                x=[1, 1],
                extra_logic_function=problem_logic_cstr,
                dict_extvar=reformulation_dict,
                tee=globaltee,
            )
            m_solved = solve_subproblem(
                m=m_fixed, subproblem_solver='baron', timelimit=100, tee=True
            )
            init_path = generate_initialization(
                m=m_solved, starting_initialization=True, model_name='cstr_' + str(NT)
            )

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
        #         new_result = {'Method': 'MINLP', 'Approach': transformation, 'Solver': solver, 'Objective': pyo.value(
        #             m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition, 'User_time': 'NA', 'NT': NT}
        #         dict_data.append(new_result)
        #         print(new_result)

        # GDPopt
        for solver in nlps:
            for strategy in strategies:
                new_result = {}
                m = build_cstrs(NT)
                m_init = initialize_model(m, json_path=init_path)
                m_solved = solve_with_gdpopt(
                    m=m_init,
                    mip='cplex',
                    nlp=solver,
                    nlp_options=nlp_opts[solver],
                    timelimit=timelimit,
                    strategy=strategy,
                    tee=globaltee,
                )
                new_result = {
                    'Method': 'GDPopt',
                    'Approach': strategy,
                    'Solver': solver,
                    'Objective': pyo.value(m_solved.obj),
                    'Time': m_solved.results.solver.user_time,
                    'Status': m_solved.results.solver.termination_condition,
                    'User_time': 'NA',
                    'NT': NT,
                }
                dict_data.append(new_result)
                print(new_result)

        visualize_cstr_superstructure(m, NT)

        # # D-SDA
        # m = build_cstrs(NT)
        # ext_ref = {m.YF: m.N, m.YR: m.N}
        # get_external_information(m, ext_ref, tee=False)

        # for solver in nlps:
        #     for k in ks:
        #         for transformation in transformations:
        #             new_result = {}
        #             m_solved, _, _ = solve_with_dsda(
        #                 model_function=build_cstrs,
        #                 model_args={'NT': NT},
        #                 starting_point=starting_point,
        #                 ext_dict=ext_ref,
        #                 ext_logic=problem_logic_cstr,
        #                 mip_transformation=True,
        #                 transformation=transformation,
        #                 k=k,
        #                 provide_starting_initialization=True,
        #                 feasible_model='cstr_' + str(NT),
        #                 subproblem_solver=solver,
        #                 subproblem_solver_options=nlp_opts[solver],
        #                 iter_timelimit=timelimit,
        #                 timelimit=timelimit,
        #                 gams_output=False,
        #                 tee=False,
        #                 global_tee=False,
        #             )
        #             new_result = {'Method': str('D-SDA_MIP_'+transformation), 'Approach': str('k='+k), 'Solver': solver, 'Objective': pyo.value(
        #                 m_solved.obj), 'Time': m_solved.dsda_time, 'Status': m_solved.dsda_status, 'User_time': m_solved.dsda_usertime, 'NT': NT}
        #             dict_data.append(new_result)
        #             print(new_result)

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
