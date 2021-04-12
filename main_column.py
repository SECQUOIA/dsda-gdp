"""Distillation column model for 2018 PSE conference"""

from __future__ import division

import csv
import logging
import os
from math import ceil, fabs

import pyomo.environ as pe
from pyomo.environ import (Block, BooleanVar, ConcreteModel, Constraint,
                           NonNegativeReals, Objective, Param, RangeSet, Set,
                           SolverFactory, Suffix, TransformationFactory, Var,
                           exactly, land, log, lor, minimize, value)
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints

from gdp.column.gdp_column import build_column
from gdp.dsda.dsda_functions import (external_ref, generate_initialization,
                                     get_external_information,
                                     initialize_model,
                                     solve_complete_external_enumeration,
                                     solve_subproblem, solve_with_dsda,
                                     solve_with_gdpopt, solve_with_minlp,
                                     visualize_dsda)


def problem_logic_column(m):
    logic_expr = []
    for n in m.intTrays:
        logic_expr.append([pe.land(~m.YR[n] for n in range(
            m.reboil_tray+1, m.feed_tray)), m.YR_is_down])
        logic_expr.append([pe.land(~m.YB[n]
                                   for n in range(m.feed_tray+1, m.max_trays)), m.YB_is_up])
    for n in m.conditional_trays:
        logic_expr.append([pe.land(pe.lor(m.YR[j] for j in range(n, m.max_trays)), pe.lor(
            pe.land(~m.YB[j] for j in range(n, m.max_trays)), m.YB[n])), m.tray[n].indicator_var])
        logic_expr.append([~pe.land(pe.lor(m.YR[j] for j in range(n, m.max_trays)), pe.lor(
            pe.land(~m.YB[j] for j in range(n, m.max_trays)), m.YB[n])), m.no_tray[n].indicator_var])
    return logic_expr


if __name__ == "__main__":
    # Inputs
    NT = 17
    timelimit = 900
    model_args = {'min_trays': 8, 'max_trays': NT, 'xD': 0.95, 'xB': 0.95}
    # Initializing at column with all trays, reboil in bottom tray and reflux in top-most tray
    starting_point = [NT-2, 1]
    globaltee = True
    # Setting logging level to ERROR to avoid printing FBBT warning of some constraints not implemented
    logging.basicConfig(level=logging.ERROR)

    csv_columns = ['Method', 'Approach', 'Solver',
                   'Objective', 'Time', 'Status', 'User_time']
    dict_data = []
    dir_path = os.path.dirname(os.path.abspath(__file__))
    csv_file = os.path.join(
        dir_path, "results", "column_results.csv")

    nlps = ['msnlp', 'knitro', 'baron']

    nlp_opts = dict((nlp, {}) for nlp in nlps)
    nlp_opts['msnlp']['add_options'] = [
        'GAMS_MODEL.optfile = 1;'
        '\n'
        '$onecho > msnlp.opt \n'
        'nlpsolver knitro \n'
        '$offecho \n'
    ]

    minlps = ['antigone', 'baron', 'scip', 'dicopt', 'sbb', 'knitro']

    minlps_opts = dict((minlp, {}) for minlp in minlps)
    minlps_opts['dicopt']['add_options'] = [
        'GAMS_MODEL.optfile = 1;'
        '\n'
        '$onecho > dicopt.opt \n'
        'stop 0 \n'
        'relaxed 2 \n'
        'maxcycles 10000 \n'
        'nlpsolver knitro \n'
        '$offecho \n'
    ]

    # NOTE DICOPT with Hull reformulation will fail in reporting the right results. See lst file (it tends to be the initialization point)

    minlps_opts['sbb']['add_options'] = [
        'GAMS_MODEL.optfile = 1;'
        '\n'
        '$onecho > sbb.opt \n'
        'rootsolver knitro \n'
        'subsolver knitro \n'
        '$offecho \n'
    ]
    transformations = ['bigm', 'hull']
    ks = ['Infinity', '2']
    strategies = ['LOA', 'GLOA', 'LBB']

    # Initializations
    json_file = os.path.join(
        dir_path, 'gdp/dsda/', 'column_' + str(NT) + '_initialization.json')
    if os.path.exists(json_file):
        init_path = json_file
    else:
        m = build_column(**model_args)
        ext_ref = {m.YB: m.intTrays, m.YR: m.intTrays}
        reformulation_dict, number_of_external_variables, lower_bounds, upper_bounds = get_external_information(
            m, ext_ref, tee=globaltee)
        m_fixed = external_ref(m=m, x=starting_point, extra_logic_function=problem_logic_column,
                               dict_extvar=reformulation_dict, tee=globaltee)
        m_solved = solve_subproblem(
            m=m_fixed, subproblem_solver='baron', timelimit=100, tee=globaltee)
        init_path = generate_initialization(
            m=m_solved, starting_initialization=True, model_name='column_'+str(NT))

    # MINLP
    # for solver in minlps:
    #     for transformation in transformations:
    #         new_result = {}
    #         m = build_column(**model_args)
    #         m_init = initialize_model(m, json_path=init_path)
    #         m_solved = solve_with_minlp(
    #             m_init,
    #             transformation=transformation,
    #             minlp=solver,
    #             minlp_options=minlps_opts[solver],
    #             timelimit=timelimit,
    #             gams_output=False,
    #             tee=globaltee,
    #         )
    #         new_result = {'Method': 'MINLP', 'Approach': transformation, 'Solver': solver, 'Objective': pe.value(
    #             m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition, 'User_time': 'NA'}
    #         dict_data.append(new_result)
    #         print(new_result)

    # # GDPopt
    # for solver in nlps:
    #     for strategy in strategies:
    #         new_result = {}
    #         m = build_column(**model_args)
    #         m_init = initialize_model(m, json_path=init_path)
    #         m_solved = solve_with_gdpopt(
    #             m_init,
    #             mip='cplex',
    #             nlp=solver,
    #             nlp_options=nlp_opts[solver],
    #             timelimit=timelimit,
    #             strategy=strategy,
    #             tee=globaltee,
    #         )
    #         new_result = {'Method': 'GDPopt', 'Approach': strategy, 'Solver': solver, 'Objective': pe.value(
    #             m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition, 'User_time': 'NA'}
    #         dict_data.append(new_result)
    #         print(new_result)

    # # D-SDA
    # m = build_column(**model_args)
    # ext_ref = {m.YB: m.intTrays, m.YR: m.intTrays}
    # get_external_information(m, ext_ref, tee=globaltee)

    # for solver in nlps:
    #     for k in ks:
    #         new_result = {}
    #         m_solved, _, _ = solve_with_dsda(
    #             model_function=build_column,
    #             model_args=model_args,
    #             starting_point=starting_point,
    #             ext_dict=ext_ref,
    #             ext_logic=problem_logic_column,
    #             k=k,
    #             provide_starting_initialization=True,
    #             feasible_model='column_' + str(NT),
    #             subproblem_solver=solver,
    #             subproblem_solver_options=nlp_opts[solver],
    #             iter_timelimit=timelimit,
    #             timelimit=timelimit,
    #             gams_output=False,
    #             tee=False,
    #             global_tee=globaltee,
    #         )
    #         new_result = {'Method': 'D-SDA', 'Approach': str('k='+k), 'Solver': solver, 'Objective': pe.value(
    #             m_solved.obj), 'Time': m_solved.dsda_time, 'Status': m_solved.dsda_status, 'User_time': m_solved.dsda_usertime}
    #         dict_data.append(new_result)
    #         print(new_result)

    # try:
    #     with open(csv_file, 'w') as csvfile:
    #         writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
    #         writer.writeheader()
    #         for data in dict_data:
    #             writer.writerow(data)
    # except IOError:
    #     print("I/O error")

    m = build_column(**model_args)
    ext_ref = {m.YB: m.intTrays, m.YR: m.intTrays}
    get_external_information(m, ext_ref, tee=False)
    iterlim = 3600
    # points = [(13, 4)]
    # points = [(14, 7), (15, 7), (15, 8), (15, 9), (7, 1),
    #           (8, 1), (9, 1), (9, 2), (10, 3), ]
    # points = [(12, 4), (12, 5), (13, 1), (13, 2), (13, 3),
    #           (13, 5), (13, 6), (13, 7), (14, 1), (14, 2),
    #           (14, 3), (14, 4), (14, 5), (14, 6), (14, 7),
    #           (15, 1), (15, 2), (15, 3), (15, 4), (15, 5),
    #           (15, 6), (15, 7), (15, 8), (15, 9), (7, 1),
    #           (8, 1), (9, 1), (9, 2), (10, 3), ]
    
    for transformation in ['hull']:
        for solver in ['knitro']:
            m_solved = solve_complete_external_enumeration(
                model_function=build_column,
                model_args=model_args,
                ext_dict=ext_ref,
                ext_logic=problem_logic_column,
                feasible_model='column_'+str(NT)+'_optimal',
                # points=points,
                subproblem_solver=solver,
                subproblem_solver_options=nlp_opts[solver],
                iter_timelimit=iterlim,
                mip_transformation=True,
                transformation=transformation,
                gams_output=False,
                tee=globaltee,
                global_tee=globaltee,
                export_csv=True,
            )
