import csv
from math import ceil, fabs

import pyomo.environ as pe
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints

from gdp.dsda.dsda_functions import (generate_initialization, initialize_model,
                                     solve_subproblem, solve_with_dsda,
                                     solve_with_gdpopt, solve_with_minlp,
                                     visualize_dsda,external_ref,get_external_information)
from gdp.small_batch.gdp_small_batch import build_small_batch

if __name__ == "__main__":
    # Inputs
    timelimit = 3600
    model_args = {}

    # D-SDA
    m = build_small_batch()
    Ext_Ref = {m.Y: m.k}

    get_external_information(m,Ext_Ref,tee=True)

    def problem_logic_batch(m): 
        logic_expr = []
        for k in m.k:
            for j in m.j:
                logic_expr.append([m.Y[k, j], m.Y_exists[k, j].indicator_var])
                logic_expr.append([~m.Y[k, j], m.Y_not_exists[k, j].indicator_var])
        return logic_expr

    ks = ['Infinity']
    starting_point = [3, 3, 3]

    nlps=['conopt']
    for solver in nlps:
        for k in ks:
            new_result = {}
            m_solved, route = solve_with_dsda(model_function=build_small_batch, model_args=model_args, starting_point=starting_point, reformulation_function=external_ref, ext_dict=Ext_Ref, ext_logic=problem_logic_batch, k=k,
                                              provide_starting_initialization=True, feasible_model='small_batch', subproblem_solver=solver, iter_timelimit=timelimit, timelimit=timelimit)
            new_result = {'Method': 'D-SDA', 'Approach': str('k = '+k), 'Solver': solver, 'Objective': pe.value(
                m_solved.obj), 'Time': m_solved.dsda_time, 'Status': m_solved.dsda_status}
    print(new_result)


    # csv_small_batchs = ['Method', 'Approach',
    #                     'Solver', 'Objective', 'Time', 'Status']
    # dict_data = []
    # csv_file = "smallbatch_results.csv"

    # # MINLPS
    # minlps = ['antigone', 'scip', 'baron',
    #           'sbb', 'dicopt', 'alphaecp', 'bonminh']
    # transformations = ['bigm', 'hull']

    # for solver in minlps:
    #     for transformation in transformations:
    #         new_result = {}
    #         m = build_small_batch(**model_args)
    #         m_init = initialize_model(
    #             m, from_feasible=True, feasible_model='small_batch')
    #         m_solved = solve_with_minlp(
    #             m_init, transformation=transformation, minlp=solver, timelimit=timelimit, gams_output=True)
    #         new_result = {'Method': 'MINLP', 'Approach': transformation, 'Solver': solver, 'Objective': pe.value(
    #             m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition}
    #         dict_data.append(new_result)

    # # GDPopt
    # nlps = ['msnlp', 'conopt']
    # strategies = ['LOA', 'GLOA']

    # for solver in nlps:
    #     for strategy in strategies:
    #         new_result = {}
    #         m = build_small_batch(**model_args)
    #         m_init = initialize_model(
    #             m, from_feasible=True, feasible_model='small_batch')
    #         m_solved = solve_with_gdpopt(
    #             m_init, mip='cplex', nlp=solver, timelimit=timelimit, strategy=strategy, nlp_output=True)
    #         new_result = {'Method': 'GDPopt', 'Approach': strategy, 'Solver': solver, 'Objective': pe.value(
    #             m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition}
    #         dict_data.append(new_result)

    # # GDPopt LBB
    # minlps = ['baron', 'scip']
    # strategies = ['LBB']

    # for solver in minlps:
    #     for strategy in strategies:
    #         new_result = {}
    #         m = build_small_batch(**model_args)
    #         m_init = initialize_model(
    #             m, from_feasible=True, feasible_model='small_batch')
    #         m_solved = solve_with_gdpopt(
    #             m_init, mip='cplex', minlp=solver, timelimit=timelimit, strategy=strategy, minlp_output=True)
    #         new_result = {'Method': 'GDPopt', 'Approach': strategy, 'Solver': solver, 'Objective': pe.value(
    #             m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition}
    #         dict_data.append(new_result)

    # # D-SDA
    # ks = ['Infinity', '2']
    # starting_point = [3, 3, 3]
    # min_allowed = {i: 1 for i in range(1, len(starting_point)+1)}
    # max_allowed = {i: 3 for i in range(1, len(starting_point)+1)}

    # for solver in nlps:
    #     for k in ks:
    #         new_result = {}
    #         m_solved, route = solve_with_dsda(model_function=build_small_batch, model_args=model_args, starting_point=starting_point, reformulation_function=external_ref, k=k,
    #                                           provide_starting_initialization=True, feasible_model='cstr', subproblem_solver=solver, min_allowed=min_allowed, max_allowed=max_allowed, iter_timelimit=timelimit, timelimit=timelimit)
    #         new_result = {'Method': 'D-SDA', 'Approach': str('k = '+k), 'Solver': solver, 'Objective': pe.value(
    #             m_solved.obj), 'Time': m_solved.dsda_time, 'Status': m_solved.dsda_status}
    #         dict_data.append(new_result)

    # print(dict_data)

    # try:
    #     with open(csv_file, 'w') as csvfile:
    #         writer = csv.DictWriter(csvfile, fieldnames=csv_small_batchs)
    #         writer.writeheader()
    #         for data in dict_data:
    #             writer.writerow(data)
    # except IOError:
    #     print("I/O error")
