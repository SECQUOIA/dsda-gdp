"""Distillation column model for 2018 PSE conference"""

from __future__ import division

import csv
from math import ceil, fabs

import pyomo.environ as pe
from pyomo.environ import (Block, BooleanVar, ConcreteModel, Constraint,
                           NonNegativeReals, Objective, Param, RangeSet, Set,
                           SolverFactory, Suffix, TransformationFactory, Var,
                           exactly, land, log, lor, minimize, value)
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints

from gdp.column.gdp_column import build_column
from gdp.dsda.dsda_functions import (generate_initialization, initialize_model,
                                     solve_subproblem, solve_with_dsda,
                                     solve_with_gdpopt, solve_with_minlp,
                                     visualize_dsda,external_ref,get_external_information)





def complete_enumeration_external(model_function=build_column, model_args={'min_trays': 8, 'max_trays': 17, 'xD': 0.95, 'xB': 0.95}, reformulation_function=external_ref, subproblem_solver='knitro', timelimit=10):
    NT = model_args['max_trays']
    X1, X2, aux, aux2, x = [], [], [], 2, {}

    for i in range(2, NT):
        X1.append(i)
        aux.append(i)
        X2.append(aux2)

    for i in range(NT-2):
        aux.pop(0)
        aux2 += 1
        for j in aux:
            X1.append(j)
            X2.append(aux2)

    print()
    feas_x, feas_y, objs = [], [], []

    print('=============================')
    print('%6s %6s %12s' % ('x1', 'x2', 'Objective'))
    print('-----------------------------')

    # Loop over all external variables and then loop over its values
    for i in range(len(X1)):
        x = [X1[i], X2[i]]
        m = model_function(**model_args)
        m_init = initialize_model(
            m, from_feasible=True, feasible_model='column')
        m_fixed = reformulation_function(m_init, x)
        m_solved = solve_subproblem(
            m_fixed, subproblem_solver=subproblem_solver, timelimit=timelimit)

        if m_solved.dsda_status == 'Optimal':
            print('%6s %6s %12s' %
                  (X1[i], X2[i], round(pe.value(m_solved.obj), 2)))
            feas_x.append(X1[i])
            feas_y.append(X2[i])
            objs.append(round(pe.value(m_solved.obj), 2))
        else:
            print('%6s %6s %12s' % (X1[i], X2[i], 'Infeasible'))

    print('=============================')
    return feas_x, feas_y, objs

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
    
    timelimit = 60
    model_args = {'min_trays': 8, 'max_trays': NT, 'xD': 0.95, 'xB': 0.95}
    starting_point = [15, 1]

    csv_columns = ['Method','Approach','Solver','Objective', 'Time', 'Status', 'User time']
    dict_data = []
    csv_file = "column_results.csv"

    nlps = ['baron', 'msnlp', 'knitro', 'ipopth']
    minlps = ['antigone', 'scip', 'baron']
    transformations = ['bigm','hull']
    ks = ['Infinity','2']
    strategies = ['LOA','LBB']

    #MINLP
    for solver in minlps:
        for transformation in transformations:
            new_result = {}
            m = build_column(**model_args)
            m_init = initialize_model(m, from_feasible=True, feasible_model='column')
            m_solved = solve_with_minlp(m_init, transformation=transformation, minlp=solver, timelimit=timelimit, gams_output=False, tee=False)
            new_result = {'Method':'MINLP', 'Approach':transformation, 'Solver':solver, 'Objective':pe.value(m_solved.obj), 'Time':m_solved.results.solver.user_time, 'Status':m_solved.results.solver.termination_condition, 'User time':'NA'}
            dict_data.append(new_result)
            print(new_result)

    #GDPopt
    for solver in nlps:
        for strategy in strategies:
            new_result = {}
            m = build_column(**model_args)
            m_init = initialize_model(m, from_feasible=True, feasible_model='column')
            m_solved = solve_with_gdpopt(m_init, mip='cplex', nlp=solver, timelimit=timelimit, strategy=strategy, tee=False)
            new_result = {'Method':'GDPopt','Approach':strategy, 'Solver':solver, 'Objective':pe.value(m_solved.obj), 'Time':m_solved.results.solver.user_time, 'Status':m_solved.results.solver.termination_condition, 'User time':'NA'}
            dict_data.append(new_result)
            print(new_result)



    # # D-SDA
    m = build_column(**model_args)
    Ext_Ref = {m.YB: m.intTrays, m.YR: m.intTrays}
    get_external_information(m,Ext_Ref,tee=False)

    for solver in nlps:
        for k in ks:
            new_result = {}
            m_solved, _ = solve_with_dsda(model_function=build_column, model_args=model_args, starting_point=starting_point, ext_dict=Ext_Ref, ext_logic=problem_logic_column,
                                    k=k, provide_starting_initialization=True, feasible_model='column', subproblem_solver=solver, iter_timelimit=timelimit, timelimit=timelimit, gams_output=False, tee=False, global_tee=False)
            new_result = {'Method':'D-SDA', 'Approach':str('k = '+k), 'Solver':solver,'Objective':pe.value(m_solved.obj), 'Time':m_solved.dsda_time, 'Status':m_solved.dsda_status, 'User time':m_solved.dsda_usertime, 'NT':NT}
            dict_data.append(new_result)
            print(new_result)

    # try:
    #     with open(csv_file, 'w') as csvfile:
    #         writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
    #         writer.writeheader()
    #         for data in dict_data:
    #             writer.writerow(data)
    # except IOError:
    #     print("I/O error")

    



    # k = 'Infinity'
    # starting_point = [15, 1]

    # m_solved, route = solve_with_dsda(model_function=build_column, model_args=model_args, starting_point=starting_point, ext_dict=Ext_Ref, ext_logic=problem_logic_column,
    #                                   k=k, provide_starting_initialization=True, feasible_model='column', subproblem_solver='knitro', 
    #                                   iter_timelimit=10, timelimit=60, gams_output=False, tee=False, global_tee=True)
    # visualize_dsda(route=route, feas_x=x, feas_y=y, objs=objs, k=k, ext1_name='YR (Reflux position)', ext2_name='YB (Boil-up position)')
    # TODO This visualization code does not work
    # print(m_solved.results)

    # Results

#     NT = 17
#     timelimit = 3600
#     model_args = {'min_trays': 8, 'max_trays': NT, 'xD': 0.95, 'xB': 0.95}

#     csv_columns = ['Method', 'Approach',
#                    'Solver', 'Objective', 'Time', 'Status']
#     dict_data = []
#     csv_file = "column_results.csv"

#     # MINLPS
#     minlps = ['antigone', 'scip', 'baron',
#               'sbb', 'dicopt', 'alphaecp', 'bonminh']
#     transformations = ['bigm', 'hull']

#     # m = build_column(**model_args)
#     # m_fixed = external_ref(m, [16, 2])
#     # m_solved = solve_subproblem(m_fixed)
#     # m_solved.YR.pprint()
#     # generate_initialization(m_solved, True, 'column')

#     # for solver in minlps:
#     #     for transformation in transformations:
#     #         new_result = {}
#     #         m = build_column(**model_args)
#     #         m_init = initialize_model(m, from_feasible=True, feasible_model='column')
#     #         m_solved = solve_with_minlp(m_init, transformation=transformation, minlp=solver, timelimit=timelimit, gams_output=True)
#     #         new_result = {'Method':'MINLP', 'Approach':transformation, 'Solver':solver, 'Objective':pe.value(m_solved.obj), 'Time':m_solved.results.solver.user_time, 'Status':m_solved.results.solver.termination_condition}
#     #         dict_data.append(new_result)

#     # GDPopt
#     # nlps = ['msnlp', 'knitro']
#     # strategies = ['LOA','GLOA']

#     # for solver in nlps:
#     #     for strategy in strategies:
#     #         new_result = {}
#     #         m = build_column(**model_args)
#     #         m_init = initialize_model(m, from_feasible=True, feasible_model='column')
#     #         m_solved = solve_with_gdpopt(m_init, mip='cplex', nlp=solver, timelimit=timelimit, strategy=strategy, nlp_output=True)
#     #         new_result = {'Method':'GDPopt','Approach':strategy, 'Solver':solver, 'Objective':pe.value(m_solved.obj), 'Time':m_solved.results.solver.user_time, 'Status':m_solved.results.solver.termination_condition}
#     #         dict_data.append(new_result)

#     # GDPopt LBB
#     #minlps = ['baron', 'scip']
#     #strategies = ['LBB']

#    # for solver in minlps:
#    #     for strategy in strategies:
#    #         new_result = {}
#    #         m = build_column(**model_args)
#    #         m_init = initialize_model(m, from_feasible=True, feasible_model='column')
#    #         m_solved = solve_with_gdpopt(m_init, mip='cplex', minlp=solver, timelimit=timelimit, strategy=strategy, minlp_output=True)
#    #         new_result = {'Method':'GDPopt','Approach':strategy, 'Solver':solver, 'Objective':pe.value(m_solved.obj), 'Time':m_solved.results.solver.user_time, 'Status':m_solved.results.solver.termination_condition}
#    #         dict_data.append(new_result)

#     # D-SDA
#     # ks = ['Infinity','2']
#     # starting_point = [16, 2]
#     # min_allowed = {i: 2 for i in range(1, len(starting_point)+1)}
#     # max_allowed = {i: NT-1 for i in range(1, len(starting_point)+1)}

#     # for solver in nlps:
#     #    for k in ks:
#     #        new_result = {}
#     #        m_solved, route = solve_with_dsda(model_function=build_column, model_args=model_args, starting_point=starting_point, reformulation_function=external_ref, k=k,
#     #                    provide_starting_initialization=True, feasible_model='cstr', subproblem_solver=solver, min_allowed=min_allowed, max_allowed=max_allowed, iter_timelimit=timelimit, timelimit=timelimit)
#     #        new_result = {'Method':'D-SDA', 'Approach':str('k = '+k), 'Solver':solver,'Objective':pe.value(m_solved.obj), 'Time':m_solved.dsda_time, 'Status':m_solved.dsda_status}
#     #        dict_data.append(new_result)

#     # print(dict_data)

#   #  try:
#   #      with open(csv_file, 'w') as csvfile:
#   #          writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
#   #          writer.writeheader()
#   #          for data in dict_data:
#   #              writer.writerow(data)
#   #  except IOError:
#   #      print("I/O error")
