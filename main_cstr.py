from __future__ import division

import csv
import json
import os
from math import ceil, fabs

import matplotlib.pyplot as plt
import networkx as nx
import pyomo.environ as pe
from pyomo.environ import SolverFactory, Suffix, value
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints

from gdp.cstr.gdp_reactor import build_cstrs
from gdp.dsda.dsda_functions import (external_ref, generate_initialization,
                                     get_external_information,
                                     initialize_model, solve_subproblem,
                                     solve_with_dsda, solve_with_gdpopt,
                                     solve_with_minlp, visualize_dsda, solve_complete_external_enumeration)


def visualize_cstr_superstructure(m, NT):
    x = list((range(1, NT+1)))

    # Initialize bypasses (b), reactors(r) and recycle
    xb = []
    xr = []
    recycle = 0

    yp = {}
    yr = {}

    # Use information from solved model
    for n in m.N:
        yp[n] = pe.value(pe.value(m.YP[n].get_associated_binary()))
        yr[n] = pe.value(pe.value(m.YR[n].get_associated_binary()))

    # Classify units in bypasses (b) or reactors(r) and determine recycle
    for i in x:
        if yp[i] > 0.5:
            xr.append(i)
        else:
            xb.append(i)
        if yr[i] > 0.5:
            recycle = i

    # Create labels for bypasses (b), reactors(r), input/output (f) and recycle(recy)
    blabels = dict(zip(range(1, len(xb)+1), xb[::-1]))
    rlabels = dict(zip(range(len(xb)+1, len(xb)+1+len(xr)), xr[::-1]))
    flabels = {0: '', NT+1: ''}
    recylabels = {'r1': '', 'r2': '', 'r3': '', 'r4': ''}

    # Create posicions (pos) for bypasses (b), reactors(r), input/output(f) and recycle(recy)
    posb = {}
    posr = {}
    posf = {0: (0.2, 0), NT+1: (NT+1, 0)}
    posrecy = {'r1': (NT+0.5, -0.0009), 'r2': (NT+0.5, 0.008),
               'r3': (NT-recycle+0.5, 0.007), 'r4': (NT-recycle+0.5, -0.0009)}

    for i in range(1, len(xb)+1):
        posb[i] = (i, 0)

    for i in range(len(xb)+1, len(xb)+1+len(xr)):
        posr[i] = (i, 0)

    # Create flow arrow from input to output
    arcsf = [(0, NT+1)]

    # Declare graph
    graph = nx.DiGraph()

    # Graph input/out(f)
    nx.draw_networkx_labels(graph, posf, flabels)
    nx.draw_networkx_edges(graph, posf, arcsf, width=1, arrowsize=10)
    nx.draw_networkx(graph, posf, node_size=1, node_color='black',
                     nodelist=flabels, with_labels=True, node_shape='', edgecolors='black')

    # Graph bypasses(b)
    nx.draw_networkx_labels(graph, posb, blabels)
    nx.draw_networkx(graph, posb, node_size=900, node_color='whitesmoke', width=1.5,
                     nodelist=blabels, with_labels=True, node_shape='s', edgecolors='black', linewidths=0.2)

    # Graph reactors(r)
    nx.draw_networkx_labels(graph, posr, rlabels)
    nx.draw_networkx(graph, posr, node_size=900, node_color='lightslategray', width=1.5,
                     nodelist=rlabels, with_labels=True, node_shape='s', edgecolors='black', linewidths=1.5)

    # Graph recycle(recy) if it exists
    if recycle != 0:
        arcsrecy = [('r1', 'r2'), ('r3', 'r4')]
        pairs = list(zip(list(arcsrecy), ['R', 'R']))
        edgelabels = dict(pairs)
        nx.draw_networkx_labels(graph, posrecy, recylabels)
        nx.draw_networkx_edges(
            graph, posrecy, arcsrecy, width=1, arrowsize=10)
        nx.draw_networkx(graph, posrecy, node_size=0, node_color='white',
                         nodelist=recylabels, node_shape='', edgecolors='black')
        nx.draw_networkx_edge_labels(
            graph, posrecy, edge_labels=edgelabels)

    plt.show()


def problem_logic_cstr(m):
    logic_expr = []
    for n in m.N:
        logic_expr.append([m.YR[n], m.YR_is_recycle[n].indicator_var])
        logic_expr.append([~m.YR[n], m.YR_is_not_recycle[n].indicator_var])
        logic_expr.append([pe.lor(pe.land(~m.YF[n2] for n2 in range(
            1, n)), m.YF[n]), m.YP_is_cstr[n].indicator_var])
        logic_expr.append([~pe.lor(pe.land(~m.YF[n2] for n2 in range(
            1, n)), m.YF[n]), m.YP_is_bypass[n].indicator_var])
        logic_expr.append([pe.lor(pe.land(~m.YF[n2] for n2 in range(
            1, n)), m.YF[n]), m.YP[n]])
    return logic_expr


if __name__ == "__main__":
    
    # Results
    NTs = range(5, 26, 1)
    # NTs = [10]
    timelimit = 900
    starting_point = [1, 1]

    globaltee = True

    csv_columns = ['Method', 'Approach', 'Solver',
                   'Objective', 'Time', 'Status', 'User_time', 'NT']
    dict_data = []
    csv_file = "cstr_results.csv"

    nlps = ['msnlp', 'knitro', 'baron']

    nlp_opts = dict((nlp, {}) for nlp in nlps)
    nlp_opts['msnlp']['add_options'] = [
        'GAMS_MODEL.optfile = 1;'
        '\n'
        '$onecho > msnlp.opt \n'
        'nlpsolver knitro \n'
        '$offecho \n'
    ]

    minlps = ['antigone', 'baron', 'scip', 'dicopt', 'sbb']

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
    transformations = ['bigm', 'hull']
    ks = ['Infinity', '2']
    strategies = ['LOA', 'GLOA', 'LBB']

    for NT in NTs:
        # Create initialization for all methods starting with a single reactor

        dir_path = os.path.dirname(os.path.abspath(__file__))
        json_file = os.path.join(
            dir_path, 'gdp/dsda/', 'cstr_' + str(NT) + '_initialization.json')
        if os.path.exists(json_file):
            init_path = json_file
        else:
            m = build_cstrs(NT)
            ext_ref = {m.YF: m.N, m.YR: m.N}
            reformulation_dict, number_of_external_variables, lower_bounds, upper_bounds = get_external_information(
                m, ext_ref, tee=globaltee)
            m_fixed = external_ref(m=m, x=[
                1, 1], other_function=problem_logic_cstr, dict_extvar=reformulation_dict, tee=True)
            m_solved = solve_subproblem(
                m=m_fixed, subproblem_solver='baron', timelimit=100, tee=True)
            init_path = generate_initialization(
                m=m_solved, starting_initialization=True, model_name='cstr_'+str(NT))

        # MINLP
        for solver in minlps:
            for transformation in transformations:
                new_result = {}
                m = build_cstrs(NT)
                m_init = initialize_model(m, json_path=init_path)
                m_solved = solve_with_minlp(
                    m=m_init,
                    transformation=transformation,
                    minlp=solver,
                    minlp_options=minlps_opts[solver],
                    timelimit=timelimit,
                    gams_output=False,
                    tee=globaltee,
                )
                new_result = {'Method': 'MINLP', 'Approach': transformation, 'Solver': solver, 'Objective': pe.value(
                    m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition, 'User_time': 'NA', 'NT': NT}
                dict_data.append(new_result)
                print(new_result)

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
                new_result = {'Method': 'GDPopt', 'Approach': strategy, 'Solver': solver, 'Objective': pe.value(
                    m_solved.obj), 'Time': m_solved.results.solver.user_time, 'Status': m_solved.results.solver.termination_condition, 'User_time': 'NA', 'NT': NT}
                dict_data.append(new_result)
                print(new_result)

        # D-SDA
        m = build_cstrs(NT)
        ext_ref = {m.YF: m.N, m.YR: m.N}
        get_external_information(m, ext_ref, tee=globaltee)

        for solver in nlps:
            for k in ks:
                new_result = {}
                m_solved, _, _ = solve_with_dsda(
                    model_function=build_cstrs,
                    model_args={'NT': NT},
                    starting_point=starting_point,
                    ext_dict=ext_ref,
                    ext_logic=problem_logic_cstr,
                    k=k,
                    provide_starting_initialization=True,
                    feasible_model='cstr_' + str(NT),
                    subproblem_solver=solver,
                    subproblem_solver_options=nlp_opts[solver],
                    iter_timelimit=timelimit,
                    timelimit=timelimit,
                    gams_output=False,
                    tee=globaltee,
                    global_tee=globaltee,
                )
                new_result = {'Method': 'D-SDA', 'Approach': str('k='+k), 'Solver': solver, 'Objective': pe.value(
                    m_solved.obj), 'Time': m_solved.dsda_time, 'Status': m_solved.dsda_status, 'User_time': m_solved.dsda_usertime, 'NT': NT}
                dict_data.append(new_result)
                print(new_result)

        try:
            with open(csv_file, 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for data in dict_data:
                    writer.writerow(data)
        except IOError:
            print("I/O error")

    # Complete enumeration
    NT = 25
    m = build_cstrs(NT)
    ext_ref = {m.YF: m.N, m.YR: m.N}
    get_external_information(m, ext_ref, tee=True)

    solve_complete_external_enumeration(build_cstrs, 
                                        model_args={'NT': NT}, 
                                        ext_dict=ext_ref, 
                                        ext_logic=problem_logic_cstr, 
                                        feasible_model='cstr_' + str(NT),
                                        subproblem_solver='knitro',
                                        iter_timelimit=30,
                                        tee=False,
                                        global_tee=True,
                                        export_csv=True)


