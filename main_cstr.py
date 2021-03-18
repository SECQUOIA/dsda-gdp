from __future__ import division

from math import ceil, fabs

import matplotlib.pyplot as plt
import networkx as nx
import pyomo.environ as pe
from pyomo.environ import SolverFactory, Suffix, value
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints

from gdp.dsda.dsda_functions import (generate_initialization, initialize_model,
                                     solve_subproblem, solve_with_dsda, solve_with_gdpopt,
                                     solve_with_minlp, visualize_dsda)
from gdp.cstr.gdp_reactor import build_cstrs


def external_ref(m, x, logic_expr=None):
    # External variable fix
    ext_var_1 = x[0]
    ext_var_2 = x[1]
    for n in m.N:
        if n == ext_var_1:
            m.YF[n].fix(True)
        else:
            m.YF[n].fix(False)

        if n == ext_var_2:
            m.YR_is_recycle[n].indicator_var.fix(True)
            m.YR_is_not_recycle[n].indicator_var.fix(False)
        else:
            m.YR_is_recycle[n].indicator_var.fix(False)
            m.YR_is_not_recycle[n].indicator_var.fix(True)

        temp = pe.value(pe.lor(pe.land(~m.YF[n2]
                                       for n2 in range(1, n)), m.YF[n]))

        if temp == True:
            m.YP_is_cstr[n].indicator_var.fix(True)
            m.YP_is_bypass[n].indicator_var.fix(False)
        else:
            m.YP_is_cstr[n].indicator_var.fix(False)
            m.YP_is_bypass[n].indicator_var.fix(True)

    pe.TransformationFactory('core.logical_to_linear').apply_to(m)
    pe.TransformationFactory('gdp.fix_disjuncts').apply_to(m)
    pe.TransformationFactory('contrib.deactivate_trivial_constraints').apply_to(
        m, tmp=False, ignore_infeasible=True)

    return m


def complete_enumeration_external(model_function=build_cstrs, model_args={'NT': 5}, reformulation_function=external_ref, subproblem_solver='msnlp', timelimit=10):
    X1 = list(range(1, NT+1))
    # TODO how to generalize for N external variables?
    X2 = list(range(1, NT+1))
    # Input variable should be dictionary of the external variables with lower and upper bounds
    print()
    feas_x, feas_y, objs = [], [], []

    print('=============================')
    print('%6s %6s %12s' % ('x1', 'x2', 'Objective'))
    print('-----------------------------')

    # Loop over all external variables and then loop over its values
    for Xi in X1:
        for Xj in X2:
            m = model_function(**model_args)
            x = [Xi, Xj]
            m_init = initialize_model(
                m, from_feasible=True,  feasible_model='cstr')
            m_fixed = reformulation_function(m_init, x)
            m_solved = solve_subproblem(
                m_fixed, subproblem_solver=subproblem_solver, timelimit=timelimit)

            if m_solved.dsda_status == 'Optimal':
                print('%6s %6s %12s' %
                      (Xi, Xj, round(pe.value(m_solved.obj), 5)))
                feas_x.append(Xi)
                feas_y.append(Xj)
                objs.append(round(pe.value(m_solved.obj), 5))
            else:
                print('%6s %6s %12s' % (Xi, Xj, 'Infeasible'))

    print('=============================')
    return feas_x, feas_y, objs


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


if __name__ == "__main__":
    # Inputs
    NT = 5
    timelimit = 10

    # # Complete enumeration
    # x, y, objs = complete_enumeration_external(
    #     model_function=build_cstrs, model_args={'NT': NT}, subproblem_solver='msnlp', timelimit=10)

    # # MINLP solution
    m = build_cstrs(NT)
    m_init = initialize_model(m, from_feasible=True, feasible_model='cstr')
    m_solved = solve_with_minlp(
        m_init, transformation='bigm', minlp='baron', timelimit=timelimit, gams_output=False)
    print(m_solved.results)
    # visualize_cstr_superstructure(m_solved, NT)

    # GDPopt method
    m = build_cstrs(NT)
    m_init = initialize_model(m, from_feasible=True, feasible_model='cstr')
    m_solved = solve_with_gdpopt(m_init, mip='cplex', nlp='conopt',
                                 timelimit=timelimit, strategy='LOA', mip_output=False, nlp_output=False)
    print(m_solved.results)
    #visualize_cstr_superstructure(m_solved, NT)

    # # D-SDA
    k = 'Infinity'
    starting_point = [1, 1]
    min_allowed = {i: 1 for i in range(1, len(starting_point)+1)}
    max_allowed = {i: NT for i in range(1, len(starting_point)+1)}

    m_solved, route = solve_with_dsda(model_function=build_cstrs, model_args={'NT': NT}, starting_point=starting_point, reformulation_function=external_ref, k=k,
                                      provide_starting_initialization=True, feasible_model='cstr', subproblem_solver='msnlp', min_allowed=min_allowed, max_allowed=max_allowed, iter_timelimit=10)
    #visualize_dsda(route=route, feas_x=x, feas_y=y, objs=objs, k=k,
    #               ext1_name='YF (Number of reactors)', ext2_name='YR (Reflux position)')
    print(m_solved.results)
    # visualize_cstr_superstructure(m_solved, NT)
