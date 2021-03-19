"""Distillation column model for 2018 PSE conference"""

from __future__ import division

from math import ceil, fabs

import pyomo.environ as pe
from pyomo.environ import (Block, BooleanVar, ConcreteModel, Constraint,
                           NonNegativeReals, Objective, Param, RangeSet, Set,
                           SolverFactory, Suffix, TransformationFactory, Var,
                           exactly, land, log, lor, minimize, value)
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints

from gdp.column.gdp_column import build_column
from gdp.dsda.dsda_functions import (
    generate_initialization, initialize_model, solve_subproblem, solve_with_dsda,
    solve_with_gdpopt, solve_with_minlp, visualize_dsda)


def external_ref(m, x, logic_expr=None):

    # Boolean variables and intTrays set definition
    m.intTrays = Set(initialize=m.trays -
                     [m.condens_tray, m.reboil_tray], doc='Interior trays of the column')
    m.YB = BooleanVar(m.intTrays, doc='Existence of boil-up flow in stage n')
    m.YR = BooleanVar(m.intTrays, doc='Existence of reflux flow in stage n')
    m.YP = BooleanVar(
        m.intTrays, doc='Boolean var associated with tray and no_tray')
    m.YB_is_up = BooleanVar()
    m.YR_is_down = BooleanVar()

    # Logical constraints

    @m.LogicalConstraint()
    def one_reflux(m):
        return exactly(1, m.YR)

    @m.LogicalConstraint()
    def one_boilup(m):
        return exactly(1, m.YB)

    @m.LogicalConstraint()
    def boilup_fix(m):
        return exactly(1, m.YB_is_up)

    @m.LogicalConstraint()
    def reflux_fix(m):
        return exactly(1, m.YR_is_down)

    @m.LogicalConstraint()
    def no_reflux_down(m):
        return m.YR_is_down.equivalent_to(land(~m.YR[n] for n in range(m.reboil_tray+1, m.feed_tray)))

    @m.LogicalConstraint()
    def no_boilup_up(m):
        return m.YB_is_up.equivalent_to(land(~m.YB[n] for n in range(m.feed_tray+1, m.max_trays)))

    @m.LogicalConstraint(m.conditional_trays)
    def YP_or_notYP(m, n):
        return m.YP[n].equivalent_to(land(lor(m.YR[j] for j in range(n, m.max_trays)), lor(land(~m.YB[j] for j in range(n, m.max_trays)), m.YB[n])))

    # Associate Boolean variables with with disjunctions
    for n in m.conditional_trays:
        m.YP[n].associate_binary_var(m.tray[n].indicator_var)

    # Fix externals

    ext_var_1 = x[0]
    ext_var_2 = x[1]

    for n in m.intTrays:
        if n == ext_var_1:
            m.YR[n].fix(True)
        else:
            m.YR[n].fix(False)

        if n == ext_var_2:
            m.YB[n].fix(True)
        else:
            m.YB[n].fix(False)

    temp = value(land(~m.YR[n]
                      for n in range(m.reboil_tray+1, m.feed_tray)))
    if temp == True:
        m.YR_is_down.fix(True)

    temp = value(land(~m.YB[n] for n in range(m.feed_tray+1, m.max_trays)))
    if temp == True:
        m.YB_is_up.fix(True)

    for n in m.conditional_trays:
        temp = value(land(lor(m.YR[j] for j in range(n, m.max_trays)), lor(
            land(~m.YB[j] for j in range(n, m.max_trays)), m.YB[n])))

        if temp == True:
            m.tray[n].indicator_var.fix(True)
            m.no_tray[n].indicator_var.fix(False)
        else:
            m.tray[n].indicator_var.fix(False)
            m.no_tray[n].indicator_var.fix(True)

    pe.TransformationFactory('core.logical_to_linear').apply_to(m)
    pe.TransformationFactory('gdp.fix_disjuncts').apply_to(m)
    pe.TransformationFactory('contrib.deactivate_trivial_constraints').apply_to(
        m, tmp=False, ignore_infeasible=True)

    return m


def complete_enumeration_external(model_function=build_column, model_args={'min_trays': 8, 'max_trays': 17, 'xD': 0.95, 'xB': 0.95}, reformulation_function=external_ref, subproblem_solver='conopt', timelimit=10):
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
        m_init = initialize_model(m, from_feasible=True, feasible_model='column')
        m_fixed = reformulation_function(m_init, x)
        m_solved = solve_subproblem(m_fixed, subproblem_solver=subproblem_solver, timelimit=timelimit)

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


if __name__ == "__main__":
    # Inputs
    NT = 17
    timelimit = 30
    model_args = {'min_trays': 8, 'max_trays': NT, 'xD': 0.95, 'xB': 0.95}

    # Complete enumeration
    # x, y, objs = complete_enumeration_external(model_function=build_column, model_args=model_args, subproblem_solver='conopt', timelimit=20)
    
    # Generate initialization
    # m = build_column(**model_args)
    # m_fixed = external_ref(m, [16, 2])
    # m_solved = solve_subproblem(
    #     m_fixed, subproblem_solver='baron', timelimit=30, tee=True)
    
    # MINLP methods
    m = build_column(**model_args)
    m_init = initialize_model(m, from_feasible=True, feasible_model='column')
    m_solved = solve_with_minlp(
       m_init, transformation='hull', minlp='antigone', timelimit=timelimit, gams_output=False)
    print(m_solved.results)

    # GDPopt methods
    # m = build_column(**model_args)
    # m_init = initialize_model(m, from_feasible=True, feasible_model='column')
    # m_solved = solve_with_gdpopt(m_init, mip='cplex',nlp='conopt', timelimit=timelimit, strategy='LOA', mip_output=False, nlp_output=False)
    # print(m_solved.results)

    # # D-SDA
    k = 'Infinity'
    starting_point = [16, 2]
    min_allowed = {i: 2 for i in range(1, len(starting_point)+1)}
    max_allowed = {i: NT-1 for i in range(1, len(starting_point)+1)}

    # m_solved, route = solve_with_dsda(model_function=build_column, model_args=model_args, starting_point=starting_point, reformulation_function=external_ref,
    #                                   k=k, provide_starting_initialization=True, feasible_model='column', subproblem_solver='conopt', min_allowed=min_allowed, max_allowed=max_allowed, iter_timelimit=10, timelimit=3600)
    # visualize_dsda(route=route, feas_x=x, feas_y=y, objs=objs, k=k, ext1_name='YR (Reflux position)', ext2_name='YB (Boil-up position)')
    # TODO This visualization code does not work
    # print(m_solved.results)
