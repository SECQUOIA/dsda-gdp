import pyomo.environ as pe
from pyomo.environ import (
    Block, ConcreteModel, Constraint, Param, log, minimize, NonNegativeReals, Objective, RangeSet, Set, Var, TransformationFactory, SolverFactory, value, BooleanVar, exactly, land, lor)
from pyomo.gdp import (Disjunct, Disjunction)
import networkx as nx
import matplotlib.pyplot as plt
import copy
import time
import numpy as np
import itertools as it
import math
from pyomo.core.base.misc import display
from pyomo.opt.base.solvers import SolverFactory
from pyomo.common.errors import InfeasibleConstraintException
from pyomo.contrib.fbbt.fbbt import fbbt
from pyomo.contrib.gdpopt.data_class import MasterProblemResult
from pyomo.core.plugins.transform.logical_to_linear import update_boolean_vars_from_binary
from pyomo.opt import SolutionStatus
from pyomo.opt import TerminationCondition as tc, SolverResults
import os
from model_serializer import to_json, from_json, StoreSpec

def build_column(min_trays, max_trays, xD, xB):
    """Builds the column model."""
    m = ConcreteModel('benzene-toluene column')
    m.comps = Set(initialize=['benzene', 'toluene'])
    min_T, max_T = 300, 400
    m.T_ref = 298.15
    max_flow = 500
    m.max_trays = max_trays
    m.condens_tray = max_trays
    m.feed_tray = math.ceil((max_trays / 2))
    m.reboil_tray = 1
    m.distillate_purity = xD
    m.bottoms_purity = xB
    m.pvap_const = {
        'benzene': {'A': -6.98273, 'B': 1.33213, 'C': -2.62863,
                    'D': -3.33399, 'Tc': 562.2, 'Pc': 48.9},
        'toluene': {'A': -7.28607, 'B': 1.38091, 'C': -2.83433,
                    'D': -2.79168, 'Tc': 591.8, 'Pc': 41.0}}
    m.vap_Cp_const = {
        'benzene': {'A': -3.392E1, 'B': 4.739E-1, 'C': -3.017E-4,
                    'D': 7.130E-8, 'E': 0},
        'toluene': {'A': -2.435E1, 'B': 5.125E-1, 'C': -2.765E-4,
                    'D': 4.911E-8, 'E': 0}}
    m.liq_Cp_const = {
        'benzene': {'A': 1.29E5, 'B': -1.7E2, 'C': 6.48E-1,
                    'D': 0, 'E': 0},
        'toluene': {'A': 1.40E5, 'B': -1.52E2, 'C': 6.95E-1,
                    'D': 0, 'E': 0}}
    m.dH_vap = {'benzene': 33.770E3, 'toluene': 38.262E3}  # J/mol

    m.trays = RangeSet(max_trays, doc='Set of potential trays')
    m.conditional_trays = Set(
        initialize=m.trays - [m.condens_tray, m.feed_tray, m.reboil_tray],
        doc="Trays that may be turned on and off.")
    m.tray = Disjunct(m.conditional_trays, doc='Disjunct for tray existence')
    m.no_tray = Disjunct(m.conditional_trays, doc='Disjunct for tray absence')

    @m.Disjunction(m.conditional_trays, doc='Tray exists or does not')
    def tray_no_tray(b, t):
        return [b.tray[t], b.no_tray[t]]
    m.minimum_num_trays = Constraint(
        expr=sum(m.tray[t].indicator_var
                 for t in m.conditional_trays) + 1  # for feed tray
        >= min_trays)

    m.T_feed = Var(
        doc='Feed temperature [K]', domain=NonNegativeReals,
        bounds=(min_T, max_T), initialize=368)
    m.feed_vap_frac = Var(
        doc='Vapor fraction of feed',
        initialize=0, bounds=(0, 1))
    m.feed = Var(
        m.comps, doc='Total component feed flow [mol/s]', initialize=50)
    m.x = Var(m.comps, m.trays, doc='Liquid mole fraction',
                bounds=(0, 1), domain=NonNegativeReals, initialize=0.5)
    m.y = Var(m.comps, m.trays, doc='Vapor mole fraction',
                bounds=(0, 1), domain=NonNegativeReals, initialize=0.5)
    m.L = Var(m.comps, m.trays,
                doc='component liquid flows from tray in kmol',
                domain=NonNegativeReals, bounds=(0, max_flow),
                initialize=50)
    m.V = Var(m.comps, m.trays,
                doc='component vapor flows from tray in kmol',
                domain=NonNegativeReals, bounds=(0, max_flow),
                initialize=50)
    m.liq = Var(m.trays, domain=NonNegativeReals,
                doc='liquid flows from tray in kmol', initialize=100,
                bounds=(0, max_flow))
    m.vap = Var(m.trays, domain=NonNegativeReals,
                doc='vapor flows from tray in kmol', initialize=100,
                bounds=(0, max_flow))
    m.B = Var(m.comps, domain=NonNegativeReals,
                doc='bottoms component flows in kmol',
                bounds=(0, max_flow), initialize=50)
    m.D = Var(m.comps, domain=NonNegativeReals,
                doc='distillate component flows in kmol',
                bounds=(0, max_flow), initialize=50)
    m.bot = Var(domain=NonNegativeReals, initialize=50, bounds=(0, 100),
                doc='bottoms flow in kmol')
    m.dis = Var(domain=NonNegativeReals, initialize=50,
                doc='distillate flow in kmol', bounds=(0, 100))
    m.reflux_ratio = Var(domain=NonNegativeReals, bounds=(0.5, 4),
                            doc='reflux ratio', initialize=1.4)
    m.reboil_ratio = Var(domain=NonNegativeReals, bounds=(1.3, 4),
                            doc='reboil ratio', initialize=0.9527)
    m.reflux_frac = Var(domain=NonNegativeReals, bounds=(0, 1 - 1E-6),
                        doc='reflux fractions')
    m.boilup_frac = Var(domain=NonNegativeReals, bounds=(0, 1 - 1E-6),
                        doc='boilup fraction')
    m.Kc = Var(
        m.comps, m.trays, doc='Phase equilibrium constant',
        domain=NonNegativeReals, initialize=1, bounds=(0, 1000))
    m.T = Var(m.trays, doc='Temperature [K]',
                domain=NonNegativeReals,
                bounds=(min_T, max_T))
    m.P = Var(doc='Pressure [bar]',
                bounds=(0, 5))
    m.gamma = Var(
        m.comps, m.trays,
        doc='liquid activity coefficent of component on tray',
        domain=NonNegativeReals, bounds=(0, 10), initialize=1)
    m.Pvap = Var(
        m.comps, m.trays,
        doc='pure component vapor pressure of component on tray in bar',
        domain=NonNegativeReals, bounds=(1E-3, 5), initialize=0.4)
    m.Pvap_X = Var(
        m.comps, m.trays,
        doc='Related to fraction of critical temperature (1 - T/Tc)',
        bounds=(0.25, 0.5), initialize=0.4)
    m.H_L = Var(
        m.comps, m.trays, bounds=(0.1, 16),
        doc='Liquid molar enthalpy of component in tray (kJ/mol)')
    m.H_V = Var(
        m.comps, m.trays, bounds=(30, 16 + 40),
        doc='Vapor molar enthalpy of component in tray (kJ/mol)')
    m.H_L_spec_feed = Var(
        m.comps, doc='Component liquid molar enthalpy in feed [kJ/mol]',
        initialize=0, bounds=(0.1, 16))
    m.H_V_spec_feed = Var(
        m.comps, doc='Component vapor molar enthalpy in feed [kJ/mol]',
        initialize=0, bounds=(30, 16 + 40))
    m.Qb = Var(domain=NonNegativeReals, doc='reboiler duty (MJ/s)',
                initialize=1, bounds=(0, 8))
    m.Qc = Var(domain=NonNegativeReals, doc='condenser duty (MJ/s)',
                initialize=1, bounds=(0, 8))

    m.partial_cond = Disjunct()
    m.total_cond = Disjunct()
    m.condenser_choice = Disjunction(expr=[m.partial_cond, m.total_cond])

    for t in m.conditional_trays:
        _build_conditional_tray_mass_balance(m, t, m.tray[t], m.no_tray[t])
    _build_feed_tray_mass_balance(m)
    _build_condenser_mass_balance(m)
    _build_reboiler_mass_balance(m)

    @m.Constraint(m.comps,
                  doc="Bottoms flow is equal to liquid leaving reboiler.")
    def bottoms_mass_balance(m, c):
        return m.B[c] == m.L[c, m.reboil_tray]

    @m.Constraint()
    def boilup_frac_defn(m):
        return m.bot == (1 - m.boilup_frac) * m.liq[m.reboil_tray + 1]

    @m.Constraint()
    def reflux_frac_defn(m):
        return m.dis == (1 - m.reflux_frac) * (
            m.vap[m.condens_tray - 1] - m.vap[m.condens_tray])

    @m.Constraint(m.trays)
    def liquid_sum(m, t):
        return sum(m.L[c, t] for c in m.comps) == m.liq[t]

    @m.Constraint(m.trays)
    def vapor_sum(m, t):
        return sum(m.V[c, t] for c in m.comps) == m.vap[t]

    m.bottoms_sum = Constraint(
        expr=sum(m.B[c] for c in m.comps) == m.bot)
    m.distil_sum = Constraint(
        expr=sum(m.D[c] for c in m.comps) == m.dis)

    @m.Constraint(m.trays)
    def monotonoic_temperature(_, t):
        return m.T[t] >= m.T[t + 1] if t < max_trays else Constraint.Skip

    for t in m.conditional_trays:
        _build_tray_phase_equilibrium(m, t, m.tray[t])
    m.feed_tray_phase_eq = Block()
    m.reboiler_phase_eq = Block()
    m.condenser_phase_eq = Block()
    _build_tray_phase_equilibrium(m, m.feed_tray, m.feed_tray_phase_eq)
    _build_tray_phase_equilibrium(m, m.reboil_tray, m.reboiler_phase_eq)
    _build_tray_phase_equilibrium(m, m.condens_tray, m.condenser_phase_eq)
    _build_column_heat_relations(m)

    @m.Constraint()
    def distillate_req(m):
        return m.D['benzene'] >= m.distillate_purity * m.dis

    @m.Constraint()
    def bottoms_req(m):
        return m.B['toluene'] >= m.bottoms_purity * m.bot

    m.obj = Objective(
        expr=(m.Qc + m.Qb) * 1E3 + 1E3 * (
            sum(m.tray[t].indicator_var for t in m.conditional_trays) + 1),
        sense=minimize)

    @m.Constraint()
    def reflux_ratio_calc(m):
        return m.reflux_frac * (m.reflux_ratio + 1) == m.reflux_ratio

    @m.Constraint()
    def reboil_ratio_calc(m):
        return m.boilup_frac * (m.reboil_ratio + 1) == m.reboil_ratio

    @m.Constraint(m.conditional_trays)
    def tray_ordering(m, t):
        """Trays close to the feed should be activated first."""
        if t + 1 < m.condens_tray and t > m.feed_tray:
            return m.tray[t].indicator_var >= m.tray[t + 1].indicator_var
        elif t > m.reboil_tray and t + 1 < m.feed_tray:
            return m.tray[t + 1].indicator_var >= m.tray[t].indicator_var
        else:
            return Constraint.NoConstraint

    # Fix feed conditions
    m.feed['benzene'].fix(50)
    m.feed['toluene'].fix(50)
    m.T_feed.fix(368)
    m.feed_vap_frac.fix(0.40395)
    m.P.fix(1.01)
    # Fix to be total condenser
    m.partial_cond.deactivate()
    m.total_cond.indicator_var.fix(1)

    m.dsda_status = 'Initialized'

    return m

def external_ref(m, x, logic_expr = None):
    
     # Boolean variables and intTrays set definition
    m.intTrays = Set(initialize=m.trays - [m.condens_tray, m.reboil_tray], doc='Interior trays of the column')
    m.YB = BooleanVar(m.intTrays, doc='Existence of boil-up flow in stage n')
    m.YR = BooleanVar(m.intTrays, doc='Existence of reflux flow in stage n')
    m.YP = BooleanVar(m.intTrays, doc='Boolean var associated with tray and no_tray')
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
    pe.TransformationFactory('contrib.deactivate_trivial_constraints').apply_to(m, tmp=False, ignore_infeasible=True)

    return m

def solve_nlp(m, nlp='msnlp', timelimit=10):
    try:
        fbbt(m)
        # SOLVE
        dir_path = os.path.dirname(os.path.abspath(__file__))
        gams_path = os.path.join(dir_path, "gamsfiles/")
        if not(os.path.exists(gams_path)):
            print('Directory for automatically generated files ' +
                  gams_path + ' does not exist. We will create it')
            os.makedirs(gams_path)

        solvername = 'gams'
        opt = SolverFactory(solvername, solver=nlp)
        m.results = opt.solve(m, tee=False,
                              # Uncomment the following lines if you want to save GAMS models
                              #keepfiles=True,
                              #tmpdir=gams_path,
                              #symbolic_solver_labels=True,
                              skip_trivial_constraints=True,
                              add_options=[
                                  'option reslim = ' + str(timelimit) + ';'
                                  'option optcr = 0.0;'
                                  # Uncomment the following lines to setup IIS computation of BARON through option file
                                  # 'GAMS_MODEL.optfile = 1;'
                                  # '\n'
                                  # '$onecho > baron.opt \n'
                                  # 'CompIIS 1 \n'
                                  # '$offecho'
                                  # 'display(execError);'
                              ])

        if m.results.solver.termination_condition == 'locallyOptimal' or m.results.solver.termination_condition == 'optimal':
            m.dsda_status = 'Optimal'
        elif m.results.solver.termination_condition == 'infeasible':
            m.dsda_status = 'Evaluated_Infeasible'

    except InfeasibleConstraintException:
        nlp_result = MasterProblemResult()
        nlp_result.feasible = False
        nlp_result.pyomo_results = SolverResults()
        nlp_result.pyomo_results.solver.termination_condition = tc.error
        m.dsda_status = 'FBBT_Infeasible'
        #print('Try an infeasible')

    return m



def list_generator(NT):
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
    return X1, X2


def complete_enumeration_external(model_function=build_column, model_args={'min_trays':8, 'max_trays':17, 'xD':0.95, 'xB':0.95}, reformulation_function=external_ref, nlp='knitro', timelimit = 10):
    NT = model_args['max_trays']
    X1, X2 = list_generator(NT)

    print()
    feas_x, feas_y, objs = [], [], []

    print('=============================')
    print('%6s %6s %12s' % ('x1', 'x2', 'Objective'))
    print('-----------------------------')

    # Loop over all external variables and then loop over its values
    for i in range(len(X1)):
        x = [X1[i], X2[i]]
        m = model_function(**model_args)
        m_init = initialize_model(m, from_feasible=True)
        m_fixed = reformulation_function(m_init, x)
        m_solved = solve_nlp(m_fixed, nlp=nlp, timelimit=timelimit)

        if m_solved.dsda_status == 'Optimal':
            print('%6s %6s %12s' % (X1[i], X2[i], round(pe.value(m_solved.obj), 2)))
            feas_x.append(X1[i])
            feas_y.append(X2[i])
            objs.append(round(pe.value(m_solved.obj), 2))
        else:
            print('%6s %6s %12s' % (X1[i], X2[i], 'Infeasible'))

    print('=============================')
    return feas_x, feas_y, objs


def _build_conditional_tray_mass_balance(m, t, tray, no_tray):
    """
    t = tray number
    tray = tray exists disjunct
    no_tray = tray absent disjunct
    """
    @tray.Constraint(m.comps)
    def mass_balance(_, c):
        return (
            # Feed in if feed tray
            (m.feed[c] if t == m.feed_tray else 0)
            # Vapor from tray t
            - m.V[c, t]
            # Loss to distillate if condenser
            - (m.D[c] if t == m.condens_tray else 0)
            # Liquid from tray above if not condenser
            + (m.L[c, t + 1] if t < m.condens_tray else 0)
            # Loss to bottoms if reboiler
            - (m.B[c] if t == m.reboil_tray else 0)
            # Liquid to tray below if not reboiler
            - (m.L[c, t] if t > m.reboil_tray else 0)
            # Vapor from tray below if not reboiler
            + (m.V[c, t - 1] if t > m.reboil_tray else 0) == 0)

    @tray.Constraint(m.comps)
    def tray_liquid_composition(_, c):
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    @tray.Constraint(m.comps)
    def tray_vapor_compositions(_, c):
        return m.V[c, t] == m.vap[t] * m.y[c, t]

    @no_tray.Constraint(m.comps)
    def liq_comp_pass_through(_, c):
        return m.x[c, t] == m.x[c, t + 1]

    @no_tray.Constraint(m.comps)
    def liq_flow_pass_through(_, c):
        return m.L[c, t] == m.L[c, t + 1]

    @no_tray.Constraint(m.comps)
    def vap_comp_pass_through(_, c):
        return m.y[c, t] == m.y[c, t - 1]

    @no_tray.Constraint(m.comps)
    def vap_flow_pass_through(_, c):
        return m.V[c, t] == m.V[c, t - 1]


def _build_feed_tray_mass_balance(m):
    t = m.feed_tray

    @m.Constraint(m.comps)
    def feed_mass_balance(_, c):
        return (
            m.feed[c]        # Feed in
            - m.V[c, t]      # Vapor from tray t
            + m.L[c, t + 1]  # Liquid from tray above
            - m.L[c, t]      # Liquid to tray below
            + m.V[c, t - 1]  # Vapor from tray below
            == 0)

    @m.Constraint(m.comps)
    def feed_tray_liquid_composition(_, c):
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    @m.Constraint(m.comps)
    def feed_tray_vapor_composition(_, c):
        return m.V[c, t] == m.vap[t] * m.y[c, t]


def _build_condenser_mass_balance(m):
    t = m.condens_tray

    @m.Constraint(m.comps)
    def condenser_mass_balance(_, c):
        return (
            - m.V[c, t]      # Vapor from tray t
            - m.D[c]         # Loss to distillate
            - m.L[c, t]      # Liquid to tray below
            + m.V[c, t - 1]  # Vapor from tray below
            == 0)

    @m.partial_cond.Constraint(m.comps)
    def condenser_liquid_composition(_, c):
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    @m.partial_cond.Constraint(m.comps)
    def condenser_vapor_composition(_, c):
        return m.V[c, t] == m.vap[t] * m.y[c, t]

    @m.total_cond.Constraint(m.comps)
    def no_vapor_flow(_, c):
        return m.V[c, t] == 0

    @m.total_cond.Constraint()
    def no_total_vapor_flow(_):
        return m.vap[t] == 0

    @m.total_cond.Constraint(m.comps)
    def liquid_fraction_pass_through(_, c):
        return m.x[c, t] == m.y[c, t - 1]

    @m.Constraint(m.comps)
    def condenser_distillate_composition(_, c):
        return m.D[c] == m.dis * m.x[c, t]


def _build_reboiler_mass_balance(m):
    t = m.reboil_tray

    @m.Constraint(m.comps)
    def reboiler_mass_balance(_, c):
        t = m.reboil_tray
        return (
            - m.V[c, t]      # Vapor from tray t
            + m.L[c, t + 1]  # Liquid from tray above
            - m.B[c]         # Loss to bottoms
            == 0)

    @m.Constraint(m.comps)
    def reboiler_liquid_composition(_, c):
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    @m.Constraint(m.comps)
    def reboiler_vapor_composition(_, c):
        return m.V[c, t] == m.vap[t] * m.y[c, t]


def _build_tray_phase_equilibrium(m, t, tray):
    @tray.Constraint(m.comps)
    def raoults_law(_, c):
        return m.y[c, t] == m.x[c, t] * m.Kc[c, t]

    @tray.Constraint(m.comps)
    def phase_equil_const(_, c):
        return m.Kc[c, t] * m.P == (
            m.gamma[c, t] * m.Pvap[c, t])

    @tray.Constraint(m.comps)
    def Pvap_relation(_, c):
        k = m.pvap_const[c]
        x = m.Pvap_X[c, t]
        return (log(m.Pvap[c, t]) - log(k['Pc'])) * (1 - x) == (
            k['A'] * x +
            k['B'] * x ** 1.5 +
            k['C'] * x ** 3 +
            k['D'] * x ** 6)

    @tray.Constraint(m.comps)
    def Pvap_X_defn(_, c):
        k = m.pvap_const[c]
        return m.Pvap_X[c, t] == 1 - m.T[t] / k['Tc']

    @tray.Constraint(m.comps)
    def gamma_calc(_, c):
        return m.gamma[c, t] == 1


def _build_column_heat_relations(m):
    @m.Expression(m.trays, m.comps)
    def liq_enthalpy_expr(_, t, c):
        k = m.liq_Cp_const[c]
        return (
            k['A'] * (m.T[t] - m.T_ref) +
            k['B'] * (m.T[t] ** 2 - m.T_ref ** 2) / 2 +
            k['C'] * (m.T[t] ** 3 - m.T_ref ** 3) / 3 +
            k['D'] * (m.T[t] ** 4 - m.T_ref ** 4) / 4 +
            k['E'] * (m.T[t] ** 5 - m.T_ref ** 5) / 5) * 1E-6

    @m.Expression(m.trays, m.comps)
    def vap_enthalpy_expr(_, t, c):
        k = m.vap_Cp_const[c]
        return (
            m.dH_vap[c] +
            k['A'] * (m.T[t] - m.T_ref) +
            k['B'] * (m.T[t] ** 2 - m.T_ref ** 2) / 2 +
            k['C'] * (m.T[t] ** 3 - m.T_ref ** 3) / 3 +
            k['D'] * (m.T[t] ** 4 - m.T_ref ** 4) / 4 +
            k['E'] * (m.T[t] ** 5 - m.T_ref ** 5) / 5) * 1E-3

    for t in m.conditional_trays:
        _build_conditional_tray_energy_balance(m, t, m.tray[t], m.no_tray[t])
    _build_feed_tray_energy_balance(m)
    _build_condenser_energy_balance(m)
    _build_reboiler_energy_balance(m)


def _build_conditional_tray_energy_balance(m, t, tray, no_tray):
    @tray.Constraint()
    def energy_balance(_):
        return sum(
            m.L[c, t + 1] * m.H_L[c, t + 1]  # heat of liquid from tray above
            - m.L[c, t] * m.H_L[c, t]  # heat of liquid to tray below
            + m.V[c, t - 1] * m.H_V[c, t - 1]  # heat of vapor from tray below
            - m.V[c, t] * m.H_V[c, t]  # heat of vapor to tray above
            for c in m.comps) * 1E-3 == 0

    @tray.Constraint(m.comps)
    def liq_enthalpy_calc(_, c):
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]

    @tray.Constraint(m.comps)
    def vap_enthalpy_calc(_, c):
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]

    @no_tray.Constraint(m.comps)
    def liq_enthalpy_pass_through(_, c):
        return m.H_L[c, t] == m.H_L[c, t + 1]

    @no_tray.Constraint(m.comps)
    def vap_enthalpy_pass_through(_, c):
        return m.H_V[c, t] == m.H_V[c, t - 1]


def _build_feed_tray_energy_balance(m):
    t = m.feed_tray

    @m.Constraint()
    def feed_tray_energy_balance(_):
        return (
            sum(m.feed[c] * (
                m.H_L_spec_feed[c] * (1 - m.feed_vap_frac) +
                m.H_V_spec_feed[c] * m.feed_vap_frac)
                for c in m.comps) +
            sum(
                # Heat of liquid from tray above
                m.L[c, t + 1] * m.H_L[c, t + 1]
                # heat of liquid to tray below
                - m.L[c, t] * m.H_L[c, t]
                # heat of vapor from tray below
                + m.V[c, t - 1] * m.H_V[c, t - 1]
                # heat of vapor to tray above
                - m.V[c, t] * m.H_V[c, t]
                for c in m.comps)) * 1E-3 == 0

    @m.Constraint(m.comps)
    def feed_tray_liq_enthalpy_calc(_, c):
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]

    @m.Constraint(m.comps)
    def feed_tray_vap_enthalpy_calc(_, c):
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]

    @m.Expression(m.comps)
    def feed_liq_enthalpy_expr(_, c):
        k = m.liq_Cp_const[c]
        return (
            k['A'] * (m.T_feed - m.T_ref) +
            k['B'] * (m.T_feed ** 2 - m.T_ref ** 2) / 2 +
            k['C'] * (m.T_feed ** 3 - m.T_ref ** 3) / 3 +
            k['D'] * (m.T_feed ** 4 - m.T_ref ** 4) / 4 +
            k['E'] * (m.T_feed ** 5 - m.T_ref ** 5) / 5) * 1E-6

    @m.Constraint(m.comps)
    def feed_liq_enthalpy_calc(_, c):
        return m.H_L_spec_feed[c] == m.feed_liq_enthalpy_expr[c]

    @m.Expression(m.comps)
    def feed_vap_enthalpy_expr(_, c):
        k = m.vap_Cp_const[c]
        return (
            m.dH_vap[c] +
            k['A'] * (m.T_feed - m.T_ref) +
            k['B'] * (m.T_feed ** 2 - m.T_ref ** 2) / 2 +
            k['C'] * (m.T_feed ** 3 - m.T_ref ** 3) / 3 +
            k['D'] * (m.T_feed ** 4 - m.T_ref ** 4) / 4 +
            k['E'] * (m.T_feed ** 5 - m.T_ref ** 5) / 5) * 1E-3

    @m.Constraint(m.comps)
    def feed_vap_enthalpy_calc(_, c):
        return m.H_V_spec_feed[c] == m.feed_vap_enthalpy_expr[c]


def _build_condenser_energy_balance(m):
    t = m.condens_tray

    @m.partial_cond.Constraint()
    def partial_condenser_energy_balance(_):
        return -m.Qc + sum(
            - m.D[c] * m.H_L[c, t]  # heat of liquid distillate
            - m.L[c, t] * m.H_L[c, t]  # heat of liquid to tray below
            + m.V[c, t - 1] * m.H_V[c, t - 1]  # heat of vapor from tray below
            - m.V[c, t] * m.H_V[c, t]  # heat of vapor from partial condenser
            for c in m.comps) * 1E-3 == 0

    @m.total_cond.Constraint()
    def total_condenser_energy_balance(_):
        return -m.Qc + sum(
            - m.D[c] * m.H_L[c, t]  # heat of liquid distillate
            - m.L[c, t] * m.H_L[c, t]  # heat of liquid to tray below
            + m.V[c, t - 1] * m.H_V[c, t - 1]  # heat of vapor from tray below
            for c in m.comps) * 1E-3 == 0

    @m.Constraint(m.comps)
    def condenser_liq_enthalpy_calc(_, c):
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]

    @m.partial_cond.Constraint(m.comps)
    def vap_enthalpy_calc(_, c):
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]


def _build_reboiler_energy_balance(m):
    t = m.reboil_tray

    @m.Constraint()
    def reboiler_energy_balance(_):
        return m.Qb + sum(
            m.L[c, t + 1] * m.H_L[c, t + 1]  # Heat of liquid from tray above
            - m.B[c] * m.H_L[c, t]  # heat of liquid bottoms if reboiler
            - m.V[c, t] * m.H_V[c, t]  # heat of vapor to tray above
            for c in m.comps) * 1E-3 == 0

    @m.Constraint(m.comps)
    def reboiler_liq_enthalpy_calc(_, c):
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]

    @m.Constraint(m.comps)
    def reboiler_vap_enthalpy_calc(_, c):
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]

def visualize_dsda(points=[], feas_x=[], feas_y=[], objs=[], k='?'):

    X1, X2 = feas_x, feas_y
    cm = plt.cm.get_cmap('viridis_r')

    def drawArrow(A, B):
        plt.arrow(A[0], A[1], B[0] - A[0], B[1] - A[1], width=0.00005,
                  head_width=0.15, head_length=0.08, color='black', shape='full')

    for i in range(len(points)-1):
        drawArrow(points[i], points[i+1])

    sc = plt.scatter(X1, X2, s=80, c=objs, cmap=cm)
    cbar = plt.colorbar(sc)
    cbar.set_label('Objective function', rotation=90)
    title_string = 'D-SDA with k = '+k
    plt.title(title_string)
    plt.xlabel("YR (Reflux position)")
    plt.ylabel("YB (Boil-up position)")
    plt.show()


def solve_with_minlp(m, transformation='bigm', minlp='baron', timelimit=10, gams_output=False):

    # Transformation step
    pe.TransformationFactory('core.logical_to_linear').apply_to(m)
    transformation_string = 'gdp.' + transformation
    pe.TransformationFactory(transformation_string).apply_to(m)

    # Solution step
    output_options = {}
    if gams_output:
        dir_path = os.path.dirname(os.path.abspath(__file__))
        gams_path = os.path.join(dir_path, "gamsfiles/")
        if not(os.path.exists(gams_path)):
            print('Directory for automatically generated files ' +
                gams_path + ' does not exist. We will create it')
            os.makedirs(gams_path)
        output_options = {'keepfiles':True,
                        'tmpdir':gams_path,
                        'symbolic_solver_labels':True}

    solvername = 'gams'
    opt = SolverFactory(solvername, solver=minlp)
    m.results = opt.solve(m, tee=True,
                          **output_options,
                          # Uncomment the following lines if you want to save GAMS models
                          add_options=[
                              'option reslim = ' + str(timelimit) + ';'
                              'option optcr = 0.0;'
                              # Uncomment the following lines to setup IIS computation of BARON through option file
                              # 'GAMS_MODEL.optfile = 1;'
                              # '\n'
                              # '$onecho > baron.opt \n'
                              # 'CompIIS 1 \n'
                              # '$offecho'
                              # 'display(execError);'
                          ])
    #update_boolean_vars_from_binary(m)
    return m


def solve_with_gdpopt(m, mip='cplex', nlp='knitro', timelimit=10, strategy='LOA', mip_output=False, nlp_output=False):
    """
    Function documentation
    """

    # Transformation step
    pe.TransformationFactory('core.logical_to_linear').apply_to(m)

    # Solution step
    mip_output_options = {}
    nlp_output_options = {}
    if mip_output:
        dir_path = os.path.dirname(os.path.abspath(__file__))
        gams_path = os.path.join(dir_path, "gamsfiles/")
        if not(os.path.exists(gams_path)):
            print('Directory for automatically generated files ' +
                gams_path + ' does not exist. We will create it')
            os.makedirs(gams_path)
        mip_output_options = {'keepfiles':True,
                        'tmpdir':gams_path,
                        'symbolic_solver_labels':True}
    
    if nlp_output:
        dir_path = os.path.dirname(os.path.abspath(__file__))
        gams_path = os.path.join(dir_path, "gamsfiles/")
        if not(os.path.exists(gams_path)):
            print('Directory for automatically generated files ' +
                gams_path + ' does not exist. We will create it')
            os.makedirs(gams_path)
        nlp_output_options = {'keepfiles':True,
                        'tmpdir':gams_path,
                        'symbolic_solver_labels':True}

    solvername = 'gdpopt'
    opt = SolverFactory(solvername)
    m.results = opt.solve(m, tee=True,
                          strategy=strategy,
                          time_limit=timelimit,
                          mip_solver='gams',
                          mip_solver_args=dict(solver=mip, warmstart=True, **mip_output_options),
                          nlp_solver='gams',
                          nlp_solver_args=dict(solver=nlp, warmstart=True, tee=True, **nlp_output_options),
                        #   mip_presolve=True,
                        #   init_strategy='fix_disjuncts',
                        #   set_cover_iterlim=0,
                          iterlim=20,
                          force_subproblem_nlp=True,
                          subproblem_presolve=False,
                        #   calc_disjunctive_bounds=True
                          )
    #update_boolean_vars_from_binary(m)
    return m

def neighborhood_k_eq_2(num_ext):
    num_neigh = 2*num_ext
    neighbors = np.concatenate((np.eye(num_ext), -np.eye(num_ext)), axis=1)
    directions = {}
    for i in range(num_neigh):
        direct = []
        directions[i+1] = direct
        for j in range(num_ext):
            direct.append(neighbors[j, i])
    return directions


def neighborhood_k_eq_inf(num_ext):  # number
    num_neigh = 3*num_ext-1
    neighbors = list(it.product([-1, 0, 1], repeat=num_ext))
    directions = {}
    for i in range(len(neighbors)):
        directions[i+1] = list(neighbors[i])
    temp = directions.copy()
    for i in directions.keys():
        if temp[i] == [0]*num_ext:
            temp.pop(i, None)
    return temp


# Creates neighbor of a given point
# Optimize option will discard out of bounds points given by min and max allowed
# Cheating will discard infeasible points for CSTR problem (X2 - X1 > 0)
# Cheating will only be used while solving GAMS issue

# start is type list and stands for the actual point
# neighborhood is type dict and is the output of a k-Neighborhood function
# newbors or new_newbors is type  dict starting in 0 with neighbor of a given point
# Neighbor 0 is the actual point

def find_actual_neighbors(start, neighborhood, optimize=True, min_allowed={}, max_allowed={}):
    neighbors = {0: start}
    for i in neighborhood.keys():
        neighbors[i] = list(map(sum, zip(start, list(neighborhood[i]))))

    if optimize:
        new_neighbors = {}
        num_vars = len(neighbors[0])
        for i in neighbors.keys():
            checked = 0
            for j in range(num_vars):
                if neighbors[i][j] >= min_allowed[j+1] and neighbors[i][j] <= max_allowed[j+1]:
                    checked += 1
            if checked == num_vars:
                new_neighbors[i] = neighbors[i]

        return new_neighbors
    return neighbors

# Evaluates a group of given points and returns the best
# ext_vars is dict with given points where 0 is actual point
# init is type dict and contains solved variables for the actual point
# fmin is type int and stands for objective at actual point
# tol is type int and stands for numerical tolereance for equallity
# fmin as return is type int and gives the best neighbor's objective
# best_var is type list and gives the best neighbor
# best_dir is type int and is the steepest direction (key in neighborhood)
# best_init is type dict and contains solved variables for the best point
# improve is type bool and shows if an improvement was made while looking for neighbors
def evaluate_neighbors(ext_vars, fmin, model_function=build_column, model_args={'min_trays':8, 'max_trays':17, 'xD':0.95, 'xB':0.95}, reformulation_function=external_ref, nlp='knitro', iter_timelimit=10, tol=0.000001):
    improve = False
    best_var = ext_vars[0]
    here = ext_vars[0]
    best_dir = 0
    temp = ext_vars
    temp.pop(0, None)
    objectives = {}
    feasibles = {}

    for i in temp.keys():

        m = model_function(**model_args)
        m_init = initialize_model(m)
        m_fixed = reformulation_function(m_init, temp[i])
        m_solved = solve_nlp(m_fixed, nlp=nlp, timelimit=iter_timelimit)
        #print(temp[i],m_solved.dsda_status)
        
        if m_solved.dsda_status == 'Optimal':
            
            objectives[i] = pe.value(m_solved.obj)
            feasibles[i] = temp[i]

    key_min = min(objectives.keys(), key=(lambda k: objectives[k]))
    min_obj = objectives[key_min]
    mins = 0
    min_points = {}

    for i in objectives.keys():
        if abs(objectives[i] - min_obj) < tol:
            min_points[i] = feasibles[i]
            mins += 1

    if mins > 1:
        ssums = {}
        for i in min_points.keys():
            ssum = 0
            for j in range(len(best_var)):
                ssum += (min_points[i][j] - here[j])**2
            ssums[i] = ssum

        key_max = max(ssums.keys(), key=(lambda k: ssums[k]))

        if objectives[key_max] + tol < fmin:
            fmin = objectives[key_max]
            best_var = ext_vars[key_max]
            best_dir = key_max
            improve = True
    else:
        if objectives[key_min] + tol < fmin:
            fmin = objectives[key_min]
            best_var = ext_vars[key_min]
            best_dir = key_min
            improve = True

    if improve == True:
        m2 = model_function(**model_args)
        m2_init = initialize_model(m2)
        m2_fixed = reformulation_function(m2_init, best_var)
        m2_solved = solve_nlp(m2_fixed, nlp=nlp, timelimit=iter_timelimit)
        generate_initialization(m2_solved)

    return fmin, best_var, best_dir, improve


# Moves from a certain start in a given direction and evaluates it
# start is type list with the actual point
# init is type dict and contains solved variables for the actual point
# fmin is type int and stands for objective at actual point
# direction is type int and is the moving direction direction (key in neighborhood)
# optimize option will discard out of bounds points given by min and max allowed
# tol is type int and stands for numerical tolereance for equallity
# fmin as return is type int and gives the best point objective (between moved and actual)
# best_var is type list and gives the best point (between moved and actual)
# move is type bool and shows if an improvement was made while looking for neighbors
# best_init is type dict and contains solved variables for the best point
def do_line_search(start, fmin, direction, model_function=build_column, model_args={'min_trays':8, 'max_trays':17, 'xD':0.95, 'xB':0.95}, reformulation_function=external_ref, nlp='knitro', optimize=True, min_allowed={}, max_allowed={}, iter_timelimit=10, tol=0.000001):
    best_var = start
    moved = False

    moved_point = list(map(sum, zip(list(start), list(direction))))

    if optimize:
        checked = 0
        for j in range(len(moved_point)):
            if moved_point[j] >= min_allowed[j+1] and moved_point[j] <= max_allowed[j+1]:
                checked += 1

        if checked == len(moved_point):
            m = model_function(**model_args)
            m_init = initialize_model(m)
            m_fixed = reformulation_function(m_init, moved_point)
            m_solved = solve_nlp(m_fixed, nlp=nlp, timelimit=iter_timelimit)

            if m_solved.dsda_status == 'Optimal':
                act_obj = pe.value(m_solved.obj)
                if act_obj + tol < fmin:
                    fmin = act_obj
                    best_var = moved_point
                    moved = True
    else:
        m = model_function(**model_args)
        m_init = initialize_model(m)
        m_fixed = reformulation_function(m_init, moved_point)
        m_solved = solve_nlp(m_fixed, nlp=nlp, timelimit=iter_timelimit)

        if m_solved.dsda_status == 'Optimal':
            act_obj = pe.value(m_solved.obj)
            if act_obj + tol < fmin:
                fmin = act_obj
                best_var = moved_point
                moved = True
    
    if moved == True:
        m2 = model_function(**model_args)
        m2_init = initialize_model(m2)
        m2_fixed = reformulation_function(m2_init, best_var)
        m2_solved = solve_nlp(m2_fixed, nlp=nlp, timelimit=iter_timelimit)
        generate_initialization(m2_solved)

    return fmin, best_var, moved



def initialize_model(m, from_feasible=False):
    wts = StoreSpec.value()

    if from_feasible:
        from_json(m, fname = 'dsda_starting_initialization.json', wts=wts)
    else:
        from_json(m, fname = 'dsda_initialization.json', wts=wts)
    return m

def generate_initialization(m, starting_initialization=False):
    wts = StoreSpec.value()

    if starting_initialization:
        to_json(m, fname = 'dsda_starting_initialization.json', human_read=True, wts=wts)
    else:
        to_json(m, fname='dsda_initialization.json', human_read=True, wts=wts)


def solve_with_dsda(k='Infinity', model_function=build_column, model_args={'min_trays':8, 'max_trays':17, 'xD':0.95, 'xB':0.95}, starting_point=[16,2], reformulation_function=external_ref, provide_starting_initialization=True, nlp='knitro', optimize=True, min_allowed={}, max_allowed={}, iter_timelimit=10, tol=0.000001):

    print('\nStarting D-SDA with k =', k)

    # Initialize
    t_start = time.process_time()
    route = []
    ext_var = starting_point
    route.append(ext_var)

    m = model_function(**model_args)
    if provide_starting_initialization:
        m_init = initialize_model(m, from_feasible=True)
        m_fixed = external_ref(m_init, ext_var)
    else:
        m_fixed = external_ref(m, ext_var)

    m_solved = solve_nlp(m_fixed, nlp=nlp, timelimit=iter_timelimit)

    fmin = pe.value(m_solved.obj)
    generate_initialization(m_solved)
    
    # Define neighborhood
    if k == '2':
        neighborhood = neighborhood_k_eq_2(len(ext_var))
    elif k == 'Infinity':
        neighborhood = neighborhood_k_eq_inf(len(ext_var))
    elif k == 'l_flat':
        neighborhood = {1: [1, 1], 2: [-1, -1],
                        3: [1, 0], 4: [-1, 0], 5: [0, 1], 6: [0, -1]}
    elif k == 'm_flat':
        neighborhood = {1: [1, -1], 2: [1, 0],
                        3: [-1, 1], 4: [0, 1], 5: [-1, 0], 6: [0, -1]}
    else:
        return "Enter a valid neighborhood ('Infinity', '2', 'l_flat' or 'm_flat')"

    looking_in_neighbors = True

    # Look in neighbors (outter cycle)
    while looking_in_neighbors:

        # Find neighbors of the actual point
        neighbors = find_actual_neighbors(ext_var, neighborhood, optimize=True,
                                 min_allowed=min_allowed, max_allowed=max_allowed)

        fmin, best_var, best_dir, improve = evaluate_neighbors(neighbors, fmin, model_function=model_function, model_args=model_args, reformulation_function=external_ref, nlp=nlp, iter_timelimit=iter_timelimit, tol=tol)

        # Stopping condition in case there is no improvement amongst neighbors
        if improve == True:
            line_searching = True
            route.append(best_var)

            # If improvement was made start line search (inner cycle)
            while line_searching:
                fmin, best_var, moved = do_line_search(best_var, fmin, neighborhood[best_dir], model_function=model_function, model_args=model_args, reformulation_function=external_ref, nlp=nlp, optimize=optimize, min_allowed=min_allowed, max_allowed=max_allowed, iter_timelimit=iter_timelimit, tol=tol)

                # Stopping condition in case no movement was done
                if moved == True:
                    route.append(best_var)
                else:
                    ext_var = best_var
                    line_searching = False

        else:
            looking_in_neighbors = False

    t_end = round(time.process_time() - t_start, 2)

    print('Objective:',round(fmin, 5))
    print('External variables:', route[-1])
    print('Execution time [s]:', t_end)

    m2 = model_function(**model_args)
    m2_init = initialize_model(m2)
    m2_fixed = external_ref(m2_init, route[-1])
    m2_solved = solve_nlp(m2_fixed, nlp=nlp, timelimit=iter_timelimit)
    m2_solved.dsda_time = t_end

    # Return visited points / final point / objective at that point / execution time
    return m2_solved, route


if __name__ == "__main__":
    # Inputs
    NT = 17
    timelimit = 10
    model_args = {'min_trays':8, 'max_trays':NT, 'xD':0.95, 'xB':0.95}

    #Complete enumeration
    x, y, objs = complete_enumeration_external(model_function=build_column, model_args=model_args, nlp='knitro', timelimit=20)

    # MINLP and GDPopt methods
    m = build_column(**model_args)
    m_init = initialize_model(m, from_feasible=True)
    m_solved = solve_with_minlp(m_init, transformation='bigm', minlp='baron', timelimit=timelimit, gams_output=False)
    # m_solved = solve_with_gdpopt(m_init, mip='cplex',nlp='knitro', timelimit=timelimit, strategy='LOA', mip_output=False, nlp_output=False)
    print(m_solved.results)
    
    # D-SDA
    k = 'Infinity'
    starting_point = [16,2]
    min_allowed = {i: 2 for i in range(1, len(starting_point)+1)}
    max_allowed = {i: NT-1 for i in range(1, len(starting_point)+1)}

    m_solved, route = solve_with_dsda(k=k, model_function=build_column, model_args=model_args, starting_point=starting_point, reformulation_function=external_ref, provide_starting_initialization=True, nlp='knitro', min_allowed=min_allowed, max_allowed=max_allowed, iter_timelimit=10)
    visualize_dsda(points=route, feas_x=x, feas_y=y, objs=objs, k=k)
    print(m_solved.results)
