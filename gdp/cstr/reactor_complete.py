import pyomo.environ as pe
from pyomo.gdp import (Disjunct, Disjunction)
import networkx as nx
import matplotlib.pyplot as plt
import copy
import time
import numpy as np
import itertools as it
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


def build_cstrs(NT: int = 5) -> pe.ConcreteModel():
    """
    Function that builds CSTR superstructure model of size NT.
    The CSTRs have a single 1st order reaction A -> B and minimizes (TODO Check)
    total reactor volume. The optimal solution should yield NT reactors with a recycle before reactor NT.
    Reference: Paper Linhan 1.

    Args:
        NT: int. Positive Integer defining the maximum number of CSTRs
    Returns:
        m = Pyomo GDP model
    """

    # PYOMO MODEL
    m = pe.ConcreteModel(name='gdp_reactors')

    # SETS
    m.I = pe.Set(initialize=['A', 'B'], doc='Set of components')
    m.N = pe.RangeSet(1, NT, doc='Set of units in the superstructure')

    # PARAMETERS
    m.k = pe.Param(initialize=2)  # Kinetic constant [L/(mol*s)]
    m.order1 = pe.Param(initialize=1)  # Partial order of reacton 1
    m.order2 = pe.Param(initialize=1)  # Partial order of reaction 2
    m.QF0 = pe.Param(initialize=1)  # Inlet volumetric flow [L/s]
    C0_Def = {'A': 0.99, 'B': 0.01}
    # Initial concentration of reagents [mol/L]
    m.C0 = pe.Param(m.I, initialize=C0_Def)

    # Inlet molar flow [mol/s]

    def F0_Def(m, i):
        return m.C0[i]*m.QF0
    m.F0 = pe.Param(m.I, initialize=F0_Def)

    # BOOLEAN VARIABLES

    # Unreacted feed in reactor n
    m.YF = pe.BooleanVar(m.N)

    # Existence of recycle flow in unit n
    m.YR = pe.BooleanVar(m.N)

    # Unit operation in n (True if unit n is a CSTR, False if unit n is a bypass)
    m.YP = pe.BooleanVar(m.N)

    # REAL VARIABLES

    # Network Variables
    # Outlet flow rate of the superstructure unit [L/s]
    m.Q = pe.Var(m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

    # Outlet flow rate recycle activation of the superstructure unit [L/s]
    m.QFR = pe.Var(m.N, initialize=0,
                   within=pe.NonNegativeReals, bounds=(0, 10))

    # Molar flow [mol/s]
    m.F = pe.Var(m.I, m.N, initialize=0,
                 within=pe.NonNegativeReals, bounds=(0, 10))

    # Molar flow  recycle activation [mol/s]
    m.FR = pe.Var(m.I, m.N, initialize=0,
                  within=pe.NonNegativeReals, bounds=(0, 10))

    # Reaction rate [mol/(L*s)]
    m.rate = pe.Var(m.I, m.N, initialize=0, within=pe.Reals, bounds=(-10, 10))

    # Reactor volume [L]
    m.V = pe.Var(m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

    # Volume activation [L]
    m.c = pe.Var(m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

    # Splitter Variables
    # Recycle flow rate  [L/s]
    m.QR = pe.Var(initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

    # Product flow rate  [L/s]
    m.QP = pe.Var(initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

    # Recycle molar flow [mol/s]
    m.R = pe.Var(m.I, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

    # Product molar flow [mol/s]
    m.P = pe.Var(m.I, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

    # CONSTRAINTS

    # Unreacted Feed Balances
    # Unreacted feed unit mole balance

    def unreact_mole_rule(m, i, n):
        if n == NT:
            return m.F0[i] + m.FR[i, n] - m.F[i, n] + m.rate[i, n]*m.V[n] == 0
        else:
            return pe.Constraint.Skip

    m.unreact_mole = pe.Constraint(m.I, m.N, rule=unreact_mole_rule)

    # Unreacted feed unit continuity

    def unreact_cont_rule(m, n):
        if n == NT:
            return m.QF0 + m.QFR[n] - m.Q[n] == 0
        else:
            return pe.Constraint.Skip

    m.unreact_cont = pe.Constraint(m.N, rule=unreact_cont_rule)

    # Reactor Balances
    # Reactor mole balance

    def react_mole_rule(m, i, n):
        if n != NT:
            return m.F[i, n+1] + m.FR[i, n] - m.F[i, n] + m.rate[i, n]*m.V[n] == 0
        else:
            return pe.Constraint.Skip

    m.react_mole = pe.Constraint(m.I, m.N, rule=react_mole_rule)

    # Reactor continuity

    def react_cont_rule(m, n):
        if n != NT:
            return m.Q[n+1] + m.QFR[n] - m.Q[n] == 0
        else:
            return pe.Constraint.Skip

    m.react_cont = pe.Constraint(m.N, rule=react_cont_rule)

    # Splitting Point Balances
    # Splitting point mole balance

    def split_mole_rule(m, i):
        return m.F[i, 1] - m.P[i] - m.R[i] == 0

    m.split_mole = pe.Constraint(m.I, rule=split_mole_rule)

    # Splitting point continuity

    def split_cont_rule(m):
        return m.Q[1] - m.QP - m.QR == 0

    m.split_cont = pe.Constraint(rule=split_cont_rule)

    # Splitting point additional constraints

    def split_add_rule(m, i):
        return m.P[i]*m.Q[1] - m.F[i, 1]*m.QP == 0

    m.split_add = pe.Constraint(m.I, rule=split_add_rule)

    # Product Specification

    def prod_spec_rule(m):
        return m.QP*0.95 - m.P['B'] == 0

    m.prod_spec = pe.Constraint(rule=prod_spec_rule)

    # Volume Constraint

    def vol_cons_rule(m, n):
        if n != 1:
            return m.V[n] - m.V[n-1] == 0
        else:
            return pe.Constraint.Skip

    m.vol_cons = pe.Constraint(m.N, rule=vol_cons_rule)

    # YD Disjunction block equation definition

    def build_cstr_equations(disjunct, n):
        m = disjunct.model()

        # Reaction rates calculation
        @disjunct.Constraint()
        def YPD_rate_calc(disjunct):
            return m.rate['A', n]*((m.Q[n])**m.order1)*((m.Q[n])**m.order2)+m.k*((m.F['A', n])**m.order1)*((m.F['B', n])**m.order2) == 0

        # Reaction rate relation
        @disjunct.Constraint()
        def YPD_rate_rel(disjunct):
            return m.rate['B', n] + m.rate['A', n] == 0

        # Volume activation
        @disjunct.Constraint()
        def YPD_vol_act(disjunct):
            return m.c[n] - m.V[n] == 0

    def build_bypass_equations(disjunct, n):
        m = disjunct.model()

        # FR desactivation
        @disjunct.Constraint(m.I)
        def neg_YPD_FR_desact(disjunct, i):
            return m.FR[i, n] == 0

        # Rate desactivation
        @disjunct.Constraint(m.I)
        def neg_YPD_rate_desact(disjunct, i):
            return m.rate[i, n] == 0

        # QFR desactivation
        @disjunct.Constraint()
        def neg_YPD_QFR_desact(disjunct):
            return m.QFR[n] == 0

        @disjunct.Constraint()
        def neg_YPD_vol_desact(disjunct):
            '''
            Volume desactivation function for defining pyomo model
            args:
                disjunct: pyomo block with disjunct to include the constraint
                n: pyomo set with reactor index
            return: 
                return constraint
            '''
            return m.c[n] == 0

    # YR Disjuction block equation definition

    def build_recycle_equations(disjunct, n):
        m = disjunct.model()

        # FR activation
        @disjunct.Constraint(m.I)
        def YRD_FR_act(disjunct, i):
            return m.FR[i, n] - m.R[i] == 0

        # QFR activation
        @disjunct.Constraint()
        def YRD_QFR_act(disjunct):
            return m.QFR[n] - m.QR == 0

    def build_no_recycle_equations(disjunct, n):
        m = disjunct.model()

        # FR desactivation
        @disjunct.Constraint(m.I)
        def neg_YRD_FR_desact(disjunct, i):
            return m.FR[i, n] == 0

        # QFR desactivation
        @disjunct.Constraint()
        def neg_YRD_QFR_desact(disjunct):
            return m.QFR[n] == 0

    # Create disjunction blocks
    m.YR_is_recycle = Disjunct(m.N, rule=build_recycle_equations)
    m.YR_is_not_recycle = Disjunct(m.N, rule=build_no_recycle_equations)

    m.YP_is_cstr = Disjunct(m.N, rule=build_cstr_equations)
    m.YP_is_bypass = Disjunct(m.N, rule=build_bypass_equations)

    # Create disjunctions

    @m.Disjunction(m.N)
    def YP_is_cstr_or_bypass(m, n):
        return [m.YP_is_cstr[n], m.YP_is_bypass[n]]

    @m.Disjunction(m.N)
    def YR_is_recycle_or_not(m, n):
        return [m.YR_is_recycle[n], m.YR_is_not_recycle[n]]

    # Associate Boolean variables with with disjunctions
    for n in m.N:
        m.YP[n].associate_binary_var(m.YP_is_cstr[n].indicator_var)
        m.YR[n].associate_binary_var(m.YR_is_recycle[n].indicator_var)

    # Logic Constraints
    # Unit must be a CSTR to include a recycle

    def cstr_if_recycle_rule(m, n):
        return m.YR[n].implies(m.YP[n])

    m.cstr_if_recycle = pe.LogicalConstraint(m.N, rule=cstr_if_recycle_rule)

    # There is only one unreacted feed

    def one_unreacted_feed_rule(m):
        return pe.exactly(1, m.YF)

    m.one_unreacted_feed = pe.LogicalConstraint(rule=one_unreacted_feed_rule)

    # There is only one recycle stream

    def one_recycle_rule(m):
        return pe.exactly(1, m.YR)

    m.one_recycle = pe.LogicalConstraint(rule=one_recycle_rule)

    # Unit operation in n constraint

    def unit_in_n_rule(m, n):
        if n == 1:
            return m.YP[n].equivalent_to(True)
        else:
            return m.YP[n].equivalent_to(pe.lor(pe.land(~m.YF[n2] for n2 in range(1, n)), m.YF[n]))

    m.unit_in_n = pe.LogicalConstraint(m.N, rule=unit_in_n_rule)

    # OBJECTIVE

    def obj_rule(m):
        return sum(m.c[n] for n in m.N)

    m.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

    m.dsda_status = 'Initialized'

    return m

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
    update_boolean_vars_from_binary(m)
    return m


def solve_with_gdpopt(m, mip='cplex', nlp='conopt', timelimit=10, strategy='LOA', mip_output=False, nlp_output=False):
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
    update_boolean_vars_from_binary(m)
    return m


def external_ref(m, x, logic_expr = None):
    # External variable fix
    ext_var_1 = x[0]
    ext_var_2 = x[1] # TODO generalize for arbitrary number of external variables
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
        # TODO Consider passing temp as argument which could be a model constraint(s) RHS
        # see logic_list in scratch.py

        if temp == True:
            m.YP_is_cstr[n].indicator_var.fix(True)
            m.YP_is_bypass[n].indicator_var.fix(False)
        else:
            m.YP_is_cstr[n].indicator_var.fix(False)
            m.YP_is_bypass[n].indicator_var.fix(True)

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
                              # keepfiles=True,
                              # tmpdir=gams_path,
                              # symbolic_solver_labels=True,
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


def complete_enumeration_external(model_function=build_cstrs, model_args={'NT':5}, reformulation_function=external_ref, nlp='msnlp', timelimit = 10):
    X1 = list(range(1, NT+1))
    X2 = list(range(1, NT+1)) # TODO how to generalize for N external variables?
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
            m_init = initialize_model(m,from_feasible=True)
            m_fixed = reformulation_function(m_init, x)
            m_solved = solve_nlp(m_fixed, nlp=nlp, timelimit=timelimit)

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
    plt.xlabel("YF (Number of reactors)")
    plt.ylabel("YR (Reflux position)")
    plt.show()


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
def evaluate_neighbors(ext_vars, fmin, model_function=build_cstrs, model_args={'NT':5}, reformulation_function=external_ref, nlp='conopt', iter_timelimit=10, tol=0.000001):
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
def do_line_search(start, fmin, direction, model_function=build_cstrs, model_args={'NT':5}, reformulation_function=external_ref, nlp='conopt', optimize=True, min_allowed={}, max_allowed={}, iter_timelimit=10, tol=0.000001):
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

def solve_with_dsda(k='Infinity', model_function=build_cstrs, model_args={'NT':5}, starting_point=[1,1], reformulation_function=external_ref, provide_starting_initialization=True, nlp='conopt', optimize=True, min_allowed={}, max_allowed={}, iter_timelimit=10, tol=0.000001):

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
    NT = 5
    timelimit = 10

    #Complete enumeration
    x, y, objs = complete_enumeration_external(model_function=build_cstrs, model_args={'NT':NT}, nlp='msnlp', timelimit=10)

    # MINLP and GDPopt methods
    m = build_cstrs(NT)
    m_init = initialize_model(m, from_feasible=True)
    m_solved = solve_with_minlp(m_init, transformation='bigm', minlp='baron', timelimit=timelimit, gams_output=False)
    # m_solved = solve_with_gdpopt(m_init, mip='cplex',nlp='conopt', timelimit=timelimit, strategy='LOA', mip_output=False, nlp_output=False)
    print(m_solved.results)
    visualize_cstr_superstructure(m_solved, NT)

    # D-SDA
    k = 'Infinity'
    starting_point = [1,1]
    min_allowed = {i: 1 for i in range(1, len(starting_point)+1)}
    max_allowed = {i: NT for i in range(1, len(starting_point)+1)}

    m_solved, route = solve_with_dsda(k=k, model_function=build_cstrs, model_args={'NT':NT}, starting_point=starting_point, reformulation_function=external_ref, provide_starting_initialization=True, nlp='msnlp', min_allowed=min_allowed, max_allowed=max_allowed, iter_timelimit=10)
    visualize_dsda(points=route, feas_x=x, feas_y=y, objs=objs, k=k)
    print(m_solved.results)
    visualize_cstr_superstructure(m_solved, NT)



    
