import pyomo.environ as pe
from pyomo.gdp import (Disjunct, Disjunction)
import networkx as nx
import matplotlib.pyplot as plt
import copy
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


def build_cstrs(NT=5):
    # INPUTS
    # NT = 5  # Size of the superstructure (This is an input parameter)

    # PYOMO MODEL
    m = pe.ConcreteModel(name='gdp_reactors')

    # SETS
    m.I = pe.Set(initialize=['A', 'B'])  # Set of components
    m.N = pe.RangeSet(1, NT)  # Set of units in the superstructure

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

    return m


def solve_with_minlp(m, transformation='bigm', minlp='baron', timelimit=10):

    # Transformation step
    pe.TransformationFactory('core.logical_to_linear').apply_to(m)
    transformation_string = 'gdp.' + transformation
    pe.TransformationFactory(transformation_string).apply_to(m)

    # Solution step
    dir_path = os.path.dirname(os.path.abspath(__file__))
    gams_path = os.path.join(dir_path, "gamsfiles/")
    if not(os.path.exists(gams_path)):
        print('Directory for automatically generated files ' +
              gams_path + ' does not exist. We will create it')
        os.makedirs(gams_path)

    solvername = 'gams'
    opt = SolverFactory(solvername, solver=minlp)
    m.results = opt.solve(m, tee=True,
                          # Uncomment the following lines if you want to save GAMS models
                          # keepfiles=True,
                          # tmpdir=gams_path,
                          # symbolic_solver_labels=True,
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


def solve_with_gdpopt(m, mip='cplex', nlp='ipopth', minlp='bonmin', timelimit=10):
    """
    Function documentation
    """

    # Transformation step
    pe.TransformationFactory('core.logical_to_linear').apply_to(m)

    # Solution step
    dir_path = os.path.dirname(os.path.abspath(__file__))
    gams_path = os.path.join(dir_path, "gamsfiles/")
    if not(os.path.exists(gams_path)):
        print('Directory for automatically generated files ' +
              gams_path + ' does not exist. We will create it')
        os.makedirs(gams_path)
    solvername = 'gdpopt'
    opt = SolverFactory(solvername)
    m.results = opt.solve(m, tee=True,
                          strategy='LOA',
                          # strategy='GLOA',
                          time_limit=timelimit,
                          mip_solver='gams',
                          mip_solver_args=dict(solver=mip, warmstart=True,
                                               keepfiles=True,
                                               tmpdir=gams_path,
                                               symbolic_solver_labels=True
                                               ),
                          nlp_solver='gams',
                          nlp_solver_args=dict(solver=nlp, warmstart=True,
                                            #    keepfiles=True,
                                            #    tmpdir=gams_path,
                                            #    symbolic_solver_labels=True
                                               ),
                          minlp_solver='gams',
                          minlp_solver_args=dict(solver=minlp, warmstart=True, tee=True,
                                                #  keepfiles=True,
                                                #  tmpdir=gams_path,
                                                #  symbolic_solver_labels=True
                                                 ),
                          subproblem_presolve=False,
                          # init_strategy='no_init',
                          set_cover_iterlim=1,
                          iterlim=1
                          # calc_disjunctive_bounds=True
                          )
    update_boolean_vars_from_binary(m)
    return m


def external_ref(m, x):
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

    return m


def solve_nlp(m, nlp='msnlp'):
    m.status = 'ok'
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
                              add_options=[
                                  'option reslim = 10;'
                                  'option optcr = 0.0;'
                                  # Uncomment the following lines to setup IIS computation of BARON through option file
                                  # 'GAMS_MODEL.optfile = 1;'
                                  # '\n'
                                  # '$onecho > baron.opt \n'
                                  # 'CompIIS 1 \n'
                                  # '$offecho'
                                  # 'display(execError);'
                              ])

        wts = StoreSpec.value()
        to_json(m, fname="test.json", human_read=True, wts=wts)
        Q_init, QFR_init, F_init, FR_init, rate_init, V_init, c_init,  R_init, P_init = {
        }, {}, {}, {}, {}, {}, {}, {}, {}

        QR_init = pe.value(m.QR)
        QP_init = pe.value(m.QP)
        for n in m.N:
            Q_init[n] = pe.value(m.Q[n])
            QFR_init[n] = pe.value(m.QFR[n])
            V_init[n] = pe.value(m.V[n])
            c_init[n] = pe.value(m.c[n])

            for i in m.I:
                F_init[i, n] = pe.value(m.F[i, n])
                FR_init[i, n] = pe.value(m.FR[i, n])
                rate_init[i, n] = pe.value(m.rate[i, n])

        for i in m.I:
            R_init[i] = pe.value(m.R[i])
            P_init[i] = pe.value(m.P[i])

        m.initialization = {'Q': Q_init, 'QFR': QFR_init, 'F': F_init, 'FR': FR_init, 'rate': rate_init,
                            'V': V_init, 'c': c_init, 'QR': QR_init, 'QP': QP_init, 'R': R_init, 'P': P_init}
        return m

    except InfeasibleConstraintException:
        nlp_result = MasterProblemResult()
        nlp_result.feasible = False
        nlp_result.pyomo_results = SolverResults()
        nlp_result.pyomo_results.solver.termination_condition = tc.error
        m.status = 'inf'
        #print('Try an infeasible')

        return m


def complete_enumeration(m, NT=5, nlp='msnlp'):
    X1 = list(range(1, NT+1))
    print()

    print('=============================')
    print('%6s %6s %12s' % ('x1', 'x2', 'Objective'))
    print('-----------------------------')

    #m2 = copy.copy(m)
    for i in range(len(X1)):
        for j in range(len(X1)):
            #m = build_cstrs(NT)
            #m2 = copy.deepcopy(m)
            # m2=m.copy()
            x = [X1[i], X1[j]]
            m_fixed = external_ref(m, x)
            m_solved = solve_nlp(m_fixed, nlp=nlp)

            if m_solved.status == 'ok':
                print('%6s %6s %12s' %
                      (X1[i], X1[j], round(pe.value(m_solved.obj), 5)))
            else:
                print('%6s %6s %12s' % (X1[i], X1[j], 'Infeasible'))

    print('=============================')


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

# start is type list and stands for the actual point
# neighborhood is type dict and is the output of a k-Neighborhood function
# newbors or new_newbors is type  dict starting in 0 with neighbor of a given point
# Neighbor 0 is the actual point


def my_neighbors(start, neighborhood, optimize=True, min_allowed={}, max_allowed={}):
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


def initialize_cstr(m):
    wts = StoreSpec.value()
    from_json(m, fname = 'test.json', wts=wts)
    # for v in m.component_objects(pe.Var, active=True):
    #     print(v.lb)
    #     m.del_component(v)
        #m.add_component(str(v), pe.Var(within=pe.NonNegativeReals, ))

    #m.Q = pe.Var(m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))
    #m.QFR = pe.Var(m.N, initialize=0,within=pe.NonNegativeReals, bounds=(0, 10))
    #m.F = pe.Var(m.I, m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))
    #m.FR = pe.Var(m.I, m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))
    #m.rate = pe.Var(m.I, m.N, initialize=0, within=pe.Reals, bounds=(-10, 10))
    #m.V = pe.Var(m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))
    #m.c = pe.Var(m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))
    #m.QR = pe.Var(initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))
    #m.QP = pe.Var(initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))
    #m.R = pe.Var(m.I, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))
    #m.P = pe.Var(m.I, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))


def visualize_solution(m, NT):
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
    # NT = 5
    # timelimit = 10
    # m = build_cstrs(NT)
    # # initialize_cstr(m)
    # m = external_ref(m,[5,5])
    # solve_nlp(m)

    # m2 = build_cstrs(NT)
    # initialize_cstr(m2)

    # m2.pprint()
    neighborhood_k_eq_inf(10)
    # complete_enumeration(m, NT, nlp='msnlp')
    # m_solved = solve_with_minlp(m, transformation='bigm', minlp='baron', timelimit=timelimit)
    # m_solved = solve_with_gdpopt(m, mip='cplex',nlp='ipopth',minlp='dicopt', timelimit=timelimit)
    # print(m_solved.results)
    # visualize_solution(m_solved,NT)
