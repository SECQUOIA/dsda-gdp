import pyomo.environ as pe
import networkx as nx
import matplotlib.pyplot as plt
from pyomo.core.base.misc import display
from pyomo.opt.base.solvers import SolverFactory
import os


def master_iplc(NT=5):
    m = pe.ConcreteModel(name='external_initialization')
    m.ext_var_1 = pe.Var(within=pe.NonNegativeIntegers, bounds=(0, NT))
    m.ext_var_2 = pe.Var(within=pe.NonNegativeIntegers, bounds=(0, NT))

    @m.Constraint()
    def feas_constraint(m):
        return m.ext_var_2 - m.ext_var_1 <= 0

    m.z_ext = pe.Var()

    @m.Constraint()
    def z_cons(m):
        return m.z_ext == 1

    def obj_rule(m):
        return m.z_ext

    m.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

    solvername = 'gams'
    opt = SolverFactory(solvername, solver='cplex')
    results = opt.solve(m, tee=True,
                        add_options=[
                            'option reslim = 20;'
                            'option optcr = 0.0;'
                        ])

    return {1: pe.value(m.ext_var_1), 2: pe.value(m.ext_var_2)}


def cstr_model(x, NT=5):
    # INPUTS
    # NT = 5  # Size of the superstructure (This is an input parameter)

    # PYOMO MODEL
    m = pe.ConcreteModel(name='cstr_model')

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

    # Compute y = y(x)

    ext_var_1 = x[1]
    ext_var_2 = x[2]

    # Existence of an unreacted feed in unit n

    def yf_Init(m, n):  # Initialization
        if n == ext_var_1:
            return 1
        else:
            return 0

    m.yf = pe.Param(m.N, initialize=yf_Init)

    # Existence of recycle flow in unit n

    def yr_Init(m, n):  # Initialization
        if n == ext_var_2:
            return 1
        else:
            return 0

    m.yr = pe.Param(m.N, initialize=yr_Init)

    # Unit operation in n (1 if unit n is a CSTR, 0 otherwise)

    def yp_Init(m, n):  # Initialization
        return 1 - (sum(m.yf[n2] for n2 in m.N if n2 <= n) - m.yf[n])

    m.yp = pe.Param(m.N, initialize=yp_Init)

    # REAL VARIABLES

    # Network Variables
    # Outlet flow rate of the superstructure unit [L/s]
    m.Q = pe.Var(m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

    # Molar flow [mol/s]
    m.F = pe.Var(m.I, m.N, initialize=0,
                 within=pe.NonNegativeReals, bounds=(0, 10))

    # Reaction rate [mol/(L*s)]
    m.rate = pe.Var(m.I, m.N, initialize=0, within=pe.Reals, bounds=(-10, 10))

    # Reactor volume [L]
    m.V = pe.Var(m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

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
            return m.F0[i] + m.yr[n]*m.R[i] - m.F[i, n] + m.yp[n]*m.rate[i, n]*m.V[n] == 0
        else:
            return pe.Constraint.Skip

    m.unreact_mole = pe.Constraint(m.I, m.N, rule=unreact_mole_rule)

    # Unreacted feed unit continuity

    def unreact_cont_rule(m, n):
        if n == NT:
            return m.QF0 + m.yr[n]*m.QR - m.Q[n] == 0
        else:
            return pe.Constraint.Skip

    m.unreact_cont = pe.Constraint(m.N, rule=unreact_cont_rule)

    # Rates
    # Reaction rates calculation

    def rate_calc_rule(m, n):
        return m.rate['A', n]*((m.Q[n])**m.order1)*((m.Q[n])**m.order2)+m.k*((m.F['A', n])**m.order1)*((m.F['B', n])**m.order2) == 0

    m.rate_calc = pe.Constraint(m.N, rule=rate_calc_rule)

    # Reaction rate relation

    def rate_rel_rule(m, n):
        return m.rate['B', n] + m.rate['A', n] == 0

    m.rate_rel = pe.Constraint(m.N, rule=rate_rel_rule)

    # Reactor Balances
    # Reactor mole balance

    def react_mole_rule(m, i, n):
        if n != NT:
            return m.F[i, n+1] + m.yr[n]*m.R[i] - m.F[i, n] + m.yp[n]*m.rate[i, n]*m.V[n] == 0
        else:
            return pe.Constraint.Skip

    m.react_mole = pe.Constraint(m.I, m.N, rule=react_mole_rule)

    # Reactor continuity

    def react_cont_rule(m, n):
        if n != NT:
            return m.Q[n+1] + m.yr[n]*m.QR - m.Q[n] == 0
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

    # OBJECTIVE

    def obj_rule(m):
        return sum(m.yp[n]*m.V[n] for n in m.N)

    m.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

    dir_path = os.path.dirname(os.path.abspath(__file__))
    gams_path = os.path.join(dir_path, "gamsfiles/")
    if not(os.path.exists(gams_path)):
        print('Directory for automatically generated files ' +
              gams_path + ' does not exist. We will create it')
        os.makedirs(gams_path)

    # SOLVE
    solvername = 'gams'
    opt = SolverFactory(solvername, solver='msnlp')
    results = opt.solve(m, tee=False,
                        # Uncomment the following lines if you want to save GAMS models
                        # keepfiles=True,
                        # tmpdir=gams_path,
                        # symbolic_solver_labels=True,

                        add_options=[
                            'option reslim = 100;'
                            'option optcr = 0.0;'
                            # Uncomment the following lines to setup IIS computation of BARON through option file
                            # 'GAMS_MODEL.optfile = 1;'
                            # '\n'
                            # '$onecho > baron.opt \n'
                            # 'CompIIS 1 \n'
                            # '$offecho'
                        ])

    return round(pe.value(m.obj), 5)


if __name__ == "__main__":
    NT = 5
    X1, X2, x, aux, aux2 = [], [], {}, [], 1

    for i in range(1, NT+1):
        X1.append(i)
        aux.append(i)
        X2.append(aux2)

    for i in range(NT-1):
        aux.pop(0)
        aux2 += 1
        for j in aux:
            X1.append(j)
            X2.append(aux2)
    print('=============================')
    print('%6s %6s %12s' % ('x1', 'x2', 'Objective'))
    print('-----------------------------')

    for i in range(len(X1)):
        x = {1: X1[i], 2: X2[i]}
        print('%6s %6s %12s' % (X1[i], X2[i], cstr_model(x, NT)))

    print('=============================')
