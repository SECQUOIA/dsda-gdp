import pyomo.environ as pe
from pyomo.gdp import (Disjunct, Disjunction)
import networkx as nx
import matplotlib.pyplot as plt
from pyomo.core.base.misc import display
from pyomo.opt.base.solvers import SolverFactory
import os


def gdp_reactors(NT=5, visualize=False):
    # INPUTS
    # NT = 5  # Size of the superstructure (This is an input parameter)
    Initial_Number_Of_Reactors = 1  # Initialization for yf
    Initial_Location_Of_Recycle = 1  # Initialization for yr

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

    # Existence of an unreacted feed in unit n
    def YF_Init(m, n):  # Initialization
        if n == Initial_Number_Of_Reactors:
            return True
        else:
            return False

    m.YF = pe.BooleanVar(m.N, initialize=YF_Init)

    # Existence of recycle flow in unit n
    def YR_Init(m, n):  # Initialization
        if n == Initial_Location_Of_Recycle:
            return True
        else:
            return False
    m.YR = pe.BooleanVar(m.N, initialize=YR_Init)

    # Unit operation in n (True if unit n is a CSTR, False if unit n is a bypass)
    def YP_Init(m, n):  # Initialization
        temp = pe.land(~m.YF[j] for j in range(1, n+1))
        return pe.lor(temp, m.YF[n])
    m.YP = pe.BooleanVar(m.N, initialize=YP_Init)

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
    # Unit must be a CSTR to include a recycle (YR -> YP)

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
        temp = pe.land(~m.YF[j] for j in range(1, n+1))
        return m.YP[n].equivalent_to(pe.lor(temp, m.YF[n]))

    m.unit_in_n = pe.LogicalConstraint(m.N, rule=unit_in_n_rule)

    # OBJECTIVE

    def obj_rule(m):
        return sum(m.c[n] for n in m.N)

    m.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

    # Transform the model using the BigM relaxation

    pe.TransformationFactory('core.logical_to_linear').apply_to(m)
    pe.TransformationFactory('gdp.bigm').apply_to(m)

    dir_path = os.path.dirname(os.path.abspath(__file__))
    gams_path = os.path.join(dir_path, "gamsfiles/")
    if not(os.path.exists(gams_path)):
        print('Directory for automatically generated files ' +
              gams_path + ' does not exist. We will create it')
        os.makedirs(gams_path)

    # SOLVE
    solvername = 'gams'
    opt = SolverFactory(solvername, solver='baron')
    results = opt.solve(m, tee=True,
                        # Uncomment the following lines if you want to save GAMS models
                        # keepfiles=True,
                        # tmpdir=gams_path,
                        # symbolic_solver_labels=True,

                        add_options=[
                            'option reslim = 20;'
                            'option optcr = 0.0;'
                            # Uncomment the following lines to setup IIS computation of BARON through option file
                            # 'GAMS_MODEL.optfile = 1;'
                            # '\n'
                            # '$onecho > baron.opt \n'
                            # 'CompIIS 1 \n'
                            # '$offecho'
                        ])

    print('Objective:', round(pe.value(m.obj), 5))

    if visualize:
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

    return results


if __name__ == "__main__":
    NT = 5
    # Visualization works best (aesthetically) for NT=5
    gdp_reactors(NT, visualize=False)
