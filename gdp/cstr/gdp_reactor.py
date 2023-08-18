import os  # Provides operating system-dependent functionality

import matplotlib.pyplot as plt  # Library for creating plots and visualizations
import networkx as nx  # Python package for working with complex networks
import pyomo.environ as pe  # Main entry point for using the Pyomo optimization modeling language
from pyomo.core.base.misc import display  # Function for printing Pyomo model components
from pyomo.gdp import (
    Disjunct,
    Disjunction,
)  # Classes for modeling disjunctive constraints in Pyomo's GDP framework
from pyomo.opt.base.solvers import (
    SolverFactory,
)  # Class for creating solver objects for optimization


def build_cstrs(NT: int = 5) -> pe.ConcreteModel():
    """
    Function that builds CSTR superstructure model of size NT.
    The CSTRs have a single 1st order reaction A -> B and minimizes (TODO Check)
    total reactor volume. The optimal solution should yield NT reactors with a recycle before reactor NT.
    Reference: Paper Linhan 1. TODO Correct reference

    Args:
        NT: int. Positive Integer defining the maximum number of CSTRs
    Returns:
        m = Pyomo GDP model
    """

    # PYOMO MODEL
    m = pe.ConcreteModel(name='gdp_reactors')

    # SETS
    m.I = pe.Set(
        initialize=['A', 'B'], doc='Set of components'
    )  # Set of components 'A' and 'B'
    m.N = pe.RangeSet(
        1, NT, doc='Set of units in the superstructure'
    )  # Set of units in the superstructure, ranging from 1 to NT

    # PARAMETERS
    m.k = pe.Param(initialize=2)  # Kinetic constant [L/(mol*s)]
    m.order1 = pe.Param(initialize=1)  # Partial order of reaction 1
    m.order2 = pe.Param(initialize=1)  # Partial order of reaction 2
    m.QF0 = pe.Param(initialize=1)  # Inlet volumetric flow [L/s]
    C0_Def = {'A': 0.99, 'B': 0.01}
    # Initial concentration of reagents [mol/L]
    m.C0 = pe.Param(m.I, initialize=C0_Def)

    # Inlet molar flow [mol/s]
    def F0_Def(m, i):
        # Function to calculate the inlet molar flow for each component i in the set m.I.
        # The calculation is done by multiplying the initial concentration of the component m.C0[i]
        # with the inlet volumetric flow m.QF0. The result is the inlet molar flow for component i.

        # Args:
        #     m: Pyomo model object
        #     i: Component index

        # Returns:
        #     Inlet molar flow for component i

        return m.C0[i] * m.QF0

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
    m.QFR = pe.Var(m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

    # Molar flow [mol/s]
    m.F = pe.Var(m.I, m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

    # Molar flow  recycle activation [mol/s]
    m.FR = pe.Var(m.I, m.N, initialize=0, within=pe.NonNegativeReals, bounds=(0, 10))

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
    # Unreacted feed unit mole balance: Ensures mole balance for the unreacted feed unit

    def unreact_mole_rule(m, i, n):
        if n == NT:
            return m.F0[i] + m.FR[i, n] - m.F[i, n] + m.rate[i, n] * m.V[n] == 0
        else:
            return pe.Constraint.Skip

    m.unreact_mole = pe.Constraint(m.I, m.N, rule=unreact_mole_rule)

    # Unreacted feed unit continuity: Ensures continuity of the unreacted feed unit

    def unreact_cont_rule(m, n):
        if n == NT:
            return m.QF0 + m.QFR[n] - m.Q[n] == 0
        else:
            return pe.Constraint.Skip

    m.unreact_cont = pe.Constraint(m.N, rule=unreact_cont_rule)

    # Reactor Balances
    # Reactor mole balance: Ensures mass balance for each component 'i' in each superstructure unit 'n'
    def react_mole_rule(m, i, n):
        if n != NT:
            return m.F[i, n + 1] + m.FR[i, n] - m.F[i, n] + m.rate[i, n] * m.V[n] == 0
        else:
            return pe.Constraint.Skip

    m.react_mole = pe.Constraint(m.I, m.N, rule=react_mole_rule)

    # Reactor continuity: Ensures continuity of the reactors

    def react_cont_rule(m, n):
        if n != NT:
            return m.Q[n + 1] + m.QFR[n] - m.Q[n] == 0
        else:
            return pe.Constraint.Skip  # Skip if n = NT

    m.react_cont = pe.Constraint(m.N, rule=react_cont_rule)

    # Splitting Point Balances
    # Splitting point mole balance: Ensures mole balance at the splitting point

    def split_mole_rule(m, i):
        return m.F[i, 1] - m.P[i] - m.R[i] == 0

    m.split_mole = pe.Constraint(m.I, rule=split_mole_rule)

    # Splitting point continuity

    def split_cont_rule(m):
        return m.Q[1] - m.QP - m.QR == 0

    m.split_cont = pe.Constraint(rule=split_cont_rule)

    # Splitting point additional constraints: Additional constraints at the splitting point

    def split_add_rule(m, i):
        return m.P[i] * m.Q[1] - m.F[i, 1] * m.QP == 0

    m.split_add = pe.Constraint(m.I, rule=split_add_rule)

    # Product Specification

    def prod_spec_rule(m):
        # This function defines the product specification rule.
        # It sets the product 'B' to be equal to 95% of QP.
        return m.QP * 0.95 - m.P['B'] == 0

    m.prod_spec = pe.Constraint(
        rule=prod_spec_rule
    )  # Sets the product specification rule as a constraint in the model.

    # Volume Constraint

    def vol_cons_rule(m, n):
        # This function defines the volume constraint rule.
        # It ensures the volume at each time step is equal to the volume at the previous time step.
        # However, for the first step, the constraint does not apply.
        if n != 1:
            return m.V[n] - m.V[n - 1] == 0
        else:
            return (
                pe.Constraint.Skip
            )  # If n equals 1, the function returns 'Constraint.Skip' which means no constraint is added for this case.

    m.vol_cons = pe.Constraint(
        m.N, rule=vol_cons_rule
    )  # Sets the volume constraint rule as a constraint in the model.

    # YD Disjunction block equation definition

    def build_cstr_equations(disjunct, n):
        # This function builds the equations for the YD disjunction block.
        m = disjunct.model()  # Accesses the underlying model from the disjunct.

        # Reaction rates calculation
        @disjunct.Constraint()
        def YPD_rate_calc(disjunct):
            # Defines the reaction rate equation based on the power-law rate equation.
            return (
                m.rate['A', n] * ((m.Q[n]) ** m.order1) * ((m.Q[n]) ** m.order2)
                + m.k * ((m.F['A', n]) ** m.order1) * ((m.F['B', n]) ** m.order2)
                == 0
            )

        # Reaction rate relation
        @disjunct.Constraint()
        def YPD_rate_rel(disjunct):
            # Defines the reaction rate relation, stating the rate of 'A' plus the rate of 'B' should be zero.
            return m.rate['B', n] + m.rate['A', n] == 0

        # Volume activation
        @disjunct.Constraint()
        def YPD_vol_act(disjunct):
            # Defines the volume activation constraint that relates the volume 'V' to the activation variable 'c'.
            return m.c[n] - m.V[n] == 0

    def build_bypass_equations(disjunct, n):
        # This function builds the equations for the bypass disjunction block.
        m = disjunct.model()  # Accesses the underlying model from the disjunct.

        # FR deactivation
        @disjunct.Constraint(m.I)
        def neg_YPD_FR_desact(disjunct, i):
            # Defines the constraint for deactivating the flow rate 'FR' in the bypass scenario.
            return m.FR[i, n] == 0

        # Rate deactivation
        @disjunct.Constraint(m.I)
        def neg_YPD_rate_desact(disjunct, i):
            # Defines the constraint for deactivating the reaction rate in the bypass scenario.
            return m.rate[i, n] == 0

        # QFR deactivation
        @disjunct.Constraint()
        def neg_YPD_QFR_desact(disjunct):
            # Defines the constraint for deactivating the flow rate 'QFR' in the bypass scenario.
            return m.QFR[n] == 0

        @disjunct.Constraint()
        def neg_YPD_vol_desact(disjunct):
            '''
            Volume deactivation function for defining pyomo model
            args:
                disjunct: pyomo block with disjunct to include the constraint
                n: pyomo set with reactor index
            return:
                return constraint
            '''
            # Defines the constraint for deactivating the volume 'c' in the bypass scenario.
            return m.c[n] == 0

    # YR Disjuction block equation definition

    def build_recycle_equations(disjunct, n):
        # This function builds the equations for the YR disjunction block.
        m = disjunct.model()  # Accesses the underlying model from the disjunct.

        # FR activation
        @disjunct.Constraint(m.I)
        def YRD_FR_act(disjunct, i):
            # Defines the constraint for activating the flow rate 'FR' in the recycle scenario.
            return m.FR[i, n] - m.R[i] == 0

        # QFR activation
        @disjunct.Constraint()
        def YRD_QFR_act(disjunct):
            # Defines the constraint for activating the flow rate 'QFR' in the recycle scenario.
            return m.QFR[n] - m.QR == 0

    def build_no_recycle_equations(disjunct, n):
        m = disjunct.model()  # Accesses the underlying model from the disjunct.

        # FR deactivation
        @disjunct.Constraint(m.I)
        def neg_YRD_FR_desact(disjunct, i):
            # Defines the constraint for deactivating the flow rate 'FR' in the non-recycle scenario.
            return m.FR[i, n] == 0

        # QFR deactivation
        @disjunct.Constraint()
        def neg_YRD_QFR_desact(disjunct):
            # Defines the constraint for deactivating the flow rate 'QFR' in the non-recycle scenario.
            return m.QFR[n] == 0

    # Create disjunction blocks
    m.YR_is_recycle = Disjunct(
        m.N, rule=build_recycle_equations
    )  # Creates a disjunct that represents the equations for the recycle scenario.
    m.YR_is_not_recycle = Disjunct(
        m.N, rule=build_no_recycle_equations
    )  # Creates a disjunct that represents the equations for the non-recycle scenario.

    m.YP_is_cstr = Disjunct(
        m.N, rule=build_cstr_equations
    )  # Creates a disjunct that represents the equations for the reactor operation.
    m.YP_is_bypass = Disjunct(
        m.N, rule=build_bypass_equations
    )  # Creates a disjunct that represents the equations for the bypass operation.

    # Create disjunctions
    # These disjunctions model the logic that one and only one scenario can be active at each time step.

    @m.Disjunction(m.N)
    def YP_is_cstr_or_bypass(m, n):
        return [m.YP_is_cstr[n], m.YP_is_bypass[n]]

    @m.Disjunction(m.N)
    def YR_is_recycle_or_not(m, n):
        return [m.YR_is_recycle[n], m.YR_is_not_recycle[n]]

    # Associate Boolean variables with with disjunctions
    # These associations link the binary variables with their corresponding disjunctions.
    for n in m.N:
        m.YP[n].associate_binary_var(m.YP_is_cstr[n].indicator_var)
        m.YR[n].associate_binary_var(m.YR_is_recycle[n].indicator_var)

    # Logic Constraints
    # These constraints capture the additional logic requirements that cannot be expressed by the disjunctions alone.

    # Unit must be a CSTR to include a recycle
    def cstr_if_recycle_rule(m, n):
        return m.YR[n].implies(
            m.YP[n]
        )  # If the recycle scenario is active, then the reactor operation must also be active.

    m.cstr_if_recycle = pe.LogicalConstraint(m.N, rule=cstr_if_recycle_rule)

    # There is only one unreacted feed
    def one_unreacted_feed_rule(m):
        return pe.exactly(1, m.YF)  # Exactly one of the feed variables can be active.

    m.one_unreacted_feed = pe.LogicalConstraint(rule=one_unreacted_feed_rule)

    # There is only one recycle stream
    def one_recycle_rule(m):
        return pe.exactly(
            1, m.YR
        )  # Exactly one of the recycle variables can be active.

    m.one_recycle = pe.LogicalConstraint(rule=one_recycle_rule)

    # Unit operation in n constraint
    def unit_in_n_rule(m, n):
        if n == 1:
            return m.YP[n].equivalent_to(
                True
            )  # The reactor operation is always active at the first step.
        else:
            # For subsequent steps, the reactor operation is active if either of the following conditions is true:
            # - None of the feed variables in the previous steps is active.
            # - The feed variable at the current step is active.
            return m.YP[n].equivalent_to(
                pe.lor(pe.land(~m.YF[n2] for n2 in range(1, n)), m.YF[n])
            )

    m.unit_in_n = pe.LogicalConstraint(m.N, rule=unit_in_n_rule)

    # OBJECTIVE
    # The objective function is to minimize the sum of the volumes 'c' over all time steps.
    def obj_rule(m):
        return sum(m.c[n] for n in m.N)

    m.obj = pe.Objective(
        rule=obj_rule, sense=pe.minimize
    )  # Defines the objective function.

    return m  # Returns the fully defined Pyomo model.
