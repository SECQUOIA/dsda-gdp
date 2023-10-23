"""
gdp_reactor.py

The code build the series of Continuous Stirred Tank Reactor models for the Generalized Disjunctive Programming superstructure. 
The CSTRs have a autocatalytic reaction A + B -> 2B and the objective function is to minimizes total reactor network volume. 
The code contains the mass and reaction balances, the logic constraints, and the objective function.
The logical constraints are the existance of the recycle flow and bypass flow, and the existance of the CSTRs.
The complete model defines a Generalized Disjunctive Programming (GDP) problem.

This model is to be imported by main_cstr.py where it is solved via differente solution methods (MINLP reformulation, GDP algorithms, and L-DSDA).

References:
[1] Linan, David A., et al. "Optimal design of superstructures for placing units and streams with multiple and ordered available locations. Part I: A new mathematical framework." Computers & Chemical Engineering 137, (2020): 106794.

"""


import os

import matplotlib.pyplot as plt
import networkx as nx
import pyomo.environ as pe
from pyomo.core.base.misc import display
from pyomo.gdp import Disjunct, Disjunction
from pyomo.opt.base.solvers import SolverFactory


def build_cstrs(NT: int = 5) -> pe.ConcreteModel():
    """
    Function that builds CSTR superstructure model of size NT.
    The CSTRs have a autocatalytic reaction A + B -> 2B and minimizes total reactor network volume.
    The optimal solution should yield NT reactors with a recycle before reactor NT.
    Reference: Optimal design of superstructures for placing units and streams with multiple and ordered available loactions.

    Args:
        NT (int): Positive Integer defining the maximum number of CSTRs
    Returns:
        m (pyomo.ConcreteModel): Pyomo GDP model
    """

    # PYOMO MODEL
    m = pe.ConcreteModel(name='gdp_reactors')

    # SETS
    m.I = pe.Set(initialize=['A', 'B'], doc='Set of components')
    m.N = pe.RangeSet(1, NT, doc='Set of units in the superstructure')

    # PARAMETERS
    m.k = pe.Param(
        initialize=2, doc="Kinetic constant [L/(mol*s)]"
    )  # Kinetic constant [L/(mol*s)]
    m.order1 = pe.Param(
        initialize=1, doc="Partial order of reaction 1"
    )  # Partial order of reacton 1
    m.order2 = pe.Param(
        initialize=1, doc="Partial order of reaction 2"
    )  # Partial order of reaction 2
    m.QF0 = pe.Param(
        initialize=1, doc="Inlet volumetric flow [L/s]"
    )  # Inlet volumetric flow [L/s]
    C0_Def = {'A': 0.99, 'B': 0.01}
    # Initial concentration of reagents [mol/L]
    m.C0 = pe.Param(
        m.I, initialize=C0_Def, doc="Initial concentration of reagents [mol/L]"
    )

    # Inlet molar flow [mol/s]

    def F0_Def(m, i):
        """
        This function calculates the molar flow for a reagent 'i' by multiplying
        its initial concentration with the inlet volumetric flow The inlet molar flow of the reagent 'i' in [mol/s].
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model
            i (str): Reagent name

        Returns:
            F0 (float): Inlet molar flow [mol/s]
        """
        return m.C0[i] * m.QF0

    m.F0 = pe.Param(m.I, initialize=F0_Def, doc="Inlet molar flow [mol/s]")

    # BOOLEAN VARIABLES

    # Unreacted feed in reactor n
    m.YF = pe.BooleanVar(m.N, doc="Unreacted feed in reactor n")

    # Existence of recycle flow in unit n
    m.YR = pe.BooleanVar(m.N, doc="Existence of recycle flow in unit n")

    # Unit operation in n (True if unit n is a CSTR, False if unit n is a bypass)
    m.YP = pe.BooleanVar(m.N, doc="Unit operation in n")

    # REAL VARIABLES

    # Network Variables
    # Outlet flow rate of the superstructure unit [L/s]
    m.Q = pe.Var(
        m.N,
        initialize=0,
        within=pe.NonNegativeReals,
        bounds=(0, 10),
        doc="Outlet flow rate of the superstructure unit [L/s]",
    )

    # Outlet flow rate recycle activation of the superstructure unit [L/s]
    m.QFR = pe.Var(
        m.N,
        initialize=0,
        within=pe.NonNegativeReals,
        bounds=(0, 10),
        doc="Outlet flow rate recycle activation of the superstructure unit [L/s]",
    )

    # Molar flow [mol/s]
    m.F = pe.Var(
        m.I,
        m.N,
        initialize=0,
        within=pe.NonNegativeReals,
        bounds=(0, 10),
        doc="Molar flow [mol/s]",
    )

    # Molar flow  recycle activation [mol/s]
    m.FR = pe.Var(
        m.I,
        m.N,
        initialize=0,
        within=pe.NonNegativeReals,
        bounds=(0, 10),
        doc="Molar flow  recycle activation [mol/s]",
    )

    # Reaction rate [mol/(L*s)]
    m.rate = pe.Var(
        m.I,
        m.N,
        initialize=0,
        within=pe.Reals,
        bounds=(-10, 10),
        doc="Reaction rate [mol/(L*s)]",
    )

    # Reactor volume [L]
    m.V = pe.Var(
        m.N,
        initialize=0,
        within=pe.NonNegativeReals,
        bounds=(0, 10),
        doc="Reactor volume [L]",
    )

    # Volume activation [L]
    m.c = pe.Var(
        m.N,
        initialize=0,
        within=pe.NonNegativeReals,
        bounds=(0, 10),
        doc="Volume activation [L]",
    )

    # Splitter Variables
    # Recycle flow rate  [L/s]
    m.QR = pe.Var(
        initialize=0,
        within=pe.NonNegativeReals,
        bounds=(0, 10),
        doc="Recycle flow rate  [L/s]",
    )

    # Product flow rate  [L/s]
    m.QP = pe.Var(
        initialize=0,
        within=pe.NonNegativeReals,
        bounds=(0, 10),
        doc="Product flow rate  [L/s]",
    )

    # Recycle molar flow [mol/s]
    m.R = pe.Var(
        m.I,
        initialize=0,
        within=pe.NonNegativeReals,
        bounds=(0, 10),
        doc="Recycle molar flow [mol/s]",
    )

    # Product molar flow [mol/s]
    m.P = pe.Var(
        m.I,
        initialize=0,
        within=pe.NonNegativeReals,
        bounds=(0, 10),
        doc="Product molar flow [mol/s]",
    )

    # CONSTRAINTS

    # Unreacted Feed Balances
    # Unreacted feed unit mole balance

    def unreact_mole_rule(m, i, n):
        """
        Unreacted feed unite: Partial mole balance, (21.D) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model
            i (str): Reagent name
            n (int): Stage number

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        if n == NT:
            return m.F0[i] + m.FR[i, n] - m.F[i, n] + m.rate[i, n] * m.V[n] == 0
        else:
            return pe.Constraint.Skip

    m.unreact_mole = pe.Constraint(
        m.I,
        m.N,
        rule=unreact_mole_rule,
        doc="Unreacted feed unite: Partial mole balance, (21.D) [1]",
    )

    # Unreacted feed unit continuity

    def unreact_cont_rule(m, n):
        if n == NT:
            """
            Unreacted feed unit: Continuity, (21.E) [1]
            
            Args:
                m (pyomo.ConcreteModel): Pyomo GDP model
                n (int): Stage number
            
            Returns:
                None. The function defines the equations for the given disjunct.
            """
            return m.QF0 + m.QFR[n] - m.Q[n] == 0
        else:
            return pe.Constraint.Skip

    m.unreact_cont = pe.Constraint(
        m.N, rule=unreact_cont_rule, doc="Unreacted feed unit: Continuity, (21.E) [1]"
    )

    # Reactor Balances
    # Reactor mole balance

    def react_mole_rule(m, i, n):
        """
        Reactor sequence: Partial Molar Balance, (21.H) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model
            i (str): Reagent name
            n (int): Stage number
        
        Returns:
            None. The function defines the equations for the given disjunct.
        """
        if n != NT:
            return m.F[i, n + 1] + m.FR[i, n] - m.F[i, n] + m.rate[i, n] * m.V[n] == 0
        else:
            return pe.Constraint.Skip

    m.react_mole = pe.Constraint(
        m.I,
        m.N,
        rule=react_mole_rule,
        doc="Reactor sequence: Partial Molar Balance, (21.H) [1]",
    )

    # Reactor continuity

    def react_cont_rule(m, n):
        """
        Reactor sequence: Continuity, (21.I) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model
            n (int): Stage number

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        if n != NT:
            return m.Q[n + 1] + m.QFR[n] - m.Q[n] == 0
        else:
            return pe.Constraint.Skip

    m.react_cont = pe.Constraint(
        m.N, rule=react_cont_rule, doc="Reactor sequence: Continuity, (21.I) [1]"
    )

    # Splitting Point Balances
    # Splitting point mole balance

    def split_mole_rule(m, i):
        """
        Splitting point: Partial mole balance, (21.L) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model
            i (str): Reagent name

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        return m.F[i, 1] - m.P[i] - m.R[i] == 0

    m.split_mole = pe.Constraint(
        m.I,
        rule=split_mole_rule,
        doc="Splitting point: Partial mole balance, (21.L) [1]",
    )

    # Splitting point continuity

    def split_cont_rule(m):
        """
        Splitting point: continuity, (21.M) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        return m.Q[1] - m.QP - m.QR == 0

    m.split_cont = pe.Constraint(
        rule=split_cont_rule, doc="Splitting point: continuity, (21.M) [1]"
    )

    # Splitting point additional constraints

    def split_add_rule(m, i):
        """
        Splitting point: additional constraints, Molarity constraints over initial and final flows,
        read as an multiplication avoid the numerical complication. (21.N) [1]
        m.P[i]/m.QP =  m.F[i,1]/m.Q[1] (molarity balance)

        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model
            i (str): Reagent name
        
        Returns:
            None. The function defines the equations for the given disjunct.
        """
        return m.P[i] * m.Q[1] - m.F[i, 1] * m.QP == 0

    m.split_add = pe.Constraint(
        m.I,
        rule=split_add_rule,
        doc="Splitting point: additional constraints, Molarity constraints over initial and final flows, (21.N) [1]",
    )

    # Product Specification

    def prod_spec_rule(m):
        """
        Product specification constraint, (21.O) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        return m.QP * 0.95 - m.P['B'] == 0

    m.prod_spec = pe.Constraint(
        rule=prod_spec_rule, doc="Product specification constraint, (21.O) [1]"
    )

    # Volume Constraint

    def vol_cons_rule(m, n):
        """
        Volume constraint, (21.P) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model
            n (int): Stage number

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        if n != 1:
            return m.V[n] - m.V[n - 1] == 0
        else:
            return pe.Constraint.Skip

    m.vol_cons = pe.Constraint(
        m.N, rule=vol_cons_rule, doc="Volume constraint, (21.P) [1]"
    )

    # YD Disjunction block equation definition

    def build_cstr_equations(disjunct, n):
        """
        Build the Continuous Stirred Tank Reactor (CSTR) equations for a given reactor/stage 'n'.
        This function defines the reaction rates, their constraints, and the volume activation
        cost for the model associated with a given disjunct.

        Args:
            m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
            n (int): The stage number for which the equations will be built.

        Return:
            None. The function defines the equations for the given disjunct.
        """
        m = disjunct.model()

        # Reaction rates calculation
        @disjunct.Constraint(doc="Reaction rates calculation")
        def YPD_rate_calc(disjunct):
            """
            Reaction rates calculation. [1]
            -r_A_n = k * (C_A_n)^order1 * (C_B_n)^order2
            -r_B_n = k * (F_A_n/Q_n)^order1 * (F_B_n/Q_n)^order2
            -r_B_n*(Q_n)^order1*(Q_n)^order2 = k * (F_A_n)^order1 * (F_B_n)^order2

            Args:
                m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
                n (int): The stage number for which the equations will be built.

            Return:
                None. The function defines the equations for the given disjunct.
            """
            return (
                m.rate['A', n] * ((m.Q[n]) ** m.order1) * ((m.Q[n]) ** m.order2)
                + m.k * ((m.F['A', n]) ** m.order1) * ((m.F['B', n]) ** m.order2)
                == 0
            )

        # Reaction rate relation
        @disjunct.Constraint(doc="Reaction rate relation")
        def YPD_rate_rel(disjunct):
            """
            The reaction rate relation between A and B.
            
            Args:
                m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
                n (int): The stage number for which the equations will be built.
            
            Return:
                None. The function defines the equations for the given disjunct.
            """
            return m.rate['B', n] + m.rate['A', n] == 0

        # Volume activation
        @disjunct.Constraint(doc="Volume activation")
        def YPD_vol_act(disjunct):
            """
            Volume activation
            
            Args:
                m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
                n (int): The stage number for which the equations will be built.
            
            Return:
                None. The function defines the equations for the given disjunct.
            """
            return m.c[n] - m.V[n] == 0

    def build_bypass_equations(disjunct, n):
        """
        Build the bypass equations for a given stage 'n' of the given CSTR superstructure.

        Args:
            m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
            n (int): The stage number for which the equations will be built.

        Return:
            None. The function defines the equations for the given disjunct.
        """
        m = disjunct.model()

        # FR deactivation
        @disjunct.Constraint(m.I, doc="Molar flow recycle deactivation")
        def neg_YPD_FR_desact(disjunct, i):
            """
            Molar flow recycle deactivation
            
            Args:
                m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
                n (int): The stage number for which the equations will be built.
            
            Return:
                None. The function defines the equations for the given disjunct.
            """
            return m.FR[i, n] == 0

        # Rate deactivation
        @disjunct.Constraint(m.I, doc="Rate deactivation")
        def neg_YPD_rate_desact(disjunct, i):
            """
            Rate deactivation
            
            Args:
                m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
                n (int): The stage number for which the equations will be built.
            
            Return:
                None. The function defines the equations for the given disjunct.
            """
            return m.rate[i, n] == 0

        # QFR deactivation
        @disjunct.Constraint(doc="Outlet flow rate recycle deactivation")
        def neg_YPD_QFR_desact(disjunct):
            """
            Outlet flow rate recycle deactivation
            
            Args:
                m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
                n (int): The stage number for which the equations will be built.
            
            Return:
                None. The function defines the equations for the given disjunct.
            """
            return m.QFR[n] == 0

        @disjunct.Constraint(doc="Volume deactivation")
        def neg_YPD_vol_desact(disjunct):
            """
            Volume deactivation
            
            Args:
                m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
                n (int): The stage number for which the equations will be built.
            
            Return:
                None. The function defines the equations for the given disjunct.
            """
            return m.c[n] == 0

    # YR Disjuction block equation definition

    def build_recycle_equations(disjunct, n):
        """
        Build the recycle equations for a given stage 'n' of the given CSTR superstructure.

        Args:
            m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
            n (int): The stage number for which the equations will be built.

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        m = disjunct.model()

        # FR activation
        @disjunct.Constraint(m.I, doc="Molar flow recycle activation")
        def YRD_FR_act(disjunct, i):
            """
            Molar flow recucle activation
            
            Args:
                m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
                n (int): The stage number for which the equations will be built.

            Return:
                None. The function defines the equations for the given disjunct.
            """
            return m.FR[i, n] - m.R[i] == 0

        # QFR activation
        @disjunct.Constraint(doc="Outlet flow rate recycle activation")
        def YRD_QFR_act(disjunct):
            """
            Outlet flow rate recycle activation via equalizing with the recycle flow rate
            
            Args:
                m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
                n (int): The stage number for which the equations will be built.

            Return:
                None. The function defines the equations for the given disjunct.
            """
            return m.QFR[n] - m.QR == 0

    def build_no_recycle_equations(disjunct, n):
        """
        Build the disjunct for non existence of recycle equations for a given stage 'n' of the given CSTR superstructure.

        arg:
            m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
            n (int): The stage number for which the equations will be built.

        return:
            None. The function defines the equations for the given disjunct.
        """
        m = disjunct.model()

        # FR deactivation
        @disjunct.Constraint(m.I, doc="Molar flow recycle deactivation")
        def neg_YRD_FR_desact(disjunct, i):
            """
            Molar flow recycle deactivation
            
            Args:
                m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
                n (int): The stage number for which the equations will be built.

            Return:
                None. The function defines the equations for the given disjunct.
            """
            return m.FR[i, n] == 0

        # QFR deactivation
        @disjunct.Constraint(doc="Outlet flow rate recycle deactivation")
        def neg_YRD_QFR_desact(disjunct):
            """
            Outlet flow rate recycle deactivation via equalizing with the recycle flow rate
            
            Args:
                m (pyomo.Disjunct): The Pyomo disjunct on which the constraints will be built.
                n (int): The stage number for which the equations will be built.
            
            Return:
                None. The function defines the equations for the given disjunct.
            """
            return m.QFR[n] == 0

    # Create disjunction blocks
    m.YR_is_recycle = Disjunct(
        m.N, rule=build_recycle_equations, doc="Disjunct for recycle equations"
    )
    m.YR_is_not_recycle = Disjunct(
        m.N,
        rule=build_no_recycle_equations,
        doc="Disjunct for non existence of recycle equations",
    )

    m.YP_is_cstr = Disjunct(
        m.N, rule=build_cstr_equations, doc="Disjunct for CSTR equations"
    )
    m.YP_is_bypass = Disjunct(
        m.N, rule=build_bypass_equations, doc="Disjunct for bypass equations"
    )

    # Create disjunctions

    @m.Disjunction(m.N, doc="Disjunction for unit operation in n")
    def YP_is_cstr_or_bypass(m, n):
        """
        Disjunction for YP, (20) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model
            n (int): Stage number

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        return [m.YP_is_cstr[n], m.YP_is_bypass[n]]

    @m.Disjunction(m.N, doc="Disjunction for existence of recycle flow in unit n")
    def YR_is_recycle_or_not(m, n):
        """
        Disjunction for YR, (19.B) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model
            n (int): Stage number
        
        Returns:
            None. The function defines the equations for the given disjunct.
        """
        return [m.YR_is_recycle[n], m.YR_is_not_recycle[n]]

    # Associate Boolean variables with with disjunctions
    for n in m.N:
        m.YP[n].associate_binary_var(m.YP_is_cstr[n].indicator_var)
        m.YR[n].associate_binary_var(m.YR_is_recycle[n].indicator_var)

    # Logic Constraints
    # Unit must be a CSTR to include a recycle

    def cstr_if_recycle_rule(m, n):
        """
        If m.YR[n] is true, then m.YP[n] must also be true.
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model
            n (int): Stage number

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        return m.YR[n].implies(m.YP[n])

    m.cstr_if_recycle = pe.LogicalConstraint(
        m.N,
        rule=cstr_if_recycle_rule,
        doc="If m.YR[n] is true, then m.YP[n] must also be true.",
    )

    # There is only one unreacted feed

    def one_unreacted_feed_rule(m):
        """
        There is only one unreacted feed, (21.B) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        return pe.exactly(1, m.YF)

    m.one_unreacted_feed = pe.LogicalConstraint(
        rule=one_unreacted_feed_rule, doc="There is only one unreacted feed, (21.B) [1]"
    )

    # There is only one recycle stream

    def one_recycle_rule(m):
        """
        There is only one recycle stream, (21.C) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        return pe.exactly(1, m.YR)

    m.one_recycle = pe.LogicalConstraint(
        rule=one_recycle_rule, doc="There is only one recycle stream, (21.C) [1]"
    )

    # Unit operation in n constraint

    def unit_in_n_rule(m, n):
        """
        P[1] is true when n=1, else YP[n] is equivalent to YF[n] or all YF[n] from 1 up to n-1 are false.
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model
            n (int): Stage number

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        if n == 1:
            return m.YP[n].equivalent_to(True)
        else:
            return m.YP[n].equivalent_to(
                pe.lor(pe.land(~m.YF[n2] for n2 in range(1, n)), m.YF[n])
            )

    m.unit_in_n = pe.LogicalConstraint(
        m.N,
        rule=unit_in_n_rule,
        doc="YP[1] is true when n=1, else YP[n] is equivalent to YF[n] or all YF[n] from 1 up to n-1 are false.",
    )

    # OBJECTIVE

    def obj_rule(m):
        """
        Objective function: Total reactor network volume, (21.Q) [1]
        
        Args:
            m (pyomo.ConcreteModel): Pyomo GDP model

        Returns:
            None. The function defines the equations for the given disjunct.
        """
        return sum(m.c[n] for n in m.N)

    m.obj = pe.Objective(
        rule=obj_rule,
        sense=pe.minimize,
        doc="Objective function: Total reactor network volume, (21.Q) [1]",
    )

    return m
