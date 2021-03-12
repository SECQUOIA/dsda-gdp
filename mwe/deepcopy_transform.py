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
from pyomo.opt import SolutionStatus
from pyomo.opt import TerminationCondition as tc, SolverResults
import os



NT = 5  # Size of the superstructure (This is an input parameter)

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
    return pe.exactly(1,m.YF)

m.one_unreacted_feed = pe.LogicalConstraint(rule=one_unreacted_feed_rule)

# There is only one recycle stream

def one_recycle_rule(m):
    return pe.exactly(1,m.YR)

m.one_recycle = pe.LogicalConstraint(rule=one_recycle_rule)

# Unit operation in n constraint

def unit_in_n_rule(m, n):
    if n == 1:
        return m.YP[n].equivalent_to(True)
    else:
        return m.YP[n].equivalent_to(pe.lor(pe.land(~m.YF[n2] for n2 in range(1,n)),m.YF[n]))

m.unit_in_n = pe.LogicalConstraint(m.N, rule=unit_in_n_rule)

# OBJECTIVE

def obj_rule(m):
    return sum(m.c[n] for n in m.N)

m.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

x = [5,5]

ext_var_1 =  x[0]
ext_var_2 =  x[1]
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
            
    temp = pe.value(pe.lor(pe.land(~m.YF[n2] for n2 in range(1,n)),m.YF[n]))
    
    if temp == True:
        m.YP_is_cstr[n].indicator_var.fix(True)
        m.YP_is_bypass[n].indicator_var.fix(False)
    else:
        m.YP_is_cstr[n].indicator_var.fix(False)
        m.YP_is_bypass[n].indicator_var.fix(True)

pe.TransformationFactory('core.logical_to_linear').apply_to(m)

m2 = copy.deepcopy(m)

pe.TransformationFactory('core.logical_to_linear').apply_to(m2)

# m3 = pe.TransformationFactory('core.logical_to_linear').create_using(m)
      

   