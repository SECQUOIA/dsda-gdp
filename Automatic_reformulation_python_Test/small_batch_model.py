import pyomo.environ as pe
from pyomo.gdp import (Disjunct, Disjunction)
from pyomo.core.base.misc import display
from pyomo.opt.base.solvers import SolverFactory
from pyomo.core.plugins.transform.logical_to_linear import update_boolean_vars_from_binary
import os

def build_small_batch():
    NK=3
    
    # Model
    m = pe.ConcreteModel()

    # Sets
    m.i = pe.Set(initialize=['a','b'])  # Set of products
    m.j = pe.Set(initialize=['mixer', 'reactor', 'centrifuge']) # Set of stages
    m.k = pe.RangeSet(NK)    # Set of potential number of parallel units

    # Parameters and Scalars

    m.h = pe.Param(initialize=6000) # Horizon time  (available time hrs)
    m.vlow = pe.Param(initialize=250) # Lower bound for size of batch unit
    m.vupp = pe.Param(initialize=2500)  # Upper bound for size of batch unit

    m.q = pe.Param(m.i, initialize={'a':200000, 'b':150000})    # Demand of product i
    m.alpha = pe.Param(m.j, initialize={'mixer':250, 'reactor':500, 'centrifuge':340})  # Cost coefficient for batch units
    m.beta = pe.Param(m.j, initialize={'mixer':0.6, 'reactor':0.6, 'centrifuge':0.6})   # Cost exponent for batch units

    def coeff_init(m,k):
        return pe.log(k)

    m.coeff = pe.Param(m.k, initialize=coeff_init)  # Represent number of parallel units

    s_init = {('a','mixer'):2, ('a','reactor'):3, ('a','centrifuge'):4, 
              ('b','mixer'):4, ('b','reactor'):6, ('b','centrifuge'):3}

    m.s = pe.Param(m.i, m.j, initialize=s_init)   # Size factor for product i in stage j (kg per l)

    t_init = {('a','mixer'):8, ('a','reactor'):20, ('a','centrifuge'):4, 
              ('b','mixer'):10, ('b','reactor'):12, ('b','centrifuge'):3}

    m.t = pe.Param(m.i, m.j, initialize=t_init)   # Processing time of product i in batch j   (hrs)

    # Variables
    m.Y = pe.BooleanVar(m.k, m.j)    # Stage existence
    m.coeffval = pe.Var(m.k, m.j,  within=pe.NonNegativeReals, bounds=(0,pe.log(NK)))  # Activation of coeff
    m.v = pe.Var(m.j, within=pe.NonNegativeReals, bounds=(pe.log(m.vlow),pe.log(m.vupp))) # Volume of stage j 
    m.b = pe.Var(m.i, within=pe.NonNegativeReals) # Batch size of product i
    m.tl = pe.Var(m.i, within=pe.NonNegativeReals)  # Cycle time of product i
    m.n = pe.Var(m.j, within=pe.NonNegativeReals)   # Number of units in parallel stage j

    # Constraints

    # Volume requirement in stage j
    @m.Constraint(m.i, m.j)
    def vol(m,i,j): 
        return m.v[j] >= pe.log(m.s[i,j]) + m.b[i]

    # Cycle time for each product i
    @m.Constraint(m.i, m.j)
    def cycle(m,i,j): 
        return m.n[j] + m.tl[i] >= pe.log(m.t[i,j])

    # Constraint for production time
    @m.Constraint()
    def time(m): 
        return sum(m.q[i]*pe.exp(m.tl[i]-m.b[i]) for i in m.i) <= m.h

    # Relating number of units to 0-1 variables
    @m.Constraint(m.j)
    def units(m,j): 
        return m.n[j] == sum(m.coeffval[k,j] for k in m.k)

    # Only one choice for parallel units is feasible
    @m.LogicalConstraint(m.j)
    def lim(m,j):
        return pe.exactly(1, m.Y[1,j], m.Y[2,j], m.Y[3,j])

    #_______ Disjunction_________

    def build_existence_equations(disjunct, k, j):
        m = disjunct.model()

        # Coeffval activation
        @disjunct.Constraint()
        def coeffval_act(disjunct):
            return m.coeffval[k,j] == m.coeff[k]


    def build_not_existence_equations(disjunct, k, j):
        m = disjunct.model()

        # Coeffval deactivation
        @disjunct.Constraint()
        def coeffval_deact(disjunct):
            return m.coeffval[k,j] == 0

    # Create disjunction block
    m.Y_exists = Disjunct(m.k, m.j, rule=build_existence_equations)
    m.Y_not_exists = Disjunct(m.k, m.j, rule=build_not_existence_equations)

    # Create disjunction

    @m.Disjunction(m.k, m.j)
    def Y_exists_or_not(m, k, j):
        return [m.Y_exists[k,j], m.Y_not_exists[k,j]]

    # Associate Boolean variables with with disjunction
    for k in m.k:
        for j in m.j:
            m.Y[k,j].associate_binary_var(m.Y_exists[k,j].indicator_var)
        
    #____________________________

    # Objective
    def obj_rule(m):
        return sum(m.alpha[j]*(pe.exp(m.n[j] + m.beta[j]*m.v[j])) for j in m.j)

    m.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

    return m